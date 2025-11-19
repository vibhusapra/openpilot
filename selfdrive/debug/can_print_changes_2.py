#!/usr/bin/env python3
import argparse
import binascii
from collections import defaultdict

from openpilot.selfdrive.debug.can_table import can_table
from openpilot.tools.lib.logreader import LogIterable, LogReader

RED = '\033[91m'
CLEAR = '\033[0m'

def collect_stable_bits(msgs, bus):
  """
  Collect bits that remain stable (don't change) throughout the message set.
  Returns: dict mapping address -> stable bit mask
  """
  dat = defaultdict(lambda: None)
  low_to_high = defaultdict(int)  # Bits ever seen as 1
  high_to_low = defaultdict(int)  # Bits ever seen as 0

  for x in msgs:
    if x.which() != 'can':
      continue

    for y in x.can:
      if y.src == bus:
        if dat[y.address] is None:
          dat[y.address] = y.dat

        i = int.from_bytes(y.dat, byteorder='big')
        low_to_high[y.address] |= i      # Accumulate bits seen as 1
        high_to_low[y.address] |= ~i     # Accumulate bits seen as 0

  # Stable bits are those that never transitioned (seen as only 0 OR only 1)
  stable_bits = {}
  stable_values = {}

  for addr in dat.keys():
    # Bits that were seen as both 0 and 1 are unstable (counters/checksums)
    unstable = low_to_high[addr] & high_to_low[addr]

    # Create mask with proper byte length
    num_bytes = len(dat[addr])
    all_bits_mask = (1 << (num_bytes * 8)) - 1

    # Stable bits are those NOT in the unstable set
    stable = all_bits_mask & ~unstable

    stable_bits[addr] = stable
    stable_values[addr] = low_to_high[addr] & stable  # The actual values of stable bits

  return stable_bits, stable_values, dat


def compare_routes(bus, init_msgs, comp_msgs, table=False):
  """
  Compare two routes by finding bits that are stable within each route
  but different between routes.
  """
  print("Analyzing route 1 (baseline)...")
  stable_bits_1, stable_values_1, dat_1 = collect_stable_bits(init_msgs, bus)

  print("Analyzing route 2 (comparison)...")
  stable_bits_2, stable_values_2, dat_2 = collect_stable_bits(comp_msgs, bus)

  print("\n" + "="*80)
  print("STABLE BIT CHANGES (filters out counters and checksums)")
  print("="*80 + "\n")

  # Find all addresses present in either route
  all_addrs = set(stable_bits_1.keys()) | set(stable_bits_2.keys())

  tables = ""
  changes_found = False

  for addr in sorted(all_addrs):
    stable_1 = stable_bits_1.get(addr, 0)
    stable_2 = stable_bits_2.get(addr, 0)
    value_1 = stable_values_1.get(addr, 0)
    value_2 = stable_values_2.get(addr, 0)

    # Only compare bits that are stable in BOTH routes
    commonly_stable = stable_1 & stable_2

    # Find bits that changed value between routes
    changed_bits = commonly_stable & (value_1 ^ value_2)

    # Also show new addresses or bits that became stable
    new_stable_bits = stable_2 & ~stable_1  # Bits stable in route2 but not route1
    new_stable_with_value = new_stable_bits & value_2

    if changed_bits == 0 and new_stable_with_value == 0:
      continue

    changes_found = True

    # Get byte length from either route
    num_bytes = len(dat_1.get(addr, dat_2.get(addr, b'\x00')))

    # Print changed bits (stable in both routes)
    if changed_bits != 0:
      b = changed_bits.to_bytes(num_bytes, byteorder='big')
      byts = ''.join([(c if c == '0' else f'{RED}{c}{CLEAR}') for c in str(binascii.hexlify(b))[2:-1]])
      header = f"{hex(addr).ljust(6)}({str(addr).ljust(4)})"
      print(f"{header} CHANGED: {byts}")

      value_1_bytes = value_1.to_bytes(num_bytes, byteorder='big')
      value_2_bytes = value_2.to_bytes(num_bytes, byteorder='big')
      print(f"       Route 1: {binascii.hexlify(value_1_bytes)}")
      print(f"       Route 2: {binascii.hexlify(value_2_bytes)}")
      print()

      tables += f"{header} CHANGED\n"
      tables += can_table(b) + "\n\n"

    # Print new stable bits
    if new_stable_with_value != 0:
      b = new_stable_with_value.to_bytes(num_bytes, byteorder='big')
      byts = ''.join([(c if c == '0' else f'{RED}{c}{CLEAR}') for c in str(binascii.hexlify(b))[2:-1]])
      header = f"{hex(addr).ljust(6)}({str(addr).ljust(4)})"
      print(f"{header} NEW:     {byts}")

      value_2_bytes = value_2.to_bytes(num_bytes, byteorder='big')
      print(f"       Route 2: {binascii.hexlify(value_2_bytes)}")
      print()

      tables += f"{header} NEW\n"
      tables += can_table(b) + "\n\n"

  if not changes_found:
    print("No stable bit changes found between routes.")

  if table and changes_found:
    print("\n" + "="*80)
    print("DETAILED BIT TABLES")
    print("="*80 + "\n")
    print(tables)


if __name__ == "__main__":
  desc = """Compares two routes and prints bits that are stable within each route but different between routes.
  This filters out counters and checksums that change within a route, revealing only meaningful signal changes.

  Usage example for finding Pilot Assist signals:
    Route 1: Drive without Pilot Assist engaged
    Route 2: Drive with Pilot Assist engaged

    python can_print_changes_2.py --bus 0 route1_segment route2_segment
  """
  parser = argparse.ArgumentParser(description=desc,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--bus", type=int, help="CAN bus to analyze", default=0)
  parser.add_argument("--table", action="store_true", help="Print detailed cabana-like tables")
  parser.add_argument("init", type=str, help="Route or segment 1 (baseline, e.g., without Pilot Assist)")
  parser.add_argument("comp", type=str, help="Route or segment 2 (comparison, e.g., with Pilot Assist)")

  args = parser.parse_args()

  init_lr: LogIterable = LogReader(args.init)
  comp_lr: LogIterable = LogReader(args.comp)

  compare_routes(args.bus, init_lr, comp_lr, table=args.table)
