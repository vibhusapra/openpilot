#!/usr/bin/env python3
"""
This script compares two routes and prints bits that are stable within each route but different between routes.
This filters out counters and checksums that change within a route, revealing only meaningful signal changes.

Difference between "NEW" and "CHANGED":

  CHANGED

  Bits that are stable in BOTH routes, but have different values

  Example:
  - Route 1 (no Pilot Assist): A bit is consistently 0 throughout the entire route
  - Route 2 (with Pilot Assist): The same bit is consistently 1 throughout the entire route
  - Result: CHANGED - this is likely a state flag (Pilot Assist on/off)

  NEW

  Bits that only became stable in route 2

  Example:
  - Route 1: A bit doesn't exist or is fluctuating/unstable
  - Route 2: The bit appears and remains stable at a specific value
  - Result: NEW - this could be a new message that only appears when Pilot Assist is active

  Real-World Example

  Imagine you're comparing:
  - Route 1: Driving normally (no Pilot Assist)
  - Route 2: Driving with Pilot Assist engaged

  CHANGED:
  PILOT_ASSIST_STATUS
  0x123 (291 ) CHANGED: 01000000000000
               Route 1: b'00000000000000'  ← bit was 0
               Route 2: b'01000000000000'  ← bit is now 1
  This is a flag that toggled when you engaged Pilot Assist.

  NEW:
  PILOT_ASSIST_CONTROL
  0x456 (1110) NEW:     ff000000000000
               Route 2: b'ff000000000000'
  This message (or these specific bits) only appeared and became stable in route 2. Maybe this message only
  exists when Pilot Assist is running, sending control commands.

  Does this make sense? The key insight is that both filter out counters/checksums (which fluctuate within
  each route), but they catch different types of signals.
  """
import argparse
import binascii
from collections import defaultdict

from opendbc.can.parser import DBC
from openpilot.selfdrive.debug.can_table import can_table
from openpilot.tools.lib.logreader import LogIterable, LogReader

RED = '\033[91m'
CLEAR = '\033[0m'

def get_msg_name(addr, dbc=None):
  """Get message name from DBC if available."""
  if dbc and addr in dbc.msgs:
    return dbc.msgs[addr].name
  return None

def format_frequency(addr, freq_1, freq_2):
  """Format frequency string, showing both if different."""
  f1 = freq_1.get(addr, 0)
  f2 = freq_2.get(addr, 0)

  if f1 == 0 and f2 == 0:
    return ""
  elif f1 == 0:
    return f"({f2} Hz)"
  elif f2 == 0:
    return f"({f1} Hz)"
  elif f1 == f2:
    return f"({f1} Hz)"
  else:
    return f"({f1}/{f2} Hz)"

def format_addr(addr):
  """Format address consistently."""
  return f"{hex(addr).ljust(6)}({str(addr).ljust(4)})"

def describe_changed_bits(changed_value, num_bytes):
  """
  Analyze changed bits and return human-readable description.
  Returns list of strings describing which bytes/nibbles/bits changed.
  """
  descriptions = []

  for byte_idx in range(num_bytes):
    # Extract byte (big-endian, so byte 0 is leftmost)
    byte_val = (changed_value >> (8 * (num_bytes - 1 - byte_idx))) & 0xFF

    if byte_val == 0:
      continue

    # Check if entire byte changed
    if byte_val == 0xFF:
      descriptions.append(f"  - Byte {byte_idx}")
      continue

    # Check nibbles
    hi_nibble = (byte_val >> 4) & 0x0F
    lo_nibble = byte_val & 0x0F

    if hi_nibble == 0x0F and lo_nibble == 0x0F:
      descriptions.append(f"  - Byte {byte_idx}")
    elif hi_nibble == 0x0F:
      descriptions.append(f"  - Byte {byte_idx} nibble HI")
    elif lo_nibble == 0x0F:
      descriptions.append(f"  - Byte {byte_idx} nibble LO")
    elif hi_nibble != 0 and lo_nibble != 0:
      # Both nibbles have some bits, list individually
      for bit_idx in range(8):
        if byte_val & (1 << (7 - bit_idx)):
          descriptions.append(f"  - Byte {byte_idx} bit {bit_idx}")
    else:
      # Only one nibble has bits, check if it's all bits or specific ones
      if hi_nibble != 0:
        if hi_nibble == 0x0F:
          descriptions.append(f"  - Byte {byte_idx} nibble HI")
        else:
          for bit_idx in range(4, 8):
            if byte_val & (1 << (7 - bit_idx)):
              descriptions.append(f"  - Byte {byte_idx} bit {bit_idx}")
      else:  # lo_nibble != 0
        if lo_nibble == 0x0F:
          descriptions.append(f"  - Byte {byte_idx} nibble LO")
        else:
          for bit_idx in range(4, 8):
            if byte_val & (1 << (7 - bit_idx)):
              descriptions.append(f"  - Byte {byte_idx} bit {bit_idx}")

  return descriptions

def collect_stable_bits(msgs, bus):
  """
  Collect bits that remain stable (don't change) throughout the message set.
  Returns: dict mapping address -> stable bit mask, stable values, data, and frequencies
  """
  dat = defaultdict(lambda: None)
  low_to_high = defaultdict(int)  # Bits ever seen as 1
  high_to_low = defaultdict(int)  # Bits ever seen as 0
  msg_count = defaultdict(int)  # Count messages per address
  first_timestamp = None
  last_timestamp = None

  for x in msgs:
    if x.which() != 'can':
      continue

    # Track timestamps for frequency calculation
    if first_timestamp is None:
      first_timestamp = x.logMonoTime
    last_timestamp = x.logMonoTime

    for y in x.can:
      if y.src == bus:
        if dat[y.address] is None:
          dat[y.address] = y.dat

        msg_count[y.address] += 1
        i = int.from_bytes(y.dat, byteorder='big')
        low_to_high[y.address] |= i      # Accumulate bits seen as 1
        high_to_low[y.address] |= ~i     # Accumulate bits seen as 0

  # Calculate frequencies in Hz
  frequencies = {}
  if first_timestamp and last_timestamp:
    duration_sec = (last_timestamp - first_timestamp) / 1e9  # Convert nanoseconds to seconds
    if duration_sec > 0:
      for addr in msg_count.keys():
        frequencies[addr] = round(msg_count[addr] / duration_sec)

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

  return stable_bits, stable_values, dat, frequencies


def compare_routes(bus, init_msgs, comp_msgs, table=False, dbc_name=None):
  """
  Compare two routes by finding bits that are stable within each route
  but different between routes.
  """
  # Load DBC if specified
  dbc = None
  if dbc_name:
    try:
      dbc = DBC(dbc_name)
      print(f"Loaded DBC: {dbc_name}")
    except Exception as e:
      print(f"Warning: Could not load DBC '{dbc_name}': {e}")
      print("Continuing without DBC...\n")

  print("Analyzing route 1 (baseline)...")
  stable_bits_1, stable_values_1, dat_1, freq_1 = collect_stable_bits(init_msgs, bus)

  print("Analyzing route 2 (comparison)...")
  stable_bits_2, stable_values_2, dat_2, freq_2 = collect_stable_bits(comp_msgs, bus)

  print("\n" + "="*80)
  print("STABLE BIT CHANGES (filters out counters and checksums)")
  print("="*80 + "\n")

  # Find all addresses present in either route
  all_addrs = set(stable_bits_1.keys()) | set(stable_bits_2.keys())

  # Sort by frequency (highest to lowest), using max of both routes
  def get_max_freq(addr):
    f1 = freq_1.get(addr, 0)
    f2 = freq_2.get(addr, 0)
    return max(f1, f2)

  sorted_addrs = sorted(all_addrs, key=get_max_freq, reverse=True)

  tables = ""
  changes_found = False

  for addr in sorted_addrs:
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

      # Print message name and frequency on separate line if available
      msg_name = get_msg_name(addr, dbc)
      freq_str = format_frequency(addr, freq_1, freq_2)
      if msg_name or freq_str:
        print(f"{msg_name} {freq_str}" if msg_name else freq_str)

      header = format_addr(addr)
      print(f"{header} CHANGED: {byts}")

      value_1_bytes = value_1.to_bytes(num_bytes, byteorder='big')
      value_2_bytes = value_2.to_bytes(num_bytes, byteorder='big')
      print(f"{'':>{len(header)}} Route 1: {binascii.hexlify(value_1_bytes)}")
      print(f"{'':>{len(header)}} Route 2: {binascii.hexlify(value_2_bytes)}")

      # Print human-readable breakdown
      descriptions = describe_changed_bits(changed_bits, num_bytes)
      if descriptions:
        print(f"{'':>{len(header)}} Changed:")
        for desc in descriptions:
          print(f"{'':>{len(header)}}{desc}")

      print()

      tables += f"{header} CHANGED\n"
      tables += can_table(b) + "\n\n"

    # Print new stable bits
    if new_stable_with_value != 0:
      b = new_stable_with_value.to_bytes(num_bytes, byteorder='big')
      byts = ''.join([(c if c == '0' else f'{RED}{c}{CLEAR}') for c in str(binascii.hexlify(b))[2:-1]])

      # Print message name and frequency on separate line if available
      msg_name = get_msg_name(addr, dbc)
      freq_str = format_frequency(addr, freq_1, freq_2)
      if msg_name or freq_str:
        print(f"{msg_name} {freq_str}" if msg_name else freq_str)

      header = format_addr(addr)
      print(f"{header} NEW:     {byts}")

      value_2_bytes = value_2.to_bytes(num_bytes, byteorder='big')
      print(f"{'':>{len(header)}} Route 2: {binascii.hexlify(value_2_bytes)}")

      # Print human-readable breakdown
      descriptions = describe_changed_bits(new_stable_with_value, num_bytes)
      if descriptions:
        print(f"{'':>{len(header)}} New:")
        for desc in descriptions:
          print(f"{'':>{len(header)}}{desc}")

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
  parser.add_argument("--dbc", type=str, help="DBC file name (e.g., 'volvo_cma') to show message names", default=None)
  parser.add_argument("init", type=str, help="Route or segment 1 (baseline, e.g., without Pilot Assist)")
  parser.add_argument("comp", type=str, help="Route or segment 2 (comparison, e.g., with Pilot Assist)")

  args = parser.parse_args()

  init_lr: LogIterable = LogReader(args.init)
  comp_lr: LogIterable = LogReader(args.comp)

  compare_routes(args.bus, init_lr, comp_lr, table=args.table, dbc_name=args.dbc)
