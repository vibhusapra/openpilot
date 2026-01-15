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
GREEN = '\033[92m'
YELLOW = '\033[93m'
CLEAR = '\033[0m'

def parse_route_with_time(route_str):
  """
  Parse route string with optional time range.
  Format: route_path:start-end
  Example: "724a74016b472e55/00000102--17b3f1d244/7:100-102"
  Returns: (route_path, time_range_tuple or None)
  """
  if ':' not in route_str:
    return route_str, None

  route_path, time_spec = route_str.rsplit(':', 1)

  # Parse time range
  if '-' not in time_spec:
    raise ValueError(f"Invalid time range format: {time_spec}. Expected start-end (e.g., 100-102)")

  try:
    start_str, end_str = time_spec.split('-', 1)
    start_sec = float(start_str)
    end_sec = float(end_str)

    if start_sec >= end_sec:
      raise ValueError(f"Start time ({start_sec}) must be less than end time ({end_sec})")

    return route_path, (start_sec, end_sec)
  except ValueError as e:
    raise ValueError(f"Invalid time range '{time_spec}': {e}")

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

def get_matching_signals(addr, changed_value, num_bytes, dbc):
  """
  Find DBC signals that match the changed bits.
  Returns list of (signal_name, signal_bits_mask, is_full_match) tuples.
  """
  matches = []

  if not dbc or addr not in dbc.msgs:
    return matches

  msg = dbc.msgs[addr]

  # For each signal, check if it overlaps with changed bits
  for sig_name, sig in msg.sigs.items():
    # Skip signals that cover all bytes (like ALL_BYTES)
    if sig.size >= num_bytes * 8:
      continue

    # Calculate which bits this signal occupies
    # For big-endian (Motorola) DBC signals:
    # - start_bit is the MSB position
    # - Bits are numbered: byte 0 has bits 0-7, byte 1 has bits 8-15, etc.
    # - Within each byte: bit 0 is MSB, bit 7 is LSB (MSB-first)
    # - For multi-byte signals, after reaching bit 7 of a byte, it continues at bit 0 of the next byte

    signal_bits = 0
    start_byte = sig.start_bit // 8
    start_bit_in_byte = sig.start_bit % 8

    remaining_bits = sig.size
    current_byte = start_byte
    current_bit_msb = start_bit_in_byte

    while remaining_bits > 0:
      # How many bits can we take from this byte?
      bits_available_in_byte = current_bit_msb + 1
      bits_to_take = min(remaining_bits, bits_available_in_byte)

      # Take bits from current_bit_msb down to current_bit_msb - bits_to_take + 1
      for i in range(bits_to_take):
        bit_in_byte_msb = current_bit_msb - i
        bit_in_byte_lsb = 7 - bit_in_byte_msb  # Convert to LSB-first
        message_bit_pos = current_byte * 8 + bit_in_byte_lsb
        signal_bits |= (1 << message_bit_pos)

      remaining_bits -= bits_to_take
      current_byte += 1
      current_bit_msb = 7  # Next byte starts from MSB (bit 7 MSB-first, which is bit 0 in DBC numbering of that byte)

    # Check if this signal overlaps with changed bits
    overlap = signal_bits & changed_value
    if overlap:
      # Check if ALL bits of the signal changed (perfect match)
      is_full_match = (overlap == signal_bits)
      matches.append((sig_name, signal_bits, is_full_match))

  return matches

def find_signal_for_bits(byte_idx, bit_mask, num_bytes, signal_matches):
  """
  Find if any DBC signal matches the given bit mask in the specified byte.
  bit_mask is relative to the specific byte (0-255), LSB-first within the byte.
  byte_idx is the byte index in big-endian order (0 = MSB byte).
  Returns (signal_name, is_full_match) or (None, False).
  """
  # Convert byte-relative bit mask to full message bit mask
  # In DBC: bit 0 = MSB of byte 0, bit 7 = LSB of byte 0, bit 8 = MSB of byte 1, etc.
  # Our byte_idx is big-endian (byte 0 = MSB), and bit_mask is LSB-first within byte
  # So we need to:
  # 1. Find the DBC bit positions for this byte
  # 2. Convert our LSB-first bit_mask to match those positions

  full_bit_mask = 0
  for bit_pos in range(8):
    if bit_mask & (1 << bit_pos):
      # This bit is set in the byte-level mask (LSB-first)
      # Convert to DBC bit position
      dbc_bit_in_byte = 7 - bit_pos  # DBC uses MSB-first within bytes
      dbc_bit_global = byte_idx * 8 + dbc_bit_in_byte

      # Now convert DBC bit to our message bit position
      dbc_byte_idx = dbc_bit_global // 8
      dbc_bit_in_byte_msb = dbc_bit_global % 8
      dbc_bit_in_byte_lsb = 7 - dbc_bit_in_byte_msb
      message_bit_pos = dbc_byte_idx * 8 + dbc_bit_in_byte_lsb

      full_bit_mask |= (1 << message_bit_pos)

  for sig_name, sig_bits, is_full_match in signal_matches:
    # Check if this signal exactly covers our bit mask
    if (sig_bits & full_bit_mask) == full_bit_mask:
      return sig_name, is_full_match

  return None, False

def describe_changed_bits(changed_value, num_bytes, addr=None, dbc=None):
  """
  Analyze changed bits and return human-readable description.
  Returns list of strings describing which bytes/nibbles/bits changed.
  """
  descriptions = []

  # Get signal matches if DBC is available
  signal_matches = []
  if dbc and addr:
    signal_matches = get_matching_signals(addr, changed_value, num_bytes, dbc)

  for byte_idx in range(num_bytes):
    # Extract byte (big-endian, so byte 0 is leftmost)
    byte_val = (changed_value >> (8 * (num_bytes - 1 - byte_idx))) & 0xFF

    if byte_val == 0:
      continue

    # Check if entire byte changed and if there's a matching signal
    if byte_val == 0xFF:
      sig_name, is_full = find_signal_for_bits(byte_idx, 0xFF, num_bytes, signal_matches)
      if sig_name:
        if is_full:
          descriptions.append(f"  - Byte {byte_idx} {GREEN}({sig_name}){CLEAR}")
        else:
          descriptions.append(f"  - Byte {byte_idx} {YELLOW}({sig_name} partial){CLEAR}")
      else:
        descriptions.append(f"  - Byte {byte_idx}")
      continue

    # Check nibbles
    hi_nibble = (byte_val >> 4) & 0x0F
    lo_nibble = byte_val & 0x0F

    if hi_nibble == 0x0F and lo_nibble == 0x0F:
      sig_name, is_full = find_signal_for_bits(byte_idx, 0xFF, num_bytes, signal_matches)
      if sig_name:
        if is_full:
          descriptions.append(f"  - Byte {byte_idx} {GREEN}({sig_name}){CLEAR}")
        else:
          descriptions.append(f"  - Byte {byte_idx} {YELLOW}({sig_name} partial){CLEAR}")
      else:
        descriptions.append(f"  - Byte {byte_idx}")
    elif hi_nibble == 0x0F:
      sig_name, is_full = find_signal_for_bits(byte_idx, 0xF0, num_bytes, signal_matches)
      if sig_name:
        if is_full:
          descriptions.append(f"  - Byte {byte_idx} nibble HI {GREEN}({sig_name}){CLEAR}")
        else:
          descriptions.append(f"  - Byte {byte_idx} nibble HI {YELLOW}({sig_name} partial){CLEAR}")
      else:
        descriptions.append(f"  - Byte {byte_idx} nibble HI")
    elif lo_nibble == 0x0F:
      sig_name, is_full = find_signal_for_bits(byte_idx, 0x0F, num_bytes, signal_matches)
      if sig_name:
        if is_full:
          descriptions.append(f"  - Byte {byte_idx} nibble LO {GREEN}({sig_name}){CLEAR}")
        else:
          descriptions.append(f"  - Byte {byte_idx} nibble LO {YELLOW}({sig_name} partial){CLEAR}")
      else:
        descriptions.append(f"  - Byte {byte_idx} nibble LO")
    elif hi_nibble != 0 and lo_nibble != 0:
      # Both nibbles have some bits, list individually
      # Use LSB-first numbering: bit 0 = LSB (rightmost), bit 7 = MSB (leftmost)
      for bit_pos in range(8):
        if byte_val & (1 << bit_pos):
          bit_mask = 1 << bit_pos
          sig_name, is_full = find_signal_for_bits(byte_idx, bit_mask, num_bytes, signal_matches)
          if sig_name:
            if is_full:
              descriptions.append(f"  - Byte {byte_idx} bit {bit_pos} {GREEN}({sig_name}){CLEAR}")
            else:
              descriptions.append(f"  - Byte {byte_idx} bit {bit_pos} {YELLOW}({sig_name} partial){CLEAR}")
          else:
            descriptions.append(f"  - Byte {byte_idx} bit {bit_pos}")
    else:
      # Only one nibble has bits, check if it's all bits or specific ones
      if hi_nibble != 0:
        if hi_nibble == 0x0F:
          sig_name, is_full = find_signal_for_bits(byte_idx, 0xF0, num_bytes, signal_matches)
          if sig_name:
            if is_full:
              descriptions.append(f"  - Byte {byte_idx} nibble HI {GREEN}({sig_name}){CLEAR}")
            else:
              descriptions.append(f"  - Byte {byte_idx} nibble HI {YELLOW}({sig_name} partial){CLEAR}")
          else:
            descriptions.append(f"  - Byte {byte_idx} nibble HI")
        else:
          # HI nibble: bits 4-7 (LSB-first numbering)
          for bit_pos in range(4, 8):
            if byte_val & (1 << bit_pos):
              bit_mask = 1 << bit_pos
              sig_name, is_full = find_signal_for_bits(byte_idx, bit_mask, num_bytes, signal_matches)
              if sig_name:
                if is_full:
                  descriptions.append(f"  - Byte {byte_idx} bit {bit_pos} {GREEN}({sig_name}){CLEAR}")
                else:
                  descriptions.append(f"  - Byte {byte_idx} bit {bit_pos} {YELLOW}({sig_name} partial){CLEAR}")
              else:
                descriptions.append(f"  - Byte {byte_idx} bit {bit_pos}")
      else:  # lo_nibble != 0
        if lo_nibble == 0x0F:
          sig_name, is_full = find_signal_for_bits(byte_idx, 0x0F, num_bytes, signal_matches)
          if sig_name:
            if is_full:
              descriptions.append(f"  - Byte {byte_idx} nibble LO {GREEN}({sig_name}){CLEAR}")
            else:
              descriptions.append(f"  - Byte {byte_idx} nibble LO {YELLOW}({sig_name} partial){CLEAR}")
          else:
            descriptions.append(f"  - Byte {byte_idx} nibble LO")
        else:
          # LO nibble: bits 0-3 (LSB-first numbering)
          for bit_pos in range(0, 4):
            if byte_val & (1 << bit_pos):
              bit_mask = 1 << bit_pos
              sig_name, is_full = find_signal_for_bits(byte_idx, bit_mask, num_bytes, signal_matches)
              if sig_name:
                if is_full:
                  descriptions.append(f"  - Byte {byte_idx} bit {bit_pos} {GREEN}({sig_name}){CLEAR}")
                else:
                  descriptions.append(f"  - Byte {byte_idx} bit {bit_pos} {YELLOW}({sig_name} partial){CLEAR}")
              else:
                descriptions.append(f"  - Byte {byte_idx} bit {bit_pos}")

  return descriptions

def collect_stable_bits(msgs, bus, time_range=None):
  """
  Collect bits that remain stable (don't change) throughout the message set.
  Returns: dict mapping address -> stable bit mask, stable values, data, and frequencies

  Args:
    msgs: Message iterator
    bus: CAN bus number
    time_range: Optional tuple (start_sec, end_sec) to filter messages by logMonoTime
  """
  dat = defaultdict(lambda: None)
  low_to_high = defaultdict(int)  # Bits ever seen as 1
  high_to_low = defaultdict(int)  # Bits ever seen as 0
  msg_count = defaultdict(int)  # Count messages per address
  first_timestamp = None
  last_timestamp = None
  route_start_time = None

  for x in msgs:
    if x.which() != 'can':
      continue

    # Track the absolute start time of the route
    if route_start_time is None:
      route_start_time = x.logMonoTime

    # Filter by time range if specified
    if time_range:
      elapsed_sec = (x.logMonoTime - route_start_time) / 1e9
      if elapsed_sec < time_range[0] or elapsed_sec > time_range[1]:
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


def compare_routes(bus, init_msgs, comp_msgs, table=False, dbc_name=None, init_time_range=None, comp_time_range=None):
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

  time_info_1 = f" (time range: {init_time_range[0]}-{init_time_range[1]}s)" if init_time_range else ""
  print(f"Analyzing route 1 (baseline){time_info_1}...")
  stable_bits_1, stable_values_1, dat_1, freq_1 = collect_stable_bits(init_msgs, bus, init_time_range)

  time_info_2 = f" (time range: {comp_time_range[0]}-{comp_time_range[1]}s)" if comp_time_range else ""
  print(f"Analyzing route 2 (comparison){time_info_2}...")
  stable_bits_2, stable_values_2, dat_2, freq_2 = collect_stable_bits(comp_msgs, bus, comp_time_range)

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
      descriptions = describe_changed_bits(changed_bits, num_bytes, addr, dbc)
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
      descriptions = describe_changed_bits(new_stable_with_value, num_bytes, addr, dbc)
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

  Time range filtering:
    You can specify time ranges (in seconds) by appending :start-end to the route path:
    python can_print_changes_2.py route1:10-20 route2:100-102 --dbc volvo_cma
    This analyzes only seconds 10-20 from route1 and seconds 100-102 from route2.
  """
  parser = argparse.ArgumentParser(description=desc,
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("--bus", type=int, help="CAN bus to analyze", default=0)
  parser.add_argument("--table", action="store_true", help="Print detailed cabana-like tables")
  parser.add_argument("--dbc", type=str, help="DBC file name (e.g., 'volvo_cma') to show message names", default=None)
  parser.add_argument("init", type=str, help="Route or segment 1 (baseline, e.g., without Pilot Assist). Optional time range: route:start-end")
  parser.add_argument("comp", type=str, help="Route or segment 2 (comparison, e.g., with Pilot Assist). Optional time range: route:start-end")

  args = parser.parse_args()

  # Parse routes and time ranges
  init_path, init_time_range = parse_route_with_time(args.init)
  comp_path, comp_time_range = parse_route_with_time(args.comp)

  init_lr: LogIterable = LogReader(init_path)
  comp_lr: LogIterable = LogReader(comp_path)

  compare_routes(args.bus, init_lr, comp_lr, table=args.table, dbc_name=args.dbc,
                 init_time_range=init_time_range, comp_time_range=comp_time_range)
