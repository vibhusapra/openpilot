#!/usr/bin/env python3
"""
CAN Message Analysis Script for LCA Reverse Engineering

Extracts vehicle signals and raw CAN message bytes from openpilot route logs
to facilitate reverse engineering of the Volvo CMA LCA (Lane Centering Assist) protocol.

Uses the current Volvo CarState class to reprocess CAN messages with the latest DBC,
ensuring signals are always parsed with up-to-date definitions even for old routes.

Usage:
  python analyze_lca_messages.py <route_id> [--output output.csv]

Example:
  python analyze_lca_messages.py 724a74016b472e55/00000045--5b72c9b386
"""

import argparse
import pandas as pd
from collections import defaultdict
import re

# Add openpilot to path
import sys
from pathlib import Path
OPENPILOT_PATH = Path(__file__).parent.parent
sys.path.insert(0, str(OPENPILOT_PATH))

from tools.lib.logreader import LogReader
from opendbc.car.volvo.carstate import CarState
from opendbc.car.volvo.values import CAR
from opendbc.car import structs, Bus


# Messages to capture raw bytes from (message_id, bus, name)
RAW_MESSAGES = [
  (0x58, 0, "LCA"),
  (0x57, 0, "LCA_3"),
  (0x69, 0, "LCA_2"),  # Contains PILOT_ASSIST_ENGAGED
  (0x16, 2, "PSCM"),
  (0x17, 2, "PSCM_RELATED"),
  (0x15, 2, "DRIVER_INPUT"),
]


def bytes_to_hex(data: bytes) -> str:
  """Convert bytes to hex string (e.g., b'\\x01\\x02' -> '0102')."""
  return data.hex().upper()


def parse_pilot_assist_engaged(lca_2_bytes: bytes) -> bool:
  """
  Parse PILOT_ASSIST_ENGAGED from LCA_2 message.
  LCA_2 message ID 0x69, PILOT_ASSIST_ENGAGED is bit 12 (LSB 0).
  """
  if len(lca_2_bytes) < 2:
    return False

  # Bit 12 is in byte 1, bit 4 (bytes are little-endian, bit 12 = byte_offset 1, bit 4)
  # bit_position = 12 -> byte_index = 12 // 8 = 1, bit_in_byte = 12 % 8 = 4
  byte_index = 12 // 8
  bit_in_byte = 12 % 8

  if len(lca_2_bytes) <= byte_index:
    return False

  return bool((lca_2_bytes[byte_index] >> bit_in_byte) & 1)


def parse_route(route_id: str, break_on_engaged: bool = False) -> pd.DataFrame:
  """
  Parse route logs and extract all relevant signals and raw message bytes.

  Uses the current Volvo CarState class to reprocess CAN with latest DBC.

  Args:
    route_id: Route identifier (e.g., '724a74016b472e55/00000045--5b72c9b386')

  Returns:
    DataFrame with timestamp-aligned signals and raw message bytes
  """
  print(f"Loading route: {route_id}")

  # Load logs
  lr = LogReader(route_id)

  # Initialize Volvo CarState with current DBC
  CP = structs.CarParams()
  CP.carFingerprint = CAR.VOLVO_XC40_RECHARGE
  car_state = CarState(CP)
  car_state.frame = 0  # Initialize frame counter (normally done by card.py)
  can_parsers = CarState.get_can_parsers(CP)

  # Storage for parsed data
  data_rows = []

  # Store latest raw message for each target message (much faster than timestamp lookup)
  latest_raw_messages = {}  # {msg_name: bytes}

  # Batch CAN messages per bus to feed to parsers
  can_batch = {Bus.main: [], Bus.pt: [], Bus.party: []}  # {bus: [(addr, data, src), ...]}
  last_can_timestamp = 0

  # Track previous values for derivative calculations
  prev_steering_angle = None
  prev_timestamp = None

  # Also track lateral accel from sensorEvents
  latest_lateral_accel = 0.0

  print("Processing logs...")
  msg_count = 0

  for msg in lr:
    msg_count += 1
    if msg_count % 10000 == 0:
      print(f"  Processed {msg_count} messages...")

    # Process CAN messages - both capture raw bytes AND batch for parsers
    if msg.which() == 'can':
      last_can_timestamp = msg.logMonoTime

      for can_msg in msg.can:
        bus = can_msg.src  # This is an integer (0, 1, 2)
        address = can_msg.address
        data = can_msg.dat

        # Capture raw bytes for target messages - just keep latest
        for msg_id, msg_bus, msg_name in RAW_MESSAGES:
          if address == msg_id and bus == msg_bus:
            latest_raw_messages[msg_name] = data

        # Batch CAN messages per bus to update parsers later
        bus_enum = [Bus.main, Bus.pt, Bus.party][bus] if bus < 3 else None
        if bus_enum and bus_enum in can_parsers:
          can_batch[bus_enum].append((address, data, bus))

    # Get lateral acceleration from sensorEvents (IMU data)
    elif msg.which() == 'sensorEvents':
      for evt in msg.sensorEvents:
        if evt.which() == 'acceleration':
          # acceleration.v[1] is lateral (y-axis)
          latest_lateral_accel = evt.acceleration.v[1]

    # Process carState messages to get timestamp sync
    # But we'll use our reprocessed CAN data, not the logged carState
    elif msg.which() == 'carState':
      timestamp = msg.logMonoTime / 1e9

      # Update all parsers with batched CAN messages
      for bus_enum in [Bus.main, Bus.pt, Bus.party]:
        if can_batch[bus_enum]:
          can_parsers[bus_enum].update([[last_can_timestamp, can_batch[bus_enum]]])
          can_batch[bus_enum] = []  # Clear batch

      # Reprocess CAN with current DBC using CarState.update()
      car_state.frame += 1  # Increment frame counter
      cs_reparsed = car_state.update(can_parsers)

      # Extract signals from reprocessed carState
      v_ego_raw = cs_reparsed.vEgoRaw  # m/s
      gas_pressed = cs_reparsed.gasPressed
      brake_pressed = cs_reparsed.brakePressed
      steering_angle_deg = cs_reparsed.steeringAngleDeg
      steering_torque = cs_reparsed.steeringTorque
      steering_pressed = cs_reparsed.steeringPressed
      gear_position = cs_reparsed.gearShifter
      cruise_enabled = cs_reparsed.cruiseState.enabled if hasattr(cs_reparsed, 'cruiseState') else False

      # Wheel speeds from reprocessed data - COMMENTED OUT (signals may not be correct)
      # if hasattr(cs_reparsed, 'wheelSpeeds'):
      #   wheel_speed_fl = cs_reparsed.wheelSpeeds.fl
      #   wheel_speed_fr = cs_reparsed.wheelSpeeds.fr
      #   wheel_speed_rl = cs_reparsed.wheelSpeeds.rl
      #   wheel_speed_rr = cs_reparsed.wheelSpeeds.rr
      # else:
      #   wheel_speed_fl = wheel_speed_fr = wheel_speed_rl = wheel_speed_rr = 0.0

      # Calculate steering angle rate
      steering_angle_rate = 0.0
      if prev_steering_angle is not None and prev_timestamp is not None:
        dt = timestamp - prev_timestamp
        if dt > 0:
          steering_angle_rate = (steering_angle_deg - prev_steering_angle) / dt

      # Get raw message bytes (latest values - no expensive timestamp search!)
      lca_bytes = bytes_to_hex(latest_raw_messages.get("LCA", b''))
      lca_3_bytes = bytes_to_hex(latest_raw_messages.get("LCA_3", b''))
      lca_2_bytes_hex = bytes_to_hex(latest_raw_messages.get("LCA_2", b''))
      pscm_bytes = bytes_to_hex(latest_raw_messages.get("PSCM", b''))
      pscm_related_bytes = bytes_to_hex(latest_raw_messages.get("PSCM_RELATED", b''))
      driver_input_bytes = bytes_to_hex(latest_raw_messages.get("DRIVER_INPUT", b''))

      # Parse PILOT_ASSIST_ENGAGED from LCA_2 raw bytes
      lca_2_raw = latest_raw_messages.get("LCA_2", b'')
      pilot_assist_engaged = parse_pilot_assist_engaged(lca_2_raw)

      # Build row
      row = {
        'timestamp': timestamp,
        'vEgoRaw_ms': v_ego_raw,
        'gasPressed': gas_pressed,
        'brakePressed': brake_pressed,
        'steeringAngleDeg': steering_angle_deg,
        'steeringTorque': steering_torque,
        'steeringPressed': steering_pressed,
        'steeringAngleRate_deg_s': steering_angle_rate,
        'lateral_accel': latest_lateral_accel,  # From IMU/sensor fusion
        'pilot_assist_engaged': pilot_assist_engaged,
        'cruise_enabled': cruise_enabled,
        'gear_position': gear_position,
        # 'wheel_speed_fl': wheel_speed_fl,
        # 'wheel_speed_fr': wheel_speed_fr,
        # 'wheel_speed_rl': wheel_speed_rl,
        # 'wheel_speed_rr': wheel_speed_rr,
        'LCA_raw_hex': lca_bytes,
        'LCA_3_raw_hex': lca_3_bytes,
        'LCA_2_raw_hex': lca_2_bytes_hex,
        'PSCM_raw_hex': pscm_bytes,
        'PSCM_RELATED_raw_hex': pscm_related_bytes,
        'DRIVER_INPUT_raw_hex': driver_input_bytes,
      }

      data_rows.append(row)

      # Break early if requested and pilot assist is engaged
      if break_on_engaged and pilot_assist_engaged:
        print(f"\n✓ Found first pilot_assist_engaged sample at {timestamp:.1f}s")
        print(f"  vEgoRaw={v_ego_raw:.2f} m/s, steeringAngleDeg={steering_angle_deg:.2f}, lateral_accel={latest_lateral_accel:.3f}")
        break

      # Update previous values
      prev_steering_angle = steering_angle_deg
      prev_timestamp = timestamp

  print(f"Total messages processed: {msg_count}")
  print(f"CarState samples extracted: {len(data_rows)}")

  # Convert to DataFrame
  df = pd.DataFrame(data_rows)

  # Sort by timestamp
  df = df.sort_values('timestamp').reset_index(drop=True)

  return df

_illegal = re.compile(r'[<>:"/\\|?*\x00-\x1F]')

def clean(name: str) -> str:
    name = _illegal.sub("", name).strip()
    return name[:150] or "file"


def main():
  parser = argparse.ArgumentParser(
    description='Extract CAN signals and raw message bytes for LCA reverse engineering',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__
  )
  parser.add_argument('route_id', help='Route ID (e.g., 724a74016b472e55/00000045--5b72c9b386)')
  parser.add_argument('--output', '-o', default='lca_analysis.csv', help='Output CSV filename')
  parser.add_argument('--break-engaged', action='store_true', help='Stop after first pilot_assist_engaged sample (for testing)')

  args = parser.parse_args()

  if args.output == 'lca_analysis.csv':
    #args.output = f'lca_analysis_{args.route_id}.csv'
    args.output = f'{clean(args.route_id)}.csv'

  # Parse route
  df = parse_route(args.route_id, break_on_engaged=args.break_engaged)

  # Export to CSV
  print(f"\nExporting to {args.output}...")
  df.to_csv(args.output, index=False)

  # Print summary statistics
  print("\n=== Summary Statistics ===")
  print(f"Total samples: {len(df)}")
  print(f"Duration: {df['timestamp'].max() - df['timestamp'].min():.1f} seconds")
  print(f"Pilot Assist engaged samples: {df['pilot_assist_engaged'].sum()} ({100*df['pilot_assist_engaged'].mean():.1f}%)")
  print(f"Speed range: {df['vEgoRaw_ms'].min():.1f} - {df['vEgoRaw_ms'].max():.1f} m/s")
  print(f"Steering angle range: {df['steeringAngleDeg'].min():.1f} - {df['steeringAngleDeg'].max():.1f} deg")
  print(f"\nFirst few rows:")
  print(df.head())
  print(f"\nCSV saved to: {args.output}")


if __name__ == '__main__':
  main()
