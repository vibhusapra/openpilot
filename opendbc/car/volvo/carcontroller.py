from opendbc.can.packer import CANPacker
from opendbc.car import Bus
from opendbc.car.interfaces import CarControllerBase
from opendbc.car.volvo.helpers import LCA3CounterSync
from opendbc.car.volvo.live_testing import LiveTestingManager
from opendbc.car.volvo.volvocan import create_lca_message, create_pscm_message, create_lca_3_message, create_lca_2_message, create_lca_4_message, create_lca_5_message, create_lca_6_message, create_lca_7_message, create_speed_message, create_speed_2_message, create_speed_3_message, create_0x1a_message, create_gear_position_message, create_egsm_message, create_pscm_related_message
from opendbc.car.volvo.values import CarControllerParams


class CarController(CarControllerBase):
  def __init__(self, dbc_names, CP):
    super().__init__(dbc_names, CP)
    self.packer = CANPacker(dbc_names[Bus.party])
    self.apply_angle_last = 0.0  # Track last applied steering angle

    # Live testing configuration manager
    self.live_testing = LiveTestingManager()
    self.liveTestingConfig = None  # Updated at 50 Hz, persists between frames

    self.gear_acc = 60
    self.lca_4_acc = 0  # Bresenham accumulator for 29 Hz

    # Counter management for LCA_2
    self.lca_2_counter_1 = None  # Will grab initial value from CarState
    self.lca_2_counter_2 = None

    # Counter management for PSCM_RELATED
    self.pscm_related_counter = None  # Will grab initial value from CarState

    # Counter management for LCA_3 (pattern-based)
    self.lca_3_counter_sync = LCA3CounterSync()

    # Counter management for LCA_5 (formerly SPEED_1)
    self.lca_5_counter = None  # Will grab initial value from CarState

    self.last_lat_active = False  # Track state

    self.lca_7_acc = 0  # Bresenham accumulator for 29 Hz
    self.lca_7_last_steer = 0 # used to calculate change in steer from last update

  def update(self, CC, CS, now_nanos):
    CS.CC_frame = self.frame
    can_sends = []
    actuators = CC.actuators

    # Detect disengagement
    if not CC.latActive and self.last_lat_active:
      #self.lca_commands.reset()  # Clear state ← IMPORTANT!
      pass

    # Explicit if condition for lat_active override (uses config loaded at 50 Hz)
    if self.liveTestingConfig:
      lat_active = self.liveTestingConfig.get('lat_active', CC.latActive)
    else:
      lat_active = CC.latActive
    #lat_active = True

    # lateral control - angle-based steering
    # NOTE: LCA message is sent every frame (even when inactive) to replace stock LCA
    # Stock LCA is permanently blocked by panda safety, so we must always send
    if self.frame % CarControllerParams.STEER_STEP == 0:  # 100 Hz
      # Get desired steering angle from controlsd (LatControlAngle)
      apply_angle = actuators.steeringAngleDeg  # degrees

      if not CC.latActive:
        apply_angle = CS.out.steeringAngleDeg  # Use current angle when inactive

      # Override apply_angle from live testing config if provided
      if self.liveTestingConfig:
        override_apply_angle = self.liveTestingConfig.get('apply_angle')
        if override_apply_angle is not None:
          apply_angle = override_apply_angle

      # No rate limiting initially - apply desired angle directly
      # TODO: Add rate limiting after basic functionality is confirmed

      # LCA - 0x58 - 100 Hz (angle-based)
      lca_overrides = self.liveTestingConfig.get('lca') if self.liveTestingConfig else None
      can_sends.append(create_lca_message(self.packer, lat_active, apply_angle, CS.msg_lca, lca_overrides))
      self.apply_angle_last = apply_angle

      # Check if PA hands-on-wheel spoof toggle is enabled (bit 7 of alternativeExperience)
      spoof_pa_hands_enabled = bool(self.CP.alternativeExperience & 128)
      spoof_pa_hands = CS.pilot_assist_engaged and spoof_pa_hands_enabled
      # PSCM (bus 2 -> 0) - 0x16 - 100 Hz
      can_sends.append(create_pscm_message(self.packer, lat_active, CS.msg_pscm, self.frame, spoof_pa_hands))
      # EGSM - 0x45 - 100 Hz
      #can_sends.append(create_egsm_message(self.packer, CS.msg_egsm))

      # PSCM_RELATED (bus 2 -> 0) - 0x17 - 100 Hz
      # Initialize counter from CarState on first run
      if self.pscm_related_counter is None:
        self.pscm_related_counter = CS.msg_pscm_related['SIG1_BYTE_1_HI_NIBBLE']

      # Increment counter by +1, wrap from 14 → 0 (modulo 15)
      self.pscm_related_counter = (self.pscm_related_counter + 1) % 15

      can_sends.append(create_pscm_related_message(self.packer, lat_active, CS.pilot_assist_engaged,
                                                     CS.msg_pscm_related, self.pscm_related_counter))

    # LCA_3 - 0x57 - avg 66.66 Hz
    #if (self.frame * 67) % 100 < 67: # if (self.frame % 3) < 2:
    # 0x57 at ~66.67 Hz: send on 2 out of every 3 frames
    # Pattern: send on frame % 3 == 0 or 2, skip when frame % 3 == 1
    if self.frame % 3 != 1:  # → 2/3 * 100 Hz = 66.67 Hz
      # Update counter with observed value, get counter to send
      counter, is_synced = self.lca_3_counter_sync.update(CS.msg_lca_3['COUNTER_1'])
      can_sends.append(create_lca_3_message(self.packer, lat_active, apply_angle, CS.msg_lca_3, counter))
      #can_sends.append(create_0x1a_message(self.packer, CS.msg_0x1a))
      pass

    # SPEED messages - 0x60, 0x68 - 50 Hz
    if self.frame % 2 == 0: # 50 Hz
      #can_sends.append(create_speed_message(self.packer, CS.msg_speed))
      #can_sends.append(create_speed_2_message(self.packer, CS.msg_speed_2))
      pass

    # LCA_2 - 0x69 - 50 Hz
    # Spoof PILOT_ASSIST_ENGAGED to keep PSCM accepting LCA commands
    if self.frame % 2 == 0: # 50 Hz
      # Initialize counters from CarState on first run
      if self.lca_2_counter_1 is None:
        self.lca_2_counter_1 = CS.msg_lca_2['COUNTER_1']
        self.lca_2_counter_2 = CS.msg_lca_2['COUNTER_2']

      # Increment counters (COUNTER_1 by +2, COUNTER_2 by +4, both modulo 16)
      self.lca_2_counter_1 = (self.lca_2_counter_1 + 2) % 16
      self.lca_2_counter_2 = (self.lca_2_counter_2 + 4) % 16

      can_sends.append(create_lca_2_message(self.packer, lat_active, CS.msg_lca_2,
                                            self.lca_2_counter_1, self.lca_2_counter_2))
      pass

    # LCA_5 (formerly SPEED_1) - 0x67 - 50 Hz
    # Contains wheel speeds + LCA signals (LCA_TURN_BITS, LCA_5_STEER)
    if self.frame % 2 == 0: # 50 Hz
      # === LIVE TESTING: Load config at 50 Hz (matches LCA_5 message frequency) ===
      # Returns None if file doesn't exist or lat_active=False in file
      # Config persists between frames for use by other messages
      self.liveTestingConfig = self.live_testing.load_config()

      # Initialize counter from CarState on first run
      if self.lca_5_counter is None:
        self.lca_5_counter = CS.msg_lca_5['COUNTER']

      # Increment counter by +4, wrap at 15 (0xF never used)
      self.lca_5_counter = (self.lca_5_counter + 4) % 15

      lca_5_overrides = self.liveTestingConfig.get('lca_5') if self.liveTestingConfig else None
      can_sends.append(create_lca_5_message(self.packer, lat_active, apply_angle,
                                            CS.msg_lca_5, self.lca_5_counter, lca_5_overrides))

    # LCA_4 - 0x90 - 29 Hz
    # Spoof LCA_ENABLE bits to maintain PA ON state when openpilot is active
    # Using Bresenham-style accumulator for precise 29 Hz
    self.lca_4_acc += 29
    if self.lca_4_acc >= 100:
      self.lca_4_acc -= 100
      lca_4_overrides = self.liveTestingConfig.get('lca_4') if self.liveTestingConfig else None
      can_sends.append(create_lca_4_message(self.packer, lat_active, CS.msg_lca_4, apply_angle, lca_4_overrides))

    # LCA_6 - 0X97 - 25 Hz
    if self.frame % 4 == 0: # 25 Hz
        lca_6_overrides = self.liveTestingConfig.get('lca_6') if self.liveTestingConfig else None
        can_sends.append(create_lca_6_message(self.packer, lat_active, CS.msg_lca_6, apply_angle, lca_6_overrides))

    # LCA_7 - 0x92 - 29 Hz
    # Using Bresenham-style accumulator for precise 29 Hz
    self.lca_7_acc += 29
    if self.lca_7_acc >= 100:
      self.lca_7_acc -= 100
      delta_steer = apply_angle - self.lca_7_last_steer
      lca_7_overrides = self.liveTestingConfig.get('lca_7') if self.liveTestingConfig else None
      can_sends.append(create_lca_7_message(self.packer, lat_active, CS.msg_lca_7, apply_angle, delta_steer, lca_7_overrides))
      self.lca_7_last_steer = apply_angle

    # GEAR_POSITION - 0x80 - 40 Hz
    #self.gear_acc += 40 # Bresenham-style approach
    #if self.gear_acc >= 100:
    #    self.gear_acc -= 100
    if self.frame % 5 == 0 or self.frame % 5 == 2:  # 2/5 * 100 Hz = 40 Hz # openpilot forward delay causes DTC in EGSM, but fixes DTC in PSCM
      #can_sends.append(create_gear_position_message(self.packer, CS.msg_gear_position))
      pass

    new_actuators = actuators.as_builder()
    new_actuators.steeringAngleDeg = self.apply_angle_last
    self.frame += 1
    self.last_lat_active = CC.latActive
    return new_actuators, can_sends
