from opendbc.car import structs, Bus
from opendbc.can.parser import CANParser
from opendbc.car.common.conversions import Conversions as CV
from opendbc.car.volvo.values import DBC, CarControllerParams
from opendbc.car.interfaces import CarStateBase

GearShifter = structs.CarState.GearShifter
TransmissionType = structs.CarParams.TransmissionType


class CarState(CarStateBase):
  def __init__(self, CP):
    super().__init__(CP)
    self.cruise_enabled_prev = False
    self.cruise_last_disabled_frame = 0
    self.cruise_double_tap_active = False
    self.gas_pressed_prev = False
    self.CC_frame = 0 # CarController frame
    self.dispatch_lca_2_msg = False
    self.msg_pscm = {}
    self.msg_lca = {}
    self.msg_lca_2 = {}
    self.msg_lca_3 = {}
    self.msg_gear_position = {}
    self.pilot_assist_engaged = False
    self.msg_lca_5 = {}  # Formerly msg_speed_1
    self.msg_speed = {}
    self.msg_speed_2 = {}
    self.msg_0x1a = {}
    self.msg_egsm = {}
    self.msg_pscm_related = {}
    self.msg_lca_4 = {}
    self.msg_lca_6 = {}
    self.msg_lca_7 = {}
  def update(self, can_parsers) -> structs.CarState:
    cp_main = can_parsers[Bus.main]
    cp_pt = can_parsers[Bus.pt]
    cp_party = can_parsers[Bus.party]
    ret = structs.CarState()

    # car speed
    # Basic vehicle state from BUS1_SPEED on PT bus
    ret.vEgoRaw = cp_pt.vl["BUS1_SPEED"]["BUS1_SPEED"]
    ret.vEgo, ret.aEgo = self.update_speed_kf(ret.vEgoRaw)
    ret.standstill = ret.vEgoRaw <= 0.1 # 0.1 m/s

    # gas
    ret.gasPressed = cp_pt.vl["ECM_1"]["ACCELERATOR_PEDAL_POS"] > 20+1 # 20 baseline + 1 tolerance

    # brake
    #ret.brakePressed = bool(cp_main.vl["LCA_2"]["BRAKE_PEDAL_PRESSED_A"] or cp_main.vl["LCA_2"]["BRAKE_PEDAL_PRESSED_B"])
    # BRAKE_PEDAL_PRESSED_A goes active when user starts pressing brake pedal, but no brake light is on yet due to tolerance
    # BRAKE_PEDAL_PRESSED_B goes active when when the brake pedal is pressed above minimum threshold, brake light is on
    ret.brakePressed = cp_main.vl["LCA_2"]["BRAKE_PEDAL_PRESSED_B"] == 1
    ret.parkingBrake = False # TODO: add parking brake

    # stability control - becomes 3 when ESC intervenes (e.g., aquaplaning)
    ret.espActive = cp_main.vl["LCA_2"]["BYTE_1_MSBS_3"] == 3

    # steering wheel
    ret.steeringAngleDeg = cp_party.vl['PSCM']['PSCM_ANGLE_SENSOR'] # openpilot expects a negative value for a right turn
    #ret.steeringAngleDeg = cp_party.vl['SAS']['SAS_ANGLE_SENSOR']

    # Driver steering torque feedback (used for driver override detection)
    ret.steeringTorque = -cp_party.vl['DRIVER_INPUT']['STEERING_DRIVER_INPUT']  # Car right turn is negative, openpilot right turn is positive
    ret.steeringPressed = abs(cp_party.vl['DRIVER_INPUT']['STEERING_DRIVER_INPUT']) > 2

    # EPS status - placeholder until actual signal is found
    self.eps_active = True  # Assume EPS is active for now

    # cruise - double-tap detection (on-off-on within 500ms/50 frames / 1000ms/100 frames)
    cruise_raw = cp_pt.vl["BUS1_CRUISE_CONTROL"]["CRUISE_CONTROL_ENABLED"] == 1 or cp_pt.vl["BUS1_CRUISE_CONTROL"]["CRUISE_CONTROL_ENABLED_IDLE_TRAFFIC"] == 1

    # Check if double-tap cruise feature is enabled (bit 6 of alternativeExperience)
    use_double_tap = bool(self.CP.alternativeExperience & 64)

    if use_double_tap:
      # Detect on-off-on double-tap pattern
      if cruise_raw and not self.cruise_enabled_prev:
        # Just turned ON - check if we turned OFF recently (within 100 frames = 1000ms)
        if self.CC_frame - self.cruise_last_disabled_frame <= 100:
          self.cruise_double_tap_active = True
      elif not cruise_raw and self.cruise_enabled_prev:
        # Just turned OFF - clear the double-tap flag
        self.cruise_last_disabled_frame = self.CC_frame
        self.cruise_double_tap_active = False

    # cruiseState.enabled always reflects raw car state (must match panda safety)
    # blockPcmEnable prevents openpilot engagement until double-tap detected
    ret.cruiseState.enabled = cruise_raw
    ret.blockPcmEnable = use_double_tap and not self.cruise_double_tap_active

    self.cruise_enabled_prev = cruise_raw
    self.gas_pressed_prev = ret.gasPressed
    ret.cruiseState.available = True  # TODO: Determine actual availability
    ret.cruiseState.speed = 0  # TODO: Find cruise set speed (not required for lateral control)
    ret.cruiseState.nonAdaptive = False
    ret.cruiseState.standstill = ret.standstill # False # Todo: Find cruise control standstill signal

    # gear
    gearPosition = cp_main.vl['GEAR_POSITION']['GEAR_POSITION'] # 0: P; 1: R; 2: N; 3: D; 4: B;
    if gearPosition == 0:
      ret.gearShifter = GearShifter.park
    elif gearPosition == 1:
      ret.gearShifter = GearShifter.reverse
    elif gearPosition == 2:
      ret.gearShifter = GearShifter.neutral
    elif gearPosition == 3:
      ret.gearShifter = GearShifter.drive
    elif gearPosition == 4:
      ret.gearShifter = GearShifter.drive

    # blinkers TODO FlexRay
    ret.leftBlinker = False
    ret.rightBlinker = False

    # lock info TODO FlexRay
    ret.doorOpen = False # TODO: add door open
    ret.seatbeltUnlatched = False # TODO: add seatbelt unlatched

    # Store entire message dictionaries
    self.msg_pscm = cp_party.vl['PSCM']
    self.msg_lca = cp_main.vl['LCA']
    self.msg_lca_2 = cp_main.vl['LCA_2']
    self.msg_lca_3 = cp_main.vl['LCA_3']
    self.msg_lca_4 = cp_main.vl['LCA_4']
    self.msg_lca_5 = cp_main.vl['LCA_5']
    self.msg_lca_6 = cp_main.vl['LCA_6']
    self.msg_lca_7 = cp_main.vl['LCA_7']
    self.msg_speed = cp_main.vl['SPEED']
    self.msg_speed_2 = cp_main.vl['SPEED_2']
    self.msg_gear_position = cp_main.vl['GEAR_POSITION']
    self.msg_egsm = cp_party.vl['EGSM']
    self.msg_pscm_related = cp_party.vl['PSCM_RELATED']

    self.pilot_assist_engaged = cp_main.vl['LCA_2']['PILOT_ASSIST_ENGAGED'] == 1
    return ret

  @staticmethod
  def get_can_parsers(CP):
    return {
      Bus.main: CANParser(DBC[CP.carFingerprint][Bus.main], [], 0),
      Bus.pt: CANParser(DBC[CP.carFingerprint][Bus.pt], [], 1),
      Bus.party: CANParser(DBC[CP.carFingerprint][Bus.party], [], 2),
    }
