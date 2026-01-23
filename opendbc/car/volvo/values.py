from dataclasses import dataclass, field

from opendbc.car.structs import CarParams
from opendbc.car import Bus, CarSpecs, DbcDict, PlatformConfig, Platforms
from opendbc.car.lateral import AngleSteeringLimits
from opendbc.car.docs_definitions import CarDocs, CarHarness, CarParts
from opendbc.car.fw_query_definitions import FwQueryConfig, Request, StdQueries

Ecu = CarParams.Ecu


class CarControllerParams:
  STEER_STEP = 1  # 100 Hz LCA command frequency (controlsd runs at 100 Hz)

  # Angle limits for rate limiting
  ANGLE_LIMITS: AngleSteeringLimits = AngleSteeringLimits(
    390, # deg - max steering angle
    ([0., 5., 25.], [2.5, 1.5, .2]),  # rate up limits at different speeds
    ([0., 5., 25.], [5., 2., .3]),    # rate down limits at different speeds
  )


@dataclass
class VolvoCarDocs(CarDocs):
  package: str = "Pilot Assist & Adaptive Cruise Control"
  car_parts: CarParts = field(default_factory=CarParts.common([CarHarness.custom]))


@dataclass
class VolvoCMAPlatformConfig(PlatformConfig):
  dbc_dict: DbcDict = field(default_factory=lambda: {
    Bus.pt: 'volvo_front_1_cma',
    Bus.main: 'volvo_mid_1',
    Bus.party: 'volvo_mid_1',
  })
  
@dataclass
class VolvoSPAPlatformConfig(PlatformConfig):
  dbc_dict: DbcDict = field(default_factory=lambda: {
    Bus.pt: 'volvo_front_1_spa',
    Bus.main: 'volvo_mid_1',
    Bus.party: 'volvo_mid_1',
  })


class CAR(Platforms):
  VOLVO_XC40_RECHARGE = VolvoCMAPlatformConfig(
    [VolvoCarDocs("Volvo XC40 Recharge 2021-2023")],
    CarSpecs(
      mass=2170,
      wheelbase=2.702,
      steerRatio=15.8,
      centerToFrontRatio=0.52,
    ),
  )
  POLESTAR_2 = VolvoCMAPlatformConfig(
    [VolvoCarDocs("Polestar 2 2020-2024")],
    CarSpecs(
      mass=2123,  # Long Range Dual Motor variant
      wheelbase=2.735,
      steerRatio=15.8,  # Same as XC40 (CMA platform)
      centerToFrontRatio=0.52,
    ),
  )
  #TODO update car specs
  VOLVO_S60_RECHARGE = VolvoSPAPlatformConfig(
    [VolvoCarDocs("Volvo S60 Recharge 2024")],
    CarSpecs(
      mass=2170,
      wheelbase=2.702,
      steerRatio=15.8,
      centerToFrontRatio=0.52,
    ),
  )


# FW Query configuration for Volvo CMA platform
# FW_QUERY_CONFIG = FwQueryConfig(
#   requests=[
#     Request(
#       [StdQueries.TESTER_PRESENT_REQUEST, StdQueries.UDS_VERSION_REQUEST],
#       [StdQueries.TESTER_PRESENT_RESPONSE, StdQueries.UDS_VERSION_RESPONSE],
#       bus=0,
#     ),
#   ],
# )
FW_QUERY_CONFIG = FwQueryConfig(
  requests=[]
)

DBC = CAR.create_dbc_map()
