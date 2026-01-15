## Testing Commands for Car Ports

To test a car port in openpilot, you have several commands and utilities available:

### 1. Car Interface Testing
```bash
pytest selfdrive/car/tests/test_car_interfaces.py -k [brand]
```
Tests car interfaces for common bugs without requiring a route. Replace `[brand]` with your target brand (e.g., `subaru`, `toyota`). [1](#0-0)

### 2. Route-Based Model Testing
```bash
python3 tools/car_porting/test_car_model.py [route_or_segment_name] --car [car_model]
python3 tools/car_porting/test_car_model.py --car VOLVO_XC40_RECHARGE [route_or_segment_name]
```
Runs comprehensive tests against a specific route to check for missing signals, blocked panda messages, and safety mismatches. [2](#0-1)

### 3. CAN Signal Analysis
```bash
tools/cabana/cabana [route]
```
Launches Cabana to view your car's CAN signals through DBC files, which openpilot uses to parse and create messages. [3](#0-2)

### 4. Fingerprint Management
```bash
python3 tools/car_porting/auto_fingerprint.py [route] [platform]
```
Automatically inserts FW fingerprints from the platform into the correct place in fingerprints.py. [4](#0-3)

### 5. Stimulus-Response Testing
For physical testing, follow the scripted test drives in `docs/car-porting/reverse-engineering.md` which include:
- Stationary ignition-only tests (doors, pedals, gearshift)
- ACC button tests
- Steering angle and torque tests
- Low speed parking lot tests
- High speed highway tests [5](#0-4)

## Notes

- The `test_car_model.py` script uses `selfdrive/car/tests/test_models.py` internally for the actual test implementation [6](#0-5)
- For route-based testing, you'll need a route with your car's data uploaded to comma connect
- The car interface tests use hypothesis-based fuzzy testing to find edge cases [7](#0-6)
- All testing commands should be run from the openpilot root directory with the virtual environment activated

Wiki pages you might want to explore:
- [Overview (commaai/openpilot)](/wiki/commaai/openpilot#1)

### Citations

**File:** tools/car_porting/README.md (L9-16)
```markdown
### [Cabana](/tools/cabana/README.md)

View your car's CAN signals through DBC files, which openpilot uses to parse and create messages that talk to the car.

Example:
```bash
> tools/cabana/cabana '1bbe6bf2d62f58a8|2022-07-14--17-11-43'
```
```

**File:** tools/car_porting/README.md (L18-26)
```markdown
### [tools/car_porting/auto_fingerprint.py](/tools/car_porting/auto_fingerprint.py)

Given a route and platform, automatically inserts FW fingerprints from the platform into the correct place in fingerprints.py

Example:
```bash
> python3 tools/car_porting/auto_fingerprint.py '1bbe6bf2d62f58a8|2022-07-14--17-11-43' 'OUTBACK'
Attempting to add fw version for:  OUTBACK
```
```

**File:** tools/car_porting/README.md (L28-40)
```markdown
### [selfdrive/car/tests/test_car_interfaces.py](/selfdrive/car/tests/test_car_interfaces.py)

Finds common bugs for car interfaces, without even requiring a route.


#### Example: Typo in signal name
```bash
> pytest selfdrive/car/tests/test_car_interfaces.py -k subaru  # replace with the brand you are working on

=====================================================================
FAILED selfdrive/car/tests/test_car_interfaces.py::TestCarInterfaces::test_car_interfaces_165_SUBARU_LEGACY_7TH_GEN - KeyError: 'CruiseControlOOPS'

```
```

**File:** tools/car_porting/test_car_model.py (L7-7)
```python
from openpilot.selfdrive.car.tests.test_models import TestCarModel
```

**File:** tools/car_porting/test_car_model.py (L21-36)
```python
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Test any route against common issues with a new car port. " +
                                               "Uses selfdrive/car/tests/test_models.py")
  parser.add_argument("route_or_segment_name", help="Specify route to run tests on")
  parser.add_argument("--car", help="Specify car model for test route")
  args = parser.parse_args()
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

  sr = SegmentRange(args.route_or_segment_name)

  test_routes = [CarTestRoute(sr.route_name, args.car, segment=seg_idx) for seg_idx in sr.seg_idxs]
  test_suite = create_test_models_suite(test_routes)

  unittest.TextTestRunner().run(test_suite)
```

**File:** docs/car-porting/reverse-engineering.md (L1-85)
```markdown
# Stimulus-Response Tests

These are example test drives that can help identify the CAN bus messaging necessary for ADAS control. Each scripted
test should be done in a separate route (ignition cycle). These tests are a guide, not necessarily exhaustive.

While testing, constant power to the comma device is highly recommended, using [comma power](https://comma.ai/shop/comma-power) if
necessary to make sure all test activity is fully captured and for ease of uploading. If constant power isn't
available, keep the ignition on for at least one minute after your test to make sure power loss doesn't result
in loss of the last minute of testing data.

## Stationary ignition-only tests, part 1

1. Ignition on, but don't start engine, remain in Park
2. Open and close each door in a defined order: driver, passenger, rear left, rear right
3. Re-enter the vehicle, close the driver's door, and fasten the driver's seatbelt
4. Slowly press and release the accelerator pedal 3 times
5. Slowly press and release the brake pedal 3 times
6. Hold the brake and move the gearshift to reverse, then neutral, then drive, then sport/eco/etc if applicable
7. Return to Park, ignition off

Brake-pressed information may show up in several messages and signals, both as on/off states and as a percentage or
pressure. It may reflect a switch on the driver's brake pedal, or a pressure-threshold state, or signals to turn on
the rear brake lights. Start by identifying all the potential signals, and confirm while driving with ACC later.

Locate signals for all four door states if possible, but some cars only expose the driver's door state on the ADAS bus.
Driver/passenger door signals may or may not change positions for LHD vs RHD cars. For cars where only the driver's
door signal is available, the same signal may follow the driver.

## Stationary ignition-only tests, part 2

1. Ignition on, but don't start engine, remain in Park
2. Press each ACC button in a defined order: main switch on/off, set, resume, cancel, accel, decel, gap adjust
3. Set the left turn signal for about five seconds
4. Operate the left turn signal one time in its touch-to-pass mode
5. Set the right turn signal for about five seconds
6. Operate the right turn signal one time in its touch-to-pass mode
7. Set the hazard / emergency indicator switch for about five seconds
8. Ignition off

Your vehicle may have a momentary-press main ACC switch or a physical toggle that remains set. Actual ACC engagement
isn't necessary for purposes of detecting the ACC button presses.

## Steering angle and steering torque tests

Power steering should be available. On ICE cars, engine RPM may be present.

1. Ignition on, start engine if applicable, remain in Park
2. Rotate the steering wheel as follows, with a few seconds pause between each step
   * Start as close to exact center as possible
   * Turn to 45 degrees right and hold
   * Turn to 90 degrees right and hold
   * Turn to 180 degrees right and hold
   * Turn to full lock right and hold, with firm pressure against lock
   * Release the wheel and allow it to bounce back slightly from lock
   * Turn to 180 degrees left and hold
   * Return to center and release
3. Ignition off

Performing the full test to the right, followed by an abbreviated test to the left, helps give additional confirmation
of signal scale, and sign/direction for both the steering wheel angle and driver input torque signals.

## Low speed / parking lot driving tests

Before this test, drive to a place like an empty parking lot where you are free to drive in a series of curves.

1. Ignition on, start engine if applicable, prepare to drive
2. Slowly (10-20mph at most) drive a figure-8 if possible, or at least one sharp left and one sharp right.
3. Come to a complete stop
4. When and where safe, drive in reverse for a short distance (10-15 feet)
5. Park the car in a safe place, ignition off

## High speed / highway driving tests

Select a place and time where you can safely set cruise control at normal travel speeds with little interference from
traffic ahead, and safely test the response of your factory lane guidance system.

1. Ignition on, start engine if applicable, prepare to drive
2. When safely able, engage adaptive cruise control below 50 mph
3. When safely able, use the ACC buttons to accelerate to 50mph, then 55mph, then 60mph
4. Disengage adaptive cruise
5. When safely able, allow your factory lane guidance to prevent lane departures, 2-3 times on both the left and right

The series of setpoints can be adjusted to local traffic regulations, and of course metric units. The specific cruise
setpoints are useful for locating the ACC HUD signals later, and confirming their precise scaling. When the car reaches
and holds the setpoint, that can also provide additional confirmation of wheel speed scaling.
```

**File:** selfdrive/car/tests/test_car_interfaces.py (L22-27)
```python
  # FIXME: Due to the lists used in carParams, Phase.target is very slow and will cause
  #  many generated examples to overrun when max_examples > ~20, don't use it
  @parameterized.expand([(car,) for car in sorted(PLATFORMS)] + [MOCK.MOCK])
  @settings(max_examples=MAX_EXAMPLES, deadline=None,
            phases=(Phase.reuse, Phase.generate, Phase.shrink))
  @given(data=st.data())
```
