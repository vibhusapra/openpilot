# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Context

This is a **Polestar 2 fork of OpenPilot 0.10.0** with full steering control support via Paper's Volvo panda safety firmware. The primary purpose is enabling angle-based steering control on Volvo/Polestar CMA platform vehicles through a Comma 3/3X device.

**Key Characteristics:**
- Base: OpenPilot v0.10.0 (Comma 3 compatible, not 0.10.3)
- opendbc submodule: Paper's fork (`github.com/paper5590/opendbc.git`, branch `master-cma`)
- Custom addition: Polestar 2 support (missing from Paper's latest)
- Safety model: Volvo (ID 35) with angle-based steering via LCA_STEER (0x58) on party bus
- Status: Validated for car testing, not production-ready

## Project Structure

### Core Architecture

**openpilot** operates as a robotics OS with modular car-specific implementations:

- **`selfdrive/car/`** - OpenPilot 0.10.0 base car interface code (minimal in 0.10.0)
- **`opendbc_repo/opendbc/car/`** - Primary car interface implementations (moved here in newer versions)
  - `volvo/` - Volvo/Polestar specific implementation
    - `interface.py` - Car detection, parameter initialization, SafetyModel.volvo setup
    - `carstate.py` - CAN message parsing, vehicle state extraction
    - `carcontroller.py` - Control command generation
    - `values.py` - Vehicle specs (mass, wheelbase, steer ratio, Polestar 2 definition)
    - `fingerprints.py` - CAN message patterns for car detection
  - Safety model mappings, DBC database configurations per brand

- **`panda/`** - Panda firmware (safety-critical C code)
  - `board/safety/modes/volvo.h` - Volvo safety model with LCA steering validation
  - Contains TX/RX hooks for message validation and relay control

- **`selfdrive/controls/`** - Planning and control loops
- **`selfdrive/locationd/torqued.py` - Lateral tuning parameters (must have 'volvo' in ALLOWED_CARS)**
- **`selfdrive/pandad/`** - Interface to panda device
- **`selfdrive/modeld/`** - Neural network inference for perception

### Three-Bus CAN Architecture (Volvo/Polestar)

```
Bus 0 (Main):   VCU1 side - body/multiplex messages
Bus 1 (PT):     Powertrain/ECM side - engine control
Bus 2 (Party):  PSCM/BCM2 side - steering/brake system
                ↑ CRITICAL: LCA_STEER (0x58) sent here
```

## Development Setup

### Initial Setup

```bash
# Clone with submodules
git clone --recursive https://github.com/vibhusapra/openpilot.git
cd openpilot
git checkout polestar2-c3-v010

# Setup environment (Ubuntu 24.04 / macOS)
tools/op.sh setup
source .venv/bin/activate

# Build openpilot
scons -u -j$(nproc)  # Linux/Mac - uses available CPU cores
```

### Build Commands

| Task | Command |
|------|---------|
| **Full build** | `scons -u -j$(nproc)` |
| **Panda firmware only** (CRITICAL for car deployment) | `cd panda && scons -j4 && cd ..` |
| **Clean build** | `scons -c && scons -u -j$(nproc)` |
| **Single module** | `scons -u selfdrive/car/volvo/` |

### Running Tests

```bash
# All tests (excluding slow tests)
pytest

# Specific test file
pytest selfdrive/car/tests/test_car_interfaces.py -v

# Volvo/Polestar tests only
pytest selfdrive/car/tests/ -k volvo -v

# With specific markers
pytest -m "not slow" --tb=short

# Single test function
pytest path/to/test_file.py::test_function_name -v
```

### Validation and Linting

```bash
# Validate Polestar 2 integration (custom script)
python3 validate_polestar2.py

# Type checking (mypy)
mypy selfdrive/car/volvo/ --ignore-missing-imports

# Code style (ruff)
ruff check selfdrive/car/volvo/
ruff format selfdrive/car/volvo/
```

## Polestar 2 Implementation Details

### Vehicle Definition

**File**: `opendbc_repo/opendbc/car/volvo/values.py`

```python
class CAR(Platforms):
  POLESTAR_2 = VolvoCMAPlatformConfig(
    [VolvoCarDocs("Polestar 2 2020-2024")],
    CarSpecs(
      mass=2123,              # kg - Long Range Dual Motor
      wheelbase=2.735,        # meters
      steerRatio=15.8,        # same as XC40 (CMA platform)
      centerToFrontRatio=0.52,
    ),
  )
```

### Car Detection Flow

1. **Fingerprinting** (`opendbc/car/volvo/fingerprints.py`)
   - 35 CAN messages define Polestar 2
   - Key message: 0x58 (LCA_STEER steering command)
   - Comma 3 CAN logger matches messages against fingerprints

2. **Safety Model Loading** (`opendbc/car/volvo/interface.py`)
   - `SafetyModel.volvo` (ID 35) loaded into panda
   - **CRITICAL**: Must build panda firmware for safety model to work

3. **Control Initialization**
   - CarController generates angle-based steering commands
   - LCA_STEER messages sent at 100Hz to party bus (0x58, 8 bytes)
   - PSCM translates to steering actuator movement

### Critical Integration Points

**OpenPilot modifications** (must be present):

1. `selfdrive/car/card.py` - Add Volvo feature flags
   ```python
   if self.CP.carFingerprint.startswith("VOLVO") or self.CP.carFingerprint.startswith("POLESTAR"):
       if self.params.get_bool("VolvoDoubleTapCruise"):
           self.CP.alternativeExperience |= 64
       if self.params.get_bool("VolvoSpoofPAHandsOnWheel"):
           self.CP.alternativeExperience |= 128
   ```

2. `selfdrive/locationd/torqued.py` - Add volvo to ALLOWED_CARS
   ```python
   ALLOWED_CARS = ['toyota', 'hyundai', 'rivian', 'honda', 'volvo']
   ```

3. `launch_openpilot.sh` - Enable environment variable support
   ```bash
   if [ -f env.sh ]; then source env.sh; fi
   ```

## Submodule Management

The opendbc submodule is **critical** and points to Paper's fork:

```bash
# Check submodule status
git submodule status
# Output: 1f38689f... opendbc_repo [master-cma]

# Update submodule to latest
git submodule update --remote

# View opendbc commit history
git -C opendbc_repo log --oneline -5

# Change opendbc branch (if needed)
git -C opendbc_repo checkout <branch>
git add opendbc_repo
git commit -m "Bump opendbc to <branch>"
```

**Important**: The opendbc submodule MUST have the Volvo safety firmware (`opendbc/safety/modes/volvo.h`) - this is why Paper's fork is required.

## Common Tasks

### Adding a New Message to Polestar 2 Fingerprint

1. **Edit fingerprints.py** - Add CAN message ID (decimal, not hex)
   ```python
   CAR.POLESTAR_2: [{
       0x58: 8,  # LCA_STEER
       # ... add new message ID here
   }]
   ```

2. **Verify in DBC** - Check `opendbc_repo/opendbc/car/volvo/*.dbc` for message definition

3. **Test detection** - Run Polestar 2 validation script
   ```bash
   python3 validate_polestar2.py
   ```

### Modifying Safety Model

**File**: `opendbc_repo/opendbc/safety/modes/volvo.h`

The safety model is written in C and runs on panda (embedded device):

1. **Edit volvo.h** for message validation rules
2. **Rebuild panda firmware**:
   ```bash
   cd panda && scons -j4
   ```
3. **Deploy to car** - The new firmware loads on next reboot

### Debugging Car Detection Issues

```bash
# On Comma 3, force Polestar 2 detection
export FINGERPRINT="POLESTAR_2"
export SKIP_FW_QUERY=1

# Monitor logs for detection messages
tail -f /data/log/swaglog.* | grep -i "volvo\|polestar\|safety"

# Check what safety model was loaded
grep safetyModel /data/log/swaglog.* | tail -5
```

## Deployment to Comma 3

### Quick Deploy (Recommended for Testing)

```bash
ssh comma@<COMMA_IP>

# Direct swap
sudo systemctl stop openpilot
cd /data
rm -rf openpilot
git clone --recursive https://github.com/vibhusapra/openpilot.git openpilot
cd openpilot
git checkout polestar2-c3-v010

# BUILD PANDA FIRMWARE (CRITICAL!)
cd panda
scons -j4

# Reboot
sudo reboot
```

### Force Polestar 2 Detection (if fingerprinting fails)

```bash
# On Comma 3, add to env.sh
echo 'export FINGERPRINT="POLESTAR_2"' >> /data/openpilot/env.sh
echo 'export SKIP_FW_QUERY=1' >> /data/openpilot/env.sh
sudo reboot
```

## Testing Strategy

**Laptop testing** (before car deployment):
- ✅ Python imports and code structure
- ✅ Fingerprint configuration
- ✅ Car interface instantiation
- ✅ Safety firmware presence
- ❌ CAN communication (requires car + panda)
- ❌ Actual steering control (requires car)

**Use `validate_polestar2.py`** to verify all laptop-testable components pass.

## Safety Considerations

**Paper's Volvo safety model is "very relaxed"** - intended for development only:

- ✅ Basic frame ID validation
- ✅ Relay control (blocks stock LCA when openpilot active)
- ❌ No angle rate limiting (unlike Toyota LTA)
- ❌ No torque validation
- ⚠️ Minimal safety enforcement

**Never deploy without:**
1. Testing in safe environment (empty parking lot)
2. Keeping hands near steering wheel
3. Being ready to take manual control
4. Starting at low speeds (< 25 mph)
5. Monitoring logs for safety errors

## Key Files Reference

| File | Purpose |
|------|---------|
| `opendbc_repo/opendbc/car/volvo/interface.py` | Car initialization, safety model setup |
| `opendbc_repo/opendbc/car/volvo/carcontroller.py` | Generates steering commands |
| `opendbc_repo/opendbc/car/volvo/carstate.py` | Parses CAN messages |
| `opendbc_repo/opendbc/car/volvo/values.py` | Vehicle specs, Polestar 2 config |
| `opendbc_repo/opendbc/car/volvo/fingerprints.py` | CAN patterns for detection |
| `opendbc_repo/opendbc/safety/modes/volvo.h` | Panda safety firmware (C) |
| `selfdrive/car/card.py` | Volvo feature flags integration |
| `selfdrive/locationd/torqued.py` | Must include 'volvo' in ALLOWED_CARS |

## Documentation

- **FINAL_CHECKLIST.md** - Pre-deployment verification checklist
- **LAPTOP_TEST_RESULTS.md** - Full test report from laptop validation
- **DEPLOYMENT_FINAL.md** - Step-by-step Comma 3 deployment guide
- **validate_polestar2.py** - Automated integration validation script

## Important URLs

- Fork: https://github.com/vibhusapra/openpilot
- Branch: `polestar2-c3-v010`
- Paper's opendbc: https://github.com/paper5590/opendbc.git (master-cma)
- OpenPilot docs: https://docs.comma.ai

## Common Pitfalls

1. **Steering doesn't work** - Panda firmware wasn't built (`cd panda && scons -j4`)
2. **Car not detected** - Use `FINGERPRINT="POLESTAR_2"` env variable
3. **Relay errors** - Check LCA_STEER message is being sent to correct bus (party bus)
4. **Safety model errors** - Verify volvo.h is in opendbc_repo and panda built
5. **Git submodule issues** - Always `git clone --recursive` and keep opendbc_repo on master-cma

---

**Status**: Validated for car testing - all pre-deployment checks passed (8/8)
**Last Updated**: 2026-01-22
**Risk Level**: Moderate (relaxed safety model, test carefully)
