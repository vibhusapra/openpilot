# Polestar 2 Development Tools

This folder contains all custom tools and documentation for the Polestar 2 OpenPilot port.

## 🚀 Quick Start

```bash
# Quick validation (5 seconds)
./quick_test.sh

# Full CI validation before deployment
./ci_test.sh

# Validate Polestar 2 integration specifically
python3 validate_polestar2.py
```

## 📁 Files Overview

### Testing Scripts
- **`ci_test.sh`** - Comprehensive CI validation suite (15 steps, ~3-5 min)
- **`quick_test.sh`** - Fast validation check (< 5 seconds)
- **`test_build.sh`** - Build validation script
- **`test_launch.py`** - Python import and launch tests
- **`validate_polestar2.py`** - Polestar 2 specific validation
- **`test_polestar2.py`** - Initial Polestar 2 tests
- **`test_structure.py`** - Code structure validation

### Documentation
- **`CLAUDE.md`** - Guide for Claude AI to understand this codebase
- **`FINAL_CHECKLIST.md`** - Pre-deployment checklist
- **`LAPTOP_TEST_RESULTS.md`** - Local test results documentation
- **`COMMA3_DEPLOYMENT.md`** - Comma 3 deployment instructions

## 🎯 Key Changes Made

### OpenPilot Core Changes
1. **`selfdrive/car/card.py`** - Added Volvo feature flags
2. **`selfdrive/locationd/torqued.py`** - Added 'volvo' to ALLOWED_CARS
3. **`selfdrive/pandad/panda.h`** - Fixed PANDA_BUS_CNT and C++ linkage
4. **`launch_openpilot.sh`** - Added env.sh support
5. **`env.sh`** - Environment variable loader

### opendbc Submodule (Paper's Fork)
- Using: `https://github.com/paper5590/opendbc.git`
- Branch: `master-cma`
- Contains Volvo panda safety firmware (11KB volvo.h)
- Added Polestar 2 to CAR enum and fingerprints

### Panda Fixes
- Added dlc_to_len include guards
- Fixed C++ extern "C" wrappers
- Fixed SConscript Dir/File issues

## ✅ What's Working

- **Build**: Core OpenPilot builds successfully
- **Polestar 2**: In CAR enum with 35 CAN messages
- **Safety**: Volvo safety model (ID 35) present
- **Steering**: Angle-based control via LCA_STEER (0x58)
- **Integration**: All required files in place

## ⚠️ Known Limitations

- Paper's safety model is "very relaxed" (development only)
- No angle rate limiting
- Minimal safety validation
- Test in safe environment first!

## 🚗 Deployment to Comma 3

```bash
ssh comma@<COMMA_IP>
cd /data
rm -rf openpilot
git clone --recursive https://github.com/vibhusapra/openpilot.git openpilot
cd openpilot
git checkout polestar2-c3-v010

# CRITICAL: Update opendbc to Paper's fork
cd opendbc_repo
git remote set-url origin https://github.com/paper5590/opendbc.git
git fetch origin
git checkout master-cma
git pull origin master-cma
cd ..

# Build panda firmware (CRITICAL!)
cd panda
scons -j4
cd ..

# Build OpenPilot
scons -u -j4

# Force Polestar 2 detection if needed
export FINGERPRINT="POLESTAR_2"
export SKIP_FW_QUERY=1
```

## 📊 Test Results

Run `./ci_test.sh` to see:
- Environment detection
- Repository validation
- Dependency checks
- Submodule verification
- Build validation
- Python import tests
- Safety model checks
- Performance metrics

All tests should pass before deploying to car!

---

**Repository**: https://github.com/vibhusapra/openpilot
**Branch**: polestar2-c3-v010
**Status**: ✅ Build validated, ready for car testing