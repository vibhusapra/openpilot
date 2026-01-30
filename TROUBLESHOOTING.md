# Polestar 2 OpenPilot Troubleshooting Guide

This guide documents known issues and their solutions for the Polestar 2 OpenPilot fork.

## Quick Diagnostics

Run the test script to check your setup:
```bash
python scripts/test_polestar2.py
```

## Known Issues and Solutions

### Issue 1: Python Version Incompatibility

**Symptom:**
```
ImportError: cannot import name 'ReprEnum' from 'enum'
```

**Cause:** The code requires Python 3.11+ but you have an older version.

**Solution:**

On **macOS**:
```bash
brew install python@3.11
python3.11 --version  # Verify installation
```

On **Ubuntu/Debian**:
```bash
sudo apt update
sudo apt install python3.11 python3.11-venv python3.11-dev
```

On **Fedora/RHEL**:
```bash
sudo dnf install python3.11 python3.11-devel
```

Then create a virtual environment:
```bash
python3.11 -m venv .venv
source .venv/bin/activate
```

---

### Issue 2: Git LFS Assets Not Downloaded

**Symptom:**
- Binary files (fonts, images, models) appear corrupted or are tiny text files
- Build fails with missing asset errors
- `git status` shows many modified binary files

**Cause:** Git LFS files weren't properly pulled.

**Solution:**
```bash
# Install git-lfs if not present
brew install git-lfs  # macOS
# or
sudo apt install git-lfs  # Ubuntu

# Initialize and pull LFS files
git lfs install
git lfs pull

# If that doesn't work, try:
git lfs fetch --all
git lfs checkout
```

---

### Issue 3: Submodule Not Initialized

**Symptom:**
```
ModuleNotFoundError: No module named 'opendbc'
```
or empty directories under `opendbc_repo/`, `panda/`, etc.

**Cause:** Git submodules weren't initialized.

**Solution:**
```bash
git submodule update --init --recursive

# Verify submodules are present
ls -la opendbc_repo/opendbc/car/volvo/
```

---

### Issue 4: controlsd Crashes on Startup

**Symptom:**
```
AttributeError: 'CarParams.LongitudinalPIDTuning' object has no attribute 'kf'
```
or similar crash in controlsd related to longitudinal tuning.

**Cause:** Missing longitudinalTuning parameters in interface.py (v0.10.0 compatibility issue).

**Solution:**

Ensure `opendbc_repo/opendbc/car/volvo/interface.py` has these lines:
```python
# Longitudinal tuning parameters (required for controlsd in v0.10.0)
ret.longitudinalTuning.kpBP = [0., 35.]
ret.longitudinalTuning.kpV = [1.2, 0.8]
ret.longitudinalTuning.kiBP = [0., 35.]
ret.longitudinalTuning.kiV = [0.18, 0.12]
```

---

### Issue 5: Car Not Fingerprinting / "No supported car detected"

**Symptom:**
- OpenPilot starts but shows "No supported car detected"
- Fingerprinting fails even with harness connected

**Possible Causes and Solutions:**

1. **Wrong CAN bus:** Ensure harness is connected to VCU1 → PSCM (MID 1 CAN)

2. **Fingerprint mismatch:** Your car's CAN IDs might differ from the stored fingerprint. Capture your car's fingerprint:
   ```bash
   # On Comma device via SSH
   cd /data/openpilot
   python selfdrive/debug/fingerprint_from_log.py
   ```

3. **Missing fingerprint entry:** Check that POLESTAR_2 is in `fingerprints.py`:
   ```python
   # opendbc_repo/opendbc/car/volvo/fingerprints.py
   CAR.POLESTAR_2: [
     {
       21: 8, 22: 8, 23: 8, ...  # CAN message IDs
     }
   ],
   ```

---

### Issue 6: Steering Not Engaging / Dashcam Only Mode

**Symptom:**
- OpenPilot shows "Dashcam Mode"
- Steering commands are not sent to the car

**Cause:** SafetyModel set to `noOutput` or `dashcamOnly = True`.

**Solution:**

Edit `opendbc_repo/opendbc/car/volvo/interface.py`:
```python
# Use Volvo safety model (not noOutput)
ret.safetyConfigs = [get_safety_config(structs.CarParams.SafetyModel.volvo)]

# Enable steering control
ret.dashcamOnly = False
```

---

### Issue 7: Build Errors with scons

**Symptom:**
```
scons: *** No SConstruct file found.
```

**Solution:**
Make sure you're in the openpilot root directory:
```bash
cd /path/to/openpilot
ls SConstruct  # Should exist
scons -j$(nproc)
```

**Symptom:**
```
error: unknown target 'arm64-apple-darwin'
```

**Cause:** Building for wrong architecture.

**Solution:**
```bash
# On Mac, ensure you're building for host
scons -j$(nproc) --mac
```

---

### Issue 8: Panda Safety Compilation Errors

**Symptom:**
```
error: redefinition of 'dlc_to_len'
```

**Cause:** Missing include guards in safety headers.

**Solution:**

The fix should already be in your fork. If not, ensure `opendbc_repo/opendbc/safety/safety.h` has:
```c
#ifndef DLC_TO_LEN_DEFINED
#define DLC_TO_LEN_DEFINED
// dlc_to_len definition here
#endif
```

---

### Issue 9: opendbc Import Errors After Pull

**Symptom:**
```
ImportError: cannot import name 'X' from 'opendbc.car'
```

**Cause:** opendbc submodule out of sync with main repo.

**Solution:**
```bash
cd opendbc_repo
git fetch origin
git checkout master-cma  # or whatever branch you need
git pull
cd ..
git add opendbc_repo
git commit -m "Update opendbc submodule"
```

---

## Environment Setup Checklist

Before building, verify:

- [ ] Python 3.11+ installed: `python3 --version`
- [ ] Git LFS installed: `git-lfs --version`
- [ ] LFS assets pulled: `git lfs pull`
- [ ] Submodules initialized: `git submodule status`
- [ ] Virtual environment activated: `which python` shows `.venv`
- [ ] scons installed: `scons --version`
- [ ] PYTHONPATH set: `echo $PYTHONPATH`

## Testing Your Setup

```bash
# 1. Set up environment
source .venv/bin/activate
export PYTHONPATH="$PWD:$PWD/opendbc_repo:$PYTHONPATH"

# 2. Run test script
python scripts/test_polestar2.py

# 3. Try building
scons -j$(nproc)
```

## Getting Help

- Check the Discord thread mentioned in the chat logs
- Review Paper's notes on known issues (angle cap, driver override asymmetry)
- Check GitHub issues on the fork repository

## Architecture Reference

```
openpilot/
├── selfdrive/
│   ├── car/
│   │   ├── card.py          # Main car process
│   │   └── car_specific.py  # Car-specific helpers
│   └── controls/
│       └── controlsd.py     # Control loop (needs longitudinalTuning)
├── opendbc_repo/            # Submodule with car interfaces
│   └── opendbc/
│       ├── car/volvo/       # Volvo/Polestar implementation
│       │   ├── interface.py
│       │   ├── carstate.py
│       │   ├── carcontroller.py
│       │   ├── values.py
│       │   └── fingerprints.py
│       ├── dbc/             # CAN database files
│       └── safety/modes/    # Safety model headers
│           └── volvo.h
└── panda/                   # CAN hardware interface
```
