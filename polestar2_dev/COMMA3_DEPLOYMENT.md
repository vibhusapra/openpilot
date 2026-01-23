# Polestar 2 OpenPilot 0.10.0 - Comma 3 Deployment Guide

## Overview
This is a fork of OpenPilot 0.10.0 with Volvo/Polestar 2 support added from Paper's implementation.

## What's Included
- ✅ OpenPilot 0.10.0 (Comma 3 compatible)
- ✅ Volvo implementation from Paper's fork
- ✅ Polestar 2 support (added to Paper's XC40 and S60)
- ✅ All necessary DBC files
- ✅ Proper car registration and discovery

## Pre-Deployment Checklist
- [ ] All structural tests pass (run `python3 test_structure.py`)
- [ ] Git repository is ready for pushing
- [ ] Comma 3 is accessible via SSH
- [ ] Polestar 2 is ready for testing

## Step 1: Push to GitHub

First, push this fork to your GitHub:

```bash
cd /tmp/openpilot_v010_polestar

# Initialize git if needed
git init

# Add your GitHub remote
git remote add origin https://github.com/YOUR_USERNAME/openpilot-polestar2.git

# Add all files
git add .
git commit -m "OpenPilot 0.10.0 with Polestar 2 support

- Based on stock OpenPilot 0.10.0 (Comma 3 compatible)
- Added Volvo implementation from Paper's fork
- Added Polestar 2 support (CMA platform)
- Uses noOutput safety mode (dashcam mode initially)"

# Push to GitHub
git push -u origin polestar2-c3-v0.10.0
```

## Step 2: Deploy to Comma 3

### Option A: Fresh Install
```bash
# SSH into your Comma 3
ssh comma@YOUR_COMMA_IP

# Stop openpilot
sudo systemctl stop openpilot

# Backup existing openpilot (optional)
cd /data
mv openpilot openpilot_backup

# Clone your fork
git clone https://github.com/YOUR_USERNAME/openpilot-polestar2.git openpilot
cd openpilot
git checkout polestar2-c3-v0.10.0

# Initialize submodules
git submodule update --init --recursive

# Reboot
sudo reboot
```

### Option B: Test Install (Parallel)
```bash
# SSH into your Comma 3
ssh comma@YOUR_COMMA_IP

# Clone to test directory
cd /data
git clone https://github.com/YOUR_USERNAME/openpilot-polestar2.git openpilot_test
cd openpilot_test
git checkout polestar2-c3-v0.10.0

# Initialize submodules
git submodule update --init --recursive

# Stop current openpilot and test
sudo systemctl stop openpilot
cd /data/openpilot_test

# Force Polestar 2 fingerprint (since CMA cars share messages)
export FINGERPRINT="POLESTAR_2"

# Run openpilot
./launch_openpilot.sh
```

## Step 3: Verify Installation

After reboot, check that Polestar 2 is detected:

```bash
# SSH back in
ssh comma@YOUR_COMMA_IP

# Check logs
grep -i "polestar\|volvo" /data/log/swaglog.*

# Monitor startup
tmux attach -t comma
```

## Expected Behavior

### On First Start:
1. **Car Detection**: Should detect as "POLESTAR_2"
2. **Safety Mode**: Will use "noOutput" (dashcam only)
3. **UI**: Should show "openpilot available" when ready
4. **Steering**: Currently DISABLED (noOutput mode)

### What Works:
- ✅ Car detection and fingerprinting
- ✅ CAN message parsing
- ✅ Basic car state reading
- ✅ Dashcam functionality

### What Doesn't Work Yet:
- ❌ Steering control (needs Volvo safety model in panda)
- ❌ ACC control (needs further development)
- ⚠️ Some Polestar 2 specific features may need tuning

## Troubleshooting

### "Unsupported Device" Error
If you get this error, you're not on 0.10.0. Make sure:
- You're using this fork (based on 0.10.0)
- NOT Paper's 0.10.3 fork
- NOT sunnypilot

### Car Not Detected
```bash
# Force fingerprint
echo 'export FINGERPRINT="POLESTAR_2"' >> /data/params/d/LaunchEnv
sudo reboot
```

### Import Errors
Check that opendbc submodule is initialized:
```bash
cd /data/openpilot
git submodule update --init --recursive
```

### View Real-Time Logs
```bash
# In one terminal
ssh comma@YOUR_COMMA_IP
tmux attach -t comma

# In another terminal
ssh comma@YOUR_COMMA_IP
tail -f /data/log/swaglog.* | grep -i "volvo\|polestar\|error"
```

## Development Notes

### Key Files Modified:
1. **opendbc_repo/opendbc/car/volvo/** - All Volvo implementation
2. **opendbc_repo/opendbc/car/values.py** - Added Volvo registration
3. **opendbc_repo/opendbc/dbc/volvo_*.dbc** - CAN database files

### Safety Mode:
Currently using `noOutput` which means:
- Openpilot can read all CAN messages
- Cannot send steering commands
- Safe for initial testing

To enable steering in the future:
1. Volvo safety model needs to be added to panda firmware
2. Change from `noOutput` to `volvo` in interface.py
3. Test thoroughly before road use

### Platform Details:
- **Polestar 2**: CMA platform (same as XC40)
- **CAN Buses**: PT (powertrain), Main (MID 1)
- **Steering**: Angle-based control
- **Mass**: 2123 kg
- **Wheelbase**: 2.735 m

## Next Steps

1. **Test Basic Detection**: Verify car is detected correctly
2. **Monitor CAN Messages**: Check that all expected messages are received
3. **Validate CarState**: Ensure speed, steering angle, etc. are correct
4. **Develop Safety Model**: Work on panda firmware for steering control
5. **Tune Parameters**: Adjust steering delays, limits, etc.
6. **Community Testing**: Share with other Polestar 2 owners

## Support

- Report issues: Create issue with logs
- Join Discord: Share findings with community
- Contribute: Submit PRs for improvements

## ⚠️ SAFETY WARNING

This is EXPERIMENTAL software. Current configuration:
- **Dashcam mode only** (no steering control)
- Not tested on public roads
- Use at your own risk
- Always maintain control of vehicle

## Credits

- Paper (@paper5590) - Original Volvo implementation
- comma.ai - OpenPilot platform
- You - Testing and validation!

---
Generated: $(date)
Version: OpenPilot 0.10.0 + Polestar 2 Support