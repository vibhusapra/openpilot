# Polestar 2 openpilot Setup Guide

## ✅ What's Ready

Your openpilot fork is now configured for Polestar 2 testing:

- ✅ **Volvo CMA implementation**: Complete steering control via LCA messages
- ✅ **PR#1 angle cap fix**: Asymmetric left/right angle encoding (fixes -4.25° issue)
- ✅ **Polestar 2 model added**: `CAR.POLESTAR_2` with CMA platform specs
- ✅ **Placeholder fingerprint**: Using XC40 baseline (will update after first drive)
- ✅ **UI toggles**: VolvoDoubleTapCruise & VolvoSpoofPAHandsOnWheel enabled
- ✅ **Torque auto-learning**: Volvo enabled in learning system
- ✅ **Safety model**: Volvo safety implemented in opendbc
- ✅ **Analysis tools**: LCA message analyzer and CAN comparison tools included

## 📋 Pre-Flight Checklist

### Phase 1: Hardware Assembly (You Are Here)

#### Harness Wiring
- [ ] **Assemble harness** following your wiring diagram:
  - MID-1 CAN #20 (VCU) → Harness Box CAN0_H
  - MID-1 CAN #19 (VCU) → Harness Box CAN0_L
  - MID-1 CAN #20 (PSCM) → Harness Box CAN2_H
  - MID-1 CAN #19 (PSCM) → Harness Box CAN2_L
  - IGN #25 → Y-split (car + harness box)
  - 12V #26 → Harness box
  - GND #32 → Harness box

#### Pre-Installation Checks
- [ ] **Multimeter verification**:
  - CAN0_H ↔ CAN0_L ≈ 60Ω (car off)
  - CAN2_H ↔ CAN2_L ≈ 60Ω (car off)
  - GND continuity between all grounds
  - IGN = 0V (car off) / 12V (car on)
  - No shorts between power rails

- [ ] **Visual inspection**:
  - All CAN twisted pairs intact
  - No exposed copper
  - WAGO connectors fully closed
  - Proper strain relief with Tesa tape
  - Labels clear and readable

### Phase 2: Software Setup

#### Flash Comma 3 Device
```bash
# On your Mac (prepare the software)
cd /Users/vibhusapra/projects/openpilot_polestar
git status  # Verify you're on master-cma branch
git log -1  # Should show "Update CLAUDE.md and add Polestar 2 support"

# Option A: Flash via SSH (if comma 3 is already running openpilot)
ssh comma@<comma3_ip>
cd /data
mv openpilot openpilot.bak  # Backup existing
git clone https://github.com/vibhusapra/openpilot_polestar openpilot
cd openpilot
git checkout master-cma
git submodule update --init --recursive
reboot

# Option B: Fresh install via comma.ai installer
# 1. Use comma installer to flash stock openpilot
# 2. SSH in and replace with your fork (see Option A)
```

#### Verify Installation
```bash
# SSH into comma 3
ssh comma@<comma3_ip>

# Check branch
cd /data/openpilot
git branch --show-current  # Should show: master-cma
git log -1 --oneline       # Should show your Polestar 2 commit

# Check opendbc submodule
cd opendbc_repo
git branch --show-current  # Should show: master-cma
ls opendbc/car/volvo/      # Should show all Volvo files

# Verify Polestar 2 model exists
python3 -c "from opendbc.car.volvo.values import CAR; print([c.name for c in CAR])"
# Should print: ['VOLVO_XC40_RECHARGE', 'POLESTAR_2']
```

### Phase 3: First Installation (Bench Test)

#### Install Harness (Car Off)
- [ ] Locate VCU1 connector near driver footwell
- [ ] Install Volvo/Polestar intercept harness inline
- [ ] Connect developer harness to harness box
- [ ] Connect harness box to Comma 3 via Ethernet
- [ ] Secure all connections
- [ ] **Do not start car yet**

#### Power-On Test (Car Off → On)
- [ ] Turn car ON (ready to drive)
- [ ] **Expected**: No dash errors
- [ ] **Expected**: Comma 3 boots normally
- [ ] **Expected**: Stock Pilot Assist still functions
- [ ] Turn car OFF

### Phase 4: Passive Monitoring (No Control)

#### First Drive - Logging Only
```bash
# SSH into comma 3 and disable lateral control temporarily
ssh comma@<comma3_ip>
cd /data/openpilot
# Edit selfdrive/controls/controlsd.py to force latActive=False for safety
# OR just don't enable openpilot during this drive
```

- [ ] **Drive conditions**: Safe area, low traffic
- [ ] **Duration**: 10-15 minutes
- [ ] **Actions**:
  - Enable stock Pilot Assist
  - Drive normally
  - Test all speeds (0-highway)
  - Test left and right turns
  - Monitor comma 3 screen for errors

- [ ] **After drive - Check logs**:
```bash
# On comma 3
cd /data/media/0/realdata/<latest_route>
# Verify rlog files exist
# Download for analysis: scp comma@<ip>:/data/media/0/realdata/<route>/* ./
```

#### Analyze First Route
```bash
# On your Mac
cd /Users/vibhusapra/projects/openpilot_polestar
python route_analysis/analyze_lca_messages.py <route_id> --output polestar2_baseline.csv

# Review CSV for:
# - LCA messages present on CAN0
# - PSCM messages present on CAN2
# - Steering angle signals look reasonable
# - No unexpected errors
```

### Phase 5: Active Steering (CRITICAL - BE CAREFUL)

#### Enable openpilot Steering
- [ ] **Prerequisites**:
  - Passive monitoring successful
  - No dash errors
  - Logs look clean
  - Comfortable with harness installation

- [ ] **Test location**: Empty parking lot or private road
- [ ] **Test procedure**:
  1. Start car
  2. Enable openpilot (or double-tap cruise if VolvoDoubleTapCruise enabled)
  3. **Keep hands on wheel**
  4. Test at low speed first (10-20 mph)
  5. Verify steering response
  6. Test left turn
  7. Test right turn
  8. Test driver override (gently resist steering)
  9. **If anything feels wrong - DISENGAGE IMMEDIATELY**

- [ ] **Expected behavior**:
  - Smooth steering
  - No jerking or sudden movements
  - Easy to override
  - Symmetrical left/right response
  - No dash errors

#### Highway Test (if parking lot successful)
- [ ] **Prerequisites**:
  - Low-speed test successful
  - No concerning behavior
  - Confident in system

- [ ] **Test conditions**:
  - Light traffic
  - Good weather
  - Straight highway section first
  - **Hands on wheel at all times**

- [ ] **Monitor for**:
  - Steering cap issues (especially right turns)
  - Asymmetric override force
  - Checksum errors on dash
  - Speed signal accuracy

### Phase 6: Tuning & Refinement

#### Known Issues to Address (from Paper's notes)

1. **Asymmetric driver override** (harder to override right turns)
   - **Location**: `opendbc_repo/opendbc/car/volvo/carcontroller.py`
   - **Likely fix**: Change torque variable from 255 → 128
   - **Test**: Compare override force left vs right

2. **City driving** (steering too slow)
   - **Investigation**: Find "high-speed steering mode" flag
   - **Tool**: Use `analyze_lca_messages.py` to compare highway vs city routes
   - **Look for**: LCA message field differences at different speeds

3. **Speed deviation** (~0.5 km/h difference)
   - **Location**: Signal scaling in `carstate.py`
   - **Fix**: Adjust conversion factor

4. **Deadband around 0°** (already removed in latest commit)
   - **Verify**: Test straight-line stability
   - **Check**: No oscillation near center

5. **Rare checksum failures**
   - **Monitor**: Dash for DTC codes
   - **Debug**: Use `can_print_changes_2.py` if persistent

#### Capture Real Polestar 2 Fingerprint
```bash
# After successful drive
cd /data/openpilot
python selfdrive/debug/print_docs_diff.py

# This will output CAN message IDs detected
# Copy the fingerprint and update opendbc_repo/opendbc/car/volvo/fingerprints.py
# Replace the placeholder with real Polestar 2 fingerprint
```

#### Enable UI Toggles (Optional)
```bash
# On comma 3 UI
# Settings → Device → Experimental Mode (if needed)
# Find: "Engage openpilot on double-tap cruise" (VolvoDoubleTapCruise)
# Find: "Pilot Assist engaged: Spoof hands on steering wheel" (VolvoSpoofPAHandsOnWheel)
# Test each individually
```

### Phase 7: Long-Term Use

#### Torque Learning
- [ ] **Drive 50+ km** at highway speeds for learning to converge
- [ ] **Check parameters**:
```bash
# On comma 3
cat /data/params/d/LiveTorqueParameters
# Should show learned values for latAccelFactor, frictionCoefficient
```

#### Route Analysis for Improvements
```bash
# Compare multiple routes
python route_analysis/analyze_lca_messages.py route1 --output route1.csv
python route_analysis/analyze_lca_messages.py route2 --output route2.csv

# Use can_print_changes_2.py to find state-dependent changes
python selfdrive/debug/can_print_changes_2.py route1:10-20 route2:10-20
```

## 🚨 Safety Reminders

**CRITICAL - READ BEFORE TESTING**:

1. **Hands on wheel at ALL times** during initial testing
2. **Be ready to disengage** - know how to turn off openpilot instantly
3. **Start in safe environments** - empty parking lots, then quiet roads
4. **Test incrementally** - don't jump straight to highway
5. **Monitor for errors** - any dash warnings = stop and investigate
6. **Driver is always responsible** - openpilot is SAE Level 2 (driver assist)

**If you experience**:
- Unexpected steering behavior → DISENGAGE IMMEDIATELY
- Dash errors → Stop and check logs
- Unable to override → PULL OVER SAFELY and power off
- Any safety concerns → Stop testing, review setup

## 📊 Success Criteria

**Minimum Viable (before regular use)**:
- ✅ No dash errors during operation
- ✅ Smooth steering on highway
- ✅ Easy driver override (< 10 lb force)
- ✅ Symmetrical left/right behavior
- ✅ Stable for 30+ minute continuous drive
- ✅ No checksum errors in logs

**Full Feature Parity (nice to have)**:
- ✅ Works in city (tight turns)
- ✅ VolvoDoubleTapCruise functional
- ✅ VolvoSpoofPAHandsOnWheel functional
- ✅ Torque learning converged
- ✅ Speed signal accurate
- ✅ All Paper's known issues resolved

## 🛠️ Troubleshooting

### Issue: No steering control
**Check**:
1. Panda safety mode (should be `volvo`, not `SAFETY_SILENT`)
2. CAN traffic on both CAN0 and CAN2 (use logs)
3. LCA messages (0x58, 0x57) being sent
4. PSCM acknowledgment in responses

**Debug**:
```bash
# SSH to comma 3
tail -f /data/media/0/realdata/<latest>/rlog.bz2 | bunzip2 | grep -i lca
```

### Issue: Steering cap / asymmetry
**Check**:
1. LCA message byte values during max turns
2. Compare with Paper's analysis
3. Verify angle encoder scaling (0.05596 deg/count)

**Debug**:
```bash
python route_analysis/analyze_lca_messages.py <route> | grep "LCA.*255\|LCA.*128"
```

### Issue: Dash errors (DTC codes)
**Check**:
1. Checksum calculations in CAN messages
2. Counter increments correct
3. Message frequencies match stock

**Debug**:
```bash
python selfdrive/debug/can_print_changes_2.py stock_route your_route
# Look for unexpected bit changes
```

## 📞 Getting Help

1. **Discord**: Check #volvo-cma channel, tag @Paper or @evi1gasm
2. **Logs**: Always share route IDs when asking for help
3. **Comparison**: Use analysis tools to compare with known-good XC40 routes
4. **This repo**: Open GitHub issues with detailed logs and descriptions

## 🎯 Current Status

**Branch**: `master-cma`
**Last commit**: Update CLAUDE.md and add Polestar 2 support
**opendbc commit**: Add Polestar 2 support to Volvo CMA platform
**Ready for**: Hardware assembly → First test drive

## 📝 Next Steps (In Order)

1. ✅ **Complete harness assembly** (follow wiring diagram)
2. ✅ **Bench test harness** (multimeter checks)
3. ✅ **Flash Comma 3** (master-cma branch)
4. ✅ **Install in car** (car off, verify no errors)
5. ✅ **Passive monitoring** (log-only drive)
6. ✅ **Analyze first route** (verify CAN messages)
7. ✅ **Low-speed steering test** (parking lot)
8. ✅ **Highway test** (if low-speed successful)
9. ✅ **Capture Polestar 2 fingerprint** (update code)
10. ✅ **Address tuning issues** (override force, deadband, etc.)
11. ✅ **Long-term testing** (torque learning, refinement)

---

**Good luck! You're building on a solid foundation - Paper's done the hard work. Stay safe! 🚗💨**
