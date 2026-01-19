# How to Flash Your Fork onto Comma 3

## 🎯 Quick Overview

You have two options:
1. **Option A**: Fresh install (recommended if comma 3 is new or having issues)
2. **Option B**: Replace existing openpilot via SSH (faster if already running)

---

## Option A: Fresh Install (Recommended for C3)

### Step 1: Prepare Your Fork for Installation

Since Comma 3 is older hardware, we need to ensure compatibility. First, let's push your changes to GitHub:

```bash
# On your Mac
cd /Users/vibhusapra/projects/openpilot_polestar

# Push to your GitHub repo
git push origin master-cma

# Also push the opendbc submodule commit
cd opendbc_repo
git push origin master-cma
cd ..
```

### Step 2: Use Comma's Custom URL Installer

The comma 3 has a built-in installer that can clone from a custom URL.

**On your Comma 3 device**:

1. **Connect to WiFi** (via settings)

2. **Access the installer**:
   - Go to Settings → Device → Uninstall openpilot (if already installed)
   - OR on fresh device, wait for installer to appear

3. **Enter custom URL**:
   - When prompted for installer URL, enter:
   ```
   https://github.com/vibhusapra/openpilot_polestar.git
   ```

4. **Select branch**: When prompted, choose `master-cma`

5. **Wait for installation**: This will take 5-10 minutes

6. **Device will reboot** automatically when done

### Step 3: Verify Installation

After reboot, SSH into your comma 3:

```bash
# Find your comma 3's IP address from Settings → Network
ssh comma@<comma3_ip>
# Default password: comma3 has no password, just press Enter
```

Verify the installation:

```bash
cd /data/openpilot

# Check branch
git branch --show-current
# Expected: master-cma

# Check recent commit
git log -1 --oneline
# Expected: "Update CLAUDE.md and add Polestar 2 support"

# Check opendbc submodule
cd opendbc_repo
git status
# Expected: On branch master-cma

# Verify Polestar 2 model exists
python3 << 'EOF'
import sys
sys.path.insert(0, '/data/openpilot/opendbc_repo')
from opendbc.car.volvo.values import CAR
print("Available models:", [c.name for c in CAR])
EOF
# Expected: ['VOLVO_XC40_RECHARGE', 'POLESTAR_2']
```

---

## Option B: SSH Install (Faster if Already Running)

If your comma 3 already has openpilot running, you can replace it via SSH.

### Step 1: SSH into Comma 3

```bash
# From your Mac
ssh comma@<comma3_ip>
```

### Step 2: Backup Existing Installation (Optional but Recommended)

```bash
cd /data
sudo mv openpilot openpilot.backup.$(date +%Y%m%d)
```

### Step 3: Clone Your Fork

```bash
cd /data
git clone --branch master-cma https://github.com/vibhusapra/openpilot_polestar.git openpilot
cd openpilot
```

### Step 4: Initialize Submodules

```bash
# This is critical - must initialize opendbc submodule
git submodule update --init --recursive

# Verify opendbc is on correct branch
cd opendbc_repo
git checkout master-cma
git pull origin master-cma
cd ..
```

### Step 5: Set Permissions

```bash
sudo chown -R comma:comma /data/openpilot
```

### Step 6: Reboot

```bash
sudo reboot
```

### Step 7: Verify Installation

After reboot, SSH back in and verify (same commands as Option A Step 3).

---

## 🔧 Troubleshooting Installation

### Issue: "git clone fails" or "Permission denied"

**Solution 1**: Check if repo is public
```bash
# On your Mac
# Go to GitHub: https://github.com/vibhusapra/openpilot_polestar
# Settings → General → Make sure it's set to Public (not Private)
```

**Solution 2**: Use HTTPS instead of SSH
```bash
# If you used git@github.com URL, switch to https://
git clone https://github.com/vibhusapra/openpilot_polestar.git openpilot
```

### Issue: "Submodule not initialized"

```bash
# On comma 3
cd /data/openpilot
git submodule update --init --recursive

# If that fails, do it manually
cd /data/openpilot
rm -rf opendbc_repo
git clone https://github.com/paper5590/opendbc.git opendbc_repo
cd opendbc_repo
git checkout master-cma
```

### Issue: "Python import errors"

The comma 3 has Python 3.11, which should work. If you see import errors:

```bash
# On comma 3
cd /data/openpilot
python3 --version  # Should be 3.11+

# Try importing manually
python3 << 'EOF'
import sys
sys.path.insert(0, '/data/openpilot')
sys.path.insert(0, '/data/openpilot/opendbc_repo')
from opendbc.car.volvo.values import CAR
print("Success!")
EOF
```

### Issue: "Device keeps rebooting" or "Won't boot"

This usually means a critical error in the code. Check logs:

```bash
# On comma 3
tail -100 /data/community/crashes/*
# OR
journalctl -u comma --no-pager | tail -100
```

If stuck, revert to backup:
```bash
cd /data
sudo rm -rf openpilot
sudo mv openpilot.backup.* openpilot
sudo reboot
```

---

## 📱 Alternative: Use the Comma App

### Method: Custom Installer URL

Some versions of comma software allow setting a custom installer URL via the phone app.

1. **Open Comma Connect app** on your phone
2. **Connect to your comma 3**
3. **Go to Settings → Developer**
4. **Set Custom Installer URL**: `https://github.com/vibhusapra/openpilot_polestar.git`
5. **Set Branch**: `master-cma`
6. **Trigger reinstall** from device

*Note: This feature may not be available on all comma 3 firmware versions.*

---

## ✅ Post-Installation Checklist

After successful installation, verify everything is ready:

### 1. Check Software Version
```bash
ssh comma@<comma3_ip>
cd /data/openpilot
cat selfdrive/version.txt  # Note the version
git log -1 --oneline        # Verify your commit
```

### 2. Check UI Settings

On the comma 3 screen:
- Settings → Device → Experimental Mode should be available
- Settings → Toggles should show:
  - "Engage openpilot on double-tap cruise" (VolvoDoubleTapCruise)
  - "Pilot Assist engaged: Spoof hands on steering wheel" (VolvoSpoofPAHandsOnWheel)

### 3. Test Drive Recognition

```bash
# On comma 3
cd /data/openpilot
python3 << 'EOF'
import sys
sys.path.insert(0, '/data/openpilot')
from selfdrive.car.fingerprints import _get_interface

# This will attempt to detect car
# Should see Volvo XC40 or Polestar 2 as possible matches
EOF
```

### 4. Check Logs Are Working

```bash
# On comma 3
ls -lh /data/media/0/realdata/
# Should see route folders created during drives
```

---

## 🎯 Expected Results After Flash

### On Comma 3 Screen:

**Without car connected**:
- Should boot to normal openpilot UI
- Settings accessible
- No crashes or error screens

**With car connected (ignition on)**:
- Should show "openpilot Unavailable" (expected - waiting for fingerprint match)
- OR if fingerprint matches: "openpilot Ready"
- On first connection, may take 10-30 seconds to fingerprint

### What's Normal:

✅ "Waiting for controls start" - normal during first boot
✅ "Car unrecognized" on first connection - expected (placeholder fingerprint)
✅ Taking time to start - comma 3 is slower than comma 3X

### What's NOT Normal:

❌ Continuous rebooting
❌ Black screen / frozen UI
❌ Python tracebacks on screen
❌ "No data" errors that persist

---

## 🔄 Updating Your Fork Later

When you make changes and want to update the comma 3:

```bash
# On comma 3
cd /data/openpilot
git fetch origin
git checkout master-cma
git pull origin master-cma
git submodule update --init --recursive
sudo reboot
```

---

## 🆘 Emergency Recovery

If something goes wrong and you need to get back to stock openpilot:

### Method 1: Via SSH
```bash
# On comma 3
cd /data
sudo rm -rf openpilot
git clone https://github.com/commaai/openpilot.git
cd openpilot
git checkout v0.9.7  # Or latest release
git submodule update --init --recursive
sudo reboot
```

### Method 2: Factory Reset
1. Settings → Device → Reset to Factory Settings
2. This will reinstall stock openpilot
3. Then you can try installing your fork again

---

## 📞 Need Help?

If you encounter issues during installation:

1. **Check comma 3 logs**:
```bash
ssh comma@<comma3_ip>
tail -100 /data/community/crashes/*
journalctl -u comma --no-pager | tail -100
```

2. **Verify GitHub access**:
```bash
# On comma 3
git ls-remote https://github.com/vibhusapra/openpilot_polestar.git
# Should list branches
```

3. **Check disk space**:
```bash
df -h /data
# Should have several GB free
```

4. **Ask for help**: Share error logs and describe what step failed

---

## 🎉 Next Steps After Successful Flash

Once your fork is installed on comma 3:

1. ✅ **Verify installation** (checklist above)
2. ✅ **Install harness in car**
3. ✅ **Follow POLESTAR2_SETUP.md** testing phases
4. ✅ **Start with passive monitoring** (no steering)
5. ✅ **Test steering in safe environment**

Good luck! 🚀
