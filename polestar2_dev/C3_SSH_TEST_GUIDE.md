# Comma 3 SSH Testing Guide

## 🏠 Testing OpenPilot Inside (Without Car)

This guide walks through testing your Polestar 2 OpenPilot build on the Comma 3 device inside your house, without needing to be in the car.

## 📋 Prerequisites

1. **Comma 3 device** powered on (USB-C power)
2. **SSH access** to the Comma 3
3. **WiFi connection** for both your computer and Comma 3
4. **GitHub fork** updated with our changes

## 🔌 Step 1: Connect to Comma 3

```bash
# Find your Comma 3's IP address (check your router or use the Comma app)
# Default SSH password is 'comma'

ssh comma@<COMMA_IP>
# Example: ssh comma@192.168.1.100
```

## 🚀 Step 2: Deploy OpenPilot

**Option A: Quick Deploy (if you trust the scripts)**
```bash
# Download and run deployment script
cd /data
curl -L https://raw.githubusercontent.com/vibhusapra/openpilot/polestar2-steering-v0.10.0/polestar2_dev/c3_deploy.sh -o deploy.sh
chmod +x deploy.sh
./deploy.sh
```

**Option B: Manual Deploy (more control)**
```bash
# 1. Stop any running OpenPilot
sudo systemctl stop openpilot
tmux kill-server 2>/dev/null || true

# 2. Backup existing
[ -d /data/openpilot ] && mv /data/openpilot /data/openpilot_backup_$(date +%Y%m%d)

# 3. Clone your fork
cd /data
git clone --recursive -b polestar2-steering-v0.10.0 https://github.com/vibhusapra/openpilot.git

# 4. Update opendbc to Paper's fork
cd /data/openpilot/opendbc_repo
git remote set-url origin https://github.com/paper5590/opendbc.git
git fetch origin
git checkout master-cma
git pull origin master-cma

# 5. Build (THIS IS THE CRITICAL PART)
cd /data/openpilot

# Build panda firmware first
cd panda
scons -j4
cd ..

# Build OpenPilot
scons -u -j4
```

## 🧪 Step 3: Run Tests

```bash
cd /data/openpilot

# Make scripts executable
chmod +x polestar2_dev/*.sh

# Run comprehensive test suite
./polestar2_dev/c3_test.sh
```

### Expected Test Results:

✅ **Should PASS:**
- Directory structure checks
- Built libraries (cereal, msgq, common)
- Python modules (.so files)
- Panda firmware files
- Critical fixes (PANDA_BUS_CNT, etc.)
- Volvo/Polestar 2 integration

❌ **May FAIL (OK without car):**
- Some process starts (need hardware)
- CAN communication tests
- Camera tests

## 🔍 Step 4: Process Testing

```bash
# Test if processes can start
./polestar2_dev/c3_process_test.sh
```

This tests:
1. Manager.py startup
2. Pandad initialization
3. Car detection modules
4. Launch script

**Note:** Some errors are expected without hardware connected.

## 🐛 Step 5: Debugging

### Check Build Logs
```bash
# If build failed
cat /tmp/scons_log_* | grep -i error

# Check what was built
ls -la cereal/*.a msgq_repo/*.a common/*.a
find . -name "*.so" -type f | grep -E "(msgq|common)"
```

### Test Python Imports Manually
```bash
cd /data/openpilot
python3
```

```python
import sys
sys.path.insert(0, '.')
sys.path.insert(0, './opendbc_repo')

# Test critical imports
from opendbc.car.volvo.interface import CarInterface
from opendbc.car.volvo.values import CAR
print(f"POLESTAR_2: {CAR.POLESTAR_2}")

# Test manager
from selfdrive.manager.manager import manager_init
print("Manager imports OK")
```

### Monitor System Logs
```bash
# Watch for errors
journalctl -f

# In another SSH session, try to start
cd /data/openpilot
./launch_openpilot.sh

# Check dmesg for kernel issues
dmesg | tail -50
```

## ✅ Step 6: Verify Success

**Build is successful if:**
1. `c3_test.sh` shows most tests passing
2. No Python import errors
3. Manager.py starts without crashing
4. Core libraries exist (.a files)
5. Python modules built (.so files)

## 🚗 Step 7: Ready for Car Testing

If all tests pass:

1. **Set environment variables:**
```bash
cd /data/openpilot
cat > env.sh << 'EOF'
export FINGERPRINT="POLESTAR_2"
export SKIP_FW_QUERY=1
export BLOCK_INTERNET=0
export COMPLETED_TRAINING=1
EOF
```

2. **Shutdown C3:**
```bash
sudo shutdown now
```

3. **Install in car** and test with engine on (not driving)

## ⚠️ Common Issues & Fixes

### Build Fails
```bash
# Missing dependency
sudo apt update
sudo apt install -y python3-dev

# Submodule issues
git submodule update --init --recursive

# Clean and rebuild
scons -c
scons -u -j4
```

### Python Import Errors
```bash
# Check Python path
python3 -c "import sys; print(sys.path)"

# Rebuild Python modules
scons --clean
scons -u -j4
```

### Manager Won't Start
```bash
# Check logs
cat /data/openpilot/selfdrive/manager/manager.log

# Try with debug
PYTHONPATH=/data/openpilot python3 selfdrive/manager/manager.py
```

## 📝 What to Report Back

After testing, note:

1. **Build Status:** Did it complete? Any errors?
2. **Test Results:** Which tests passed/failed?
3. **Import Tests:** Can Python import Volvo modules?
4. **Process Tests:** Does manager.py start?
5. **Error Messages:** Any specific errors or crashes?

## 🔄 Updating Code

If we need to make changes:

```bash
cd /data/openpilot
git fetch origin
git pull origin polestar2-steering-v0.10.0

# Rebuild
scons -u -j4

# Re-test
./polestar2_dev/c3_test.sh
```

## 💡 Pro Tips

1. **Use tmux** for multiple sessions:
```bash
tmux new -s build    # Building
tmux new -s test     # Testing
tmux new -s monitor  # Monitoring logs
```

2. **Save logs** for debugging:
```bash
./polestar2_dev/c3_test.sh 2>&1 | tee test_results.log
```

3. **Check resource usage**:
```bash
htop  # CPU/Memory
df -h # Disk space
```

---

## Ready? Let's Test! 🚀

Start with Step 1 and work through systematically. Report any issues and we'll debug together via SSH.