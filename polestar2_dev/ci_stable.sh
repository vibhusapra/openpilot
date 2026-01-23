#!/bin/bash

# Stable CI Test - No fancy progress, just works

set +e  # Don't exit on error
set -o pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}  OpenPilot CI Test - Stable Version${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo

# 1. Check environment
echo -e "${BLUE}[1/7] Environment${NC}"
echo "Python: $(python3 --version)"
echo "Scons: $(scons --version 2>&1 | head -1)"
echo "Branch: $(git branch --show-current)"
echo

# 2. Check critical fixes
echo -e "${BLUE}[2/7] Critical Fixes${NC}"

if grep -q "PANDA_BUS_CNT" "selfdrive/pandad/panda.h"; then
    echo -e "${GREEN}✓${NC} PANDA_BUS_CNT defined"
else
    echo -e "${RED}✗${NC} PANDA_BUS_CNT missing"
fi

if grep -q 'extern "C"' "selfdrive/pandad/panda.h"; then
    echo -e "${GREEN}✓${NC} extern \"C\" wrapper present"
else
    echo -e "${RED}✗${NC} extern \"C\" wrapper missing"
fi

if grep -q "DLC_TO_LEN_DEFINED" "panda/board/can.h"; then
    echo -e "${GREEN}✓${NC} dlc_to_len guard present"
else
    echo -e "${YELLOW}⚠${NC} dlc_to_len guard missing"
fi
echo

# 3. Check Polestar 2 files
echo -e "${BLUE}[3/7] Polestar 2 Integration${NC}"

if [ -f "opendbc_repo/opendbc/car/volvo/interface.py" ]; then
    echo -e "${GREEN}✓${NC} Volvo interface.py"
else
    echo -e "${RED}✗${NC} Volvo interface.py missing"
fi

if grep -q "POLESTAR_2" "opendbc_repo/opendbc/car/volvo/values.py" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} POLESTAR_2 in values.py"
else
    echo -e "${RED}✗${NC} POLESTAR_2 not found"
fi

if [ -f "opendbc_repo/opendbc/safety/modes/volvo.h" ]; then
    SIZE=$(stat -f%z "opendbc_repo/opendbc/safety/modes/volvo.h" 2>/dev/null || stat -c%s "opendbc_repo/opendbc/safety/modes/volvo.h" 2>/dev/null || echo "0")
    echo -e "${GREEN}✓${NC} Volvo safety firmware (${SIZE} bytes)"
else
    echo -e "${RED}✗${NC} Volvo safety firmware missing"
fi
echo

# 4. Clean build
echo -e "${BLUE}[4/7] Clean Previous Build${NC}"
scons -c -s > /dev/null 2>&1
rm -rf .sconsign.dblite cereal/gen 2>/dev/null
echo -e "${GREEN}✓${NC} Cleaned"
echo

# 5. Build - SIMPLE VERSION
echo -e "${BLUE}[5/7] Building OpenPilot${NC}"
echo "This will take 2-5 minutes..."
echo "Running: scons -u -j$(sysctl -n hw.ncpu || echo 4)"
echo

# Just run scons directly - no pipes, no progress bars
BUILD_START=$(date +%s)

if scons -u -j$(sysctl -n hw.ncpu || echo 4); then
    BUILD_RESULT=0
    echo -e "\n${GREEN}✓${NC} Build completed successfully"
else
    BUILD_RESULT=$?
    echo -e "\n${RED}✗${NC} Build failed with exit code $BUILD_RESULT"
fi

BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))
echo "Build time: ${BUILD_TIME} seconds"
echo

# 6. Check artifacts
echo -e "${BLUE}[6/7] Build Artifacts${NC}"

if [ -f "cereal/libcereal.a" ]; then
    echo -e "${GREEN}✓${NC} cereal/libcereal.a ($(ls -lh cereal/libcereal.a | awk '{print $5}'))"
else
    echo -e "${RED}✗${NC} cereal/libcereal.a missing"
fi

if [ -f "msgq_repo/libmsgq.a" ]; then
    echo -e "${GREEN}✓${NC} msgq_repo/libmsgq.a ($(ls -lh msgq_repo/libmsgq.a | awk '{print $5}'))"
else
    echo -e "${RED}✗${NC} msgq_repo/libmsgq.a missing"
fi

if [ -f "common/libcommon.a" ]; then
    echo -e "${GREEN}✓${NC} common/libcommon.a ($(ls -lh common/libcommon.a | awk '{print $5}'))"
else
    echo -e "${RED}✗${NC} common/libcommon.a missing"
fi
echo

# 7. Python imports
echo -e "${BLUE}[7/7] Python Import Test${NC}"

python3 << 'EOF'
import sys
import os
os.chdir('.')
sys.path.insert(0, '.')
sys.path.insert(0, './opendbc_repo')

try:
    from opendbc.car.volvo.values import CAR
    print("✓ Volvo values imported")

    if hasattr(CAR, 'POLESTAR_2'):
        print(f"✓ POLESTAR_2 exists: {CAR.POLESTAR_2}")
    else:
        print("✗ POLESTAR_2 not found")

    from opendbc.car.volvo.fingerprints import FINGERPRINTS
    if CAR.POLESTAR_2 in FINGERPRINTS:
        fp = FINGERPRINTS[CAR.POLESTAR_2][0]
        print(f"✓ Fingerprint: {len(fp)} CAN messages")
    else:
        print("✗ Fingerprint missing")
except Exception as e:
    print(f"✗ Import failed: {e}")
EOF

echo

# Summary
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
if [ $BUILD_RESULT -eq 0 ]; then
    echo -e "${GREEN}${BOLD}✓ BUILD PASSED${NC}"
    echo "OpenPilot 0.10.0 with Polestar 2 support is ready!"
else
    echo -e "${RED}${BOLD}✗ BUILD FAILED${NC}"
    echo "Fix errors before deployment"
fi
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

exit $BUILD_RESULT