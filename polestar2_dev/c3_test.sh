#!/bin/bash

# Comma 3 Test Script - Run after deployment
# Tests if everything built correctly

set +e  # Don't exit on error - we want to see all test results

echo "=========================================="
echo "  Comma 3 OpenPilot Test Suite"
echo "  Polestar 2 Integration"
echo "=========================================="
echo

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

TESTS_PASSED=0
TESTS_FAILED=0

# Test function
run_test() {
    local test_name="$1"
    local test_cmd="$2"

    echo -e "${BLUE}Testing:${NC} $test_name"

    if eval "$test_cmd"; then
        echo -e "${GREEN}✓${NC} $test_name passed"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}✗${NC} $test_name failed"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    echo
}

# 1. Check directory structure
echo -e "${BLUE}[1/8] Directory Structure${NC}"
echo "----------------------------------------"
run_test "OpenPilot directory exists" "[ -d /data/openpilot ]"
run_test "Polestar2 dev tools exist" "[ -d /data/openpilot/polestar2_dev ]"
run_test "opendbc_repo exists" "[ -d /data/openpilot/opendbc_repo ]"
run_test "panda directory exists" "[ -d /data/openpilot/panda ]"

# 2. Check built libraries
echo -e "${BLUE}[2/8] Built Libraries${NC}"
echo "----------------------------------------"
run_test "cereal library built" "[ -f /data/openpilot/cereal/libcereal.a ]"
run_test "msgq library built" "[ -f /data/openpilot/msgq_repo/libmsgq.a ]"
run_test "common library built" "[ -f /data/openpilot/common/libcommon.a ]"

# 3. Check Python modules
echo -e "${BLUE}[3/8] Python Extension Modules${NC}"
echo "----------------------------------------"
run_test "msgq ipc module" "[ -f /data/openpilot/msgq_repo/msgq/ipc_pyx.so ]"
run_test "visionipc module" "[ -f /data/openpilot/msgq_repo/msgq/visionipc/visionipc_pyx.so ]"
run_test "params module" "[ -f /data/openpilot/common/params_pyx.so ]"
run_test "transformations module" "[ -f /data/openpilot/common/transformations/transformations.so ]"

# 4. Check Panda firmware
echo -e "${BLUE}[4/8] Panda Firmware${NC}"
echo "----------------------------------------"
run_test "Panda H7 recovery" "[ -f /data/openpilot/panda/board/obj/panda_h7_recovery.bin ]"
run_test "Panda H7 main" "[ -f /data/openpilot/panda/board/obj/panda_h7.bin ]"

# 5. Check critical fixes
echo -e "${BLUE}[5/8] Critical Fixes${NC}"
echo "----------------------------------------"
run_test "PANDA_BUS_CNT defined" "grep -q 'PANDA_BUS_CNT' /data/openpilot/selfdrive/pandad/panda.h"
run_test "extern C wrapper present" "grep -q 'extern \"C\"' /data/openpilot/selfdrive/pandad/panda.h"
run_test "DBC generation disabled" "grep -q '^# generated = env.Command' /data/openpilot/opendbc_repo/opendbc/dbc/SConscript"

# 6. Check Polestar 2 integration
echo -e "${BLUE}[6/8] Polestar 2 Integration${NC}"
echo "----------------------------------------"
run_test "Volvo interface.py exists" "[ -f /data/openpilot/opendbc_repo/opendbc/car/volvo/interface.py ]"
run_test "POLESTAR_2 in values.py" "grep -q 'POLESTAR_2' /data/openpilot/opendbc_repo/opendbc/car/volvo/values.py"
run_test "Volvo safety firmware" "[ -f /data/openpilot/opendbc_repo/opendbc/safety/modes/volvo.h ]"
run_test "Volvo allowed in torqued" "grep -q 'volvo' /data/openpilot/selfdrive/locationd/torqued.py"

# 7. Python imports test
echo -e "${BLUE}[7/8] Python Import Tests${NC}"
echo "----------------------------------------"

cd /data/openpilot

# Test Volvo module imports
python3 << 'EOF'
import sys
sys.path.insert(0, '.')
sys.path.insert(0, './opendbc_repo')

try:
    from opendbc.car.volvo.interface import CarInterface
    print("✓ Volvo interface imports")
    exit(0)
except ImportError as e:
    print(f"✗ Volvo interface import failed: {e}")
    exit(1)
EOF

if [ $? -eq 0 ]; then
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# Test POLESTAR_2 exists
python3 << 'EOF'
import sys
sys.path.insert(0, '.')
sys.path.insert(0, './opendbc_repo')

try:
    from opendbc.car.volvo.values import CAR
    if hasattr(CAR, 'POLESTAR_2'):
        print(f"✓ POLESTAR_2 exists: {CAR.POLESTAR_2}")
        exit(0)
    else:
        print("✗ POLESTAR_2 not in CAR enum")
        exit(1)
except Exception as e:
    print(f"✗ Error: {e}")
    exit(1)
EOF

if [ $? -eq 0 ]; then
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# 8. Check if manager can import
echo -e "${BLUE}[8/8] Manager Import Test${NC}"
echo "----------------------------------------"

python3 << 'EOF'
import sys
sys.path.insert(0, '.')

try:
    from selfdrive.manager.manager import manager_init
    print("✓ Manager can be imported")
    exit(0)
except ImportError as e:
    print(f"✗ Manager import failed: {e}")
    exit(1)
EOF

if [ $? -eq 0 ]; then
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

# Summary
echo
echo "=========================================="
echo "  TEST SUMMARY"
echo "=========================================="
echo -e "${GREEN}Passed:${NC} $TESTS_PASSED"
echo -e "${RED}Failed:${NC} $TESTS_FAILED"

if [ $TESTS_FAILED -eq 0 ]; then
    echo
    echo -e "${GREEN}✅ ALL TESTS PASSED!${NC}"
    echo
    echo "OpenPilot should be ready to run."
    echo "Try starting it with:"
    echo "  cd /data/openpilot && ./launch_openpilot.sh"
else
    echo
    echo -e "${RED}⚠️ SOME TESTS FAILED${NC}"
    echo
    echo "Debug the failures above before trying to run OpenPilot."
    echo "Common issues:"
    echo "  - Missing .so files: Build may have failed"
    echo "  - Import errors: Python path or dependency issues"
    echo "  - Missing panda firmware: Panda build may have failed"
fi

echo "=========================================="

exit $TESTS_FAILED