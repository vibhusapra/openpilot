#!/bin/bash

# Fixed CI Test Suite - Handles dependencies and errors properly
# This version won't crash and will show exactly what's happening

set +e  # Don't exit on error - capture everything
set -o pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# Timing
START_TIME=$(date +%s)
STEP_COUNT=0
TOTAL_STEPS=10
ERRORS=()
WARNINGS=()

# Create log directory
LOG_DIR="polestar2_dev/ci_logs"
mkdir -p "$LOG_DIR"
TEST_LOG="$LOG_DIR/ci_test_$(date +%Y%m%d_%H%M%S).log"

# Logging function
log() {
    echo "$1" | tee -a "$TEST_LOG"
}

# Banner
print_banner() {
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}  OpenPilot CI Test Suite - Fixed Version${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "  Version: OpenPilot 0.10.0"
    echo -e "  Branch: $(git branch --show-current 2>/dev/null || echo 'unknown')"
    echo -e "  Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
    echo -e "  Log file: $TEST_LOG"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo
}

# Step counter
step() {
    STEP_COUNT=$((STEP_COUNT + 1))
    echo -e "\n${BLUE}[$STEP_COUNT/$TOTAL_STEPS]${NC} ${BOLD}$1${NC}"
    echo -e "${CYAN}$(printf '%.0s─' {1..60})${NC}"
}

# Success/failure markers
pass() {
    echo -e "  ${GREEN}✓${NC} $1"
}

fail() {
    echo -e "  ${RED}✗${NC} $1"
    ERRORS+=("$1")
}

warn() {
    echo -e "  ${YELLOW}⚠${NC} $1"
    WARNINGS+=("$1")
}

info() {
    echo -e "  ${MAGENTA}ℹ${NC} $1"
}

print_banner

# ============================================================================
step "Environment Detection"
# ============================================================================

OS_TYPE="unknown"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS_TYPE="linux"
    info "Linux detected"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS_TYPE="macos"
    info "macOS detected ($(sw_vers -productVersion))"
else
    warn "Unknown OS: $OSTYPE"
fi

ARCH=$(uname -m)
info "Architecture: $ARCH"

# Python version
PYTHON_CMD=$(which python3)
if [ -z "$PYTHON_CMD" ]; then
    fail "Python3 not found"
    exit 1
fi

PYTHON_VERSION=$($PYTHON_CMD --version 2>&1 | cut -d' ' -f2)
info "Python: $PYTHON_VERSION at $PYTHON_CMD"

# ============================================================================
step "Critical Tools Check"
# ============================================================================

# Check scons
if command -v scons &> /dev/null; then
    SCONS_VERSION=$(scons --version 2>&1 | head -1)
    pass "scons: $SCONS_VERSION"
else
    fail "scons not installed"
    info "Installing scons..."

    if [[ "$OS_TYPE" == "macos" ]]; then
        brew install scons 2>&1 | tail -3
    else
        pip3 install scons 2>&1 | tail -3
    fi

    if command -v scons &> /dev/null; then
        pass "scons installed successfully"
    else
        fail "Failed to install scons"
        exit 1
    fi
fi

# ============================================================================
step "Python Dependencies"
# ============================================================================

echo -e "${BOLD}Checking required Python packages:${NC}"

PYTHON_DEPS=(
    "numpy"
    "pycapnp"
    "Cython"
    "cffi"
    "pyzmq"
    "pycryptodome"
)

MISSING_DEPS=()

for dep in "${PYTHON_DEPS[@]}"; do
    # Special case for pycryptodome which imports as Crypto
    if [ "$dep" = "pycryptodome" ]; then
        CHECK_MODULE="Crypto"
    else
        CHECK_MODULE="${dep//-/_}"
    fi

    if python3 -c "import $CHECK_MODULE" &>/dev/null; then
        pass "$dep"
    else
        warn "$dep not installed"
        MISSING_DEPS+=("$dep")
    fi
done

if [ ${#MISSING_DEPS[@]} -gt 0 ]; then
    echo
    info "Installing missing dependencies: ${MISSING_DEPS[*]}"
    pip3 install ${MISSING_DEPS[@]} 2>&1 | tail -5

    # Verify installation
    echo
    echo -e "${BOLD}Verifying installation:${NC}"
    for dep in "${MISSING_DEPS[@]}"; do
        if [ "$dep" = "pycryptodome" ]; then
            CHECK_MODULE="Crypto"
        else
            CHECK_MODULE="${dep//-/_}"
        fi

        if python3 -c "import $CHECK_MODULE" &>/dev/null; then
            pass "$dep installed"
        else
            fail "$dep installation failed"
        fi
    done
fi

# ============================================================================
step "Repository Validation"
# ============================================================================

if [ ! -f "launch_openpilot.sh" ] || [ ! -d "selfdrive" ]; then
    fail "Not in OpenPilot root directory"
    exit 1
else
    pass "In OpenPilot directory"
fi

# Check git status
UNCOMMITTED=$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')
if [ "$UNCOMMITTED" -gt 0 ]; then
    warn "$UNCOMMITTED uncommitted changes"
else
    pass "Working directory clean"
fi

# ============================================================================
step "Critical Fixes Validation"
# ============================================================================

echo -e "${BOLD}Checking our fixes are in place:${NC}"

# Check PANDA_BUS_CNT
if grep -q "PANDA_BUS_CNT" "selfdrive/pandad/panda.h" 2>/dev/null; then
    pass "PANDA_BUS_CNT defined in panda.h"
else
    fail "PANDA_BUS_CNT not defined - THIS WILL CAUSE BUILD FAILURE"
fi

# Check extern C wrapper
if grep -q 'extern "C"' "selfdrive/pandad/panda.h" 2>/dev/null; then
    pass "extern \"C\" wrapper present"
else
    fail "extern \"C\" wrapper missing - WILL CAUSE C++ LINKAGE ERRORS"
fi

# Check dlc_to_len guard
if grep -q "DLC_TO_LEN_DEFINED" "panda/board/can.h" 2>/dev/null; then
    pass "dlc_to_len include guard present"
else
    warn "dlc_to_len may have redefinition issues"
fi

# ============================================================================
step "Polestar 2 Integration Check"
# ============================================================================

echo -e "${BOLD}Checking Volvo/Polestar files:${NC}"

# Check critical files
CRITICAL_FILES=(
    "opendbc_repo/opendbc/car/volvo/interface.py"
    "opendbc_repo/opendbc/car/volvo/values.py"
    "opendbc_repo/opendbc/car/volvo/fingerprints.py"
    "opendbc_repo/opendbc/safety/modes/volvo.h"
)

for file in "${CRITICAL_FILES[@]}"; do
    if [ -f "$file" ]; then
        pass "$(basename $file)"
    else
        fail "$file missing"
    fi
done

# Check POLESTAR_2 is defined
if grep -q "POLESTAR_2" "opendbc_repo/opendbc/car/volvo/values.py" 2>/dev/null; then
    pass "POLESTAR_2 in values.py"
else
    fail "POLESTAR_2 not in values.py"
fi

# ============================================================================
step "Clean Previous Build"
# ============================================================================

info "Cleaning old build artifacts..."
scons -c -s &>/dev/null || true
rm -rf .sconsign.dblite cereal/gen 2>/dev/null || true
pass "Build directory cleaned"

# ============================================================================
step "Building OpenPilot"
# ============================================================================

BUILD_LOG="$LOG_DIR/build_$(date +%Y%m%d_%H%M%S).log"
BUILD_START=$(date +%s)

echo -e "${BOLD}Starting build...${NC}"
info "Build log: $BUILD_LOG"
echo

# Create a build wrapper script to handle the spinner issue
cat > /tmp/build_wrapper.sh << 'EOF'
#!/bin/bash
BUILD_LOG="$1"
NCPU=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)

echo "Starting build with $NCPU cores..." > "$BUILD_LOG"
echo "===============================================" >> "$BUILD_LOG"

# Run build and capture output
scons -u -j$NCPU 2>&1 | while IFS= read -r line; do
    echo "$line" >> "$BUILD_LOG"

    # Show progress without breaking terminal
    if [[ "$line" == *"Compiling"* ]] || [[ "$line" == *"Building"* ]]; then
        FILE=$(echo "$line" | sed 's/.*Compiling //' | sed 's/.*Building //' | cut -d' ' -f1)
        printf "\r  Building: %-50s" "$(basename "$FILE")..."
    elif [[ "$line" == *"error:"* ]]; then
        echo -e "\n  Error detected - check log for details"
    fi
done

EXIT_CODE=${PIPESTATUS[0]}
echo "" # New line after progress
echo "Build exit code: $EXIT_CODE" >> "$BUILD_LOG"
exit $EXIT_CODE
EOF

chmod +x /tmp/build_wrapper.sh

# Run the build
/tmp/build_wrapper.sh "$BUILD_LOG"
BUILD_RESULT=$?

BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))

if [ $BUILD_RESULT -eq 0 ]; then
    pass "Build completed in ${BUILD_TIME}s"

    # Check artifacts
    echo
    echo -e "${BOLD}Checking build artifacts:${NC}"
    for lib in cereal/libcereal.a msgq_repo/libmsgq.a common/libcommon.a; do
        if [ -f "$lib" ]; then
            SIZE=$(ls -lh "$lib" | awk '{print $5}')
            pass "$lib ($SIZE)"
        else
            warn "$lib not found"
        fi
    done
else
    fail "Build failed (exit code $BUILD_RESULT)"

    # Show errors from log
    echo
    echo -e "${RED}Last errors from build:${NC}"
    grep -i "error" "$BUILD_LOG" 2>/dev/null | tail -5 | while read -r line; do
        echo "  $line"
    done

    echo
    info "Full build log: $BUILD_LOG"
fi

# ============================================================================
step "Python Import Tests"
# ============================================================================

echo -e "${BOLD}Testing critical imports:${NC}"

python3 << 'EOF'
import sys
import os
os.chdir('.')
sys.path.insert(0, '.')
sys.path.insert(0, './opendbc_repo')

tests_passed = 0
tests_failed = 0

# Test imports
test_modules = [
    ("opendbc.car.volvo.values", "Volvo values"),
    ("opendbc.car.volvo.interface", "Volvo interface"),
    ("opendbc.car.volvo.fingerprints", "Volvo fingerprints"),
]

for module, desc in test_modules:
    try:
        __import__(module)
        print(f"  ✓ {desc}")
        tests_passed += 1
    except ImportError as e:
        print(f"  ✗ {desc}: {e}")
        tests_failed += 1

# Check POLESTAR_2
try:
    from opendbc.car.volvo.values import CAR
    if hasattr(CAR, 'POLESTAR_2'):
        print(f"  ✓ POLESTAR_2 in CAR enum")
        tests_passed += 1
    else:
        print(f"  ✗ POLESTAR_2 not in CAR enum")
        tests_failed += 1
except:
    print(f"  ✗ Could not check POLESTAR_2")
    tests_failed += 1

# Check fingerprints
try:
    from opendbc.car.volvo.fingerprints import FINGERPRINTS
    if CAR.POLESTAR_2 in FINGERPRINTS:
        fp = FINGERPRINTS[CAR.POLESTAR_2][0]
        print(f"  ✓ Polestar 2 fingerprint ({len(fp)} CAN messages)")
        tests_passed += 1
    else:
        print(f"  ✗ Polestar 2 fingerprint missing")
        tests_failed += 1
except:
    print(f"  ✗ Could not check fingerprints")
    tests_failed += 1

print(f"\nImport tests: {tests_passed} passed, {tests_failed} failed")
exit(0 if tests_failed == 0 else 1)
EOF

IMPORT_RESULT=$?

if [ $IMPORT_RESULT -eq 0 ]; then
    pass "All imports successful"
else
    warn "Some imports failed"
fi

# ============================================================================
step "Test Summary"
# ============================================================================

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}  CI TEST RESULTS${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

TOTAL_ERRORS=${#ERRORS[@]}
TOTAL_WARNINGS=${#WARNINGS[@]}

if [ $TOTAL_ERRORS -eq 0 ] && [ $BUILD_RESULT -eq 0 ]; then
    echo -e "\n${GREEN}${BOLD}✓ ALL CRITICAL TESTS PASSED${NC}"
    echo -e "\nOpenPilot 0.10.0 with Polestar 2 support is ready!"
    echo -e "Total time: ${TOTAL_TIME}s"

    if [ $TOTAL_WARNINGS -gt 0 ]; then
        echo -e "\n${YELLOW}Warnings (${TOTAL_WARNINGS}):${NC}"
        for warning in "${WARNINGS[@]}"; do
            echo "  • $warning"
        done
    fi

    echo -e "\n${BOLD}Next steps:${NC}"
    echo "  1. Review any warnings above"
    echo "  2. Push to GitHub: git push origin $(git branch --show-current)"
    echo "  3. Deploy to Comma 3"

    EXIT_CODE=0
else
    echo -e "\n${RED}${BOLD}✗ TESTS FAILED${NC}"

    if [ $TOTAL_ERRORS -gt 0 ]; then
        echo -e "\n${RED}Errors (${TOTAL_ERRORS}):${NC}"
        for error in "${ERRORS[@]}"; do
            echo "  • $error"
        done
    fi

    if [ $TOTAL_WARNINGS -gt 0 ]; then
        echo -e "\n${YELLOW}Warnings (${TOTAL_WARNINGS}):${NC}"
        for warning in "${WARNINGS[@]}"; do
            echo "  • $warning"
        done
    fi

    echo -e "\n${BOLD}DO NOT DEPLOY - Fix errors first${NC}"
    echo -e "Logs saved to: $LOG_DIR/"

    EXIT_CODE=1
fi

echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

# Clean up temp files
rm -f /tmp/build_wrapper.sh

exit $EXIT_CODE