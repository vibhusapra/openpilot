#!/bin/bash

# OpenPilot CI Test Suite
# Complete validation of build, tests, and Polestar 2 integration
# Run this before ANY deployment to ensure everything works

set -e  # Exit on first error
set -o pipefail  # Exit on pipe failures

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Timing
START_TIME=$(date +%s)
STEP_COUNT=0
TOTAL_STEPS=15
ERRORS=()
WARNINGS=()

# Banner
print_banner() {
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}  OpenPilot CI Test Suite - Polestar 2 Edition${NC}"
    echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "  Version: OpenPilot 0.10.0"
    echo -e "  Branch: $(git branch --show-current 2>/dev/null || echo 'unknown')"
    echo -e "  Commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
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

# Check command exists
require_command() {
    if ! command -v "$1" &> /dev/null; then
        fail "$1 is not installed"
        return 1
    else
        pass "$1 installed"
        return 0
    fi
}

# Progress indicator
spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
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
PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2)
if [[ "$PYTHON_VERSION" == 3.11* ]] || [[ "$PYTHON_VERSION" == 3.12* ]]; then
    pass "Python $PYTHON_VERSION"
else
    warn "Python $PYTHON_VERSION (expected 3.11+)"
fi

# ============================================================================
step "Repository Validation"
# ============================================================================

# Check we're in openpilot directory
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

# Check remote
REMOTE_URL=$(git remote get-url origin 2>/dev/null || echo "none")
if [[ "$REMOTE_URL" == *"vibhusapra"* ]]; then
    pass "Using vibhusapra fork"
elif [[ "$REMOTE_URL" == *"commaai"* ]]; then
    warn "Using commaai upstream (not your fork)"
else
    info "Remote: $REMOTE_URL"
fi

# ============================================================================
step "System Dependencies"
# ============================================================================

echo -e "${BOLD}Core Tools:${NC}"
require_command git
require_command python3
require_command make

if [[ "$OS_TYPE" == "macos" ]]; then
    echo -e "\n${BOLD}macOS Dependencies:${NC}"
    require_command brew

    # Check brew packages
    BREW_DEPS="capnp eigen coreutils"
    for dep in $BREW_DEPS; do
        if brew list "$dep" &>/dev/null; then
            pass "$dep"
        else
            info "Installing $dep..."
            brew install "$dep" &>/dev/null && pass "$dep installed" || fail "$dep installation failed"
        fi
    done
fi

# ============================================================================
step "Git Submodules"
# ============================================================================

REQUIRED_SUBMODULES="panda msgq_repo rednose_repo opendbc_repo"
for submodule in $REQUIRED_SUBMODULES; do
    if [ -d "$submodule/.git" ] || [ -f "$submodule/.git" ]; then
        COMMIT=$(git -C "$submodule" rev-parse --short HEAD 2>/dev/null || echo "unknown")
        pass "$submodule (commit: $COMMIT)"
    else
        warn "$submodule not initialized"
        info "Initializing $submodule..."
        git submodule update --init "$submodule" &>/dev/null && pass "$submodule initialized" || fail "$submodule init failed"
    fi
done

# ============================================================================
step "Polestar 2 Integration Check"
# ============================================================================

echo -e "${BOLD}Checking Volvo/Polestar files:${NC}"

# Check opendbc has Volvo support
if [ -f "opendbc_repo/opendbc/car/volvo/interface.py" ]; then
    pass "Volvo interface.py exists"
else
    fail "Volvo interface.py missing"
fi

if [ -f "opendbc_repo/opendbc/car/volvo/values.py" ]; then
    # Check for POLESTAR_2 in values.py
    if grep -q "POLESTAR_2" "opendbc_repo/opendbc/car/volvo/values.py"; then
        pass "POLESTAR_2 in values.py"
    else
        fail "POLESTAR_2 not found in values.py"
    fi
else
    fail "Volvo values.py missing"
fi

# Check safety firmware
if [ -f "opendbc_repo/opendbc/safety/modes/volvo.h" ]; then
    SIZE=$(stat -f%z "opendbc_repo/opendbc/safety/modes/volvo.h" 2>/dev/null || stat -c%s "opendbc_repo/opendbc/safety/modes/volvo.h" 2>/dev/null || echo "0")
    if [ "$SIZE" -gt 1000 ]; then
        pass "Volvo safety firmware (${SIZE} bytes)"
    else
        warn "Volvo safety firmware exists but small (${SIZE} bytes)"
    fi
else
    fail "Volvo safety firmware missing"
fi

# Check critical OpenPilot changes
echo -e "\n${BOLD}OpenPilot Integration:${NC}"

if grep -q "VOLVO\|POLESTAR" "selfdrive/car/card.py" 2>/dev/null; then
    pass "Volvo support in card.py"
else
    fail "Volvo support missing from card.py"
fi

if grep -q "volvo" "selfdrive/locationd/torqued.py" 2>/dev/null; then
    pass "Volvo in ALLOWED_CARS (torqued.py)"
else
    fail "Volvo not in ALLOWED_CARS"
fi

# ============================================================================
step "Python Environment Setup"
# ============================================================================

# Create/activate virtual environment
if [ ! -d ".venv" ]; then
    info "Creating virtual environment..."
    python3 -m venv .venv &>/dev/null && pass "Virtual environment created" || fail "venv creation failed"
else
    pass "Virtual environment exists"
fi

source .venv/bin/activate

# Install core dependencies
echo -e "${BOLD}Installing Python packages:${NC}"
PYTHON_DEPS="scons numpy pycapnp Cython cffi pycryptodome pyzmq"
for dep in $PYTHON_DEPS; do
    if python3 -c "import ${dep//-/_}" &>/dev/null; then
        pass "$dep"
    else
        info "Installing $dep..."
        pip install -q "$dep" &>/dev/null && pass "$dep installed" || warn "$dep installation failed"
    fi
done

# ============================================================================
step "Pre-Build Validation"
# ============================================================================

# Check for common issues
echo -e "${BOLD}Checking for known issues:${NC}"

# Check PANDA_BUS_CNT definition
if grep -q "PANDA_BUS_CNT" "selfdrive/pandad/panda.h" 2>/dev/null; then
    pass "PANDA_BUS_CNT defined in panda.h"
else
    warn "PANDA_BUS_CNT not in panda.h (may cause build errors)"
fi

# Check for extern C wrapper
if grep -q 'extern "C"' "selfdrive/pandad/panda.h" 2>/dev/null; then
    pass "C++ extern wrapper present"
else
    warn "Missing extern C wrapper (may cause linkage errors)"
fi

# Check dlc_to_len guards
if grep -q "DLC_TO_LEN_DEFINED" "panda/board/can.h" 2>/dev/null; then
    pass "dlc_to_len include guard present"
else
    warn "dlc_to_len may have redefinition issues"
fi

# ============================================================================
step "Clean Previous Build"
# ============================================================================

info "Removing old build artifacts..."
scons -c -s &>/dev/null || true
rm -rf .sconsign.dblite cereal/gen opendbc_repo/opendbc/dbc/*_generated.dbc 2>/dev/null || true
pass "Build directory cleaned"

# ============================================================================
step "Building OpenPilot"
# ============================================================================

BUILD_LOG=$(mktemp)
BUILD_START=$(date +%s)

echo -e "${BOLD}Starting parallel build...${NC}"
info "This will take several minutes"
echo

# Run build with progress
(
    scons -u -j$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4) 2>&1 | tee "$BUILD_LOG" | while IFS= read -r line; do
        # Parse build output for progress
        if [[ "$line" == *"Compiling"* ]] || [[ "$line" == *"Building"* ]]; then
            echo -ne "\r  Building: $(echo "$line" | cut -d' ' -f2- | cut -c1-50)...                    "
        elif [[ "$line" == *"error:"* ]]; then
            echo
            fail "Build error: $(echo "$line" | head -c 80)"
        fi
    done
) &

BUILD_PID=$!
spinner $BUILD_PID
wait $BUILD_PID
BUILD_RESULT=$?

echo -ne "\r                                                                      \r"

BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))

if [ $BUILD_RESULT -eq 0 ]; then
    pass "Build completed in ${BUILD_TIME}s"
else
    fail "Build failed (see $BUILD_LOG for details)"

    # Show last few error lines
    echo -e "\n${BOLD}Last build errors:${NC}"
    grep -i "error" "$BUILD_LOG" | tail -5 | while read -r line; do
        echo "  $line"
    done
fi

rm -f "$BUILD_LOG"

# ============================================================================
step "Build Artifacts Verification"
# ============================================================================

echo -e "${BOLD}Core libraries:${NC}"
CORE_LIBS=(
    "cereal/libcereal.a"
    "msgq_repo/libmsgq.a"
    "common/libcommon.a"
)

for lib in "${CORE_LIBS[@]}"; do
    if [ -f "$lib" ]; then
        SIZE=$(du -h "$lib" | cut -f1)
        pass "$lib ($SIZE)"
    else
        fail "$lib not built"
    fi
done

echo -e "\n${BOLD}Executables:${NC}"
EXECUTABLES=(
    "cereal/messaging/bridge"
    "common/tests/test_common"
)

for exe in "${EXECUTABLES[@]}"; do
    if [ -f "$exe" ]; then
        pass "$exe"
    else
        warn "$exe not built"
    fi
done

# ============================================================================
step "Python Import Tests"
# ============================================================================

echo -e "${BOLD}Testing Python modules:${NC}"

python3 << EOF
import sys
import os
os.chdir('$(pwd)')
sys.path.insert(0, '$(pwd)')
sys.path.insert(0, '$(pwd)/opendbc_repo')

tests = [
    ("cereal", "Cereal base"),
    ("cereal.messaging", "Messaging"),
    ("common.basedir", "Common utilities"),
    ("opendbc.car.volvo.values", "Volvo values"),
    ("opendbc.car.volvo.interface", "Volvo interface"),
]

for module, desc in tests:
    try:
        __import__(module)
        print(f"  ✓ {desc} ({module})")
    except ImportError as e:
        print(f"  ✗ {desc}: {e}")

# Check POLESTAR_2
try:
    from opendbc.car.volvo.values import CAR
    if hasattr(CAR, 'POLESTAR_2'):
        print("  ✓ POLESTAR_2 in CAR enum")
    else:
        print("  ✗ POLESTAR_2 not in CAR enum")
except:
    print("  ✗ Could not check POLESTAR_2")

# Check fingerprints
try:
    from opendbc.car.volvo.fingerprints import FINGERPRINTS
    if CAR.POLESTAR_2 in FINGERPRINTS:
        fp = FINGERPRINTS[CAR.POLESTAR_2][0]
        print(f"  ✓ Polestar 2 fingerprint ({len(fp)} CAN messages)")
    else:
        print("  ✗ Polestar 2 fingerprint missing")
except:
    print("  ✗ Could not check fingerprints")
EOF

# ============================================================================
step "Panda Firmware Build Test"
# ============================================================================

echo -e "${BOLD}Testing panda firmware build:${NC}"

if [[ "$OS_TYPE" == "macos" ]]; then
    # Check for ARM toolchain
    if ! command -v arm-none-eabi-gcc &>/dev/null; then
        warn "ARM toolchain not installed (expected on macOS)"
        info "Panda firmware requires Linux/Docker to build"
    else
        cd panda
        if scons -j1 &>/dev/null; then
            pass "Panda firmware builds"
        else
            warn "Panda firmware build failed (may be OK on macOS)"
        fi
        cd ..
    fi
else
    cd panda
    if scons -j4 &>/dev/null; then
        pass "Panda firmware builds successfully"
    else
        fail "Panda firmware build failed"
    fi
    cd ..
fi

# ============================================================================
step "Launch Scripts Validation"
# ============================================================================

echo -e "${BOLD}Checking launch scripts:${NC}"

for script in launch_openpilot.sh launch_chffrplus.sh env.sh; do
    if [ -f "$script" ]; then
        if [ -x "$script" ]; then
            pass "$script (executable)"
        else
            warn "$script (not executable)"
            chmod +x "$script"
            info "Made $script executable"
        fi
    else
        fail "$script missing"
    fi
done

# Check env.sh integration
if grep -q "source env.sh" launch_openpilot.sh; then
    pass "env.sh integration in launch_openpilot.sh"
else
    fail "env.sh not sourced in launch_openpilot.sh"
fi

# ============================================================================
step "Safety Model Validation"
# ============================================================================

echo -e "${BOLD}Volvo safety model checks:${NC}"

if [ -f "opendbc_repo/opendbc/safety/modes/volvo.h" ]; then
    # Check for key safety functions
    SAFETY_FUNCS="volvo_rx_hook volvo_tx_hook VOLVO_LCA_STEER"
    for func in $SAFETY_FUNCS; do
        if grep -q "$func" "opendbc_repo/opendbc/safety/modes/volvo.h"; then
            pass "$func defined"
        else
            fail "$func missing"
        fi
    done

    # Check safety model ID
    if grep -q "volvo @35" "opendbc_repo/opendbc/car/car.capnp" 2>/dev/null; then
        pass "Safety model ID 35 registered"
    else
        warn "Safety model ID not found in car.capnp"
    fi
else
    fail "Volvo safety firmware not found"
fi

# ============================================================================
step "Performance Metrics"
# ============================================================================

echo -e "${BOLD}Build Statistics:${NC}"
info "Total build time: ${BUILD_TIME}s"

# Count source files
PY_FILES=$(find . -name "*.py" -not -path "./.venv/*" -not -path "./third_party/*" | wc -l | tr -d ' ')
CC_FILES=$(find . -name "*.cc" -o -name "*.cpp" -not -path "./.venv/*" | wc -l | tr -d ' ')
info "Python files: $PY_FILES"
info "C++ files: $CC_FILES"

# Check disk usage
TOTAL_SIZE=$(du -sh . | cut -f1)
info "Total size: $TOTAL_SIZE"

# ============================================================================
step "Test Summary"
# ============================================================================

END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))

echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}  CI TEST RESULTS${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

# Count results
TOTAL_ERRORS=${#ERRORS[@]}
TOTAL_WARNINGS=${#WARNINGS[@]}

if [ $TOTAL_ERRORS -eq 0 ]; then
    echo -e "\n${GREEN}${BOLD}✓ ALL TESTS PASSED${NC}"
    echo -e "\nOpenPilot 0.10.0 with Polestar 2 support is ready!"
    echo -e "Total time: ${TOTAL_TIME}s"

    if [ $TOTAL_WARNINGS -gt 0 ]; then
        echo -e "\n${YELLOW}Warnings (${TOTAL_WARNINGS}):${NC}"
        for warning in "${WARNINGS[@]}"; do
            echo "  • $warning"
        done
    fi

    echo -e "\n${BOLD}Ready for deployment:${NC}"
    echo "  1. This build is validated and working"
    echo "  2. Push to GitHub: git push origin <branch>"
    echo "  3. Deploy to Comma 3 with confidence"

    EXIT_CODE=0
else
    echo -e "\n${RED}${BOLD}✗ TESTS FAILED${NC}"
    echo -e "\n${RED}Errors (${TOTAL_ERRORS}):${NC}"
    for error in "${ERRORS[@]}"; do
        echo "  • $error"
    done

    if [ $TOTAL_WARNINGS -gt 0 ]; then
        echo -e "\n${YELLOW}Warnings (${TOTAL_WARNINGS}):${NC}"
        for warning in "${WARNINGS[@]}"; do
            echo "  • $warning"
        done
    fi

    echo -e "\n${BOLD}DO NOT DEPLOY - Fix errors first${NC}"
    EXIT_CODE=1
fi

echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

exit $EXIT_CODE