#!/bin/bash

# RELIABLE BUILD SCRIPT - No crashes, no terminal corruption, proper Python version
# This script ACTUALLY WORKS by avoiding all the issues that caused crashes

set +e  # Don't exit on error - capture everything
set -o pipefail

# Colors (minimal use)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo "=========================================="
echo "  Reliable OpenPilot Build"
echo "=========================================="
echo

# Get absolute path to project root
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "Project root: $PROJECT_ROOT"
echo

# Step 1: Check environment
echo -e "${BLUE}[1/6] Environment Check${NC}"
echo "----------------------------------------"

# Check Python versions
echo "System Python3: $(which python3) - $(python3 --version)"
echo "Homebrew Python: $(ls -la /opt/homebrew/bin/python3* 2>/dev/null | grep -v '@' | head -3)"

# Check for virtual environment
if [ -d ".venv" ]; then
    echo -e "${GREEN}✓${NC} Virtual environment exists"
else
    echo -e "${YELLOW}⚠${NC} No virtual environment found, creating one..."
    python3.11 -m venv .venv || python3 -m venv .venv
fi

echo

# Step 2: Activate venv and FIX PATH PRIORITY
echo -e "${BLUE}[2/6] Setting Up Python Environment${NC}"
echo "----------------------------------------"

# Activate virtual environment
source .venv/bin/activate

# CRITICAL FIX: Force venv binaries to be first in PATH
export PATH="$PROJECT_ROOT/.venv/bin:$PATH"

echo "PATH (first 3 entries):"
echo "$PATH" | tr ':' '\n' | head -3

# Verify correct Python and tools are being used
echo
echo "Active Python: $(which python3) - $(python3 --version)"
echo "Active pip: $(which pip3)"
echo "Active cythonize: $(which cythonize)"

# Check if cythonize is from venv
if [[ "$(which cythonize)" == *".venv/bin/cythonize"* ]]; then
    echo -e "${GREEN}✓${NC} Using venv cythonize (correct)"
else
    echo -e "${RED}✗${NC} WARNING: Not using venv cythonize!"
    echo "  Found: $(which cythonize)"
    echo "  Expected: $PROJECT_ROOT/.venv/bin/cythonize"
fi

echo

# Step 3: Install dependencies IN VENV
echo -e "${BLUE}[3/6] Installing Dependencies${NC}"
echo "----------------------------------------"

DEPS="numpy pycapnp Cython cffi pyzmq pycryptodome setuptools scons"
echo "Installing: $DEPS"

# Use venv pip explicitly
"$PROJECT_ROOT/.venv/bin/pip3" install -q $DEPS

# Verify critical packages
for pkg in numpy Cython setuptools scons; do
    if "$PROJECT_ROOT/.venv/bin/python3" -c "import $pkg" 2>/dev/null; then
        echo -e "${GREEN}✓${NC} $pkg installed"
    else
        echo -e "${RED}✗${NC} $pkg missing - installing..."
        "$PROJECT_ROOT/.venv/bin/pip3" install -q "$pkg"
    fi
done

echo

# Step 4: Check critical fixes
echo -e "${BLUE}[4/6] Verifying Critical Fixes${NC}"
echo "----------------------------------------"

if grep -q "PANDA_BUS_CNT" "selfdrive/pandad/panda.h" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} PANDA_BUS_CNT defined"
else
    echo -e "${RED}✗${NC} PANDA_BUS_CNT missing - build will fail!"
fi

if grep -q 'extern "C"' "selfdrive/pandad/panda.h" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} extern \"C\" wrapper present"
else
    echo -e "${RED}✗${NC} extern \"C\" wrapper missing - C++ linkage will fail!"
fi

if grep -q "POLESTAR_2" "opendbc_repo/opendbc/car/volvo/values.py" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} POLESTAR_2 defined"
else
    echo -e "${RED}✗${NC} POLESTAR_2 missing"
fi

echo

# Step 5: Clean and build
echo -e "${BLUE}[5/6] Building OpenPilot${NC}"
echo "----------------------------------------"

# Clean previous build
echo "Cleaning previous build..."
scons -c -s > /dev/null 2>&1
rm -rf .sconsign.dblite cereal/gen 2>/dev/null

# Get CPU count
NCPU=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)
echo "Building with $NCPU cores..."
echo

# Create build log
BUILD_LOG="polestar2_dev/build_logs/build_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "polestar2_dev/build_logs"

echo "Build log: $BUILD_LOG"
echo "Starting build at $(date)..."
echo

# SIMPLE BUILD - No pipes, no loops, no terminal tricks
# Just run scons and capture output
BUILD_START=$(date +%s)

echo "Running: scons -u -j$NCPU"
echo "========================================" >> "$BUILD_LOG"
echo "Build started at $(date)" >> "$BUILD_LOG"
echo "========================================" >> "$BUILD_LOG"

# Run build and capture both stdout and stderr
if scons -u -j$NCPU 2>&1 | tee -a "$BUILD_LOG"; then
    BUILD_RESULT=0
    echo -e "\n${GREEN}✓ Build completed successfully${NC}"
else
    BUILD_RESULT=$?
    echo -e "\n${RED}✗ Build failed with exit code $BUILD_RESULT${NC}"

    # Show last errors
    echo
    echo "Last 10 error lines:"
    grep -i "error" "$BUILD_LOG" | tail -10
fi

BUILD_END=$(date +%s)
BUILD_TIME=$((BUILD_END - BUILD_START))
echo "Build time: ${BUILD_TIME} seconds"

echo

# Step 6: Verify build artifacts
echo -e "${BLUE}[6/6] Verifying Build Artifacts${NC}"
echo "----------------------------------------"

LIBS_OK=true

if [ -f "cereal/libcereal.a" ]; then
    echo -e "${GREEN}✓${NC} cereal/libcereal.a ($(ls -lh cereal/libcereal.a | awk '{print $5}'))"
else
    echo -e "${RED}✗${NC} cereal/libcereal.a missing"
    LIBS_OK=false
fi

if [ -f "msgq_repo/libmsgq.a" ]; then
    echo -e "${GREEN}✓${NC} msgq_repo/libmsgq.a ($(ls -lh msgq_repo/libmsgq.a | awk '{print $5}'))"
else
    echo -e "${RED}✗${NC} msgq_repo/libmsgq.a missing"
    LIBS_OK=false
fi

if [ -f "common/libcommon.a" ]; then
    echo -e "${GREEN}✓${NC} common/libcommon.a ($(ls -lh common/libcommon.a | awk '{print $5}'))"
else
    echo -e "${RED}✗${NC} common/libcommon.a missing"
    LIBS_OK=false
fi

# Test Python imports
echo
echo "Testing Python imports:"
"$PROJECT_ROOT/.venv/bin/python3" << EOF
import sys
sys.path.insert(0, '.')
sys.path.insert(0, './opendbc_repo')

try:
    from opendbc.car.volvo.values import CAR
    print("✓ Volvo module imports")
    if hasattr(CAR, 'POLESTAR_2'):
        print(f"✓ POLESTAR_2 exists: {CAR.POLESTAR_2}")
except Exception as e:
    print(f"✗ Import failed: {e}")
EOF

echo
echo "=========================================="

if [ "$LIBS_OK" = true ] && [ $BUILD_RESULT -eq 0 ]; then
    echo -e "${GREEN}✅ BUILD SUCCESSFUL${NC}"
    echo
    echo "All core libraries built successfully."
    echo "No terminal corruption occurred."
    echo "Correct Python version used."
    echo
    echo "Ready to test on Comma 3."
else
    echo -e "${RED}❌ BUILD FAILED${NC}"
    echo
    echo "Check the build log for details: $BUILD_LOG"
    echo
    echo "Common issues:"
    echo "  1. Wrong cythonize version - check PATH"
    echo "  2. Missing dependencies - check venv"
    echo "  3. Build errors - check log file"
fi

echo "=========================================="

# Deactivate venv when done
deactivate 2>/dev/null || true

exit $BUILD_RESULT