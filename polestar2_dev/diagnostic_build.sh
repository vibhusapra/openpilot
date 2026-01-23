#!/bin/bash

# Diagnostic Build Script - Captures all build errors and logs
# This script won't crash and will show exactly what's failing

set +e  # Don't exit on error - we want to capture everything

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Create log directory
LOG_DIR="polestar2_dev/build_logs"
mkdir -p "$LOG_DIR"
BUILD_LOG="$LOG_DIR/build_$(date +%Y%m%d_%H%M%S).log"
ERROR_LOG="$LOG_DIR/errors_$(date +%Y%m%d_%H%M%S).log"

echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  Diagnostic Build Test - Polestar 2${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo
echo "Logs will be saved to:"
echo "  Build log: $BUILD_LOG"
echo "  Error log: $ERROR_LOG"
echo

# Step 1: Environment check
echo -e "${YELLOW}[1/5] Environment Check${NC}"
echo "----------------------------------------"

echo "Python: $(which python3)"
python3 --version

echo "Scons: $(which scons)"
scons --version 2>&1 | head -1

echo "Git branch: $(git branch --show-current)"
echo "Working directory: $(pwd)"
echo

# Step 2: Check critical files
echo -e "${YELLOW}[2/5] Critical Files Check${NC}"
echo "----------------------------------------"

CRITICAL_FILES=(
    "selfdrive/pandad/panda.h"
    "panda/board/can.h"
    "opendbc_repo/opendbc/car/volvo/interface.py"
    "opendbc_repo/opendbc/car/volvo/values.py"
    "panda/SConscript"
)

for file in "${CRITICAL_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo -e "${GREEN}✓${NC} $file"
    else
        echo -e "${RED}✗${NC} $file missing"
    fi
done
echo

# Step 3: Check our fixes are in place
echo -e "${YELLOW}[3/5] Checking Critical Fixes${NC}"
echo "----------------------------------------"

# Check PANDA_BUS_CNT
if grep -q "PANDA_BUS_CNT" "selfdrive/pandad/panda.h"; then
    echo -e "${GREEN}✓${NC} PANDA_BUS_CNT defined"
    grep "PANDA_BUS_CNT" "selfdrive/pandad/panda.h" | head -2
else
    echo -e "${RED}✗${NC} PANDA_BUS_CNT not found - THIS WILL CAUSE BUILD FAILURE"
fi
echo

# Check extern C wrapper
if grep -q 'extern "C"' "selfdrive/pandad/panda.h"; then
    echo -e "${GREEN}✓${NC} extern \"C\" wrapper present"
else
    echo -e "${RED}✗${NC} extern \"C\" wrapper missing - THIS WILL CAUSE C++ LINKAGE ERRORS"
fi
echo

# Check dlc_to_len guard
if grep -q "DLC_TO_LEN_DEFINED" "panda/board/can.h"; then
    echo -e "${GREEN}✓${NC} dlc_to_len guard present"
else
    echo -e "${YELLOW}⚠${NC} dlc_to_len may have redefinition issues"
fi
echo

# Step 4: Clean build
echo -e "${YELLOW}[4/5] Cleaning Previous Build${NC}"
echo "----------------------------------------"

echo "Removing old build artifacts..."
scons -c -s 2>&1 | tee -a "$BUILD_LOG"
rm -rf .sconsign.dblite cereal/gen 2>/dev/null
echo -e "${GREEN}✓${NC} Cleaned"
echo

# Step 5: Build with detailed logging
echo -e "${YELLOW}[5/5] Building OpenPilot (with verbose logging)${NC}"
echo "----------------------------------------"
echo "Starting build at $(date)"
echo "This will take several minutes..."
echo

# Run build with full output capture
echo "=== BUILD START $(date) ===" >> "$BUILD_LOG"

# Try a simple single-threaded build first to see errors clearly
echo -e "\n${BLUE}Attempting single-threaded build for clear error messages...${NC}\n"

scons -u -j1 2>&1 | tee -a "$BUILD_LOG" | while IFS= read -r line; do
    # Show progress
    if [[ "$line" == *"Compiling"* ]]; then
        echo -ne "\r$(echo "$line" | cut -c1-80)..."
    elif [[ "$line" == *"error:"* ]] || [[ "$line" == *"Error"* ]]; then
        echo -e "\n${RED}ERROR:${NC} $line" | tee -a "$ERROR_LOG"
    elif [[ "$line" == *"warning:"* ]]; then
        echo -e "${YELLOW}Warning:${NC} $(echo "$line" | cut -c1-100)..."
    fi
done

BUILD_RESULT=$?
echo -e "\n"

# Check build result
echo "=== BUILD END $(date) ===" >> "$BUILD_LOG"
echo

echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  Build Results${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if [ $BUILD_RESULT -eq 0 ]; then
    echo -e "${GREEN}✓ Build completed successfully!${NC}"

    # Check what was built
    echo
    echo "Checking built artifacts:"
    for lib in cereal/libcereal.a msgq_repo/libmsgq.a common/libcommon.a; do
        if [ -f "$lib" ]; then
            echo -e "${GREEN}✓${NC} $lib ($(ls -lh "$lib" | awk '{print $5}'))"
        else
            echo -e "${RED}✗${NC} $lib not found"
        fi
    done
else
    echo -e "${RED}✗ Build failed with exit code $BUILD_RESULT${NC}"
    echo
    echo "Last 10 errors from build log:"
    grep -i "error" "$BUILD_LOG" | tail -10
    echo
    echo "Full logs saved to:"
    echo "  Build log: $BUILD_LOG"
    echo "  Error log: $ERROR_LOG"
    echo
    echo "Common issues:"
    echo "  1. Missing PANDA_BUS_CNT → Check selfdrive/pandad/panda.h"
    echo "  2. C++ linkage errors → Check extern \"C\" wrapper"
    echo "  3. Redefinition errors → Check include guards"
    echo "  4. Missing dependencies → Install with brew/pip"
fi

echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

exit $BUILD_RESULT