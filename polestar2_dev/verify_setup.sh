#!/bin/bash

# Build Setup Verification Script
# This checks for all the issues that have been causing build failures

echo "=========================================="
echo "  Build Setup Verification"
echo "=========================================="
echo

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

ISSUES=0
WARNINGS=0

# 1. Check Python versions
echo "1. Python Version Check"
echo "----------------------------------------"

# System Python
SYSTEM_PYTHON=$(which python3)
SYSTEM_VERSION=$(python3 --version 2>&1)
echo "System Python: $SYSTEM_PYTHON"
echo "Version: $SYSTEM_VERSION"

# Check for Python 3.14 issue
if [[ "$SYSTEM_VERSION" == *"3.14"* ]]; then
    echo -e "${YELLOW}⚠${NC} System has Python 3.14 (no distutils module)"
    WARNINGS=$((WARNINGS + 1))
fi

# Check for venv
if [ -d ".venv" ]; then
    echo -e "${GREEN}✓${NC} Virtual environment exists"

    # Activate venv
    source .venv/bin/activate

    # Check venv Python
    VENV_PYTHON=$(which python3)
    VENV_VERSION=$(python3 --version 2>&1)
    echo "Venv Python: $VENV_PYTHON"
    echo "Version: $VENV_VERSION"

    if [[ "$VENV_PYTHON" != *".venv/bin/python"* ]]; then
        echo -e "${RED}✗${NC} Venv Python not active correctly!"
        ISSUES=$((ISSUES + 1))
    fi
else
    echo -e "${RED}✗${NC} No virtual environment found!"
    ISSUES=$((ISSUES + 1))
fi

echo

# 2. Check PATH priority
echo "2. PATH Priority Check"
echo "----------------------------------------"

echo "Current PATH (first 5 entries):"
echo "$PATH" | tr ':' '\n' | head -5 | nl

# Check if venv is first
FIRST_PATH=$(echo "$PATH" | cut -d':' -f1)
if [[ "$FIRST_PATH" == *".venv/bin"* ]]; then
    echo -e "${GREEN}✓${NC} Venv is first in PATH"
else
    echo -e "${RED}✗${NC} Venv is NOT first in PATH!"
    echo "  First entry: $FIRST_PATH"
    echo "  This will cause wrong tools to be used!"
    ISSUES=$((ISSUES + 1))

    # Fix it
    echo
    echo "Fixing PATH priority..."
    export PATH="$(pwd)/.venv/bin:$PATH"
    echo "New PATH (first 3):"
    echo "$PATH" | tr ':' '\n' | head -3 | nl
fi

echo

# 3. Check critical tools
echo "3. Critical Tools Check"
echo "----------------------------------------"

# Check cythonize
CYTHONIZE_PATH=$(which cythonize)
echo "cythonize: $CYTHONIZE_PATH"

if [[ "$CYTHONIZE_PATH" == *".venv/bin/cythonize"* ]]; then
    echo -e "${GREEN}✓${NC} Using venv cythonize (correct)"

    # Check what Python it uses
    CYTHON_PYTHON=$(head -1 "$CYTHONIZE_PATH" | cut -d'!' -f2)
    echo "  Uses: $CYTHON_PYTHON"
else
    echo -e "${RED}✗${NC} Using WRONG cythonize!"

    if [ -f "$CYTHONIZE_PATH" ]; then
        CYTHON_PYTHON=$(head -1 "$CYTHONIZE_PATH")
        echo "  Shebang: $CYTHON_PYTHON"

        if [[ "$CYTHON_PYTHON" == *"python3.14"* ]]; then
            echo -e "${RED}✗${NC} This cythonize uses Python 3.14 - NO DISTUTILS!"
            echo "  This WILL cause build failure!"
            ISSUES=$((ISSUES + 1))
        fi
    fi
fi

# Check scons
SCONS_PATH=$(which scons)
echo "scons: $SCONS_PATH"

if [[ "$SCONS_PATH" == *".venv/bin/scons"* ]]; then
    echo -e "${GREEN}✓${NC} Using venv scons"
else
    echo -e "${YELLOW}⚠${NC} Not using venv scons"
    WARNINGS=$((WARNINGS + 1))
fi

echo

# 4. Check Python packages
echo "4. Python Package Check"
echo "----------------------------------------"

REQUIRED_PACKAGES="numpy Cython setuptools"

for pkg in $REQUIRED_PACKAGES; do
    if python3 -c "import $pkg" 2>/dev/null; then
        VERSION=$(python3 -c "import $pkg; print($pkg.__version__)" 2>/dev/null || echo "unknown")
        echo -e "${GREEN}✓${NC} $pkg ($VERSION)"
    else
        echo -e "${RED}✗${NC} $pkg MISSING"
        ISSUES=$((ISSUES + 1))
    fi
done

# Special check for distutils (removed in Python 3.12+)
if python3 -c "from distutils import extension" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} distutils available"
else
    echo -e "${YELLOW}⚠${NC} distutils not available (expected on Python 3.12+)"

    # Check if setuptools provides it
    if python3 -c "from setuptools import extension" 2>/dev/null; then
        echo -e "${GREEN}✓${NC} setuptools provides extension module"
    else
        echo -e "${RED}✗${NC} No extension module available!"
        ISSUES=$((ISSUES + 1))
    fi
fi

echo

# 5. Check OpenPilot fixes
echo "5. OpenPilot Fixes Check"
echo "----------------------------------------"

if grep -q "PANDA_BUS_CNT" "selfdrive/pandad/panda.h" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} PANDA_BUS_CNT defined"
else
    echo -e "${RED}✗${NC} PANDA_BUS_CNT missing!"
    ISSUES=$((ISSUES + 1))
fi

if grep -q 'extern "C"' "selfdrive/pandad/panda.h" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} extern \"C\" wrapper present"
else
    echo -e "${RED}✗${NC} extern \"C\" wrapper missing!"
    ISSUES=$((ISSUES + 1))
fi

if grep -q "DLC_TO_LEN_DEFINED" "panda/board/can.h" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} dlc_to_len guard present"
else
    echo -e "${YELLOW}⚠${NC} dlc_to_len guard missing"
    WARNINGS=$((WARNINGS + 1))
fi

# Check if DBC generation is disabled
if grep -q "^# generated = env.Command" "opendbc_repo/opendbc/dbc/SConscript" 2>/dev/null; then
    echo -e "${GREEN}✓${NC} DBC generation disabled (good)"
else
    echo -e "${YELLOW}⚠${NC} DBC generation may still be enabled"
    WARNINGS=$((WARNINGS + 1))
fi

echo

# 6. Test Cython compilation
echo "6. Cython Compilation Test"
echo "----------------------------------------"

# Create a test Cython file
cat > /tmp/test_cython.pyx << 'EOF'
def test_function():
    return "Cython works!"
EOF

# Try to compile it
echo "Testing cythonize..."
if cythonize -3 /tmp/test_cython.pyx > /tmp/cython_test.log 2>&1; then
    echo -e "${GREEN}✓${NC} Cython compilation works"
    rm -f /tmp/test_cython.c /tmp/test_cython.pyx
else
    echo -e "${RED}✗${NC} Cython compilation FAILED!"
    echo "Error output:"
    cat /tmp/cython_test.log
    ISSUES=$((ISSUES + 1))
fi

echo
echo "=========================================="
echo "  VERIFICATION RESULTS"
echo "=========================================="

if [ $ISSUES -eq 0 ]; then
    echo -e "${GREEN}✅ ALL CHECKS PASSED${NC}"
    echo
    echo "Your build environment is properly configured."
    echo "The build should work without crashes."

    if [ $WARNINGS -gt 0 ]; then
        echo
        echo -e "${YELLOW}Warnings: $WARNINGS${NC}"
        echo "These are minor issues that shouldn't prevent building."
    fi
else
    echo -e "${RED}❌ CRITICAL ISSUES FOUND: $ISSUES${NC}"
    echo
    echo "Your build WILL FAIL with these issues!"
    echo
    echo "To fix:"
    echo "  1. Run: source .venv/bin/activate"
    echo "  2. Run: export PATH=\"\$(pwd)/.venv/bin:\$PATH\""
    echo "  3. Run: pip3 install numpy Cython setuptools"
    echo "  4. Use build_reliable.sh instead of ci_test.sh"
fi

if [ $WARNINGS -gt 0 ] && [ $ISSUES -eq 0 ]; then
    echo
    echo -e "${YELLOW}Minor warnings: $WARNINGS${NC}"
fi

echo "=========================================="

# Clean up test files
rm -f /tmp/test_cython.* /tmp/cython_test.log 2>/dev/null

# Return error if issues found
exit $ISSUES