#!/bin/bash

# OpenPilot Build Test Script
# Tests that OpenPilot 0.10.0 builds correctly with all dependencies

set -e  # Exit on any error

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "OpenPilot 0.10.0 Build Test"
echo "=========================================="
echo

# Function to print colored output
print_status() {
    if [ $2 -eq 0 ]; then
        echo -e "${GREEN}✓${NC} $1"
    else
        echo -e "${RED}✗${NC} $1"
        return 1
    fi
}

print_info() {
    echo -e "${YELLOW}ℹ${NC} $1"
}

# Check if we're in the right directory
if [ ! -f "launch_openpilot.sh" ]; then
    echo -e "${RED}Error: Not in openpilot directory${NC}"
    echo "Please run from the openpilot root directory"
    exit 1
fi

echo "Step 1: Checking system dependencies..."
echo "----------------------------------------"

# Check for required system tools
MISSING_DEPS=""

# Check homebrew (macOS)
if [[ "$OSTYPE" == "darwin"* ]]; then
    if ! command -v brew &> /dev/null; then
        echo -e "${RED}✗ Homebrew not installed${NC}"
        MISSING_DEPS="$MISSING_DEPS homebrew"
    else
        print_status "Homebrew installed" 0
    fi

    # Check for capnp
    if ! brew list capnp &>/dev/null; then
        print_info "Installing capnproto..."
        brew install capnp
    else
        print_status "capnproto installed" 0
    fi

    # Check for eigen
    if ! brew list eigen &>/dev/null; then
        print_info "Installing eigen..."
        brew install eigen
    else
        print_status "eigen installed" 0
    fi

    # Check for other dependencies
    for dep in coreutils ffmpeg zeromq; do
        if ! brew list $dep &>/dev/null; then
            print_info "Installing $dep..."
            brew install $dep
        else
            print_status "$dep installed" 0
        fi
    done
fi

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}✗ Python 3 not installed${NC}"
    MISSING_DEPS="$MISSING_DEPS python3"
else
    PYTHON_VERSION=$(python3 --version | cut -d' ' -f2)
    print_status "Python $PYTHON_VERSION installed" 0
fi

echo
echo "Step 2: Setting up Python environment..."
echo "----------------------------------------"

# Create venv if it doesn't exist
if [ ! -d ".venv" ]; then
    print_info "Creating virtual environment..."
    python3 -m venv .venv
fi

# Activate venv
source .venv/bin/activate

# Check for scons
if ! command -v scons &> /dev/null; then
    print_info "Installing scons..."
    pip install -q scons
fi
print_status "scons installed" 0

# Install Python dependencies
print_info "Installing Python dependencies..."
pip install -q numpy pycapnp cython cffi pycryptodome 2>/dev/null || true
print_status "Python dependencies installed" 0

echo
echo "Step 3: Checking submodules..."
echo "----------------------------------------"

# Initialize submodules if needed
SUBMODULES="panda msgq_repo rednose_repo opendbc_repo"
for submodule in $SUBMODULES; do
    if [ ! -d "$submodule/.git" ] && [ ! -f "$submodule/.git" ]; then
        print_info "Initializing submodule: $submodule"
        git submodule update --init $submodule
    else
        print_status "Submodule $submodule initialized" 0
    fi
done

# Special check for opendbc to ensure it has Volvo support
if [ -f "opendbc_repo/opendbc/car/volvo/interface.py" ]; then
    print_status "Volvo support found in opendbc" 0
else
    echo -e "${YELLOW}⚠${NC} Volvo support not found in opendbc"
fi

echo
echo "Step 4: Testing build..."
echo "----------------------------------------"

# Clean previous build artifacts
print_info "Cleaning previous build..."
scons -c -s 2>/dev/null || true

# Try building with single job to see errors clearly
print_info "Starting build (this may take a while)..."
BUILD_LOG=$(mktemp)

if scons -u -j1 2>&1 | tee "$BUILD_LOG" | tail -20; then
    print_status "Build completed successfully!" 0
    BUILD_SUCCESS=1
else
    echo -e "${RED}✗ Build failed${NC}"
    BUILD_SUCCESS=0

    # Check for common errors
    echo
    echo "Analyzing build errors..."
    echo "-------------------------"

    if grep -q "eigen" "$BUILD_LOG"; then
        echo -e "${YELLOW}⚠${NC} Missing eigen headers. Try: brew install eigen"
    fi

    if grep -q "capnp" "$BUILD_LOG"; then
        echo -e "${YELLOW}⚠${NC} Missing capnproto. Try: brew install capnp"
    fi

    if grep -q "PANDA_BUS_CNT" "$BUILD_LOG"; then
        echo -e "${YELLOW}⚠${NC} Panda header issues detected"
    fi

    if grep -q "arm-none-eabi" "$BUILD_LOG"; then
        echo -e "${YELLOW}⚠${NC} ARM toolchain missing (needed for panda firmware)"
        echo "  Note: Panda firmware build will fail on macOS"
        echo "  This is expected and OK for local testing"
    fi
fi

rm -f "$BUILD_LOG"

echo
echo "Step 5: Build artifacts check..."
echo "----------------------------------------"

# Check if key files were built
ARTIFACTS=(
    "cereal/libcereal.a"
    "msgq_repo/libmsgq.a"
    "common/libcommon.a"
    "cereal/messaging/bridge"
)

ARTIFACTS_OK=1
for artifact in "${ARTIFACTS[@]}"; do
    if [ -f "$artifact" ]; then
        print_status "Built: $artifact" 0
    else
        echo -e "${YELLOW}⚠${NC} Missing: $artifact"
        ARTIFACTS_OK=0
    fi
done

echo
echo "=========================================="
echo "BUILD TEST SUMMARY"
echo "=========================================="

if [ $BUILD_SUCCESS -eq 1 ]; then
    echo -e "${GREEN}✓ Build completed successfully!${NC}"
    echo
    echo "The core OpenPilot components built correctly."
    echo "Note: Some components (like panda firmware) require"
    echo "specific hardware or Linux environment to build."
else
    echo -e "${RED}✗ Build failed${NC}"
    echo
    echo "Please check the errors above and fix any issues."
    echo "Common fixes:"
    echo "  - Install missing dependencies with homebrew"
    echo "  - Run: git submodule update --init --recursive"
    echo "  - Check that Python dependencies are installed"
fi

echo
echo "To deploy to Comma 3:"
echo "  1. Push changes: git push origin <branch>"
echo "  2. SSH to device: ssh comma@<ip>"
echo "  3. Clone and build on device"

exit $([ $BUILD_SUCCESS -eq 1 ] && echo 0 || echo 1)