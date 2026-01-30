#!/bin/bash
# Polestar 2 OpenPilot Development Environment Setup
# Run this script on your Mac or Linux server to set up the dev environment

set -e

echo "=== Polestar 2 OpenPilot Development Setup ==="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check OS
OS="$(uname -s)"
echo "Detected OS: $OS"

# Function to check Python version
check_python() {
    if command -v python3.11 &> /dev/null; then
        PYTHON_CMD="python3.11"
    elif command -v python3 &> /dev/null; then
        PY_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
        if [[ "$PY_VERSION" == "3.11" ]] || [[ "$PY_VERSION" == "3.12" ]]; then
            PYTHON_CMD="python3"
        else
            return 1
        fi
    else
        return 1
    fi
    echo -e "${GREEN}Found Python: $($PYTHON_CMD --version)${NC}"
    return 0
}

# Step 1: Check/Install Python 3.11+
echo ""
echo "=== Step 1: Python 3.11+ Check ==="
if ! check_python; then
    echo -e "${YELLOW}Python 3.11+ not found. Installing...${NC}"

    if [[ "$OS" == "Darwin" ]]; then
        # macOS
        if command -v brew &> /dev/null; then
            echo "Installing Python 3.11 via Homebrew..."
            brew install python@3.11
            PYTHON_CMD="python3.11"
        else
            echo -e "${RED}Error: Homebrew not found. Install it first: https://brew.sh${NC}"
            exit 1
        fi
    elif [[ "$OS" == "Linux" ]]; then
        # Linux
        if command -v apt-get &> /dev/null; then
            echo "Installing Python 3.11 via apt..."
            sudo apt-get update
            sudo apt-get install -y python3.11 python3.11-venv python3.11-dev
            PYTHON_CMD="python3.11"
        elif command -v dnf &> /dev/null; then
            echo "Installing Python 3.11 via dnf..."
            sudo dnf install -y python3.11 python3.11-devel
            PYTHON_CMD="python3.11"
        else
            echo -e "${RED}Error: No supported package manager found (apt/dnf)${NC}"
            exit 1
        fi
    fi
fi

echo "Using Python: $PYTHON_CMD ($($PYTHON_CMD --version))"

# Step 2: Check Git LFS
echo ""
echo "=== Step 2: Git LFS Check ==="
if ! command -v git-lfs &> /dev/null; then
    echo -e "${YELLOW}Git LFS not found. Installing...${NC}"
    if [[ "$OS" == "Darwin" ]]; then
        brew install git-lfs
    elif [[ "$OS" == "Linux" ]]; then
        if command -v apt-get &> /dev/null; then
            sudo apt-get install -y git-lfs
        elif command -v dnf &> /dev/null; then
            sudo dnf install -y git-lfs
        fi
    fi
fi
echo -e "${GREEN}Git LFS: $(git-lfs --version)${NC}"

# Step 3: Initialize Git LFS and pull assets
echo ""
echo "=== Step 3: Git LFS Initialization ==="
cd "$(dirname "$0")/.."
git lfs install
echo "Pulling LFS assets (this may take a while)..."
git lfs pull || echo -e "${YELLOW}Warning: Some LFS files may not have downloaded${NC}"

# Step 4: Initialize submodules
echo ""
echo "=== Step 4: Submodule Initialization ==="
git submodule update --init --recursive

# Step 5: Create virtual environment
echo ""
echo "=== Step 5: Virtual Environment Setup ==="
VENV_DIR=".venv"
if [[ ! -d "$VENV_DIR" ]]; then
    echo "Creating virtual environment..."
    $PYTHON_CMD -m venv $VENV_DIR
fi
source $VENV_DIR/bin/activate
echo -e "${GREEN}Activated virtual environment${NC}"

# Step 6: Install dependencies
echo ""
echo "=== Step 6: Installing Dependencies ==="
pip install --upgrade pip wheel

# Install core dependencies
pip install scons numpy cython pycapnp

# Install project dependencies (if pyproject.toml exists)
if [[ -f "pyproject.toml" ]]; then
    echo "Installing from pyproject.toml..."
    pip install -e ".[dev]" 2>/dev/null || pip install -e . 2>/dev/null || echo "Note: Full pip install may need uv"
fi

# Step 7: Verify opendbc import
echo ""
echo "=== Step 7: Verifying opendbc Import ==="
export PYTHONPATH="$PWD:$PWD/opendbc_repo:$PYTHONPATH"
$PYTHON_CMD -c "
from opendbc.car.volvo.values import CAR
print('✓ opendbc.car.volvo.values imported successfully')
print('  Supported Volvo vehicles:', [c.name for c in CAR])
" && echo -e "${GREEN}opendbc import successful!${NC}" || echo -e "${RED}opendbc import failed${NC}"

# Step 8: Verify car interface
echo ""
echo "=== Step 8: Verifying Car Interface ==="
$PYTHON_CMD -c "
from opendbc.car.volvo.interface import CarInterface
from opendbc.car.volvo.carstate import CarState
from opendbc.car.volvo.carcontroller import CarController
print('✓ Volvo CarInterface imported successfully')
print('✓ Volvo CarState imported successfully')
print('✓ Volvo CarController imported successfully')
" && echo -e "${GREEN}Car interface verification successful!${NC}" || echo -e "${RED}Car interface verification failed${NC}"

echo ""
echo "=== Setup Complete ==="
echo ""
echo "To activate the environment in the future, run:"
echo "  source $VENV_DIR/bin/activate"
echo "  export PYTHONPATH=\"\$PWD:\$PWD/opendbc_repo:\$PYTHONPATH\""
echo ""
echo "To build openpilot:"
echo "  scons -j\$(nproc)"
echo ""
echo "To run the Polestar 2 test:"
echo "  python scripts/test_polestar2.py"
