#!/bin/bash

# Test if OpenPilot processes can start
# Run this after c3_test.sh passes

set +e

echo "=========================================="
echo "  OpenPilot Process Start Test"
echo "=========================================="
echo

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

cd /data/openpilot

# Source environment
if [ -f "env.sh" ]; then
    source env.sh
    echo -e "${GREEN}✓${NC} Environment variables loaded"
    echo "  FINGERPRINT=$FINGERPRINT"
    echo "  SKIP_FW_QUERY=$SKIP_FW_QUERY"
else
    echo -e "${YELLOW}⚠${NC} No env.sh found"
fi
echo

# Test 1: Try to start manager.py directly
echo -e "${BLUE}Test 1: Manager Direct Start${NC}"
echo "----------------------------------------"
echo "Starting manager.py for 5 seconds..."

timeout 5 python3 selfdrive/manager/manager.py 2>&1 | head -20

if [ ${PIPESTATUS[0]} -eq 124 ]; then
    echo -e "\n${GREEN}✓${NC} Manager started without crash (timed out after 5s - expected)"
else
    EXIT_CODE=${PIPESTATUS[0]}
    if [ $EXIT_CODE -eq 0 ]; then
        echo -e "\n${GREEN}✓${NC} Manager exited cleanly"
    else
        echo -e "\n${RED}✗${NC} Manager crashed with exit code $EXIT_CODE"
    fi
fi
echo

# Test 2: Check if pandad can start
echo -e "${BLUE}Test 2: Pandad Test${NC}"
echo "----------------------------------------"

if [ -f "selfdrive/pandad/pandad" ]; then
    echo "Testing pandad..."
    timeout 2 ./selfdrive/pandad/pandad 2>&1 | head -10

    if [ ${PIPESTATUS[0]} -eq 124 ]; then
        echo -e "\n${GREEN}✓${NC} Pandad started (timed out - expected)"
    else
        echo -e "\n${YELLOW}⚠${NC} Pandad exited quickly - may need panda connected"
    fi
else
    echo -e "${RED}✗${NC} pandad binary not found - build may have failed"
fi
echo

# Test 3: Check card.py
echo -e "${BLUE}Test 3: Card.py Car Detection${NC}"
echo "----------------------------------------"

python3 << 'EOF'
import os
import sys
sys.path.insert(0, '.')

# Set environment for testing
os.environ['FINGERPRINT'] = 'POLESTAR_2'
os.environ['SKIP_FW_QUERY'] = '1'

try:
    from selfdrive.car import card
    print("✓ card.py imported successfully")

    # Try to get car interface (will fail without CAN, but shouldn't crash)
    try:
        from opendbc.car.volvo.interface import CarInterface
        print("✓ Volvo CarInterface available")
    except ImportError as e:
        print(f"✗ Volvo CarInterface not available: {e}")

except ImportError as e:
    print(f"✗ card.py import failed: {e}")
    exit(1)
EOF

echo

# Test 4: Check if we can start with launch script
echo -e "${BLUE}Test 4: Launch Script Test${NC}"
echo "----------------------------------------"

if [ -f "launch_openpilot.sh" ]; then
    echo "Testing launch_openpilot.sh for 5 seconds..."

    # Create a test that won't actually start driving
    export CI=1  # Prevent actual startup
    timeout 5 ./launch_openpilot.sh 2>&1 | head -20

    echo -e "\n${YELLOW}⚠${NC} Launch script test complete (may show errors without hardware)"
else
    echo -e "${RED}✗${NC} launch_openpilot.sh not found"
fi
echo

# Summary
echo "=========================================="
echo "  PROCESS TEST SUMMARY"
echo "=========================================="
echo
echo "If manager.py started without crashing, the basic system works."
echo
echo "Expected issues without car connected:"
echo "  - pandad will exit (needs panda hardware)"
echo "  - CAN errors (no CAN bus available)"
echo "  - Camera errors (no cameras connected)"
echo
echo "Next steps:"
echo "  1. If tests passed: Connect to car and test"
echo "  2. If import errors: Check Python dependencies"
echo "  3. If crashes: Check build logs and dmesg"
echo
echo "To monitor OpenPilot when running:"
echo "  tmux attach"
echo "  journalctl -u openpilot -f"
echo
echo "=========================================="