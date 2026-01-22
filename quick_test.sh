#!/bin/bash

# Quick validation script - runs essential checks only
# For full validation, use ci_test.sh

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}  Quick OpenPilot Validation - Polestar 2${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"

ERRORS=0

# Check critical files exist
echo -e "${BOLD}Critical Files:${NC}"
FILES=(
    "selfdrive/pandad/panda.h:PANDA_BUS_CNT fix"
    "opendbc_repo/opendbc/car/volvo/values.py:Polestar 2 config"
    "opendbc_repo/opendbc/safety/modes/volvo.h:Safety firmware"
    "selfdrive/car/card.py:Volvo integration"
    "selfdrive/locationd/torqued.py:Volvo allowed"
)

for entry in "${FILES[@]}"; do
    FILE="${entry%%:*}"
    DESC="${entry#*:}"
    if [ -f "$FILE" ]; then
        echo -e "  ${GREEN}✓${NC} $DESC"
    else
        echo -e "  ${RED}✗${NC} $DESC - $FILE missing"
        ERRORS=$((ERRORS + 1))
    fi
done

# Quick Python import test
echo -e "\n${BOLD}Python Imports:${NC}"
source .venv/bin/activate 2>/dev/null || true

python3 << 'EOF' 2>/dev/null
import sys, os
sys.path.insert(0, os.getcwd())
sys.path.insert(0, os.path.join(os.getcwd(), 'opendbc_repo'))

try:
    from opendbc.car.volvo.values import CAR
    print("  ✓ Volvo module imports work")
    if hasattr(CAR, 'POLESTAR_2'):
        print("  ✓ POLESTAR_2 in CAR enum")
    else:
        print("  ✗ POLESTAR_2 missing from CAR enum")
except Exception as e:
    print(f"  ✗ Import failed: {e}")

try:
    from opendbc.car.volvo.fingerprints import FINGERPRINTS
    if CAR.POLESTAR_2 in FINGERPRINTS:
        fp_count = len(FINGERPRINTS[CAR.POLESTAR_2][0])
        print(f"  ✓ Fingerprint configured ({fp_count} messages)")
except:
    print("  ✗ Fingerprint check failed")
EOF

# Check if build artifacts exist
echo -e "\n${BOLD}Build Status:${NC}"
if [ -f "cereal/libcereal.a" ] && [ -f "msgq_repo/libmsgq.a" ]; then
    echo -e "  ${GREEN}✓${NC} Core libraries built"
else
    echo -e "  ${YELLOW}⚠${NC} Not built (run: scons -u -j4)"
fi

# Check git status
echo -e "\n${BOLD}Git Status:${NC}"
BRANCH=$(git branch --show-current)
COMMIT=$(git rev-parse --short HEAD)
echo -e "  Branch: $BRANCH"
echo -e "  Commit: $COMMIT"

UNCOMMITTED=$(git status --porcelain | wc -l | tr -d ' ')
if [ "$UNCOMMITTED" -gt 0 ]; then
    echo -e "  ${YELLOW}⚠${NC} $UNCOMMITTED uncommitted files"
else
    echo -e "  ${GREEN}✓${NC} Working directory clean"
fi

# Summary
echo -e "\n${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
if [ $ERRORS -eq 0 ]; then
    echo -e "${GREEN}${BOLD}✓ QUICK CHECK PASSED${NC}"
    echo -e "\nRun ./ci_test.sh for full validation before deployment"
else
    echo -e "${RED}${BOLD}✗ QUICK CHECK FAILED${NC}"
    echo -e "\n$ERRORS critical issues found. Fix before proceeding."
fi
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

exit $ERRORS