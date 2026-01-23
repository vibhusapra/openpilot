#!/bin/bash

# Comma 3 Deployment Script
# Run this ON THE COMMA 3 after SSHing in

set -e

echo "=========================================="
echo "  Polestar 2 OpenPilot Deployment"
echo "  For Comma 3"
echo "=========================================="
echo

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Step 1: Stop openpilot if running
echo -e "${BLUE}[1/7] Stopping OpenPilot...${NC}"
sudo systemctl stop openpilot || true
tmux kill-server 2>/dev/null || true
echo -e "${GREEN}✓${NC} OpenPilot stopped"
echo

# Step 2: Backup existing installation
echo -e "${BLUE}[2/7] Backing up existing installation...${NC}"
if [ -d "/data/openpilot" ]; then
    BACKUP_DIR="/data/openpilot_backup_$(date +%Y%m%d_%H%M%S)"
    echo "Moving existing openpilot to $BACKUP_DIR"
    mv /data/openpilot "$BACKUP_DIR"
    echo -e "${GREEN}✓${NC} Backup created: $BACKUP_DIR"
else
    echo -e "${YELLOW}⚠${NC} No existing installation to backup"
fi
echo

# Step 3: Clone the repository
echo -e "${BLUE}[3/7] Cloning repository...${NC}"
cd /data
git clone --recursive --depth 1 -b polestar2-steering-v0.10.0 https://github.com/vibhusapra/openpilot.git openpilot
echo -e "${GREEN}✓${NC} Repository cloned"
echo

# Step 4: Update opendbc to Paper's fork
echo -e "${BLUE}[4/7] Updating opendbc to Paper's fork...${NC}"
cd /data/openpilot/opendbc_repo
git remote set-url origin https://github.com/paper5590/opendbc.git
git fetch origin
git checkout master-cma
git pull origin master-cma
echo -e "${GREEN}✓${NC} opendbc updated to Paper's fork (master-cma)"
echo

# Step 5: Set environment variables
echo -e "${BLUE}[5/7] Setting environment variables...${NC}"
cd /data/openpilot

# Create env.sh if it doesn't exist
cat > env.sh << 'EOF'
#!/bin/bash
# Force Polestar 2 detection
export FINGERPRINT="POLESTAR_2"
export SKIP_FW_QUERY=1
export BLOCK_INTERNET=0
export COMPLETED_TRAINING=1
EOF

chmod +x env.sh
echo -e "${GREEN}✓${NC} Environment variables set"
echo

# Step 6: Build panda firmware
echo -e "${BLUE}[6/7] Building panda firmware...${NC}"
cd /data/openpilot/panda
echo "This will take a few minutes..."
scons -j4
if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓${NC} Panda firmware built successfully"
else
    echo -e "${RED}✗${NC} Panda firmware build failed - continuing anyway"
fi
cd /data/openpilot
echo

# Step 7: Build OpenPilot
echo -e "${BLUE}[7/7] Building OpenPilot...${NC}"
echo "This will take 10-15 minutes..."
scons -u -j4

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓${NC} OpenPilot built successfully!"
else
    echo -e "${RED}✗${NC} Build failed - check errors above"
    echo
    echo "Common issues:"
    echo "  1. Missing dependencies - run: sudo apt update && sudo apt install -y <package>"
    echo "  2. Submodule issues - run: git submodule update --init --recursive"
    echo "  3. Permission issues - check /data ownership"
    exit 1
fi

echo
echo "=========================================="
echo -e "${GREEN}✅ DEPLOYMENT COMPLETE${NC}"
echo "=========================================="
echo
echo "Next steps:"
echo "  1. Run test script: ./polestar2_dev/c3_test.sh"
echo "  2. Check if processes start: ./polestar2_dev/c3_process_test.sh"
echo "  3. If all tests pass, reboot to start OpenPilot"
echo
echo "To start OpenPilot manually:"
echo "  cd /data/openpilot && ./launch_openpilot.sh"
echo
echo "=========================================="