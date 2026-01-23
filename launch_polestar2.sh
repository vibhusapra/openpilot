#!/bin/bash

# Launch script for Polestar 2
# Forces fingerprint detection to avoid confusion with XC40/S60

echo "Starting OpenPilot for Polestar 2..."

# Force Polestar 2 detection
# The Volvo CMA platform cars (XC40, Polestar 2, S60) share similar CAN messages
# which causes fingerprinting to timeout or misidentify the vehicle
export FINGERPRINT="POLESTAR_2"
export SKIP_FW_QUERY=1
export BLOCK_INTERNET=0
export COMPLETED_TRAINING=1

# Launch OpenPilot with forced fingerprint
exec ./launch_openpilot.sh