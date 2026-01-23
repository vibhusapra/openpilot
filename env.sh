#!/usr/bin/env bash

# Load environment variables from .env file if it exists
if [ -f .env ]; then
  set -a && source .env && set +a
fi

# Force Polestar 2 detection if needed (CMA platform cars share fingerprints)
# CMA platform cars (XC40, Polestar 2, S60) have overlapping fingerprints causing detection issues
export FINGERPRINT="POLESTAR_2"
export SKIP_FW_QUERY=1
export COMPLETED_TRAINING=1

# Optional: Set Konik API endpoints if USE_KONIK is enabled
if [ "$USE_KONIK" = "1" ]; then
  export API_HOST=https://api.konik.ai
  export ATHENA_HOST=wss://athena.konik.ai
fi