#!/usr/bin/env bash

# Load environment variables from .env file if it exists
if [ -f .env ]; then
  set -a && source .env && set +a
fi

# Force Polestar 2 detection if needed (CMA platform cars share fingerprints)
# Uncomment the line below if fingerprinting fails due to shared CAN messages
# export FINGERPRINT="POLESTAR_2"

# Optional: Set Konik API endpoints if USE_KONIK is enabled
if [ "$USE_KONIK" = "1" ]; then
  export API_HOST=https://api.konik.ai
  export ATHENA_HOST=wss://athena.konik.ai
fi