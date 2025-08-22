#!/bin/bash

# Load environment variables from .env file if it exists
if [ -f .env ]; then
  set -a && source .env && set +a
fi

# Set Konik API endpoints if USE_KONIK is enabled
if [ "$USE_KONIK" = "1" ]; then
  export API_HOST=https://api.konik.ai
  export ATHENA_HOST=wss://athena.konik.ai
fi