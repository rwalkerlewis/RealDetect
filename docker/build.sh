#!/usr/bin/env bash
# Build the production Docker image.
set -euo pipefail

cd "$(dirname "$0")/.."

echo "Building RealDetect production image..."
docker compose build app
echo ""
echo "Done. Run with:  docker compose up app"
