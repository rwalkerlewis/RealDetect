#!/usr/bin/env bash
# Run the test suite inside a Docker container.
set -euo pipefail

cd "$(dirname "$0")/.."

echo "Running RealDetect tests in Docker..."
echo ""

docker compose run --rm test
