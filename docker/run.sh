#!/usr/bin/env bash
# Run a RealDetect example inside the dev container.
#
# Usage:
#   ./docker/run.sh ridgecrest
#   ./docker/run.sh dprk
#   ./docker/run.sh ensemble_ridgecrest
#   ./docker/run.sh ensemble_dprk
#
# The binary is rebuilt if the build volume is stale (cmake + ninja).
# Output is written to ./output/<example>/ on the host.
#
set -euo pipefail

cd "$(dirname "$0")/.."

EXAMPLE="${1:-ridgecrest}"
VALID=("ridgecrest" "dprk" "ensemble_ridgecrest" "ensemble_dprk")

# Validate
OK=false
for v in "${VALID[@]}"; do [ "$EXAMPLE" = "$v" ] && OK=true; done
if ! $OK; then
  echo "Unknown example: '${EXAMPLE}'"
  echo "Valid options: ${VALID[*]}"
  exit 1
fi

echo "Running example: ${EXAMPLE}"
echo ""

EXAMPLE="${EXAMPLE}" docker compose run --rm examples
