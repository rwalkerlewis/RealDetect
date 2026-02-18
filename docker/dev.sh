#!/usr/bin/env bash
# Quick-start script for the development environment.
# Builds the dev container and starts the web GUI at http://localhost:8080
set -euo pipefail

cd "$(dirname "$0")/.."

echo "Starting RealDetect development environment..."
echo "  Web GUI will be available at http://localhost:8080"
echo ""

docker compose up --build dev
