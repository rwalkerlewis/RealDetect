#!/bin/bash
set -e

echo "========================================"
echo " RealDetect Development Environment"
echo "========================================"
echo ""

# Build C++ application
echo "[1/3] Building C++ application..."
mkdir -p /app/build
cd /app/build
cmake /app/src_mount -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON 2>&1
make -j$(nproc) 2>&1
echo "  -> Build complete."
echo ""

# Run tests
echo "[2/3] Running tests..."
./realdetect_tests 2>&1 || echo "  -> Some tests failed (continuing anyway)."
echo ""

# Start web GUI
echo "[3/3] Starting web GUI on port 8080..."
cd /app
exec python3 /app/src_mount/gui/server.py \
    --stations /app/src_mount/config/stations.txt \
    --config /app/src_mount/config/realdetect.conf \
    --db /app/data/realdetect_catalog.db \
    --build-dir /app/build \
    --host 0.0.0.0 \
    --port 8080
