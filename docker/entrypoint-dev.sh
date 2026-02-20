#!/bin/bash
# RealDetect development/test entrypoint
# Controlled by environment variables set in docker-compose.yml:
#   BUILD_TYPE : Debug | Release            (default: Debug)
#   SANITIZER  : asan | ubsan | tsan | ""  (default: "")
#   RUN_MODE   : dev | test | examples      (default: dev)
#   EXAMPLE    : ridgecrest | dprk | ensemble_ridgecrest | ensemble_dprk
set -e

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()    { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERR]${NC}   $*"; }

# ── Defaults ──────────────────────────────────────────────────────────────────
BUILD_TYPE="${BUILD_TYPE:-Debug}"
SANITIZER="${SANITIZER:-}"
RUN_MODE="${RUN_MODE:-dev}"
EXAMPLE="${EXAMPLE:-ridgecrest}"
SRC="/app/src_mount"
BUILD="/app/build"

echo "========================================================"
echo "  RealDetect Container"
echo "  Mode      : ${RUN_MODE}"
echo "  Build type: ${BUILD_TYPE}"
echo "  Sanitizer : ${SANITIZER:-none}"
echo "========================================================"
echo ""

# ── CMake extra flags based on sanitizer ─────────────────────────────────────
SANITIZER_FLAGS=""
case "${SANITIZER}" in
  asan)
    info "AddressSanitizer enabled"
    SANITIZER_FLAGS="-DCMAKE_CXX_FLAGS='-fsanitize=address,leak -fno-omit-frame-pointer' \
                     -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,leak' \
                     -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
    ;;
  ubsan)
    info "UndefinedBehaviorSanitizer enabled"
    SANITIZER_FLAGS="-DCMAKE_CXX_FLAGS='-fsanitize=undefined -fno-omit-frame-pointer' \
                     -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=undefined' \
                     -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
    ;;
  tsan)
    info "ThreadSanitizer enabled"
    SANITIZER_FLAGS="-DCMAKE_CXX_FLAGS='-fsanitize=thread -fno-omit-frame-pointer' \
                     -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=thread' \
                     -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
    ;;
  "")
    : # no sanitizer
    ;;
  *)
    warn "Unknown SANITIZER='${SANITIZER}' – skipping sanitizer flags"
    ;;
esac

# ── Build ─────────────────────────────────────────────────────────────────────
build_project() {
    info "Building (${BUILD_TYPE}${SANITIZER:+, ${SANITIZER}})…"
    mkdir -p "${BUILD}" && cd "${BUILD}"

    # eval is needed to expand SANITIZER_FLAGS as individual cmake args
    eval cmake "${SRC}" \
        -G Ninja \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DBUILD_TESTS=ON \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        ${SANITIZER_FLAGS} \
        2>&1

    # Symlink compile_commands.json to source root for IDE/clangd
    ln -sf "${BUILD}/compile_commands.json" "${SRC}/compile_commands.json" 2>/dev/null || true

    ninja -j"$(nproc)" 2>&1
    info "Build complete."
    echo ""
}

# ── Run tests ─────────────────────────────────────────────────────────────────
run_tests() {
    info "Running test suite…"
    cd "${BUILD}"
    ./realdetect_tests 2>&1
    STATUS=$?
    if [ "${STATUS}" -eq 0 ]; then
        info "All tests passed."
    else
        warn "Some tests failed (exit code ${STATUS})."
    fi
    return "${STATUS}"
}

# ── Run an example ────────────────────────────────────────────────────────────
run_example() {
    local ex="${EXAMPLE}"
    info "Running example: ${ex}"
    cd "${BUILD}"

    case "${ex}" in
      ridgecrest)
        ./ridgecrest_example \
            --stations "${SRC}/data/ridgecrest/stations.txt" \
            --waveforms "${SRC}/data/ridgecrest/waveforms.mseed" \
            --output /app/output/ridgecrest \
            --config  "${SRC}/config/realdetect.conf"
        ;;
      dprk)
        ./dprk_example \
            --stations "${SRC}/data/dprk/stations.txt" \
            --waveforms "${SRC}/data/dprk/waveforms.mseed" \
            --output /app/output/dprk \
            --config  "${SRC}/config/realdetect.conf"
        ;;
      ensemble_ridgecrest)
        ./ensemble_ridgecrest \
            --stations "${SRC}/data/ridgecrest/stations.txt" \
            --waveforms "${SRC}/data/ridgecrest/waveforms.mseed" \
            --output /app/output/ensemble_ridgecrest \
            --config  "${SRC}/config/realdetect.conf"
        ;;
      ensemble_dprk)
        ./ensemble_dprk \
            --stations "${SRC}/data/dprk/stations.txt" \
            --waveforms "${SRC}/data/dprk/waveforms.mseed" \
            --output /app/output/ensemble_dprk \
            --config  "${SRC}/config/realdetect.conf"
        ;;
      *)
        error "Unknown EXAMPLE='${ex}'. Options: ridgecrest | dprk | ensemble_ridgecrest | ensemble_dprk"
        exit 1
        ;;
    esac
}

# ── Main dispatch ─────────────────────────────────────────────────────────────
case "${RUN_MODE}" in
  dev)
    build_project
    run_tests || true   # don't abort dev server on test failure
    info "Starting web GUI on port 8080…"
    exec python3 "${SRC}/gui/server.py" \
        --stations  "${SRC}/config/stations.txt" \
        --config    "${SRC}/config/realdetect.conf" \
        --db        /app/data/realdetect_catalog.db \
        --build-dir "${BUILD}" \
        --host 0.0.0.0 \
        --port 8080
    ;;

  test)
    build_project
    run_tests
    ;;

  examples)
    build_project
    run_example
    ;;

  *)
    error "Unknown RUN_MODE='${RUN_MODE}'. Options: dev | test | examples"
    exit 1
    ;;
esac
