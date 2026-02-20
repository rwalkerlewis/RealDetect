#!/usr/bin/env bash
# Build multi-architecture Docker images using docker buildx.
#
# Supported platforms: linux/amd64  linux/arm64
#
# Usage:
#   ./docker/cross.sh                        # build both arches, load to local daemon
#   ./docker/cross.sh --push                 # build & push to registry
#   ./docker/cross.sh --platform linux/arm64 # single arch
#   IMAGE_TAG=myrepo/realdetect:v1 ./docker/cross.sh --push
#
set -euo pipefail

cd "$(dirname "$0")/.."

# ── Config ────────────────────────────────────────────────────────────────────
IMAGE_TAG="${IMAGE_TAG:-realdetect:latest}"
PLATFORMS="${PLATFORMS:-linux/amd64,linux/arm64}"
BUILDER_NAME="realdetect-builder"
DOCKERFILE="Dockerfile"
PUSH=false
LOAD=false

# ── Parse args ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --push)              PUSH=true ;;
    --load)              LOAD=true ;;
    --platform)          PLATFORMS="$2"; shift ;;
    --platform=*)        PLATFORMS="${1#*=}" ;;
    --tag)               IMAGE_TAG="$2";    shift ;;
    --tag=*)             IMAGE_TAG="${1#*=}" ;;
    --dev)               DOCKERFILE="Dockerfile.dev" ;;
    -h|--help)
      echo "Usage: $0 [--push] [--load] [--platform P] [--tag T] [--dev]"
      exit 0
      ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
  shift
done

# --load only works for a single platform
if $LOAD && [[ "$PLATFORMS" == *","* ]]; then
  echo "⚠  --load is only supported for a single platform."
  echo "   Defaulting to linux/amd64.  Use --platform linux/arm64 to override."
  PLATFORMS="linux/amd64"
fi

# ── Ensure buildx builder with multi-platform support ─────────────────────────
if ! docker buildx inspect "${BUILDER_NAME}" &>/dev/null; then
  echo "Creating buildx builder '${BUILDER_NAME}'…"
  docker buildx create --name "${BUILDER_NAME}" --driver docker-container --bootstrap
fi
docker buildx use "${BUILDER_NAME}"

echo "========================================================"
echo "  RealDetect cross-platform build"
echo "  Image     : ${IMAGE_TAG}"
echo "  Platforms : ${PLATFORMS}"
echo "  Dockerfile: ${DOCKERFILE}"
echo "========================================================"
echo ""

# ── Build ─────────────────────────────────────────────────────────────────────
BUILD_ARGS=(
  buildx build
  --builder "${BUILDER_NAME}"
  --platform "${PLATFORMS}"
  --file "${DOCKERFILE}"
  --tag "${IMAGE_TAG}"
)

$PUSH && BUILD_ARGS+=( --push )
$LOAD && BUILD_ARGS+=( --load )

# If neither --push nor --load, cache-only build (CI smoke test)
if ! $PUSH && ! $LOAD; then
  echo "ℹ  Neither --push nor --load specified; building to cache only."
  BUILD_ARGS+=( --output type=cacheonly )
fi

docker "${BUILD_ARGS[@]}" .

echo ""
echo "Done."
if $LOAD; then
  echo "Image loaded: ${IMAGE_TAG}"
  echo "Run with:  docker run --rm -p 8080:8080 ${IMAGE_TAG}"
fi
if $PUSH; then
  echo "Image pushed: ${IMAGE_TAG}"
fi
