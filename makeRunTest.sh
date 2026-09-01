#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-${PROJECT_DIR}/build-release}"

if command -v nproc >/dev/null 2>&1; then
    DEFAULT_BUILD_JOBS="$(nproc)"
else
    DEFAULT_BUILD_JOBS=1
fi
BUILD_JOBS="${BUILD_JOBS:-${DEFAULT_BUILD_JOBS}}"

cmake \
    -S "${PROJECT_DIR}" \
    -B "${BUILD_DIR}" \
    -DCMAKE_BUILD_TYPE=Release

cmake \
    --build "${BUILD_DIR}" \
    --config Release \
    --parallel "${BUILD_JOBS}"

VWSS_EXECUTABLE="${BUILD_DIR}/test"
if [[ ! -x "${VWSS_EXECUTABLE}" ]]; then
    VWSS_EXECUTABLE="${BUILD_DIR}/Release/test"
fi

if [[ ! -x "${VWSS_EXECUTABLE}" ]]; then
    echo "Release executable was not produced in ${BUILD_DIR}." >&2
    exit 1
fi

echo "Running VWSS release benchmark: ${VWSS_EXECUTABLE}"
exec "${VWSS_EXECUTABLE}" "$@"
