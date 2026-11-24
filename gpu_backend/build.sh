#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

NVCC="${NVCC:-nvcc}"
CUDA_ARCH="${CUDA_ARCH:-sm_89}"
CXX="${CXX:-g++}"
CUDA_LIB_DIR="${CUDA_LIB_DIR:-/usr/local/cuda/lib64}"

"${NVCC}" -O3 -std=c++17 -arch="${CUDA_ARCH}" --shared \
  -Xcompiler "-O3 -DNDEBUG -fPIC" \
  "${SCRIPT_DIR}/c2_gpu_engine.cu" \
  -o "${SCRIPT_DIR}/libibamr_c2_gpu.so" \
  -L"${CUDA_LIB_DIR}" -lcudart -ldl \
  -Xlinker -rpath,"${CUDA_LIB_DIR}"

"${CXX}" -O3 -std=c++17 -DNDEBUG \
  -I"${SCRIPT_DIR}" \
  "${SCRIPT_DIR}/test_c2_gpu_lib.cpp" \
  -o "${SCRIPT_DIR}/test_c2_gpu_lib" \
  -ldl

ls -lh "${SCRIPT_DIR}/libibamr_c2_gpu.so" "${SCRIPT_DIR}/test_c2_gpu_lib"
