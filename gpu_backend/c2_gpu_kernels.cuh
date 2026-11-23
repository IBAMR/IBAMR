// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef IBAMR_C2_GPU_KERNELS_CUH
#define IBAMR_C2_GPU_KERNELS_CUH

#include <cuda_runtime.h>

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace c2_gpu
{

// Runtime grid/fluid constants uploaded via cudaMemcpyToSymbol in c2_gpu_create
__constant__ int c_nx0;
__constant__ int c_ny0;
__constant__ int c_nz0;
__constant__ size_t c_n_l0_cells;
__constant__ double c_dx0;
__constant__ double c_dx1;
__constant__ double c_rho0;
__constant__ double c_mu0;
__constant__ double c_d_vel_const;

// Indexers for patch arrays
__host__ __device__ inline int cidx_g(int i, int j, int k, int nx, int ny, int g)
{
    return (i + g) + (nx + 2 * g) * ((j + g) + (ny + 2 * g) * (k + g));
}
__host__ __device__ inline int cidx_g0(int i, int j, int k, int nx, int ny) { return i + nx * (j + ny * k); }
__host__ __device__ inline int cidx_g1(int i, int j, int k, int nx, int ny) { return cidx_g(i, j, k, nx, ny, 1); }

__host__ __device__ inline int sidx_g(int axis, int i, int j, int k, int nx, int ny, int g)
{
    const int sx = nx + (axis == 0 ? 1 : 0) + 2 * g;
    const int sy = ny + (axis == 1 ? 1 : 0) + 2 * g;
    return (i + g) + sx * ((j + g) + sy * (k + g));
}
__host__ __device__ inline int sidx_g0(int axis, int i, int j, int k, int nx, int ny)
{
    return sidx_g(axis, i, j, k, nx, ny, 0);
}
__host__ __device__ inline int sidx_g1(int axis, int i, int j, int k, int nx, int ny)
{
    return sidx_g(axis, i, j, k, nx, ny, 1);
}
__host__ __device__ inline int sidx_g3(int axis, int i, int j, int k, int nx, int ny)
{
    return sidx_g(axis, i, j, k, nx, ny, 3);
}

struct Box3i
{
    int lo[3];
    int hi[3];
    int nx, ny, nz;
    size_t off_c_g0;
    size_t off_c_g1;
    size_t off_s_g0[3];
    size_t off_s_g1[3];
    size_t off_s_g2[3];
    size_t off_s_g3[3];
};

__device__ inline int lookup_l1_patch(int gi, int gj, int gk, const int16_t* l1_map)
{
    if (gi < 0 || gi >= 2 * c_nx0 || gj < 0 || gj >= 2 * c_ny0 || gk < 0 || gk >= 2 * c_nz0) return -1;
    int i0 = gi >> 1, j0 = gj >> 1, k0 = gk >> 1;
    return static_cast<int>(l1_map[i0 + c_nx0 * (j0 + c_ny0 * k0)]);
}

// Return owning L1 patch and local indices for an L1 side index (gi, gj, gk) along `axis`.
// Returns -1 if not on L1 (either outside domain or in pure coarse L0 region).
__device__ inline bool lookup_l1_side_owner(
    int axis, int gi, int gj, int gk,
    const int16_t* l1_map, const Box3i* boxes,
    int& out_pat, int& out_li, int& out_lj, int& out_lk)
{
    if (axis == 0)
    {
        if (gi < 0 || gi > 2 * c_nx0 || gj < 0 || gj >= 2 * c_ny0 || gk < 0 || gk >= 2 * c_nz0) return false;
        int q = (gi < 2 * c_nx0) ? lookup_l1_patch(gi, gj, gk, l1_map) : -1;
        if (q < 0 && gi > 0) q = lookup_l1_patch(gi - 1, gj, gk, l1_map);
        if (q < 0) return false;
        const Box3i& qb = boxes[q];
        out_pat = q;
        out_li = gi - qb.lo[0];
        out_lj = gj - qb.lo[1];
        out_lk = gk - qb.lo[2];
        return true;
    }
    if (axis == 1)
    {
        if (gi < 0 || gi >= 2 * c_nx0 || gj < 0 || gj > 2 * c_ny0 || gk < 0 || gk >= 2 * c_nz0) return false;
        int q = (gj < 2 * c_ny0) ? lookup_l1_patch(gi, gj, gk, l1_map) : -1;
        if (q < 0 && gj > 0) q = lookup_l1_patch(gi, gj - 1, gk, l1_map);
        if (q < 0) return false;
        const Box3i& qb = boxes[q];
        out_pat = q;
        out_li = gi - qb.lo[0];
        out_lj = gj - qb.lo[1];
        out_lk = gk - qb.lo[2];
        return true;
    }
    if (gi < 0 || gi >= 2 * c_nx0 || gj < 0 || gj >= 2 * c_ny0 || gk < 0 || gk > 2 * c_nz0) return false;
    int q = (gk < 2 * c_nz0) ? lookup_l1_patch(gi, gj, gk, l1_map) : -1;
    if (q < 0 && gk > 0) q = lookup_l1_patch(gi, gj, gk - 1, l1_map);
    if (q < 0) return false;
    const Box3i& qb = boxes[q];
    out_pat = q;
    out_li = gi - qb.lo[0];
    out_lj = gj - qb.lo[1];
    out_lk = gk - qb.lo[2];
    return true;
}

// ============================================================================
// Vector & Reduction Kernels
// ============================================================================
__global__ void wgt_dot1_kernel(const double* x, const double* w, size_t n, double* out)
{
    __shared__ double sh[256];
    double s = 0.0;
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
         i += static_cast<size_t>(gridDim.x) * blockDim.x)
        s += w[i] * x[i] * x[i];
    sh[threadIdx.x] = s;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1)
    {
        if (threadIdx.x < st) sh[threadIdx.x] += sh[threadIdx.x + st];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = sh[0];
}

__global__ void wgt_dot_ab_kernel(
    const double* a0, const double* a1, const double* a2, const double* ap,
    const double* b0, const double* b1, const double* b2, const double* bp,
    const double* w0, const double* w1, const double* w2, const double* wp,
    size_t n0, size_t n1, size_t n2, size_t np, double* out)
{
    __shared__ double sh[256];
    double s = 0.0;
    size_t tid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;
    size_t tot = n0 + n1 + n2 + np;
    for (size_t i = tid; i < tot; i += stride)
    {
        if (i < n0) s += w0[i] * a0[i] * b0[i];
        else if (i < n0 + n1) { size_t j = i - n0; s += w1[j] * a1[j] * b1[j]; }
        else if (i < n0 + n1 + n2) { size_t j = i - (n0 + n1); s += w2[j] * a2[j] * b2[j]; }
        else { size_t j = i - (n0 + n1 + n2); s += wp[j] * ap[j] * bp[j]; }
    }
    sh[threadIdx.x] = s;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1)
    {
        if (threadIdx.x < st) sh[threadIdx.x] += sh[threadIdx.x + st];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = sh[0];
}

__global__ void max_abs_diff_kernel(const double* a, const double* b, size_t n, double* out_sq, double* out_max)
{
    __shared__ double sh_sq[256];
    __shared__ double sh_mx[256];
    double sq = 0.0, mx = 0.0;
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
         i += static_cast<size_t>(gridDim.x) * blockDim.x)
    {
        double d = fabs(a[i] - b[i]);
        sq += d * d;
        if (!(d <= mx)) mx = d; // catches NaN too
    }
    sh_sq[threadIdx.x] = sq;
    sh_mx[threadIdx.x] = mx;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1)
    {
        if (threadIdx.x < st)
        {
            sh_sq[threadIdx.x] += sh_sq[threadIdx.x + st];
            double bval = sh_mx[threadIdx.x + st];
            if (!(bval <= sh_mx[threadIdx.x])) sh_mx[threadIdx.x] = bval;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        out_sq[blockIdx.x] = sh_sq[0];
        out_max[blockIdx.x] = sh_mx[0];
    }
}

__global__ void max_abs_1_kernel(const double* a, size_t n, double* out_max)
{
    __shared__ double sh_mx[256];
    double mx = 0.0;
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
         i += static_cast<size_t>(gridDim.x) * blockDim.x)
    {
        double v = fabs(a[i]);
        if (!(v <= mx)) mx = v;
    }
    sh_mx[threadIdx.x] = mx;
    __syncthreads();
    for (int st = blockDim.x >> 1; st > 0; st >>= 1)
    {
        if (threadIdx.x < st)
        {
            double bval = sh_mx[threadIdx.x + st];
            if (!(bval <= sh_mx[threadIdx.x])) sh_mx[threadIdx.x] = bval;
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) out_max[blockIdx.x] = sh_mx[0];
}

template <typename T>
__global__ void fill_const_kernel(T* dst, T val, size_t n)
{
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
         i += static_cast<size_t>(gridDim.x) * blockDim.x)
        dst[i] = val;
}

__global__ void axpy4_kernel(
    double alpha,
    const double* x0, const double* x1, const double* x2, const double* xp,
    double* y0, double* y1, double* y2, double* yp,
    size_t n0, size_t n1, size_t n2, size_t np)
{
    size_t tid = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t stride = static_cast<size_t>(gridDim.x) * blockDim.x;
    size_t tot = n0 + n1 + n2 + np;
    for (size_t i = tid; i < tot; i += stride)
    {
        if (i < n0) y0[i] += alpha * x0[i];
        else if (i < n0 + n1) { size_t j = i - n0; y1[j] += alpha * x1[j]; }
        else if (i < n0 + n1 + n2) { size_t j = i - (n0 + n1); y2[j] += alpha * x2[j]; }
        else { size_t j = i - (n0 + n1 + n2); yp[j] += alpha * xp[j]; }
    }
}

__global__ void linear_sum_side_kernel(
    double* dst, double w0, const double* x0, double w1, const double* x1, size_t n)
{
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
         i += static_cast<size_t>(gridDim.x) * blockDim.x)
        dst[i] = w0 * x0[i] + w1 * x1[i];
}

// Host helpers
template <typename RealT>
__global__ void pack_g0_to_g_side_l0_kernel(const double* src, RealT* dst, int axis, int nx, int ny, int nz, int g)
{
    int nxa = nx + (axis == 0), nya = ny + (axis == 1), nza = nz + (axis == 2);
    int tot = nxa * nya * nza;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
        dst[sidx_g(axis, i, j, k, nx, ny, g)] = static_cast<RealT>(src[idx]);
    }
}

template <typename RealT>
__global__ void pack_g0_to_g_side_l1_kernel(
    const double* src, RealT* dst, const Box3i* boxes, int n_l1, int axis, int g)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        int tot = nxa * nya * nza;
        size_t o0 = b.off_s_g0[axis];
        size_t og = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
            dst[og + sidx_g(axis, i, j, k, b.nx, b.ny, g)] = static_cast<RealT>(src[o0 + idx]);
        }
    }
}

__global__ void unpack_g_to_g0_side_l0_kernel(const double* src, double* dst, int axis, int nx, int ny, int nz, int g)
{
    int nxa = nx + (axis == 0), nya = ny + (axis == 1), nza = nz + (axis == 2);
    int tot = nxa * nya * nza;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
        dst[idx] = src[sidx_g(axis, i, j, k, nx, ny, g)];
    }
}

__global__ void unpack_g_to_g0_side_l1_kernel(
    const double* src, double* dst, const Box3i* boxes, int n_l1, int axis, int g)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        int tot = nxa * nya * nza;
        size_t o0 = b.off_s_g0[axis];
        size_t og = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
            dst[o0 + idx] = src[og + sidx_g(axis, i, j, k, b.nx, b.ny, g)];
        }
    }
}

template <typename RealT>
__global__ void pack_g0_to_g1_cell_l0_kernel(const double* src, RealT* dst, int nx, int ny, int nz)
{
    int tot = nx * ny * nz;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, t = idx / nx, j = t % ny, k = t / ny;
        dst[cidx_g1(i, j, k, nx, ny)] = static_cast<RealT>(src[idx]);
    }
}

template <typename RealT>
__global__ void pack_g0_to_g1_cell_l1_kernel(const double* src, RealT* dst, const Box3i* boxes, int n_l1)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        int tot = b.nx * b.ny * b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, t = idx / b.nx, j = t % b.ny, k = t / b.ny;
            dst[b.off_c_g1 + cidx_g1(i, j, k, b.nx, b.ny)] = static_cast<RealT>(src[b.off_c_g0 + idx]);
        }
    }
}

__global__ void pack_g1_to_g0_cell_l0_kernel(const double* src, double* dst, int nx, int ny, int nz)
{
    int tot = nx * ny * nz;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, t = idx / nx, j = t % ny, k = t / ny;
        dst[idx] = src[cidx_g1(i, j, k, nx, ny)];
    }
}

__global__ void pack_g1_to_g0_cell_l1_kernel(const double* src, double* dst, const Box3i* boxes, int n_l1)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        int tot = b.nx * b.ny * b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, t = idx / b.nx, j = t % b.ny, k = t / b.ny;
            dst[b.off_c_g0 + idx] = src[b.off_c_g1 + cidx_g1(i, j, k, b.nx, b.ny)];
        }
    }
}

// Copy L1 interior -> overlapping L1 patch-interior ghosts (pure fine-fine intra-level exchange) for g=1 or g=3
// Also synchronizes shared fine patch-boundary faces (note: SideData keeps lower-patch canonical copy on interior, or both already equal).
template <typename RealT>
__global__ void sync_l1_intralevel_side_kernel(
    RealT* u_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis, int g)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const int sx = b.nx + (axis == 0) + 2 * g;
        const int sy = b.ny + (axis == 1) + 2 * g;
        const int sz = b.nz + (axis == 2) + 2 * g;
        const int tot = sx * sy * sz;
        const size_t my_off = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % sx) - g;
            int t = idx / sx;
            int lj = (t % sy) - g;
            int lk = (t / sy) - g;
            // Skip strictly interior points owned by this patch; wait: also copy shared boundary face from lower owner if needed?
            // All physical interior points on patch boundaries: copy from canonical owner so patches always agree bit-for-bit.
            int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
            bool in_box = (li >= 0 && li < nxa && lj >= 0 && lj < nya && lk >= 0 && lk < nza);
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            int q = -1, ql = 0, qj = 0, qk = 0;
            if (lookup_l1_side_owner(axis, gi, gj, gk, l1_map, boxes, q, ql, qj, qk))
            {
                if (q != pat && (!in_box || q < pat))
                {
                    const Box3i& qb = boxes[q];
                    size_t q_off = (g == 1) ? qb.off_s_g1[axis] : qb.off_s_g3[axis];
                    u_l1[my_off + sidx_g(axis, li, lj, lk, b.nx, b.ny, g)] =
                        u_l1[q_off + sidx_g(axis, ql, qj, qk, qb.nx, qb.ny, g)];
                }
            }
        }
    }
}
template <typename RealT>
__global__ void sync_l1_intralevel_cell_kernel(
    RealT* c_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const int sx = b.nx + 2, sy = b.ny + 2, sz = b.nz + 2;
        const int tot = sx * sy * sz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % sx) - 1;
            int t = idx / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;
            if (li >= 0 && li < b.nx && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz) continue;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            int q = lookup_l1_patch(gi, gj, gk, l1_map);
            if (q >= 0 && q != pat)
            {
                const Box3i& qb = boxes[q];
                c_l1[b.off_c_g1 + cidx_g1(li, lj, lk, b.nx, b.ny)] =
                    c_l1[qb.off_c_g1 + cidx_g1(gi - qb.lo[0], gj - qb.lo[1], gk - qb.lo[2], qb.nx, qb.ny)];
            }
        }
    }
}


// Wait: in SSPRK3, does L1->L1 intra-level ghost copy only copy where ANOTHER L1 patch owns the cell inside `q`'s interior?
// Note look closely: if (gi,gj,gk) lies on the shared face between pat and q (e.g. li=b.nx, ql=0), copying from q < pat guarantees unique value;
// what if a ghost cell (li < 0) sits OUTSIDE pat's interior and is INSIDE q's interior? It gets q's interior value!

// Coarsen L1 -> L0 (CONSERVATIVE_COARSEN = unweighted 4-face average, matching CartSideDoubleWeightedAverage and CartSideDoubleCubicCoarsen for ratio=2<4)
template <typename RealT>
__global__ void restrict_l1_to_l0_side_kernel(
    const RealT* u_l1, RealT* u_l0, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis, int g)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const size_t o1 = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        const int cnx = b.nx >> 1, cny = b.ny >> 1, cnz = b.nz >> 1;
        const int cnxa = cnx + (axis == 0), cnya = cny + (axis == 1), cnza = cnz + (axis == 2);
        const int tot = cnxa * cnya * cnza;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int ci = idx % cnxa, t = idx / cnxa, cj = t % cnya, ck = t / cnya;
            int ic = (b.lo[0] >> 1) + ci, jc = (b.lo[1] >> 1) + cj, kc = (b.lo[2] >> 1) + ck;
            // Avoid double write on upper patch face if covered by next patch (shared face between patches has same values after intra sync)
            int li = ci * 2, lj = cj * 2, lk = ck * 2;
            double sum = 0.0;
            if (axis == 0)
            {
                if (ci == cnx && ic < c_nx0 && lookup_l1_patch(ic * 2, jc * 2, kc * 2, l1_map) >= 0) continue;
                for (int dk = 0; dk < 2; ++dk)
                    for (int dj = 0; dj < 2; ++dj)
                        sum += static_cast<double>(u_l1[o1 + sidx_g(0, li, lj + dj, lk + dk, b.nx, b.ny, g)]);
            }
            else if (axis == 1)
            {
                if (cj == cny && jc < c_ny0 && lookup_l1_patch(ic * 2, jc * 2, kc * 2, l1_map) >= 0) continue;
                for (int dk = 0; dk < 2; ++dk)
                    for (int di = 0; di < 2; ++di)
                        sum += static_cast<double>(u_l1[o1 + sidx_g(1, li + di, lj, lk + dk, b.nx, b.ny, g)]);
            }
            else
            {
                if (ck == cnz && kc < c_nz0 && lookup_l1_patch(ic * 2, jc * 2, kc * 2, l1_map) >= 0) continue;
                for (int dj = 0; dj < 2; ++dj)
                    for (int di = 0; di < 2; ++di)
                        sum += static_cast<double>(u_l1[o1 + sidx_g(2, li + di, lj + dj, lk, b.nx, b.ny, g)]);
            }
            u_l0[sidx_g(axis, ic, jc, kc, c_nx0, c_ny0, g)] = static_cast<RealT>(0.25 * sum);
        }
    }
}

template <typename RealT>
__global__ void restrict_l1_to_l0_cell_kernel(
    const RealT* c_l1, RealT* c_l0, const Box3i* boxes, int n_l1)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const int cnx = b.nx >> 1, cny = b.ny >> 1, cnz = b.nz >> 1;
        const int tot = cnx * cny * cnz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int ci = idx % cnx, t = idx / cnx, cj = t % cny, ck = t / cny;
            int ic = (b.lo[0] >> 1) + ci, jc = (b.lo[1] >> 1) + cj, kc = (b.lo[2] >> 1) + ck;
            int fi = ci * 2, fj = cj * 2, fk = ck * 2;
            double s = 0.0;
            for (int dk = 0; dk < 2; ++dk)
                for (int dj = 0; dj < 2; ++dj)
                    for (int di = 0; di < 2; ++di)
                        s += static_cast<double>(c_l1[b.off_c_g1 + cidx_g1(fi + di, fj + dj, fk + dk, b.nx, b.ny)]);
            c_l0[cidx_g1(ic, jc, kc, c_nx0, c_ny0)] = static_cast<RealT>(0.125 * s);
        }
    }
}

// Also support g=0 restriction directly into g=0 L0 for coefficient syncs & diagnostics
template <typename RealT>
__global__ void restrict_g0_l1_to_l0_side_kernel(
    const RealT* u_l1, RealT* u_l0, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const size_t o1 = b.off_s_g0[axis];
        const int cnx = b.nx >> 1, cny = b.ny >> 1, cnz = b.nz >> 1;
        const int cnxa = cnx + (axis == 0), cnya = cny + (axis == 1), cnza = cnz + (axis == 2);
        const int tot = cnxa * cnya * cnza;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int ci = idx % cnxa, t = idx / cnxa, cj = t % cnya, ck = t / cnya;
            int ic = (b.lo[0] >> 1) + ci, jc = (b.lo[1] >> 1) + cj, kc = (b.lo[2] >> 1) + ck;
            int li = ci * 2, lj = cj * 2, lk = ck * 2;
            double sum = 0.0;
            if (axis == 0)
            {
                if (ci == cnx && ic < c_nx0 && lookup_l1_patch(ic * 2, jc * 2, kc * 2, l1_map) >= 0) continue;
                for (int dk = 0; dk < 2; ++dk)
                    for (int dj = 0; dj < 2; ++dj)
                        sum += static_cast<double>(u_l1[o1 + sidx_g0(0, li, lj + dj, lk + dk, b.nx, b.ny)]);
            }
            else if (axis == 1)
            {
                if (cj == cny && jc < c_ny0 && lookup_l1_patch(ic * 2, jc * 2, kc * 2, l1_map) >= 0) continue;
                for (int dk = 0; dk < 2; ++dk)
                    for (int di = 0; di < 2; ++di)
                        sum += static_cast<double>(u_l1[o1 + sidx_g0(1, li + di, lj, lk + dk, b.nx, b.ny)]);
            }
            else
            {
                if (ck == cnz && kc < c_nz0 && lookup_l1_patch(ic * 2, jc * 2, kc * 2, l1_map) >= 0) continue;
                for (int dj = 0; dj < 2; ++dj)
                    for (int di = 0; di < 2; ++di)
                        sum += static_cast<double>(u_l1[o1 + sidx_g0(2, li + di, lj + dj, lk, b.nx, b.ny)]);
            }
            u_l0[sidx_g0(axis, ic, jc, kc, c_nx0, c_ny0)] = static_cast<RealT>(0.25 * sum);
        }
    }
}

// Shared L1 upper patch face average sync on g=0 (matching SideDataSynchronization where adjacent L1 patches meet)
template <typename RealT>
__global__ void sync_l1_shared_faces_g0_kernel(
    RealT* u_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        if (axis == 0)
        {
            int li = b.nx;
            int gi = b.lo[0] + li;
            if (gi >= 2 * c_nx0) continue;
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
            {
                int lj = idx % b.ny, lk = idx / b.ny;
                int q = lookup_l1_patch(gi, b.lo[1] + lj, b.lo[2] + lk, l1_map);
                if (q >= 0 && q != pat)
                {
                    const Box3i& qb = boxes[q];
                    u_l1[b.off_s_g0[0] + sidx_g0(0, li, lj, lk, b.nx, b.ny)] =
                        u_l1[qb.off_s_g0[0] + sidx_g0(0, 0, (b.lo[1] + lj) - qb.lo[1], (b.lo[2] + lk) - qb.lo[2], qb.nx, qb.ny)];
                }
            }
        }
        else if (axis == 1)
        {
            int lj = b.ny;
            int gj = b.lo[1] + lj;
            if (gj >= 2 * c_ny0) continue;
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            {
                int li = idx % b.nx, lk = idx / b.nx;
                int q = lookup_l1_patch(b.lo[0] + li, gj, b.lo[2] + lk, l1_map);
                if (q >= 0 && q != pat)
                {
                    const Box3i& qb = boxes[q];
                    u_l1[b.off_s_g0[1] + sidx_g0(1, li, lj, lk, b.nx, b.ny)] =
                        u_l1[qb.off_s_g0[1] + sidx_g0(1, (b.lo[0] + li) - qb.lo[0], 0, (b.lo[2] + lk) - qb.lo[2], qb.nx, qb.ny)];
                }
            }
        }
        else
        {
            int lk = b.nz;
            int gk = b.lo[2] + lk;
            if (gk >= 2 * c_nz0) continue;
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            {
                int li = idx % b.nx, lj = idx / b.nx;
                int q = lookup_l1_patch(b.lo[0] + li, b.lo[1] + lj, gk, l1_map);
                if (q >= 0 && q != pat)
                {
                    const Box3i& qb = boxes[q];
                    u_l1[b.off_s_g0[2] + sidx_g0(2, li, lj, lk, b.nx, b.ny)] =
                        u_l1[qb.off_s_g0[2] + sidx_g0(2, (b.lo[0] + li) - qb.lo[0], (b.lo[1] + lj) - qb.lo[1], 0, qb.nx, qb.ny)];
                }
            }
        }
    }
}


// ============================================================================
// Physical Boundary Extrapolation (CartExtrapPhysBdryOp + CartSideRobinPhysBdryOp + DivFree)
// ============================================================================
// Apply CartExtrapPhysBdryOp along codim 1,2,3 (extrap_mode: 0 = CONSTANT, 1 = LINEAR)
__global__ void cart_extrap_phys_side_kernel(
    double* u, int axis, int nx, int ny, int nz,
    int dom_lo0, int dom_hi0, int dom_lo1, int dom_hi1, int dom_lo2, int dom_hi2,
    int b_lo0, int b_lo1, int b_lo2, int bnx, int bny, int bnz,
    int g, int extrap_mode, size_t u_off)
{
    int snx = bnx + (axis == 0) + 2 * g;
    int sny = bny + (axis == 1) + 2 * g;
    int snz = bnz + (axis == 2) + 2 * g;
    int tot = snx * sny * snz;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int li = (idx % snx) - g;
        int t = idx / snx;
        int lj = (t % sny) - g;
        int lk = (t / sny) - g;
        int nxa = bnx + (axis == 0), nya = bny + (axis == 1), nza = bnz + (axis == 2);
        if (li >= 0 && li < nxa && lj >= 0 && lj < nya && lk >= 0 && lk < nza) continue;
        int gi = b_lo0 + li, gj = b_lo1 + lj, gk = b_lo2 + lk;
        int d_max0 = dom_hi0 + (axis == 0 ? 1 : 0);
        int d_max1 = dom_hi1 + (axis == 1 ? 1 : 0);
        int d_max2 = dom_hi2 + (axis == 2 ? 1 : 0);
        if (gi >= dom_lo0 && gi <= d_max0 && gj >= dom_lo1 && gj <= d_max1 && gk >= dom_lo2 && gk <= d_max2) continue;

        // Project onto nearest physical patch boundary face inside this patch (only touches physical boundaries intersecting this patch!)
        int bi = li, bj = lj, bk = lk;
        int si = 0, sj = 0, sk = 0;
        bool phys_bdry = false;
        if (gi < dom_lo0 && b_lo0 == dom_lo0) { bi = 0; si = 1; phys_bdry = true; }
        else if (gi > d_max0 && (b_lo0 + bnx - 1) == dom_hi0) { bi = nxa - 1; si = -1; phys_bdry = true; }
        if (gj < dom_lo1 && b_lo1 == dom_lo1) { bj = 0; sj = 1; phys_bdry = true; }
        else if (gj > d_max1 && (b_lo1 + bny - 1) == dom_hi1) { bj = nya - 1; sj = -1; phys_bdry = true; }
        if (gk < dom_lo2 && b_lo2 == dom_lo2) { bk = 0; sk = 1; phys_bdry = true; }
        else if (gk > d_max2 && (b_lo2 + bnz - 1) == dom_hi2) { bk = nza - 1; sk = -1; phys_bdry = true; }
        if (!phys_bdry) continue;

        double val = u[u_off + sidx_g(axis, bi, bj, bk, bnx, bny, g)];
        if (extrap_mode == 1)
        {
            if (si != 0)
            {
                double f1 = u[u_off + sidx_g(axis, bi + si, bj, bk, bnx, bny, g)];
                val += (u[u_off + sidx_g(axis, bi, bj, bk, bnx, bny, g)] - f1) * abs(li - bi);
            }
            if (sj != 0)
            {
                double f1 = u[u_off + sidx_g(axis, bi, bj + sj, bk, bnx, bny, g)];
                val += (u[u_off + sidx_g(axis, bi, bj, bk, bnx, bny, g)] - f1) * abs(lj - bj);
            }
            if (sk != 0)
            {
                double f1 = u[u_off + sidx_g(axis, bi, bj, bk + sk, bnx, bny, g)];
                val += (u[u_off + sidx_g(axis, bi, bj, bk, bnx, bny, g)] - f1) * abs(lk - bk);
            }
        }
        u[u_off + sidx_g(axis, li, lj, lk, bnx, bny, g)] = val;
    }
}

__global__ void cart_extrap_phys_side_l1_kernel(
    double* u, int axis, const Box3i* boxes, int n_l1, int g, int extrap_mode)
{
    int dom_hi0 = 2 * c_nx0 - 1, dom_hi1 = 2 * c_ny0 - 1, dom_hi2 = 2 * c_nz0 - 1;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t u_off = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        int snx = nxa + 2 * g, sny = nya + 2 * g, snz = nza + 2 * g;
        int tot = snx * sny * snz;
        int d_max0 = dom_hi0 + (axis == 0 ? 1 : 0);
        int d_max1 = dom_hi1 + (axis == 1 ? 1 : 0);
        int d_max2 = dom_hi2 + (axis == 2 ? 1 : 0);
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % snx) - g;
            int t = idx / snx;
            int lj = (t % sny) - g;
            int lk = (t / sny) - g;
            if (li >= 0 && li < nxa && lj >= 0 && lj < nya && lk >= 0 && lk < nza) continue;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            if (gi >= 0 && gi <= d_max0 && gj >= 0 && gj <= d_max1 && gk >= 0 && gk <= d_max2) continue;
            int bi = li, bj = lj, bk = lk;
            int si = 0, sj = 0, sk = 0;
            bool phys_bdry = false;
            if (gi < 0 && b.lo[0] == 0) { bi = 0; si = 1; phys_bdry = true; }
            else if (gi > d_max0 && b.hi[0] == dom_hi0) { bi = nxa - 1; si = -1; phys_bdry = true; }
            if (gj < 0 && b.lo[1] == 0) { bj = 0; sj = 1; phys_bdry = true; }
            else if (gj > d_max1 && b.hi[1] == dom_hi1) { bj = nya - 1; sj = -1; phys_bdry = true; }
            if (gk < 0 && b.lo[2] == 0) { bk = 0; sk = 1; phys_bdry = true; }
            else if (gk > d_max2 && b.hi[2] == dom_hi2) { bk = nza - 1; sk = -1; phys_bdry = true; }
            if (!phys_bdry) continue;
            double val = u[u_off + sidx_g(axis, bi, bj, bk, b.nx, b.ny, g)];
            if (extrap_mode == 1)
            {
                if (si != 0)
                    val += (u[u_off + sidx_g(axis, bi, bj, bk, b.nx, b.ny, g)] -
                            u[u_off + sidx_g(axis, bi + si, bj, bk, b.nx, b.ny, g)]) * abs(li - bi);
                if (sj != 0)
                    val += (u[u_off + sidx_g(axis, bi, bj, bk, b.nx, b.ny, g)] -
                            u[u_off + sidx_g(axis, bi, bj + sj, bk, b.nx, b.ny, g)]) * abs(lj - bj);
                if (sk != 0)
                    val += (u[u_off + sidx_g(axis, bi, bj, bk, b.nx, b.ny, g)] -
                            u[u_off + sidx_g(axis, bi, bj, bk + sk, b.nx, b.ny, g)]) * abs(lk - bk);
            }
            u[u_off + sidx_g(axis, li, lj, lk, b.nx, b.ny, g)] = val;
        }
    }
}

// CartSideRobinPhysBdryOp (codim 1 + codim 2 + codim 3) for 3D channel + z-midplane symmetry velocity (`is_rho == false`) or homogeneous Neumann `rho` (`is_rho == true`):
// Note:
// - If is_rho: Neumann (a=0, b=1, g=0) on ALL 6 physical faces:
//     normal component: u(i_b +/- k) = u(i_b -/+ k); tangential: u(i_g) = u(i_i)
// - If !is_rho (3D channel velocity with z-midplane symmetry):
//     x_lo, x_up (axis=0): u0 Neumann 0 (until DivFree overwrites normal ghosts); u1, u2 Dirichlet 0 (u_g = -u_i)
//     y_lo, y_up (axis=1): u1 Dirichlet 0 (u_b=0, u_g=-u_i); u0, u2 Dirichlet 0 (u_g=-u_i)
//     z_lo (symmetry plane, loc=4): u2 Dirichlet 0 (u_b=0, u_g=-u_i); u0, u1 Neumann 0 (u_g=+u_i)
//     z_up (loc=5): u2 Dirichlet 0 (u_b=0, u_g=-u_i); u0, u1 Dirichlet 0 (u_g=-u_i)
// Followed on codim2 & codim3 corners by linear extrapolation (`ccrobinphysbdryop2`, `scrobinphysbdryop2`, `scrobinphysbdryop3`)
// Note in CartSideRobinPhysBdryOp::fillGhostCellValuesCodim1Normal/Transverse:
//   codim1 sweeps in physical box trimmed to ghost width g, so `i in [1..g]` mirrors into interior!
template <typename RealT>
__global__ void apply_channel_robin_side_phys_kernel(
    RealT* u, int axis, int b_lo0, int b_lo1, int b_lo2, int bnx, int bny, int bnz,
    int dom_hi0, int dom_hi1, int dom_hi2, int g, bool is_rho, size_t u_off)
{
    // Pass 1: codim 1 physical faces touching this patch
    int nxa = bnx + (axis == 0), nya = bny + (axis == 1), nza = bnz + (axis == 2);
    // x_lo / x_up faces (for indices with lj in [0..nya-1], lk in [0..nza-1])
    if (b_lo0 == 0)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nya * nza; idx += gridDim.x * blockDim.x)
        {
            int j = idx % nya, k = idx / nya;
            for (int r = 1; r <= g; ++r)
            {
                if (axis == 0) u[u_off + sidx_g(0, -r, j, k, bnx, bny, g)] = u[u_off + sidx_g(0, r, j, k, bnx, bny, g)];
                else
                {
                    RealT v = u[u_off + sidx_g(axis, r - 1, j, k, bnx, bny, g)];
                    u[u_off + sidx_g(axis, -r, j, k, bnx, bny, g)] = is_rho ? v : -v;
                }
            }
        }
    }
    if (b_lo0 + bnx - 1 == dom_hi0)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nya * nza; idx += gridDim.x * blockDim.x)
        {
            int j = idx % nya, k = idx / nya;
            for (int r = 1; r <= g; ++r)
            {
                if (axis == 0) u[u_off + sidx_g(0, bnx + r, j, k, bnx, bny, g)] = u[u_off + sidx_g(0, bnx - r, j, k, bnx, bny, g)];
                else
                {
                    RealT v = u[u_off + sidx_g(axis, bnx - r, j, k, bnx, bny, g)];
                    u[u_off + sidx_g(axis, bnx - 1 + r, j, k, bnx, bny, g)] = is_rho ? v : -v;
                }
            }
        }
    }
    // y_lo / y_up faces
    if (b_lo1 == 0)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nza; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, k = idx / nxa;
            if (!is_rho && axis == 1) u[u_off + sidx_g(1, i, 0, k, bnx, bny, g)] = RealT(0.0);
            for (int r = 1; r <= g; ++r)
            {
                if (axis == 1)
                {
                    RealT v = u[u_off + sidx_g(1, i, r, k, bnx, bny, g)];
                    u[u_off + sidx_g(1, i, -r, k, bnx, bny, g)] = is_rho ? v : -v;
                }
                else
                {
                    RealT v = u[u_off + sidx_g(axis, i, r - 1, k, bnx, bny, g)];
                    u[u_off + sidx_g(axis, i, -r, k, bnx, bny, g)] = is_rho ? v : -v;
                }
            }
        }
    }
    if (b_lo1 + bny - 1 == dom_hi1)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nza; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, k = idx / nxa;
            if (!is_rho && axis == 1) u[u_off + sidx_g(1, i, bny, k, bnx, bny, g)] = RealT(0.0);
            for (int r = 1; r <= g; ++r)
            {
                if (axis == 1)
                {
                    RealT v = u[u_off + sidx_g(1, i, bny - r, k, bnx, bny, g)];
                    u[u_off + sidx_g(1, i, bny + r, k, bnx, bny, g)] = is_rho ? v : -v;
                }
                else
                {
                    RealT v = u[u_off + sidx_g(axis, i, bny - r, k, bnx, bny, g)];
                    u[u_off + sidx_g(axis, i, bny - 1 + r, k, bnx, bny, g)] = is_rho ? v : -v;
                }
            }
        }
    }
    // z_lo (symmetry plane) / z_up faces
    if (b_lo2 == 0)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nya; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, j = idx / nxa;
            if (!is_rho && axis == 2) u[u_off + sidx_g(2, i, j, 0, bnx, bny, g)] = RealT(0.0);
            for (int r = 1; r <= g; ++r)
            {
                if (axis == 2)
                {
                    RealT v = u[u_off + sidx_g(2, i, j, r, bnx, bny, g)];
                    u[u_off + sidx_g(2, i, j, -r, bnx, bny, g)] = is_rho ? v : -v;
                }
                else
                {
                    // z-midplane symmetry: tangential u0, u1 mirror with +1 (du/dz=0).
                    RealT v = u[u_off + sidx_g(axis, i, j, r - 1, bnx, bny, g)];
                    u[u_off + sidx_g(axis, i, j, -r, bnx, bny, g)] = v;
                }
            }
        }
    }
    if (b_lo2 + bnz - 1 == dom_hi2)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nya; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, j = idx / nxa;
            if (!is_rho && axis == 2) u[u_off + sidx_g(2, i, j, bnz, bnx, bny, g)] = RealT(0.0);
            for (int r = 1; r <= g; ++r)
            {
                if (axis == 2)
                {
                    RealT v = u[u_off + sidx_g(2, i, j, bnz - r, bnx, bny, g)];
                    u[u_off + sidx_g(2, i, j, bnz + r, bnx, bny, g)] = is_rho ? v : -v;
                }
                else
                {
                    RealT v = u[u_off + sidx_g(axis, i, j, bnz - r, bnx, bny, g)];
                    u[u_off + sidx_g(axis, i, j, bnz - 1 + r, bnx, bny, g)] = is_rho ? v : -v;
                }
            }
        }
    }
}


template <typename RealT>
__global__ void apply_channel_robin_side_phys_l1_kernel(
    RealT* u, const Box3i* boxes, int n_l1, int axis, int g, bool is_rho)
{
    int dom_hi0 = 2 * c_nx0 - 1, dom_hi1 = 2 * c_ny0 - 1, dom_hi2 = 2 * c_nz0 - 1;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t u_off = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        // L1 boxes in the 3D channel configuration do not touch x=0 or x=2*c_nx0-1; check y and z physical faces.
        if (b.lo[1] == 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nza; idx += gridDim.x * blockDim.x)
            {
                int i = idx % nxa, k = idx / nxa;
                if (!is_rho && axis == 1) u[u_off + sidx_g(1, i, 0, k, b.nx, b.ny, g)] = RealT(0.0);
                for (int r = 1; r <= g; ++r)
                {
                    if (axis == 1)
                    {
                        RealT v = u[u_off + sidx_g(1, i, r, k, b.nx, b.ny, g)];
                        u[u_off + sidx_g(1, i, -r, k, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                    else
                    {
                        RealT v = u[u_off + sidx_g(axis, i, r - 1, k, b.nx, b.ny, g)];
                        u[u_off + sidx_g(axis, i, -r, k, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                }
            }
        }
        if (b.hi[1] == dom_hi1)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nza; idx += gridDim.x * blockDim.x)
            {
                int i = idx % nxa, k = idx / nxa;
                if (!is_rho && axis == 1) u[u_off + sidx_g(1, i, b.ny, k, b.nx, b.ny, g)] = RealT(0.0);
                for (int r = 1; r <= g; ++r)
                {
                    if (axis == 1)
                    {
                        RealT v = u[u_off + sidx_g(1, i, b.ny - r, k, b.nx, b.ny, g)];
                        u[u_off + sidx_g(1, i, b.ny + r, k, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                    else
                    {
                        RealT v = u[u_off + sidx_g(axis, i, b.ny - r, k, b.nx, b.ny, g)];
                        u[u_off + sidx_g(axis, i, b.ny - 1 + r, k, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                }
            }
        }
        if (b.lo[2] == 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nya; idx += gridDim.x * blockDim.x)
            {
                int i = idx % nxa, j = idx / nxa;
                if (!is_rho && axis == 2) u[u_off + sidx_g(2, i, j, 0, b.nx, b.ny, g)] = RealT(0.0);
                for (int r = 1; r <= g; ++r)
                {
                    if (axis == 2)
                    {
                        RealT v = u[u_off + sidx_g(2, i, j, r, b.nx, b.ny, g)];
                        u[u_off + sidx_g(2, i, j, -r, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                    else
                    {
                        RealT v = u[u_off + sidx_g(axis, i, j, r - 1, b.nx, b.ny, g)];
                        u[u_off + sidx_g(axis, i, j, -r, b.nx, b.ny, g)] = v;
                    }
                }
            }
        }
        if (b.hi[2] == dom_hi2)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nya; idx += gridDim.x * blockDim.x)
            {
                int i = idx % nxa, j = idx / nxa;
                if (!is_rho && axis == 2) u[u_off + sidx_g(2, i, j, b.nz, b.nx, b.ny, g)] = RealT(0.0);
                for (int r = 1; r <= g; ++r)
                {
                    if (axis == 2)
                    {
                        RealT v = u[u_off + sidx_g(2, i, j, b.nz - r, b.nx, b.ny, g)];
                        u[u_off + sidx_g(2, i, j, b.nz + r, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                    else
                    {
                        RealT v = u[u_off + sidx_g(axis, i, j, b.nz - r, b.nx, b.ny, g)];
                        u[u_off + sidx_g(axis, i, j, b.nz - 1 + r, b.nx, b.ny, g)] = is_rho ? v : -v;
                    }
                }
            }
        }
    }
}


// Physical corner linear extrap (codim2 & codim3 matching CartSideRobinPhysBdryOp)
// Actually CartSideRobinPhysBdryOp uses linear extrap: u(g_a, g_b) = u(bdry_a, g_b) + (u(bdry_a, g_b) - u(bdry_a + shift, g_b)) * dist
// We compute codim2 after codim1, then codim3 after codim2.
__global__ void apply_channel_robin_side_codim23_kernel(
    double* u, int axis, int nx, int ny, int nz,
    int b_lo0, int b_lo1, int b_lo2, int dom_hi0, int dom_hi1, int dom_hi2, int g, size_t u_off)
{
    int nxa = nx + (axis == 0), nya = ny + (axis == 1), nza = nz + (axis == 2);
    int snx = nxa + 2 * g, sny = nya + 2 * g, snz = nza + 2 * g;
    int tot = snx * sny * snz;
    for (int pass = 2; pass <= 3; ++pass)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % snx) - g;
            int t = idx / snx;
            int lj = (t % sny) - g;
            int lk = (t / sny) - g;
            bool ox = (li < 0 && b_lo0 == 0) || (li >= nxa && b_lo0 + nx - 1 == dom_hi0);
            bool oy = (lj < 0 && b_lo1 == 0) || (lj >= nya && b_lo1 + ny - 1 == dom_hi1);
            bool oz = (lk < 0 && b_lo2 == 0) || (lk >= nza && b_lo2 + nz - 1 == dom_hi2);
            int oob = (ox ? 1 : 0) + (oy ? 1 : 0) + (oz ? 1 : 0);
            if (oob != pass) continue;
            if (!((li >= 0 && li < nxa) || ox)) continue;
            if (!((lj >= 0 && lj < nya) || oy)) continue;
            if (!((lk >= 0 && lk < nza) || oz)) continue;
            int bi = li, si = 0;
            if (li < 0) { bi = 0; si = 1; } else if (li >= nxa) { bi = nxa - 1; si = -1; }
            int bj = lj, sj = 0;
            if (lj < 0) { bj = 0; sj = 1; } else if (lj >= nya) { bj = nya - 1; sj = -1; }
            int bk = lk, sk = 0;
            if (lk < 0) { bk = 0; sk = 1; } else if (lk >= nza) { bk = nza - 1; sk = -1; }
            // Wait: if pass==2 & scrobinphysbdryop23 uses average of two extrapolations or linear along both shifts:
            // CartSideRobinPhysBdryOp `scrobinphysbdryop23` (let's use standard IBTK codim2 linear extrap:
            // val = u(bi, bj, bk) + sum_d (u(bi,bj,bk) - u(bi+sd, ...)) or u(g along one)? Notice in SSPRK3 all corners are outside channel walls where u=0 or symmetry!)
            double v0 = u[u_off + sidx_g(axis, bi, bj, bk, nx, ny, g)];
            double val = v0;
            if (si != 0) val += (v0 - u[u_off + sidx_g(axis, bi + si, bj, bk, nx, ny, g)]) * abs(li - bi);
            if (sj != 0) val += (v0 - u[u_off + sidx_g(axis, bi, bj + sj, bk, nx, ny, g)]) * abs(lj - bj);
            if (sk != 0) val += (v0 - u[u_off + sidx_g(axis, bi, bj, bk + sk, nx, ny, g)]) * abs(lk - bk);
            u[u_off + sidx_g(axis, li, lj, lk, nx, ny, g)] = val;
        }
    }
}

__global__ void apply_channel_robin_side_codim23_l1_kernel(
    double* u, const Box3i* boxes, int n_l1, int axis, int g)
{
    int dom_hi0 = 2 * c_nx0 - 1, dom_hi1 = 2 * c_ny0 - 1, dom_hi2 = 2 * c_nz0 - 1;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t u_off = (g == 1) ? b.off_s_g1[axis] : b.off_s_g3[axis];
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        int snx = nxa + 2 * g, sny = nya + 2 * g, snz = nza + 2 * g;
        int tot = snx * sny * snz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % snx) - g;
            int t = idx / snx;
            int lj = (t % sny) - g;
            int lk = (t / sny) - g;
            bool ox = (li < 0 && b.lo[0] == 0) || (li >= nxa && b.hi[0] == dom_hi0);
            bool oy = (lj < 0 && b.lo[1] == 0) || (lj >= nya && b.hi[1] == dom_hi1);
            bool oz = (lk < 0 && b.lo[2] == 0) || (lk >= nza && b.hi[2] == dom_hi2);
            int oob = (ox ? 1 : 0) + (oy ? 1 : 0) + (oz ? 1 : 0);
            if (oob < 2) continue;
            if (!((li >= 0 && li < nxa) || ox)) continue;
            if (!((lj >= 0 && lj < nya) || oy)) continue;
            if (!((lk >= 0 && lk < nza) || oz)) continue;
            int bi = li, si = 0; if (li < 0) { bi = 0; si = 1; } else if (li >= nxa) { bi = nxa - 1; si = -1; }
            int bj = lj, sj = 0; if (lj < 0) { bj = 0; sj = 1; } else if (lj >= nya) { bj = nya - 1; sj = -1; }
            int bk = lk, sk = 0; if (lk < 0) { bk = 0; sk = 1; } else if (lk >= nza) { bk = nza - 1; sk = -1; }
            double v0 = u[u_off + sidx_g(axis, bi, bj, bk, b.nx, b.ny, g)];
            double val = v0;
            if (si != 0) val += (v0 - u[u_off + sidx_g(axis, bi + si, bj, bk, b.nx, b.ny, g)]) * abs(li - bi);
            if (sj != 0) val += (v0 - u[u_off + sidx_g(axis, bi, bj + sj, bk, b.nx, b.ny, g)]) * abs(lj - bj);
            if (sk != 0) val += (v0 - u[u_off + sidx_g(axis, bi, bj, bk + sk, b.nx, b.ny, g)]) * abs(lk - bk);
            u[u_off + sidx_g(axis, li, lj, lk, b.nx, b.ny, g)] = val;
        }
    }
}

// StaggeredStokesPhysicalBoundaryHelper::enforceDivergenceFreeConditionAtBoundary
// Sweeps k = 0..g-1 outward from EVERY physical boundary face (x_lo, x_up, y_lo, y_up, z_lo, z_up; dx_0=dx_1=dx_2 so factor dx_norm/dx_axis=1)
__global__ void enforce_div_free_phys_bdry_l0_kernel(
    double* u0, double* u1, double* u2, int nx, int ny, int nz, int g, bool all_bdry)
{
    // Pass 0: x_lo (i_g = -1..-g, j in 0..ny-1, k in 0..nz-1) and x_up (i_g = nx..nx+g-1)
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < ny * nz; idx += gridDim.x * blockDim.x)
    {
        int j = idx % ny, k = idx / ny;
        for (int step = 0; step < g; ++step)
        {
            int ig = -1 - step;
            double dv = u0[sidx_g(0, ig + 1, j, k, nx, ny, g)] +
                        (u1[sidx_g(1, ig, j + 1, k, nx, ny, g)] - u1[sidx_g(1, ig, j, k, nx, ny, g)]) +
                        (u2[sidx_g(2, ig, j, k + 1, nx, ny, g)] - u2[sidx_g(2, ig, j, k, nx, ny, g)]);
            u0[sidx_g(0, ig, j, k, nx, ny, g)] = dv;

            ig = nx + step;
            dv = u0[sidx_g(0, ig, j, k, nx, ny, g)] -
                 (u1[sidx_g(1, ig, j + 1, k, nx, ny, g)] - u1[sidx_g(1, ig, j, k, nx, ny, g)]) -
                 (u2[sidx_g(2, ig, j, k + 1, nx, ny, g)] - u2[sidx_g(2, ig, j, k, nx, ny, g)]);
            u0[sidx_g(0, ig + 1, j, k, nx, ny, g)] = dv;
        }
    }
    if (!all_bdry) return;
    // Sweep y_lo/y_up and z_lo/z_up only when ALL_BDRY (SSPRK3 / viscous_rhs), NOT for NORMAL_TRACTION_BDRY (StokesOperator::apply)
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx * nz; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, k = idx / nx;
        for (int step = 0; step < g; ++step)
        {
            int jg = -1 - step;
            double dv = u1[sidx_g(1, i, jg + 1, k, nx, ny, g)] +
                        (u0[sidx_g(0, i + 1, jg, k, nx, ny, g)] - u0[sidx_g(0, i, jg, k, nx, ny, g)]) +
                        (u2[sidx_g(2, i, jg, k + 1, nx, ny, g)] - u2[sidx_g(2, i, jg, k, nx, ny, g)]);
            u1[sidx_g(1, i, jg, k, nx, ny, g)] = dv;

            jg = ny + step;
            dv = u1[sidx_g(1, i, jg, k, nx, ny, g)] -
                 (u0[sidx_g(0, i + 1, jg, k, nx, ny, g)] - u0[sidx_g(0, i, jg, k, nx, ny, g)]) -
                 (u2[sidx_g(2, i, jg, k + 1, nx, ny, g)] - u2[sidx_g(2, i, jg, k, nx, ny, g)]);
            u1[sidx_g(1, i, jg + 1, k, nx, ny, g)] = dv;
        }
    }
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx * ny; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, j = idx / nx;
        for (int step = 0; step < g; ++step)
        {
            int kg = -1 - step;
            double dv = u2[sidx_g(2, i, j, kg + 1, nx, ny, g)] +
                        (u0[sidx_g(0, i + 1, j, kg, nx, ny, g)] - u0[sidx_g(0, i, j, kg, nx, ny, g)]) +
                        (u1[sidx_g(1, i, j + 1, kg, nx, ny, g)] - u1[sidx_g(1, i, j, kg, nx, ny, g)]);
            u2[sidx_g(2, i, j, kg, nx, ny, g)] = dv;

            kg = nz + step;
            dv = u2[sidx_g(2, i, j, kg, nx, ny, g)] -
                 (u0[sidx_g(0, i + 1, j, kg, nx, ny, g)] - u0[sidx_g(0, i, j, kg, nx, ny, g)]) -
                 (u1[sidx_g(1, i, j + 1, kg, nx, ny, g)] - u1[sidx_g(1, i, j, kg, nx, ny, g)]);
            u2[sidx_g(2, i, j, kg + 1, nx, ny, g)] = dv;
        }
    }
}

__global__ void enforce_div_free_phys_bdry_l1_kernel(
    double* u0, double* u1, double* u2, const Box3i* boxes, int n_l1, int g)
{
    int dom_hi1 = 2 * c_ny0 - 1, dom_hi2 = 2 * c_nz0 - 1;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t o0 = (g == 1) ? b.off_s_g1[0] : b.off_s_g3[0];
        size_t o1 = (g == 1) ? b.off_s_g1[1] : b.off_s_g3[1];
        size_t o2 = (g == 1) ? b.off_s_g1[2] : b.off_s_g3[2];
        if (b.lo[1] == 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                for (int step = 0; step < g; ++step)
                {
                    int jg = -1 - step;
                    u1[o1 + sidx_g(1, i, jg, k, b.nx, b.ny, g)] =
                        u1[o1 + sidx_g(1, i, jg + 1, k, b.nx, b.ny, g)] +
                        (u0[o0 + sidx_g(0, i + 1, jg, k, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, jg, k, b.nx, b.ny, g)]) +
                        (u2[o2 + sidx_g(2, i, jg, k + 1, b.nx, b.ny, g)] - u2[o2 + sidx_g(2, i, jg, k, b.nx, b.ny, g)]);
                }
            }
        }
        if (b.hi[1] == dom_hi1)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                for (int step = 0; step < g; ++step)
                {
                    int jg = b.ny + step;
                    u1[o1 + sidx_g(1, i, jg + 1, k, b.nx, b.ny, g)] =
                        u1[o1 + sidx_g(1, i, jg, k, b.nx, b.ny, g)] -
                        (u0[o0 + sidx_g(0, i + 1, jg, k, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, jg, k, b.nx, b.ny, g)]) -
                        (u2[o2 + sidx_g(2, i, jg, k + 1, b.nx, b.ny, g)] - u2[o2 + sidx_g(2, i, jg, k, b.nx, b.ny, g)]);
                }
            }
        }
        if (b.lo[2] == 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                for (int step = 0; step < g; ++step)
                {
                    int kg = -1 - step;
                    u2[o2 + sidx_g(2, i, j, kg, b.nx, b.ny, g)] =
                        u2[o2 + sidx_g(2, i, j, kg + 1, b.nx, b.ny, g)] +
                        (u0[o0 + sidx_g(0, i + 1, j, kg, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, j, kg, b.nx, b.ny, g)]) +
                        (u1[o1 + sidx_g(1, i, j + 1, kg, b.nx, b.ny, g)] - u1[o1 + sidx_g(1, i, j, kg, b.nx, b.ny, g)]);
                }
            }
        }
        if (b.hi[2] == dom_hi2)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                for (int step = 0; step < g; ++step)
                {
                    int kg = b.nz + step;
                    u2[o2 + sidx_g(2, i, j, kg + 1, b.nx, b.ny, g)] =
                        u2[o2 + sidx_g(2, i, j, kg, b.nx, b.ny, g)] -
                        (u0[o0 + sidx_g(0, i + 1, j, kg, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, j, kg, b.nx, b.ny, g)]) -
                        (u1[o1 + sidx_g(1, i, j + 1, kg, b.nx, b.ny, g)] - u1[o1 + sidx_g(1, i, j, kg, b.nx, b.ny, g)]);
                }
            }
        }
    }
}


// ============================================================================
// SAMRAI CONSERVATIVE_LINEAR_REFINE (L0 -> L1 ghosts) for g=3 Side Data
// ============================================================================
__device__ inline double mc_lim(double dm, double dp)
{
    if (dm * dp <= 0.0) return 0.0;
    double abs_m = fabs(dm), abs_p = fabs(dp), abs_c = 0.5 * fabs(dm + dp);
    double m = fmin(2.0 * abs_m, fmin(2.0 * abs_p, abs_c));
    return (dm > 0.0) ? m : -m;
}

__global__ void refine_l0_to_l1_cons_linear_side_g3_kernel(
    const double* u_l0, double* u_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis)
{
    const int g = 3;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const size_t o1 = b.off_s_g3[axis];
        const int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        const int snx = nxa + 2 * g, sny = nya + 2 * g, snz = nza + 2 * g;
        const int tot = snx * sny * snz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % snx) - g;
            int t = idx / snx;
            int lj = (t % sny) - g;
            int lk = (t / sny) - g;
            if (li >= 0 && li < nxa && lj >= 0 && lj < nya && lk >= 0 && lk < nza) continue;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            int ql, qj, qk, q;
            // Skip physical outside domain (handled by physical boundary op) and skip fine-owned faces (handled by intra-level sync)
            int max_gi = 2 * c_nx0 - (axis == 0 ? 0 : 1);
            int max_gj = 2 * c_ny0 - (axis == 1 ? 0 : 1);
            int max_gk = 2 * c_nz0 - (axis == 2 ? 0 : 1);
            if (gi < 0 || gi > max_gi || gj < 0 || gj > max_gj || gk < 0 || gk > max_gk) continue;
            if (lookup_l1_side_owner(axis, gi, gj, gk, l1_map, boxes, q, ql, qj, qk)) continue;

            // Map fine side index to parent coarse side index + sub-offsets matching cartclinrefsidedoub3d
            int ic = gi >> 1, jc = gj >> 1, kc = gk >> 1;
            int sub_i = gi - 2 * ic, sub_j = gj - 2 * jc, sub_k = gk - 2 * kc;
            double ox = (axis == 0) ? (0.5 * sub_i) : (-0.25 + 0.5 * sub_i);
            double oy = (axis == 1) ? (0.5 * sub_j) : (-0.25 + 0.5 * sub_j);
            double oz = (axis == 2) ? (0.5 * sub_k) : (-0.25 + 0.5 * sub_k);

            double vc = u_l0[sidx_g3(axis, ic, jc, kc, c_nx0, c_ny0)];
            double sx = mc_lim(vc - u_l0[sidx_g3(axis, ic - 1, jc, kc, c_nx0, c_ny0)],
                               u_l0[sidx_g3(axis, ic + 1, jc, kc, c_nx0, c_ny0)] - vc);
            double sy = mc_lim(vc - u_l0[sidx_g3(axis, ic, jc - 1, kc, c_nx0, c_ny0)],
                               u_l0[sidx_g3(axis, ic, jc + 1, kc, c_nx0, c_ny0)] - vc);
            double sz = mc_lim(vc - u_l0[sidx_g3(axis, ic, jc, kc - 1, c_nx0, c_ny0)],
                               u_l0[sidx_g3(axis, ic, jc, kc + 1, c_nx0, c_ny0)] - vc);

            u_l1[o1 + sidx_g3(axis, li, lj, lk, b.nx, b.ny)] = vc + ox * sx + oy * sy + oz * sz;
        }
    }
}

// INSVCStaggeredConservativeMassMomentumRKIntegrator::enforceDivergenceFreeConditionAtCoarseFineInterface
// Sweeps k = 0..g-1 outward on CF_GHOST_WIDTH = 1 strip along every codim-1 C-F face of L1 patch `pat`
__global__ void enforce_div_free_cf_l1_kernel(
    double* u0, double* u1, double* u2, const Box3i* boxes, int n_l1, const int16_t* l1_map, int g)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        size_t o0 = (g == 1) ? b.off_s_g1[0] : b.off_s_g3[0];
        size_t o1 = (g == 1) ? b.off_s_g1[1] : b.off_s_g3[1];
        size_t o2 = (g == 1) ? b.off_s_g1[2] : b.off_s_g3[2];

        // loc 0 (x_lo, is_lower=true) and loc 1 (x_up, is_lower=false)
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
        {
            int j = idx % b.ny, k = idx / b.ny;
            // x_lo: bdry cell in L1 ghost adjacent to x_lo face is li = -1 (`gi = b.lo[0] - 1`). Only C-F if gi >= 0 and lookup < 0.
            if (b.lo[0] > 0 && lookup_l1_patch(b.lo[0] - 1, b.lo[1] + j, b.lo[2] + k, l1_map) < 0)
            {
                for (int step = 0; step < g; ++step)
                {
                    int is = -1 - step;
                    u0[o0 + sidx_g(0, is, j, k, b.nx, b.ny, g)] =
                        u0[o0 + sidx_g(0, is + 1, j, k, b.nx, b.ny, g)] +
                        (u1[o1 + sidx_g(1, is, j + 1, k, b.nx, b.ny, g)] - u1[o1 + sidx_g(1, is, j, k, b.nx, b.ny, g)]) +
                        (u2[o2 + sidx_g(2, is, j, k + 1, b.nx, b.ny, g)] - u2[o2 + sidx_g(2, is, j, k, b.nx, b.ny, g)]);
                }
            }
            if (b.hi[0] + 1 < 2 * c_nx0 && lookup_l1_patch(b.hi[0] + 1, b.lo[1] + j, b.lo[2] + k, l1_map) < 0)
            {
                for (int step = 0; step < g; ++step)
                {
                    int is = b.nx + step;
                    u0[o0 + sidx_g(0, is + 1, j, k, b.nx, b.ny, g)] =
                        u0[o0 + sidx_g(0, is, j, k, b.nx, b.ny, g)] -
                        (u1[o1 + sidx_g(1, is, j + 1, k, b.nx, b.ny, g)] - u1[o1 + sidx_g(1, is, j, k, b.nx, b.ny, g)]) -
                        (u2[o2 + sidx_g(2, is, j, k + 1, b.nx, b.ny, g)] - u2[o2 + sidx_g(2, is, j, k, b.nx, b.ny, g)]);
                }
            }
        }
        // loc 2 (y_lo) and loc 3 (y_up)
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, k = idx / b.nx;
            if (b.lo[1] > 0 && lookup_l1_patch(b.lo[0] + i, b.lo[1] - 1, b.lo[2] + k, l1_map) < 0)
            {
                for (int step = 0; step < g; ++step)
                {
                    int js = -1 - step;
                    u1[o1 + sidx_g(1, i, js, k, b.nx, b.ny, g)] =
                        u1[o1 + sidx_g(1, i, js + 1, k, b.nx, b.ny, g)] +
                        (u0[o0 + sidx_g(0, i + 1, js, k, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, js, k, b.nx, b.ny, g)]) +
                        (u2[o2 + sidx_g(2, i, js, k + 1, b.nx, b.ny, g)] - u2[o2 + sidx_g(2, i, js, k, b.nx, b.ny, g)]);
                }
            }
            if (b.hi[1] + 1 < 2 * c_ny0 && lookup_l1_patch(b.lo[0] + i, b.hi[1] + 1, b.lo[2] + k, l1_map) < 0)
            {
                for (int step = 0; step < g; ++step)
                {
                    int js = b.ny + step;
                    u1[o1 + sidx_g(1, i, js + 1, k, b.nx, b.ny, g)] =
                        u1[o1 + sidx_g(1, i, js, k, b.nx, b.ny, g)] -
                        (u0[o0 + sidx_g(0, i + 1, js, k, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, js, k, b.nx, b.ny, g)]) -
                        (u2[o2 + sidx_g(2, i, js, k + 1, b.nx, b.ny, g)] - u2[o2 + sidx_g(2, i, js, k, b.nx, b.ny, g)]);
                }
            }
        }
        // loc 4 (z_lo) and loc 5 (z_up)
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, j = idx / b.nx;
            if (b.lo[2] > 0 && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] - 1, l1_map) < 0)
            {
                for (int step = 0; step < g; ++step)
                {
                    int ks = -1 - step;
                    u2[o2 + sidx_g(2, i, j, ks, b.nx, b.ny, g)] =
                        u2[o2 + sidx_g(2, i, j, ks + 1, b.nx, b.ny, g)] +
                        (u0[o0 + sidx_g(0, i + 1, j, ks, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, j, ks, b.nx, b.ny, g)]) +
                        (u1[o1 + sidx_g(1, i, j + 1, ks, b.nx, b.ny, g)] - u1[o1 + sidx_g(1, i, j, ks, b.nx, b.ny, g)]);
                }
            }
            if (b.hi[2] + 1 < 2 * c_nz0 && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.hi[2] + 1, l1_map) < 0)
            {
                for (int step = 0; step < g; ++step)
                {
                    int ks = b.nz + step;
                    u2[o2 + sidx_g(2, i, j, ks + 1, b.nx, b.ny, g)] =
                        u2[o2 + sidx_g(2, i, j, ks, b.nx, b.ny, g)] -
                        (u0[o0 + sidx_g(0, i + 1, j, ks, b.nx, b.ny, g)] - u0[o0 + sidx_g(0, i, j, ks, b.nx, b.ny, g)]) -
                        (u1[o1 + sidx_g(1, i, j + 1, ks, b.nx, b.ny, g)] - u1[o1 + sidx_g(1, i, j, ks, b.nx, b.ny, g)]);
                }
            }
        }
    }
}


// ============================================================================
// Fused SSPRK3 + CUI Stencil Kernels (matching navier_stokes_staggered_helpers3d.f.m4)
// ============================================================================
__device__ inline double cui_hr(double qU, double qC, double qD)
{
    double denom = qD - qU;
    if (denom == 0.0) return qC;
    double ac = (qC - qU) / denom;
    double af;
    if (0.0 < ac && ac <= (2.0 / 13.0)) af = 3.0 * ac;
    else if ((2.0 / 13.0) < ac && ac <= 0.8) af = ac * (5.0 / 6.0) + (1.0 / 3.0);
    else if (0.8 < ac && ac <= 1.0) af = 1.0;
    else af = ac;
    return af * denom + qU;
}

// Helper: compute flux u_adv * R_half (and optionally u_adv * R_half * V_half) on face (fi, fj, fk) along `flux_dir`
// for control volumes centered at side `comp_a`!
__device__ inline double face_v_adv(
    int comp_a, int flux_dir, int fi, int fj, int fk,
    const double* u0, const double* u1, const double* u2,
    int nx, int ny, int g)
{
    const double* u_f = (flux_dir == 0) ? u0 : (flux_dir == 1 ? u1 : u2);
    int di = (comp_a == 0) ? 1 : 0;
    int dj = (comp_a == 1) ? 1 : 0;
    int dk = (comp_a == 2) ? 1 : 0;
    return 0.5 * (u_f[sidx_g(flux_dir, fi - di, fj - dj, fk - dk, nx, ny, g)] +
                  u_f[sidx_g(flux_dir, fi, fj, fk, nx, ny, g)]);
}

__device__ inline double face_cui_interp(
    double v_adv, int comp_a, int flux_dir, int fi, int fj, int fk,
    const double* q_a, int nx, int ny, int g)
{
    int si = (flux_dir == 0) ? 1 : 0;
    int sj = (flux_dir == 1) ? 1 : 0;
    int sk = (flux_dir == 2) ? 1 : 0;
    double qU, qC, qD;
    if (v_adv >= 0.0)
    {
        qC = q_a[sidx_g(comp_a, fi - si, fj - sj, fk - sk, nx, ny, g)];
        qU = q_a[sidx_g(comp_a, fi - 2 * si, fj - 2 * sj, fk - 2 * sk, nx, ny, g)];
        qD = q_a[sidx_g(comp_a, fi, fj, fk, nx, ny, g)];
    }
    else
    {
        qC = q_a[sidx_g(comp_a, fi, fj, fk, nx, ny, g)];
        qU = q_a[sidx_g(comp_a, fi + si, fj + sj, fk + sk, nx, ny, g)];
        qD = q_a[sidx_g(comp_a, fi - si, fj - sj, fk - sk, nx, ny, g)];
    }
    return cui_hr(qU, qC, qD);
}

// Compute SSPRK stage density update (step = 0, 1, 2):
// R_new = a0 * R_cur + a1 * R_pre - a2 * dt * div(V_adv * R_half)
// and when compute_final == true (step == 2):
//   P_half = R_half * V_half
//   N_out = div(V_adv * P_half) - V_cur * ((R_new - R_cur) / dt + div(V_adv * R_half))
__global__ void ssprk3_stage_update_kernel(
    int comp_a, double dt, double a0, double a1, double a2, bool compute_final,
    const double* r_cur_g0, const double* r_pre_g3,
    const double* v0_g3, const double* v1_g3, const double* v2_g3,
    const double* v_cur_g0,
    double* r_new_g0, double* n_out_g0,
    int nx, int ny, int nz, double h,
    bool is_l1, const Box3i* boxes, int n_l1,
    size_t r_cur_off0, size_t r_pre_off3, size_t v_off3_0, size_t v_off3_1, size_t v_off3_2)
{
    const int g = 3;
    const double inv_h = 1.0 / h;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o_g0 = is_l1 ? b.off_s_g0[comp_a] : 0;
        size_t o_pre3 = is_l1 ? b.off_s_g3[comp_a] : 0;
        size_t o_v0 = is_l1 ? b.off_s_g3[0] : 0;
        size_t o_v1 = is_l1 ? b.off_s_g3[1] : 0;
        size_t o_v2 = is_l1 ? b.off_s_g3[2] : 0;

        const double* rc = r_cur_g0 + o_g0;
        const double* rp = r_pre_g3 + o_pre3;
        const double* vu0 = v0_g3 + o_v0;
        const double* vu1 = v1_g3 + o_v1;
        const double* vu2 = v2_g3 + o_v2;
        const double* va = (comp_a == 0) ? vu0 : (comp_a == 1 ? vu1 : vu2);
        const double* vc = v_cur_g0 + o_g0;
        double* rn = r_new_g0 + o_g0;
        double* no = compute_final ? (n_out_g0 + o_g0) : nullptr;

        int nxa = b.nx + (comp_a == 0), nya = b.ny + (comp_a == 1), nza = b.nz + (comp_a == 2);
        int tot = nxa * nya * nza;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
            // Flux x
            double vx_p = face_v_adv(comp_a, 0, i + 1, j, k, vu0, vu1, vu2, b.nx, b.ny, g);
            double vx_m = face_v_adv(comp_a, 0, i,     j, k, vu0, vu1, vu2, b.nx, b.ny, g);
            double rx_p = face_cui_interp(vx_p, comp_a, 0, i + 1, j, k, rp, b.nx, b.ny, g);
            double rx_m = face_cui_interp(vx_m, comp_a, 0, i,     j, k, rp, b.nx, b.ny, g);
            double px0 = (rx_p * vx_p - rx_m * vx_m) * inv_h;

            // Flux y
            double vy_p = face_v_adv(comp_a, 1, i, j + 1, k, vu0, vu1, vu2, b.nx, b.ny, g);
            double vy_m = face_v_adv(comp_a, 1, i, j,     k, vu0, vu1, vu2, b.nx, b.ny, g);
            double ry_p = face_cui_interp(vy_p, comp_a, 1, i, j + 1, k, rp, b.nx, b.ny, g);
            double ry_m = face_cui_interp(vy_m, comp_a, 1, i, j,     k, rp, b.nx, b.ny, g);
            double px1 = (ry_p * vy_p - ry_m * vy_m) * inv_h;

            // Flux z
            double vz_p = face_v_adv(comp_a, 2, i, j, k + 1, vu0, vu1, vu2, b.nx, b.ny, g);
            double vz_m = face_v_adv(comp_a, 2, i, j, k,     vu0, vu1, vu2, b.nx, b.ny, g);
            double rz_p = face_cui_interp(vz_p, comp_a, 2, i, j, k + 1, rp, b.nx, b.ny, g);
            double rz_m = face_cui_interp(vz_m, comp_a, 2, i, j, k,     rp, b.nx, b.ny, g);
            double px2 = (rz_p * vz_p - rz_m * vz_m) * inv_h;

            double div_mass = px0 + px1 + px2;
            double r_old = rc[idx];
            double r_pre = rp[sidx_g3(comp_a, i, j, k, b.nx, b.ny)];
            // Exact evaluation order matching vc_update_density3d (div subtracted along x then y then z):
            double r_upd = a0 * r_old + a1 * r_pre - a2 * dt * px0;
            r_upd = r_upd - a2 * dt * px1;
            r_upd = r_upd - a2 * dt * px2;
            rn[idx] = r_upd;

            if (compute_final)
            {
                double ux_p = face_cui_interp(vx_p, comp_a, 0, i + 1, j, k, va, b.nx, b.ny, g);
                double ux_m = face_cui_interp(vx_m, comp_a, 0, i,     j, k, va, b.nx, b.ny, g);
                double uy_p = face_cui_interp(vy_p, comp_a, 1, i, j + 1, k, va, b.nx, b.ny, g);
                double uy_m = face_cui_interp(vy_m, comp_a, 1, i, j,     k, va, b.nx, b.ny, g);
                double uz_p = face_cui_interp(vz_p, comp_a, 2, i, j, k + 1, va, b.nx, b.ny, g);
                double uz_m = face_cui_interp(vz_m, comp_a, 2, i, j, k,     va, b.nx, b.ny, g);

                double n_val = (vx_p * (rx_p * ux_p) - vx_m * (rx_m * ux_m)) * inv_h;
                n_val += (vy_p * (ry_p * uy_p) - vy_m * (ry_m * uy_m)) * inv_h;
                n_val += (vz_p * (rz_p * uz_p) - vz_m * (rz_m * uz_m)) * inv_h;

                // vc_mass_conservation_error3d: E = (Rnew - Rold)/dt + px0 + px1 + px2
                double E = ((r_upd - r_old) / dt + px0) + px1 + px2;
                n_val -= vc[idx] * E;
                no[idx] = n_val;
            }
        }
    }
    (void)r_cur_off0; (void)r_pre_off3; (void)v_off3_0; (void)v_off3_1; (void)v_off3_2;
}

// ============================================================================
// Quadratic C-F Interpolation Kernels (g=1, matching quadcfinterpolation3d.f.m4)
// Used strictly by Native IBAMR VCStaggeredStokesOperator::apply!
// ============================================================================
__global__ void quad_cf_tangential_cell_l1_kernel(
    const double* p_l0, double* p_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const size_t o1 = b.off_c_g1;
        // Fill 1-cell ghost layer on each C-F face
        // loc 0 (x_lo, li = -1) & loc 1 (x_up, li = b.nx)
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
        {
            int j = idx % b.ny, k = idx / b.ny;
            int gj = b.lo[1] + j, gk = b.lo[2] + k;
            int jc = gj >> 1, kc = gk >> 1;
            int jp = gj & 1, kp = gk & 1;
            double t1 = jp + 0.5 - 1.0, t2 = kp + 0.5 - 1.0;
            double w1[3] = {0.125 * t1 * t1 - 0.25 * t1, -0.25 * t1 * t1 + 1.0, 0.125 * t1 * t1 + 0.25 * t1};
            double w2[3] = {0.125 * t2 * t2 - 0.25 * t2, -0.25 * t2 * t2 + 1.0, 0.125 * t2 * t2 + 0.25 * t2};
            if (b.lo[0] > 0 && lookup_l1_patch(b.lo[0] - 1, gj, gk, l1_map) < 0)
            {
                int ic = (b.lo[0] - 1) >> 1;
                double sum = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wj = -1; wj <= 1; ++wj)
                        sum += w1[wj + 1] * w2[wk + 1] * p_l0[cidx_g1(ic, jc + wj, kc + wk, c_nx0, c_ny0)];
                p_l1[o1 + cidx_g1(-1, j, k, b.nx, b.ny)] = sum;
            }
            if (b.hi[0] + 1 < 2 * c_nx0 && lookup_l1_patch(b.hi[0] + 1, gj, gk, l1_map) < 0)
            {
                int ic = (b.hi[0] + 1) >> 1;
                double sum = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wj = -1; wj <= 1; ++wj)
                        sum += w1[wj + 1] * w2[wk + 1] * p_l0[cidx_g1(ic, jc + wj, kc + wk, c_nx0, c_ny0)];
                p_l1[o1 + cidx_g1(b.nx, j, k, b.nx, b.ny)] = sum;
            }
        }
        // loc 2 & loc 3 (y_lo, y_up)
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, k = idx / b.nx;
            int gi = b.lo[0] + i, gk = b.lo[2] + k;
            int ic = gi >> 1, kc = gk >> 1;
            int ip = gi & 1, kp = gk & 1;
            double t0 = ip + 0.5 - 1.0, t2 = kp + 0.5 - 1.0;
            double w0[3] = {0.125 * t0 * t0 - 0.25 * t0, -0.25 * t0 * t0 + 1.0, 0.125 * t0 * t0 + 0.25 * t0};
            double w2[3] = {0.125 * t2 * t2 - 0.25 * t2, -0.25 * t2 * t2 + 1.0, 0.125 * t2 * t2 + 0.25 * t2};
            if (b.lo[1] > 0 && lookup_l1_patch(gi, b.lo[1] - 1, gk, l1_map) < 0)
            {
                int jc = (b.lo[1] - 1) >> 1;
                double sum = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wi = -1; wi <= 1; ++wi)
                        sum += w0[wi + 1] * w2[wk + 1] * p_l0[cidx_g1(ic + wi, jc, kc + wk, c_nx0, c_ny0)];
                p_l1[o1 + cidx_g1(i, -1, k, b.nx, b.ny)] = sum;
            }
            if (b.hi[1] + 1 < 2 * c_ny0 && lookup_l1_patch(gi, b.hi[1] + 1, gk, l1_map) < 0)
            {
                int jc = (b.hi[1] + 1) >> 1;
                double sum = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wi = -1; wi <= 1; ++wi)
                        sum += w0[wi + 1] * w2[wk + 1] * p_l0[cidx_g1(ic + wi, jc, kc + wk, c_nx0, c_ny0)];
                p_l1[o1 + cidx_g1(i, b.ny, k, b.nx, b.ny)] = sum;
            }
        }
        // loc 4 & loc 5 (z_lo, z_up)
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, j = idx / b.nx;
            int gi = b.lo[0] + i, gj = b.lo[1] + j;
            int ic = gi >> 1, jc = gj >> 1;
            int ip = gi & 1, jp = gj & 1;
            double t0 = ip + 0.5 - 1.0, t1 = jp + 0.5 - 1.0;
            double w0[3] = {0.125 * t0 * t0 - 0.25 * t0, -0.25 * t0 * t0 + 1.0, 0.125 * t0 * t0 + 0.25 * t0};
            double w1[3] = {0.125 * t1 * t1 - 0.25 * t1, -0.25 * t1 * t1 + 1.0, 0.125 * t1 * t1 + 0.25 * t1};
            if (b.lo[2] > 0 && lookup_l1_patch(gi, gj, b.lo[2] - 1, l1_map) < 0)
            {
                int kc = (b.lo[2] - 1) >> 1;
                double sum = 0.0;
                for (int wj = -1; wj <= 1; ++wj)
                    for (int wi = -1; wi <= 1; ++wi)
                        sum += w0[wi + 1] * w1[wj + 1] * p_l0[cidx_g1(ic + wi, jc + wj, kc, c_nx0, c_ny0)];
                p_l1[o1 + cidx_g1(i, j, -1, b.nx, b.ny)] = sum;
            }
            if (b.hi[2] + 1 < 2 * c_nz0 && lookup_l1_patch(gi, gj, b.hi[2] + 1, l1_map) < 0)
            {
                int kc = (b.hi[2] + 1) >> 1;
                double sum = 0.0;
                for (int wj = -1; wj <= 1; ++wj)
                    for (int wi = -1; wi <= 1; ++wi)
                        sum += w0[wi + 1] * w1[wj + 1] * p_l0[cidx_g1(ic + wi, jc + wj, kc, c_nx0, c_ny0)];
                p_l1[o1 + cidx_g1(i, j, b.nz, b.nx, b.ny)] = sum;
            }
        }
    }
}


__global__ void quad_cf_normal_cell_l1_kernel(
    double* p_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    // Wait: CartCellDoubleQuadraticCFInterpolation normal extension at g=1 (ratio R=2):
    // t=0.5 => wgt(-1) = -0.5*(2*0.5+2-1-1)/5 = -0.5*(2)/5 = -0.2,
    //          wgt(0)  = 0.5*(2*0.5+6-1-3)/3 = 0.5*(3)/3 = +0.5? Wait check exact lines 375-380 of quadcfinterpolation3d.f.m4:
    // Let's check: R=2, t=0.5:
    //  2*R*t + R - 4*t^2 - 2*t = 2 + 2 - 1 - 1 = 2 => wgt(-1) = -0.5 * 2 / 5 = -0.2
    //  2*R*t + 3*R - 4*t^2 - 6*t = 2 + 6 - 1 - 3 = 4 => wgt(0) = 0.5 * 4 / 3 = +2/3
    //  4*t^2 + 8*t + 3 = 1 + 4 + 3 = 8 => wgt(+1) = 8 / (3*5) = 8/15! (Note: -0.2 + 2/3 + 8/15 = (-3+10+8)/15 = 15/15 = 1.0!)
    // Exact! And note only overwrite if tangential C-F ghost exists and not covered by another L1 patch.
    const double wm = -0.2, w0 = 2.0 / 3.0, wp = 8.0 / 15.0;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const size_t o1 = b.off_c_g1;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
        {
            int j = idx % b.ny, k = idx / b.ny;
            if (b.lo[0] > 0 && lookup_l1_patch(b.lo[0] - 1, b.lo[1] + j, b.lo[2] + k, l1_map) < 0)
            {
                double wc = p_l1[o1 + cidx_g1(-1, j, k, b.nx, b.ny)];
                double wi0 = p_l1[o1 + cidx_g1(0, j, k, b.nx, b.ny)];
                double wi1 = p_l1[o1 + cidx_g1(1, j, k, b.nx, b.ny)];
                p_l1[o1 + cidx_g1(-1, j, k, b.nx, b.ny)] = wp * wc + w0 * wi0 + wm * wi1;
            }
            if (b.hi[0] + 1 < 2 * c_nx0 && lookup_l1_patch(b.hi[0] + 1, b.lo[1] + j, b.lo[2] + k, l1_map) < 0)
            {
                double wc = p_l1[o1 + cidx_g1(b.nx, j, k, b.nx, b.ny)];
                double wi0 = p_l1[o1 + cidx_g1(b.nx - 1, j, k, b.nx, b.ny)];
                double wi1 = p_l1[o1 + cidx_g1(b.nx - 2, j, k, b.nx, b.ny)];
                p_l1[o1 + cidx_g1(b.nx, j, k, b.nx, b.ny)] = wp * wc + w0 * wi0 + wm * wi1;
            }
        }
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, k = idx / b.nx;
            if (b.lo[1] > 0 && lookup_l1_patch(b.lo[0] + i, b.lo[1] - 1, b.lo[2] + k, l1_map) < 0)
            {
                double wc = p_l1[o1 + cidx_g1(i, -1, k, b.nx, b.ny)];
                double wi0 = p_l1[o1 + cidx_g1(i, 0, k, b.nx, b.ny)];
                double wi1 = p_l1[o1 + cidx_g1(i, 1, k, b.nx, b.ny)];
                p_l1[o1 + cidx_g1(i, -1, k, b.nx, b.ny)] = wp * wc + w0 * wi0 + wm * wi1;
            }
            if (b.hi[1] + 1 < 2 * c_ny0 && lookup_l1_patch(b.lo[0] + i, b.hi[1] + 1, b.lo[2] + k, l1_map) < 0)
            {
                double wc = p_l1[o1 + cidx_g1(i, b.ny, k, b.nx, b.ny)];
                double wi0 = p_l1[o1 + cidx_g1(i, b.ny - 1, k, b.nx, b.ny)];
                double wi1 = p_l1[o1 + cidx_g1(i, b.ny - 2, k, b.nx, b.ny)];
                p_l1[o1 + cidx_g1(i, b.ny, k, b.nx, b.ny)] = wp * wc + w0 * wi0 + wm * wi1;
            }
        }
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
        {
            int i = idx % b.nx, j = idx / b.nx;
            if (b.lo[2] > 0 && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] - 1, l1_map) < 0)
            {
                double wc = p_l1[o1 + cidx_g1(i, j, -1, b.nx, b.ny)];
                double wi0 = p_l1[o1 + cidx_g1(i, j, 0, b.nx, b.ny)];
                double wi1 = p_l1[o1 + cidx_g1(i, j, 1, b.nx, b.ny)];
                p_l1[o1 + cidx_g1(i, j, -1, b.nx, b.ny)] = wp * wc + w0 * wi0 + wm * wi1;
            }
            if (b.hi[2] + 1 < 2 * c_nz0 && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.hi[2] + 1, l1_map) < 0)
            {
                double wc = p_l1[o1 + cidx_g1(i, j, b.nz, b.nx, b.ny)];
                double wi0 = p_l1[o1 + cidx_g1(i, j, b.nz - 1, b.nx, b.ny)];
                double wi1 = p_l1[o1 + cidx_g1(i, j, b.nz - 2, b.nx, b.ny)];
                p_l1[o1 + cidx_g1(i, j, b.nz, b.nx, b.ny)] = wp * wc + w0 * wi0 + wm * wi1;
            }
        }
    }
}


// Helper for SC quadratic+cubic weights in quadcfinterpolation3d.f.m4 (R=2):
// quad_wgt(m in {-1,0,+1}, p in {0,1}): t = p + 0.5 (p=0 -> t=0.5, p=1 -> t=1.5)
//   quad(-1, t) = 0.125*(4*t^2 - 16*t + 12)/4 -> p=0: 5/32 = 0.15625; p=1: -3/32 = -0.09375
//   quad( 0, t) = 0.25*(-4*t^2 + 8*t + 12)/4 -> p=0: 15/16 = 0.9375; p=1: 15/16 = 0.9375
//   quad(+1, t) = -0.125*(-4*t^2 + 4)/4     -> p=0: -3/32 = -0.09375; p=1: 5/32 = 0.15625
// Wait check sum: p=0: (5+30-3)/32 = 1.0! p=1: (-3+30+5)/32 = 1.0!
// cubic_wgt(m in {-1,0,1,2}, p in {0,1,2}): t = p (R=2, denom=8):
//   p=0 -> {0, 1, 0, 0}; p=2 -> {0, 0, 1, 0};
//   p=1 (t=1):
//     m=-1: -(1)*(1-6+8)/48 = -3/48 = -1/16 = -0.0625
//     m= 0: 0.5*(1-4-4+16)/8 = 9/16 = 0.5625
//     m=+1: 0.5*1*(-1+2+8)/8 = 9/16 = 0.5625
//     m=+2: -(1)*(-1+4)/48 = -3/48 = -1/16 = -0.0625
__device__ inline void sc_quad_w(int p, double w[3])
{
    if (p == 0) { w[0] = 0.15625; w[1] = 0.9375; w[2] = -0.09375; }
    else { w[0] = -0.09375; w[1] = 0.9375; w[2] = 0.15625; }
}
__device__ inline void sc_cubic_w(int p, double w[4])
{
    if (p == 0) { w[0] = 0.0; w[1] = 1.0; w[2] = 0.0; w[3] = 0.0; }
    else if (p == 1) { w[0] = -0.0625; w[1] = 0.5625; w[2] = 0.5625; w[3] = -0.0625; }
    else { w[0] = 0.0; w[1] = 0.0; w[2] = 1.0; w[3] = 0.0; }
}

// Side-centered C-F quadratic tangential + normal interpolation (R=2, g=1), matching scquadtangentialinterpolation3d & scquadnormalinterpolation3d
// Two-step clean execution:
//  Step A: Tangential interp into W_l1 (or directly into U_l1 on faces where indicator != 1, leaving interior untouched).
//  Step B: Copy U_l1 into W_l1 (or compute W_l1 from U_l0 and U_l1), then apply normal quadratic interp into U_l1!
// Notice on normal component `loc 0/1` (`x_lo/x_up`):
//  `indicator0` == 1 on the C-F boundary face itself (`li = 0` or `li = b.nx`), so tangential sets `li = -1` (x_lo) or `li = b.nx+1` (x_up) using `i_c = (b.lo[0]-1)/2` for loc 0 and `i_c = (b.hi[0]+1)/2 + 1` for loc 1!
//  Then normal uses `t0 = 1` (`R=2`) -> for `scquadnormalinterpolation3d` normal component (`k=0`):
//  `quad_wgt(-1) = -1*(-1+2)/3 = -1/3`, `quad_wgt(0) = (2+2-1-1)/2 = 1.0`?? WAIT: check line 1238-1240:
//  `t0 = 1.0, R0 = 2.0`: `quad(-1) = -1*(-1+2)/3 = -1/3`; `quad(0) = (2+2-1-1)/2 = 1`; `quad(+1) = 1*2 / 2 / 3 = +1/3`! (-1/3 + 1 + 1/3 = 1.0!!)
//  And look at indices in line 1248 vs 1261:
//  `loc=0`: `bupper(0) = ilowerf(0) - 1`, `U0(bupper(0), j, k)` IS `li = -1`! Look: `wp * W0(-1) + w0 * W0(0) + wm * W0(+1)`!! (uses `li = -1`, `0`, `+1` -> where `W0(-1)` was just filled by tangential C-F interp at `x = -1.5*dx1`? No, coarse face `i_c = b.lo[0]/2 - 1` is at `x = -2*dx1` (`X = -2`), `li = 0` is at `X = 0`, `li = +1` is at `X = +1`, and target `li = -1` is at `X = -1`!! Check Lagrange polynomial through `(-2, 0, +1)` at `-1`:
//  At `X = -1`: `w(-2) = (-1-0)(-1-1)/((-2-0)(-2-1)) = 2/6 = +1/3`! `w(0) = (-1+2)(-1-1)/((0+2)(0-1)) = -2/-2 = +1`! `w(+1) = (-1+2)(-1-0)/((1+2)(1-0)) = -1/3`!! EXACT!!)
//  And for tangential component on codim-1 `loc 0/1` (`lines 1272-1281`): `t0 = 0.5, R0 = 2.0`:
//  Same cell normal coeffs `{wp=8/15, w0=2/3, wm=-1/5}` (`U1(li=-1) = 8/15*W1(-1) + 2/3*W1(0) - 1/5*W1(+1)`)!
__global__ void quad_cf_tangential_side_l1_kernel(
    const double* u0_l0, const double* u1_l0, const double* u2_l0,
    double* w0_l1, double* w1_l1, double* w2_l1,
    const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        size_t o0 = b.off_s_g1[0], o1 = b.off_s_g1[1], o2 = b.off_s_g1[2];
        // loc 0 (x_lo, li = -1) & loc 1 (x_up, li = b.nx)
        auto do_x_face = [&](bool is_lo) {
            int li_g = is_lo ? -1 : b.nx;
            int gi_g = b.lo[0] + li_g;
            int ic_n = is_lo ? (gi_g >> 1) : ((gi_g >> 1) + 1); // shifted for normal u0!
            int ic_t = is_lo ? ic_n : (ic_n - 1);               // shifted back for tangential u1, u2!
            int li_n = is_lo ? -1 : (b.nx + 1);
            int li_t = is_lo ? -1 : b.nx;
            // u0 (n=x) at li_n
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
            {
                int j = idx % b.ny, k = idx / b.ny;
                int gj = b.lo[1] + j, gk = b.lo[2] + k;
                if (lookup_l1_patch(gi_g, gj, gk, l1_map) >= 0) continue;
                double w1[3], w2[3]; sc_quad_w(gj & 1, w1); sc_quad_w(gk & 1, w2);
                int jc = gj >> 1, kc = gk >> 1;
                double s = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wj = -1; wj <= 1; ++wj)
                    {
                        double wt = w1[wj + 1] * w2[wk + 1];
                        if (wt != 0.0) s += wt * u0_l0[sidx_g1(0, ic_n, jc + wj, kc + wk, c_nx0, c_ny0)];
                    }
                w0_l1[o0 + sidx_g1(0, li_n, j, k, b.nx, b.ny)] = s;
            }
            // u1 (cubic in y, quad in z)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.ny + 1) * b.nz; idx += gridDim.x * blockDim.x)
            {
                int j = idx % (b.ny + 1), k = idx / (b.ny + 1);
                int gj = b.lo[1] + j, gk = b.lo[2] + k;
                int ql, qj, qk, q;
                if (lookup_l1_side_owner(1, b.lo[0] + li_t, gj, gk, l1_map, boxes, q, ql, qj, qk)) continue;
                int jc = (gj & 1) ? ((gj - 1) >> 1) : (gj >> 1);
                int jp = (gj & 1) ? 1 : 0;
                int kc = gk >> 1;
                double c1[4], w2[3]; sc_cubic_w(jp, c1); sc_quad_w(gk & 1, w2);
                double s = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wj = -1; wj <= 2; ++wj)
                    {
                        double wt = c1[wj + 1] * w2[wk + 1];
                        if (wt != 0.0) s += wt * u1_l0[sidx_g1(1, ic_t, jc + wj, kc + wk, c_nx0, c_ny0)];
                    }
                w1_l1[o1 + sidx_g1(1, li_t, j, k, b.nx, b.ny)] = s;
            }
            // u2 (quad in y, cubic in z)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * (b.nz + 1); idx += gridDim.x * blockDim.x)
            {
                int j = idx % b.ny, k = idx / b.ny;
                int gj = b.lo[1] + j, gk = b.lo[2] + k;
                int ql, qj, qk, q;
                if (lookup_l1_side_owner(2, b.lo[0] + li_t, gj, gk, l1_map, boxes, q, ql, qj, qk)) continue;
                int jc = gj >> 1;
                int kc = (gk & 1) ? ((gk - 1) >> 1) : (gk >> 1);
                int kp = (gk & 1) ? 1 : 0;
                double w1[3], c2[4]; sc_quad_w(gj & 1, w1); sc_cubic_w(kp, c2);
                double s = 0.0;
                for (int wk = -1; wk <= 2; ++wk)
                    for (int wj = -1; wj <= 1; ++wj)
                    {
                        double wt = w1[wj + 1] * c2[wk + 1];
                        if (wt != 0.0) s += wt * u2_l0[sidx_g1(2, ic_t, jc + wj, kc + wk, c_nx0, c_ny0)];
                    }
                w2_l1[o2 + sidx_g1(2, li_t, j, k, b.nx, b.ny)] = s;
            }
        };
        if (b.lo[0] > 0) do_x_face(true);
        if (b.hi[0] + 1 < 2 * c_nx0) do_x_face(false);

        // loc 2 (y_lo, lj = -1) & loc 3 (y_up, lj = b.ny)
        auto do_y_face = [&](bool is_lo) {
            int lj_g = is_lo ? -1 : b.ny;
            int gj_g = b.lo[1] + lj_g;
            int jc_n = is_lo ? (gj_g >> 1) : ((gj_g >> 1) + 1);
            int jc_t = is_lo ? jc_n : (jc_n - 1);
            int lj_n = is_lo ? -1 : (b.ny + 1);
            int lj_t = is_lo ? -1 : b.ny;
            // u1 (n=y)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                int gi = b.lo[0] + i, gk = b.lo[2] + k;
                if (lookup_l1_patch(gi, gj_g, gk, l1_map) >= 0) continue;
                double w0[3], w2[3]; sc_quad_w(gi & 1, w0); sc_quad_w(gk & 1, w2);
                int ic = gi >> 1, kc = gk >> 1;
                double s = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wi = -1; wi <= 1; ++wi)
                    {
                        double wt = w0[wi + 1] * w2[wk + 1];
                        if (wt != 0.0) s += wt * u1_l0[sidx_g1(1, ic + wi, jc_n, kc + wk, c_nx0, c_ny0)];
                    }
                w1_l1[o1 + sidx_g1(1, i, lj_n, k, b.nx, b.ny)] = s;
            }
            // u0 (cubic in x, quad in z)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.nx + 1) * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % (b.nx + 1), k = idx / (b.nx + 1);
                int gi = b.lo[0] + i, gk = b.lo[2] + k;
                int ql, qj, qk, q;
                if (lookup_l1_side_owner(0, gi, b.lo[1] + lj_t, gk, l1_map, boxes, q, ql, qj, qk)) continue;
                int ic = (gi & 1) ? ((gi - 1) >> 1) : (gi >> 1);
                int ip = (gi & 1) ? 1 : 0;
                int kc = gk >> 1;
                double c0[4], w2[3]; sc_cubic_w(ip, c0); sc_quad_w(gk & 1, w2);
                double s = 0.0;
                for (int wk = -1; wk <= 1; ++wk)
                    for (int wi = -1; wi <= 2; ++wi)
                    {
                        double wt = c0[wi + 1] * w2[wk + 1];
                        if (wt != 0.0) s += wt * u0_l0[sidx_g1(0, ic + wi, jc_t, kc + wk, c_nx0, c_ny0)];
                    }
                w0_l1[o0 + sidx_g1(0, i, lj_t, k, b.nx, b.ny)] = s;
            }
            // u2 (quad in x, cubic in z)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * (b.nz + 1); idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                int gi = b.lo[0] + i, gk = b.lo[2] + k;
                int ql, qj, qk, q;
                if (lookup_l1_side_owner(2, gi, b.lo[1] + lj_t, gk, l1_map, boxes, q, ql, qj, qk)) continue;
                int ic = gi >> 1;
                int kc = (gk & 1) ? ((gk - 1) >> 1) : (gk >> 1);
                int kp = (gk & 1) ? 1 : 0;
                double w0[3], c2[4]; sc_quad_w(gi & 1, w0); sc_cubic_w(kp, c2);
                double s = 0.0;
                for (int wk = -1; wk <= 2; ++wk)
                    for (int wi = -1; wi <= 1; ++wi)
                    {
                        double wt = w0[wi + 1] * c2[wk + 1];
                        if (wt != 0.0) s += wt * u2_l0[sidx_g1(2, ic + wi, jc_t, kc + wk, c_nx0, c_ny0)];
                    }
                w2_l1[o2 + sidx_g1(2, i, lj_t, k, b.nx, b.ny)] = s;
            }
        };
        if (b.lo[1] > 0) do_y_face(true);
        if (b.hi[1] + 1 < 2 * c_ny0) do_y_face(false);

        // loc 4 (z_lo, lk = -1) & loc 5 (z_up, lk = b.nz)
        auto do_z_face = [&](bool is_lo) {
            int lk_g = is_lo ? -1 : b.nz;
            int gk_g = b.lo[2] + lk_g;
            int kc_n = is_lo ? (gk_g >> 1) : ((gk_g >> 1) + 1);
            int kc_t = is_lo ? kc_n : (kc_n - 1);
            int lk_n = is_lo ? -1 : (b.nz + 1);
            int lk_t = is_lo ? -1 : b.nz;
            // u2 (n=z)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                int gi = b.lo[0] + i, gj = b.lo[1] + j;
                if (lookup_l1_patch(gi, gj, gk_g, l1_map) >= 0) continue;
                double w0[3], w1[3]; sc_quad_w(gi & 1, w0); sc_quad_w(gj & 1, w1);
                int ic = gi >> 1, jc = gj >> 1;
                double s = 0.0;
                for (int wj = -1; wj <= 1; ++wj)
                    for (int wi = -1; wi <= 1; ++wi)
                    {
                        double wt = w0[wi + 1] * w1[wj + 1];
                        if (wt != 0.0) s += wt * u2_l0[sidx_g1(2, ic + wi, jc + wj, kc_n, c_nx0, c_ny0)];
                    }
                w2_l1[o2 + sidx_g1(2, i, j, lk_n, b.nx, b.ny)] = s;
            }
            // u0 (cubic in x, quad in y)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.nx + 1) * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % (b.nx + 1), j = idx / (b.nx + 1);
                int gi = b.lo[0] + i, gj = b.lo[1] + j;
                int ql, qj, qk, q;
                if (lookup_l1_side_owner(0, gi, gj, b.lo[2] + lk_t, l1_map, boxes, q, ql, qj, qk)) continue;
                int ic = (gi & 1) ? ((gi - 1) >> 1) : (gi >> 1);
                int ip = (gi & 1) ? 1 : 0;
                int jc = gj >> 1;
                double c0[4], w1[3]; sc_cubic_w(ip, c0); sc_quad_w(gj & 1, w1);
                double s = 0.0;
                for (int wj = -1; wj <= 1; ++wj)
                    for (int wi = -1; wi <= 2; ++wi)
                    {
                        double wt = c0[wi + 1] * w1[wj + 1];
                        if (wt != 0.0) s += wt * u0_l0[sidx_g1(0, ic + wi, jc + wj, kc_t, c_nx0, c_ny0)];
                    }
                w0_l1[o0 + sidx_g1(0, i, j, lk_t, b.nx, b.ny)] = s;
            }
            // u1 (quad in x, cubic in y)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * (b.ny + 1); idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                int gi = b.lo[0] + i, gj = b.lo[1] + j;
                int ql, qj, qk, q;
                if (lookup_l1_side_owner(1, gi, gj, b.lo[2] + lk_t, l1_map, boxes, q, ql, qj, qk)) continue;
                int ic = gi >> 1;
                int jc = (gj & 1) ? ((gj - 1) >> 1) : (gj >> 1);
                int jp = (gj & 1) ? 1 : 0;
                double w0[3], c1[4]; sc_quad_w(gi & 1, w0); sc_cubic_w(jp, c1);
                double s = 0.0;
                for (int wj = -1; wj <= 2; ++wj)
                    for (int wi = -1; wi <= 1; ++wi)
                    {
                        double wt = w0[wi + 1] * c1[wj + 1];
                        if (wt != 0.0) s += wt * u1_l0[sidx_g1(1, ic + wi, jc + wj, kc_t, c_nx0, c_ny0)];
                    }
                w1_l1[o1 + sidx_g1(1, i, j, lk_t, b.nx, b.ny)] = s;
            }
        };
        if (b.lo[2] > 0) do_z_face(true);
        if (b.hi[2] + 1 < 2 * c_nz0) do_z_face(false);
    }
}

__global__ void quad_cf_normal_side_l1_kernel(
    const double* w0_l1, const double* w1_l1, const double* w2_l1,
    double* u0_l1, double* u1_l1, double* u2_l1,
    const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    const double wm = -0.2, w0 = 2.0 / 3.0, wp = 8.0 / 15.0; // tangential side
    const double nm = -1.0 / 3.0, n0 = 1.0, np = 1.0 / 3.0;  // normal side
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        size_t o0 = b.off_s_g1[0], o1 = b.off_s_g1[1], o2 = b.off_s_g1[2];
        // loc 0 (x_lo) & loc 1 (x_up)
        if (b.lo[0] > 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
            {
                int j = idx % b.ny, k = idx / b.ny;
                if (lookup_l1_patch(b.lo[0] - 1, b.lo[1] + j, b.lo[2] + k, l1_map) < 0)
                    u0_l1[o0 + sidx_g1(0, -1, j, k, b.nx, b.ny)] =
                        np * w0_l1[o0 + sidx_g1(0, -1, j, k, b.nx, b.ny)] +
                        n0 * w0_l1[o0 + sidx_g1(0,  0, j, k, b.nx, b.ny)] +
                        nm * w0_l1[o0 + sidx_g1(0,  1, j, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.ny + 1) * b.nz; idx += gridDim.x * blockDim.x)
            {
                int j = idx % (b.ny + 1), k = idx / (b.ny + 1);
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(1, b.lo[0] - 1, b.lo[1] + j, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u1_l1[o1 + sidx_g1(1, -1, j, k, b.nx, b.ny)] =
                        wp * w1_l1[o1 + sidx_g1(1, -1, j, k, b.nx, b.ny)] +
                        w0 * w1_l1[o1 + sidx_g1(1,  0, j, k, b.nx, b.ny)] +
                        wm * w1_l1[o1 + sidx_g1(1,  1, j, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * (b.nz + 1); idx += gridDim.x * blockDim.x)
            {
                int j = idx % b.ny, k = idx / b.ny;
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(2, b.lo[0] - 1, b.lo[1] + j, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u2_l1[o2 + sidx_g1(2, -1, j, k, b.nx, b.ny)] =
                        wp * w2_l1[o2 + sidx_g1(2, -1, j, k, b.nx, b.ny)] +
                        w0 * w2_l1[o2 + sidx_g1(2,  0, j, k, b.nx, b.ny)] +
                        wm * w2_l1[o2 + sidx_g1(2,  1, j, k, b.nx, b.ny)];
            }
        }
        if (b.hi[0] + 1 < 2 * c_nx0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
            {
                int j = idx % b.ny, k = idx / b.ny;
                if (lookup_l1_patch(b.hi[0] + 1, b.lo[1] + j, b.lo[2] + k, l1_map) < 0)
                    u0_l1[o0 + sidx_g1(0, b.nx + 1, j, k, b.nx, b.ny)] =
                        np * w0_l1[o0 + sidx_g1(0, b.nx + 1, j, k, b.nx, b.ny)] +
                        n0 * w0_l1[o0 + sidx_g1(0, b.nx,     j, k, b.nx, b.ny)] +
                        nm * w0_l1[o0 + sidx_g1(0, b.nx - 1, j, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.ny + 1) * b.nz; idx += gridDim.x * blockDim.x)
            {
                int j = idx % (b.ny + 1), k = idx / (b.ny + 1);
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(1, b.hi[0] + 1, b.lo[1] + j, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u1_l1[o1 + sidx_g1(1, b.nx, j, k, b.nx, b.ny)] =
                        wp * w1_l1[o1 + sidx_g1(1, b.nx,     j, k, b.nx, b.ny)] +
                        w0 * w1_l1[o1 + sidx_g1(1, b.nx - 1, j, k, b.nx, b.ny)] +
                        wm * w1_l1[o1 + sidx_g1(1, b.nx - 2, j, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * (b.nz + 1); idx += gridDim.x * blockDim.x)
            {
                int j = idx % b.ny, k = idx / b.ny;
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(2, b.hi[0] + 1, b.lo[1] + j, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u2_l1[o2 + sidx_g1(2, b.nx, j, k, b.nx, b.ny)] =
                        wp * w2_l1[o2 + sidx_g1(2, b.nx,     j, k, b.nx, b.ny)] +
                        w0 * w2_l1[o2 + sidx_g1(2, b.nx - 1, j, k, b.nx, b.ny)] +
                        wm * w2_l1[o2 + sidx_g1(2, b.nx - 2, j, k, b.nx, b.ny)];
            }
        }
        // loc 2 (y_lo) & loc 3 (y_up)
        if (b.lo[1] > 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                if (lookup_l1_patch(b.lo[0] + i, b.lo[1] - 1, b.lo[2] + k, l1_map) < 0)
                    u1_l1[o1 + sidx_g1(1, i, -1, k, b.nx, b.ny)] =
                        np * w1_l1[o1 + sidx_g1(1, i, -1, k, b.nx, b.ny)] +
                        n0 * w1_l1[o1 + sidx_g1(1, i,  0, k, b.nx, b.ny)] +
                        nm * w1_l1[o1 + sidx_g1(1, i,  1, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.nx + 1) * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % (b.nx + 1), k = idx / (b.nx + 1);
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(0, b.lo[0] + i, b.lo[1] - 1, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u0_l1[o0 + sidx_g1(0, i, -1, k, b.nx, b.ny)] =
                        wp * w0_l1[o0 + sidx_g1(0, i, -1, k, b.nx, b.ny)] +
                        w0 * w0_l1[o0 + sidx_g1(0, i,  0, k, b.nx, b.ny)] +
                        wm * w0_l1[o0 + sidx_g1(0, i,  1, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * (b.nz + 1); idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(2, b.lo[0] + i, b.lo[1] - 1, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u2_l1[o2 + sidx_g1(2, i, -1, k, b.nx, b.ny)] =
                        wp * w2_l1[o2 + sidx_g1(2, i, -1, k, b.nx, b.ny)] +
                        w0 * w2_l1[o2 + sidx_g1(2, i,  0, k, b.nx, b.ny)] +
                        wm * w2_l1[o2 + sidx_g1(2, i,  1, k, b.nx, b.ny)];
            }
        }
        if (b.hi[1] + 1 < 2 * c_ny0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                if (lookup_l1_patch(b.lo[0] + i, b.hi[1] + 1, b.lo[2] + k, l1_map) < 0)
                    u1_l1[o1 + sidx_g1(1, i, b.ny + 1, k, b.nx, b.ny)] =
                        np * w1_l1[o1 + sidx_g1(1, i, b.ny + 1, k, b.nx, b.ny)] +
                        n0 * w1_l1[o1 + sidx_g1(1, i, b.ny,     k, b.nx, b.ny)] +
                        nm * w1_l1[o1 + sidx_g1(1, i, b.ny - 1, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.nx + 1) * b.nz; idx += gridDim.x * blockDim.x)
            {
                int i = idx % (b.nx + 1), k = idx / (b.nx + 1);
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(0, b.lo[0] + i, b.hi[1] + 1, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u0_l1[o0 + sidx_g1(0, i, b.ny, k, b.nx, b.ny)] =
                        wp * w0_l1[o0 + sidx_g1(0, i, b.ny,     k, b.nx, b.ny)] +
                        w0 * w0_l1[o0 + sidx_g1(0, i, b.ny - 1, k, b.nx, b.ny)] +
                        wm * w0_l1[o0 + sidx_g1(0, i, b.ny - 2, k, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * (b.nz + 1); idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, k = idx / b.nx;
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(2, b.lo[0] + i, b.hi[1] + 1, b.lo[2] + k, l1_map, boxes, q, ql, qj, qk))
                    u2_l1[o2 + sidx_g1(2, i, b.ny, k, b.nx, b.ny)] =
                        wp * w2_l1[o2 + sidx_g1(2, i, b.ny,     k, b.nx, b.ny)] +
                        w0 * w2_l1[o2 + sidx_g1(2, i, b.ny - 1, k, b.nx, b.ny)] +
                        wm * w2_l1[o2 + sidx_g1(2, i, b.ny - 2, k, b.nx, b.ny)];
            }
        }
        // loc 4 (z_lo) & loc 5 (z_up)
        if (b.lo[2] > 0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                if (lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] - 1, l1_map) < 0)
                    u2_l1[o2 + sidx_g1(2, i, j, -1, b.nx, b.ny)] =
                        np * w2_l1[o2 + sidx_g1(2, i, j, -1, b.nx, b.ny)] +
                        n0 * w2_l1[o2 + sidx_g1(2, i, j,  0, b.nx, b.ny)] +
                        nm * w2_l1[o2 + sidx_g1(2, i, j,  1, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.nx + 1) * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % (b.nx + 1), j = idx / (b.nx + 1);
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(0, b.lo[0] + i, b.lo[1] + j, b.lo[2] - 1, l1_map, boxes, q, ql, qj, qk))
                    u0_l1[o0 + sidx_g1(0, i, j, -1, b.nx, b.ny)] =
                        wp * w0_l1[o0 + sidx_g1(0, i, j, -1, b.nx, b.ny)] +
                        w0 * w0_l1[o0 + sidx_g1(0, i, j,  0, b.nx, b.ny)] +
                        wm * w0_l1[o0 + sidx_g1(0, i, j,  1, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * (b.ny + 1); idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(1, b.lo[0] + i, b.lo[1] + j, b.lo[2] - 1, l1_map, boxes, q, ql, qj, qk))
                    u1_l1[o1 + sidx_g1(1, i, j, -1, b.nx, b.ny)] =
                        wp * w1_l1[o1 + sidx_g1(1, i, j, -1, b.nx, b.ny)] +
                        w0 * w1_l1[o1 + sidx_g1(1, i, j,  0, b.nx, b.ny)] +
                        wm * w1_l1[o1 + sidx_g1(1, i, j,  1, b.nx, b.ny)];
            }
        }
        if (b.hi[2] + 1 < 2 * c_nz0)
        {
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                if (lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.hi[2] + 1, l1_map) < 0)
                    u2_l1[o2 + sidx_g1(2, i, j, b.nz + 1, b.nx, b.ny)] =
                        np * w2_l1[o2 + sidx_g1(2, i, j, b.nz + 1, b.nx, b.ny)] +
                        n0 * w2_l1[o2 + sidx_g1(2, i, j, b.nz,     b.nx, b.ny)] +
                        nm * w2_l1[o2 + sidx_g1(2, i, j, b.nz - 1, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < (b.nx + 1) * b.ny; idx += gridDim.x * blockDim.x)
            {
                int i = idx % (b.nx + 1), j = idx / (b.nx + 1);
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(0, b.lo[0] + i, b.lo[1] + j, b.hi[2] + 1, l1_map, boxes, q, ql, qj, qk))
                    u0_l1[o0 + sidx_g1(0, i, j, b.nz, b.nx, b.ny)] =
                        wp * w0_l1[o0 + sidx_g1(0, i, j, b.nz,     b.nx, b.ny)] +
                        w0 * w0_l1[o0 + sidx_g1(0, i, j, b.nz - 1, b.nx, b.ny)] +
                        wm * w0_l1[o0 + sidx_g1(0, i, j, b.nz - 2, b.nx, b.ny)];
            }
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * (b.ny + 1); idx += gridDim.x * blockDim.x)
            {
                int i = idx % b.nx, j = idx / b.nx;
                int q, ql, qj, qk;
                if (!lookup_l1_side_owner(1, b.lo[0] + i, b.lo[1] + j, b.hi[2] + 1, l1_map, boxes, q, ql, qj, qk))
                    u1_l1[o1 + sidx_g1(1, i, j, b.nz, b.nx, b.ny)] =
                        wp * w1_l1[o1 + sidx_g1(1, i, j, b.nz,     b.nx, b.ny)] +
                        w0 * w1_l1[o1 + sidx_g1(1, i, j, b.nz - 1, b.nx, b.ny)] +
                        wm * w1_l1[o1 + sidx_g1(1, i, j, b.nz - 2, b.nx, b.ny)];
            }
        }
    }
}


// ============================================================================
// Native IBAMR Pressure BCs (TRACTION at x_lo/x_up, Neumann on y/z) & Velocity Operator
// ============================================================================
__global__ void apply_pressure_phys_bc_l0_kernel(
    double* p_l0, const double* u0_l0, const double* u0_cur_l0,
    int nx, int ny, int nz, double h, double mu, double p_in, bool homogeneous)
{
    // x_lo / x_up: TRACTION BC -p + 2*mu*du_n/dx_n = g -> p_ghost = 2*gamma - p_int
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < ny * nz; idx += gridDim.x * blockDim.x)
    {
        int j = idx % ny, k = idx / ny;
        double u_new_0 = u0_l0[sidx_g1(0, 0, j, k, nx, ny)];
        double u_new_1 = u0_l0[sidx_g1(0, 1, j, k, nx, ny)];
        double u_new_2 = u0_l0[sidx_g1(0, 2, j, k, nx, ny)];
        double du_new_lo = (+2.0 * u_new_1 - 1.5 * u_new_0 - 0.5 * u_new_2) / h;
        double gamma_lo = mu * du_new_lo;
        if (!homogeneous)
        {
            double uc0 = u0_cur_l0[sidx_g0(0, 0, j, k, nx, ny)];
            double uc1 = u0_cur_l0[sidx_g0(0, 1, j, k, nx, ny)];
            double uc2 = u0_cur_l0[sidx_g0(0, 2, j, k, nx, ny)];
            double du_cur_lo = (+2.0 * uc1 - 1.5 * uc0 - 0.5 * uc2) / h;
            gamma_lo += mu * du_cur_lo + p_in; // gcoef = -p_in, -gcoef = +p_in
        }
        p_l0[cidx_g1(-1, j, k, nx, ny)] = 2.0 * gamma_lo - p_l0[cidx_g1(0, j, k, nx, ny)];

        double un_0 = u0_l0[sidx_g1(0, nx,     j, k, nx, ny)];
        double un_1 = u0_l0[sidx_g1(0, nx - 1, j, k, nx, ny)];
        double un_2 = u0_l0[sidx_g1(0, nx - 2, j, k, nx, ny)];
        double du_new_up = -(2.0 * un_1 - 1.5 * un_0 - 0.5 * un_2) / h;
        double gamma_up = mu * du_new_up;
        if (!homogeneous)
        {
            double uc0 = u0_cur_l0[sidx_g0(0, nx,     j, k, nx, ny)];
            double uc1 = u0_cur_l0[sidx_g0(0, nx - 1, j, k, nx, ny)];
            double uc2 = u0_cur_l0[sidx_g0(0, nx - 2, j, k, nx, ny)];
            double du_cur_up = -(2.0 * uc1 - 1.5 * uc0 - 0.5 * uc2) / h;
            gamma_up += mu * du_cur_up;
        }
        p_l0[cidx_g1(nx, j, k, nx, ny)] = 2.0 * gamma_up - p_l0[cidx_g1(nx - 1, j, k, nx, ny)];
    }
    // y_lo / y_up & z_lo / z_up: homogeneous Neumann dp/dn = 0 -> p_g = p_i
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx * nz; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, k = idx / nx;
        p_l0[cidx_g1(i, -1, k, nx, ny)] = p_l0[cidx_g1(i, 0, k, nx, ny)];
        p_l0[cidx_g1(i, ny, k, nx, ny)] = p_l0[cidx_g1(i, ny - 1, k, nx, ny)];
    }
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx * ny; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, j = idx / nx;
        p_l0[cidx_g1(i, j, -1, nx, ny)] = p_l0[cidx_g1(i, j, 0, nx, ny)];
        p_l0[cidx_g1(i, j, nz, nx, ny)] = p_l0[cidx_g1(i, j, nz - 1, nx, ny)];
    }
}

__global__ void apply_pressure_neumann_l1_kernel(
    double* p_l1, const Box3i* boxes, int n_l1)
{
    int dom_hi1 = 2 * c_ny0 - 1, dom_hi2 = 2 * c_nz0 - 1;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t o1 = b.off_c_g1;
        if (b.lo[1] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; p_l1[o1 + cidx_g1(i, -1, k, b.nx, b.ny)] = p_l1[o1 + cidx_g1(i, 0, k, b.nx, b.ny)]; }
        if (b.hi[1] == dom_hi1)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; p_l1[o1 + cidx_g1(i, b.ny, k, b.nx, b.ny)] = p_l1[o1 + cidx_g1(i, b.ny - 1, k, b.nx, b.ny)]; }
        if (b.lo[2] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; p_l1[o1 + cidx_g1(i, j, -1, b.nx, b.ny)] = p_l1[o1 + cidx_g1(i, j, 0, b.nx, b.ny)]; }
        if (b.hi[2] == dom_hi2)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; p_l1[o1 + cidx_g1(i, j, b.nz, b.nx, b.ny)] = p_l1[o1 + cidx_g1(i, j, b.nz - 1, b.nx, b.ny)]; }
    }
}

// Physical pseudo-traction Dirichlet p_ghost = -p_int at x_lo/x_up ONLY used inside Phi projection preconditioner!
__global__ void apply_phi_phys_bc_l0_kernel(double* phi_l0, int nx, int ny, int nz)
{
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < ny * nz; idx += gridDim.x * blockDim.x)
    {
        int j = idx % ny, k = idx / ny;
        phi_l0[cidx_g1(-1, j, k, nx, ny)] = -phi_l0[cidx_g1(0, j, k, nx, ny)];
        phi_l0[cidx_g1(nx, j, k, nx, ny)] = -phi_l0[cidx_g1(nx - 1, j, k, nx, ny)];
    }
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx * nz; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, k = idx / nx;
        phi_l0[cidx_g1(i, -1, k, nx, ny)] = phi_l0[cidx_g1(i, 0, k, nx, ny)];
        phi_l0[cidx_g1(i, ny, k, nx, ny)] = phi_l0[cidx_g1(i, ny - 1, k, nx, ny)];
    }
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx * ny; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nx, j = idx / nx;
        phi_l0[cidx_g1(i, j, -1, nx, ny)] = phi_l0[cidx_g1(i, j, 0, nx, ny)];
        phi_l0[cidx_g1(i, j, nz, nx, ny)] = phi_l0[cidx_g1(i, j, nz - 1, nx, ny)];
    }
}

// CartExtrapPhysBdryOp("LINEAR") for cell-centered L0 g=1 ghosts (codim 1 + 2 + 3) before C-F quadratic interp
__global__ void cart_extrap_phys_cell_l0_kernel(double* p, int nx, int ny, int nz)
{
    int sx = nx + 2, sy = ny + 2, sz = nz + 2;
    int tot = sx * sy * sz;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int li = (idx % sx) - 1;
        int t = idx / sx;
        int lj = (t % sy) - 1;
        int lk = (t / sy) - 1;
        if (li >= 0 && li < nx && lj >= 0 && lj < ny && lk >= 0 && lk < nz) continue;
        int bi = li, si = 0;
        if (li < 0) { bi = 0; si = 1; } else if (li >= nx) { bi = nx - 1; si = -1; }
        int bj = lj, sj = 0;
        if (lj < 0) { bj = 0; sj = 1; } else if (lj >= ny) { bj = ny - 1; sj = -1; }
        int bk = lk, sk = 0;
        if (lk < 0) { bk = 0; sk = 1; } else if (lk >= nz) { bk = nz - 1; sk = -1; }
        double v0 = p[cidx_g1(bi, bj, bk, nx, ny)];
        double val = v0;
        if (si != 0) val += (v0 - p[cidx_g1(bi + si, bj, bk, nx, ny)]) * abs(li - bi);
        if (sj != 0) val += (v0 - p[cidx_g1(bi, bj + sj, bk, nx, ny)]) * abs(lj - bj);
        if (sk != 0) val += (v0 - p[cidx_g1(bi, bj, bk + sk, nx, ny)]) * abs(lk - bk);
        p[cidx_g1(li, lj, lk, nx, ny)] = val;
    }
}

// Match HierarchyMathOps::vc_laplace xeqScheduleOutersideRestriction (0.25 * sum of 4 L1 boundary faces -> L0 C-F face cr ^ cl)
__global__ void outerside_restrict_cf_l1_to_l0_side_kernel(
    const double* u_l1, double* u_l0, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const size_t o1 = b.off_s_g0[axis];
        const int cnx = b.nx >> 1, cny = b.ny >> 1, cnz = b.nz >> 1;
        if (axis == 0)
        {
            int tot = cny * cnz;
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
            {
                int jc = idx % cny, kc = idx / cny;
                int gj = b.lo[1] + 2 * jc, gk = b.lo[2] + 2 * kc;
                if (lookup_l1_patch(b.lo[0] - 1, gj, gk, l1_map) < 0)
                {
                    double s = 0.0;
                    for (int dk = 0; dk < 2; ++dk)
                        for (int dj = 0; dj < 2; ++dj)
                            s += u_l1[o1 + sidx_g0(0, 0, 2 * jc + dj, 2 * kc + dk, b.nx, b.ny)];
                    u_l0[sidx_g0(0, b.lo[0] >> 1, gj >> 1, gk >> 1, c_nx0, c_ny0)] = 0.25 * s;
                }
                int gi_up = b.lo[0] + b.nx;
                if (lookup_l1_patch(gi_up, gj, gk, l1_map) < 0)
                {
                    double s = 0.0;
                    for (int dk = 0; dk < 2; ++dk)
                        for (int dj = 0; dj < 2; ++dj)
                            s += u_l1[o1 + sidx_g0(0, b.nx, 2 * jc + dj, 2 * kc + dk, b.nx, b.ny)];
                    u_l0[sidx_g0(0, gi_up >> 1, gj >> 1, gk >> 1, c_nx0, c_ny0)] = 0.25 * s;
                }
            }
        }
        else if (axis == 1)
        {
            int tot = cnx * cnz;
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
            {
                int ic = idx % cnx, kc = idx / cnx;
                int gi = b.lo[0] + 2 * ic, gk = b.lo[2] + 2 * kc;
                if (lookup_l1_patch(gi, b.lo[1] - 1, gk, l1_map) < 0)
                {
                    double s = 0.0;
                    for (int dk = 0; dk < 2; ++dk)
                        for (int di = 0; di < 2; ++di)
                            s += u_l1[o1 + sidx_g0(1, 2 * ic + di, 0, 2 * kc + dk, b.nx, b.ny)];
                    u_l0[sidx_g0(1, gi >> 1, b.lo[1] >> 1, gk >> 1, c_nx0, c_ny0)] = 0.25 * s;
                }
                int gj_up = b.lo[1] + b.ny;
                if (lookup_l1_patch(gi, gj_up, gk, l1_map) < 0)
                {
                    double s = 0.0;
                    for (int dk = 0; dk < 2; ++dk)
                        for (int di = 0; di < 2; ++di)
                            s += u_l1[o1 + sidx_g0(1, 2 * ic + di, b.ny, 2 * kc + dk, b.nx, b.ny)];
                    u_l0[sidx_g0(1, gi >> 1, gj_up >> 1, gk >> 1, c_nx0, c_ny0)] = 0.25 * s;
                }
            }
        }
        else
        {
            int tot = cnx * cny;
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
            {
                int ic = idx % cnx, jc = idx / cnx;
                int gi = b.lo[0] + 2 * ic, gj = b.lo[1] + 2 * jc;
                if (lookup_l1_patch(gi, gj, b.lo[2] - 1, l1_map) < 0)
                {
                    double s = 0.0;
                    for (int dj = 0; dj < 2; ++dj)
                        for (int di = 0; di < 2; ++di)
                            s += u_l1[o1 + sidx_g0(2, 2 * ic + di, 2 * jc + dj, 0, b.nx, b.ny)];
                    u_l0[sidx_g0(2, gi >> 1, gj >> 1, b.lo[2] >> 1, c_nx0, c_ny0)] = 0.25 * s;
                }
                int gk_up = b.lo[2] + b.nz;
                if (lookup_l1_patch(gi, gj, gk_up, l1_map) < 0)
                {
                    double s = 0.0;
                    for (int dj = 0; dj < 2; ++dj)
                        for (int di = 0; di < 2; ++di)
                            s += u_l1[o1 + sidx_g0(2, 2 * ic + di, 2 * jc + dj, b.nz, b.nx, b.ny)];
                    u_l0[sidx_g0(2, gi >> 1, gj >> 1, gk_up >> 1, c_nx0, c_ny0)] = 0.25 * s;
                }
            }
        }
    }
}

__global__ void add_scaled_prod_side_g1_to_g0_kernel(
    double* dst_g0, const double* a_g0, const double* b_g1, double scale,
    int axis, int nx, int ny, int nz, bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o0 = is_l1 ? b.off_s_g0[axis] : 0;
        size_t o1 = is_l1 ? b.off_s_g1[axis] : 0;
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        int tot = nxa * nya * nza;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
            if (is_l1)
            {
                if (axis == 0 && i == b.nx && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
                if (axis == 1 && j == b.ny && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
                if (axis == 2 && k == b.nz && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            }
            dst_g0[o0 + idx] += scale * a_g0[o0 + idx] * b_g1[o1 + sidx_g1(axis, i, j, k, b.nx, b.ny)];
        }
    }
}

__global__ void copy_dirichlet_bc_g0_kernel(
    const double* x1, const double* x2, double* y1, double* y2,
    int nx, int ny, int nz, bool is_l1, const Box3i* boxes, int n_l1)
{
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o1 = is_l1 ? b.off_s_g0[1] : 0;
        size_t o2 = is_l1 ? b.off_s_g0[2] : 0;
        if (b.lo[1] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; size_t off = o1 + sidx_g0(1, i, 0, k, b.nx, b.ny); y1[off] = x1[off]; }
        if (b.lo[1] + b.ny == gy_max)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; size_t off = o1 + sidx_g0(1, i, b.ny, k, b.nx, b.ny); y1[off] = x1[off]; }
        if (b.lo[2] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; size_t off = o2 + sidx_g0(2, i, j, 0, b.nx, b.ny); y2[off] = x2[off]; }
        if (b.lo[2] + b.nz == gz_max)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; size_t off = o2 + sidx_g0(2, i, j, b.nz, b.nx, b.ny); y2[off] = x2[off]; }
    }
}

template <typename RealT>
__global__ void enforce_homog_phys_bcs_side_kernel(
    RealT* u0, RealT* u1, RealT* u2, int nx, int ny, int nz, bool is_l1, const Box3i* boxes, int n_patches)
{
    int gx_max = is_l1 ? 2 * c_nx0 : nx;
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;

    for (int p = blockIdx.y; p < (is_l1 ? n_patches : 1); p += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[p];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t off0 = is_l1 ? b.off_s_g1[0] : 0;
        size_t off1 = is_l1 ? b.off_s_g1[1] : 0;
        size_t off2 = is_l1 ? b.off_s_g1[2] : 0;

        int sx = b.nx + 3, sy = b.ny + 3, sz = b.nz + 3;
        int tot = sx * sy * sz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int li = (lin % sx) - 1;
            int t = lin / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;

            int gi = b.lo[0] + li;
            int gj = b.lo[1] + lj;
            int gk = b.lo[2] + lk;

            if (lj >= 0 && lj <= b.ny && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
            {
                if (gj == 0 || gj == gy_max)
                    u1[off1 + sidx_g1(1, li, lj, lk, b.nx, b.ny)] = RealT(0.0);
            }
            if (lk >= 0 && lk <= b.nz && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
            {
                if (gk == 0 || gk == gz_max)
                    u2[off2 + sidx_g1(2, li, lj, lk, b.nx, b.ny)] = RealT(0.0);
            }

            if (gi == -1 && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz && li == -1)
                u0[off0 + sidx_g1(0, -1, lj, lk, b.nx, b.ny)] = u0[off0 + sidx_g1(0, 1, lj, lk, b.nx, b.ny)];
            if (gi == gx_max + 1 && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz && li == b.nx + 1)
                u0[off0 + sidx_g1(0, b.nx + 1, lj, lk, b.nx, b.ny)] = u0[off0 + sidx_g1(0, b.nx - 1, lj, lk, b.nx, b.ny)];

            if (gj == -1 && lj == -1 && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
            {
                u0[off0 + sidx_g1(0, li, -1, lk, b.nx, b.ny)] = -u0[off0 + sidx_g1(0, li, 0, lk, b.nx, b.ny)];
                u2[off2 + sidx_g1(2, li, -1, lk, b.nx, b.ny)] = -u2[off2 + sidx_g1(2, li, 0, lk, b.nx, b.ny)];
            }
            if (gj == -1 && lj == -1 && li == b.nx && gi <= gx_max && lk >= 0 && lk < b.nz)
                u0[off0 + sidx_g1(0, b.nx, -1, lk, b.nx, b.ny)] = -u0[off0 + sidx_g1(0, b.nx, 0, lk, b.nx, b.ny)];
            if (gj == gy_max && lj == b.ny && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
            {
                u0[off0 + sidx_g1(0, li, b.ny, lk, b.nx, b.ny)] = -u0[off0 + sidx_g1(0, li, b.ny - 1, lk, b.nx, b.ny)];
                u2[off2 + sidx_g1(2, li, b.ny, lk, b.nx, b.ny)] = -u2[off2 + sidx_g1(2, li, b.ny - 1, lk, b.nx, b.ny)];
            }
            if (gj == gy_max && lj == b.ny && li == b.nx && gi <= gx_max && lk >= 0 && lk < b.nz)
                u0[off0 + sidx_g1(0, b.nx, b.ny, lk, b.nx, b.ny)] = -u0[off0 + sidx_g1(0, b.nx, b.ny - 1, lk, b.nx, b.ny)];

            if (gk == -1 && lk == -1 && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
            {
                u0[off0 + sidx_g1(0, li, lj, -1, b.nx, b.ny)] = u0[off0 + sidx_g1(0, li, lj, 0, b.nx, b.ny)];
                u1[off1 + sidx_g1(1, li, lj, -1, b.nx, b.ny)] = u1[off1 + sidx_g1(1, li, lj, 0, b.nx, b.ny)];
            }
            if (gk == -1 && lk == -1 && li == b.nx && gi <= gx_max && lj >= 0 && lj < b.ny)
                u0[off0 + sidx_g1(0, b.nx, lj, -1, b.nx, b.ny)] = u0[off0 + sidx_g1(0, b.nx, lj, 0, b.nx, b.ny)];

            if (gk == gz_max && lk == b.nz && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
            {
                u0[off0 + sidx_g1(0, li, lj, b.nz, b.nx, b.ny)] = -u0[off0 + sidx_g1(0, li, lj, b.nz - 1, b.nx, b.ny)];
                u1[off1 + sidx_g1(1, li, lj, b.nz, b.nx, b.ny)] = -u1[off1 + sidx_g1(1, li, lj, b.nz - 1, b.nx, b.ny)];
            }
            if (gk == gz_max && lk == b.nz && li == b.nx && gi <= gx_max && lj >= 0 && lj < b.ny)
                u0[off0 + sidx_g1(0, b.nx, lj, b.nz, b.nx, b.ny)] = -u0[off0 + sidx_g1(0, b.nx, lj, b.nz - 1, b.nx, b.ny)];

            if (gi == -1 && li == -1 && lj >= 0 && lj <= b.ny && lk >= 0 && lk < b.nz)
                u1[off1 + sidx_g1(1, -1, lj, lk, b.nx, b.ny)] = -u1[off1 + sidx_g1(1, 0, lj, lk, b.nx, b.ny)];
            if (gi == -1 && li == -1 && lj >= 0 && lj < b.ny && lk >= 0 && lk <= b.nz)
                u2[off2 + sidx_g1(2, -1, lj, lk, b.nx, b.ny)] = -u2[off2 + sidx_g1(2, 0, lj, lk, b.nx, b.ny)];
            if (gi == gx_max && li == b.nx && lj >= 0 && lj <= b.ny && lk >= 0 && lk < b.nz)
                u1[off1 + sidx_g1(1, b.nx, lj, lk, b.nx, b.ny)] = -u1[off1 + sidx_g1(1, b.nx - 1, lj, lk, b.nx, b.ny)];
            if (gi == gx_max && li == b.nx && lj >= 0 && lj < b.ny && lk >= 0 && lk <= b.nz)
                u2[off2 + sidx_g1(2, b.nx, lj, lk, b.nx, b.ny)] = -u2[off2 + sidx_g1(2, b.nx - 1, lj, lk, b.nx, b.ny)];

            if (gj == -1 && lj == -1 && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
                u1[off1 + sidx_g1(1, li, -1, lk, b.nx, b.ny)] = -u1[off1 + sidx_g1(1, li, 1, lk, b.nx, b.ny)];
            if (gj == gy_max + 1 && lj == b.ny + 1 && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
                u1[off1 + sidx_g1(1, li, b.ny + 1, lk, b.nx, b.ny)] = -u1[off1 + sidx_g1(1, li, b.ny - 1, lk, b.nx, b.ny)];

            if (gk == -1 && lk == -1 && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
                u2[off2 + sidx_g1(2, li, lj, -1, b.nx, b.ny)] = -u2[off2 + sidx_g1(2, li, lj, 1, b.nx, b.ny)];
            if (gk == gz_max + 1 && lk == b.nz + 1 && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
                u2[off2 + sidx_g1(2, li, lj, b.nz + 1, b.nx, b.ny)] = -u2[off2 + sidx_g1(2, li, lj, b.nz - 1, b.nx, b.ny)];
        }
    }
}

template <typename RealT>
__global__ void enforce_homog_phys_bcs_cell_kernel(
    RealT* p, int nx, int ny, int nz, bool is_l1, const Box3i* boxes, int n_patches)
{
    int gx_max = is_l1 ? 2 * c_nx0 : nx;
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;

    for (int pat = blockIdx.y; pat < (is_l1 ? n_patches : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t off = is_l1 ? b.off_c_g1 : 0;

        int sx = b.nx + 2, sy = b.ny + 2, sz = b.nz + 2;
        int tot = sx * sy * sz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int li = (lin % sx) - 1;
            int t = lin / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            if (gi == -1 && li == -1 && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz)
                p[off + cidx_g1(-1, lj, lk, b.nx, b.ny)] = -p[off + cidx_g1(0, lj, lk, b.nx, b.ny)];
            if (gi == gx_max && li == b.nx && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz)
                p[off + cidx_g1(b.nx, lj, lk, b.nx, b.ny)] = -p[off + cidx_g1(b.nx - 1, lj, lk, b.nx, b.ny)];
            if (gj == -1 && lj == -1 && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
                p[off + cidx_g1(li, -1, lk, b.nx, b.ny)] = p[off + cidx_g1(li, 0, lk, b.nx, b.ny)];
            if (gj == gy_max && lj == b.ny && li >= 0 && li < b.nx && lk >= 0 && lk < b.nz)
                p[off + cidx_g1(li, b.ny, lk, b.nx, b.ny)] = p[off + cidx_g1(li, b.ny - 1, lk, b.nx, b.ny)];
            if (gk == -1 && lk == -1 && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
                p[off + cidx_g1(li, lj, -1, b.nx, b.ny)] = p[off + cidx_g1(li, lj, 0, b.nx, b.ny)];
            if (gk == gz_max && lk == b.nz && li >= 0 && li < b.nx && lj >= 0 && lj < b.ny)
                p[off + cidx_g1(li, lj, b.nz, b.nx, b.ny)] = p[off + cidx_g1(li, lj, b.nz - 1, b.nx, b.ny)];
        }
    }
}

template <typename RealT>
__global__ void fill_l1_cell_ghosts_kernel(
    RealT* p_l1, const RealT* p_l0, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        int sx = b.nx + 2, sy = b.ny + 2, sz = b.nz + 2;
        int tot = sx * sy * sz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int li = (lin % sx) - 1;
            int t = lin / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;
            if (li >= 0 && li < b.nx && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz) continue;

            int gi = b.lo[0] + li;
            int gj = b.lo[1] + lj;
            int gk = b.lo[2] + lk;
            if (gi < 0 || gi >= 2 * c_nx0 || gj < 0 || gj >= 2 * c_ny0 || gk < 0 || gk >= 2 * c_nz0) continue;

            int q = lookup_l1_patch(gi, gj, gk, l1_map);
            RealT val = RealT(0.0);
            if (q >= 0)
            {
                const Box3i qb = boxes[q];
                val = p_l1[qb.off_c_g1 + cidx_g1(gi - qb.lo[0], gj - qb.lo[1], gk - qb.lo[2], qb.nx, qb.ny)];
            }
            else
            {
                int i0 = gi >> 1, j0 = gj >> 1, k0 = gk >> 1;
                RealT ox = (gi & 1) ? RealT(0.25) : RealT(-0.25);
                RealT oy = (gj & 1) ? RealT(0.25) : RealT(-0.25);
                RealT oz = (gk & 1) ? RealT(0.25) : RealT(-0.25);
                RealT pc = p_l0[cidx_g1(i0, j0, k0, c_nx0, c_ny0)];
                RealT gx = RealT(0.5) * (p_l0[cidx_g1(i0 + 1, j0, k0, c_nx0, c_ny0)] - p_l0[cidx_g1(i0 - 1, j0, k0, c_nx0, c_ny0)]);
                RealT gy = RealT(0.5) * (p_l0[cidx_g1(i0, j0 + 1, k0, c_nx0, c_ny0)] - p_l0[cidx_g1(i0, j0 - 1, k0, c_nx0, c_ny0)]);
                RealT gz = RealT(0.5) * (p_l0[cidx_g1(i0, j0, k0 + 1, c_nx0, c_ny0)] - p_l0[cidx_g1(i0, j0, k0 - 1, c_nx0, c_ny0)]);
                val = pc + ox * gx + oy * gy + oz * gz;
            }
            p_l1[b.off_c_g1 + cidx_g1(li, lj, lk, b.nx, b.ny)] = val;
        }
    }
}

template <typename RealT>
__global__ void fill_l1_side_ghosts_kernel(
    RealT* u_l1_a, const RealT* u_l0_a, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis)
{
    int gx_lim = 2 * c_nx0 + (axis == 0 ? 1 : 0);
    int gy_lim = 2 * c_ny0 + (axis == 1 ? 1 : 0);
    int gz_lim = 2 * c_nz0 + (axis == 2 ? 1 : 0);

    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        int nxa = b.nx + (axis == 0 ? 1 : 0);
        int nya = b.ny + (axis == 1 ? 1 : 0);
        int nza = b.nz + (axis == 2 ? 1 : 0);
        int sx = nxa + 2, sy = nya + 2, sz = nza + 2;
        int tot = sx * sy * sz;

        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int li = (lin % sx) - 1;
            int t = lin / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;

            int gi = b.lo[0] + li;
            int gj = b.lo[1] + lj;
            int gk = b.lo[2] + lk;

            if (gi < 0 || gi >= gx_lim || gj < 0 || gj >= gy_lim || gk < 0 || gk >= gz_lim) continue;

            bool is_ghost = (li < 0 || li >= nxa || lj < 0 || lj >= nya || lk < 0 || lk >= nza);
            bool is_shared_upper = false;
            if (axis == 0 && li == b.nx) is_shared_upper = (lookup_l1_patch(gi, gj, gk, l1_map) >= 0);
            if (axis == 1 && lj == b.ny) is_shared_upper = (lookup_l1_patch(gi, gj, gk, l1_map) >= 0);
            if (axis == 2 && lk == b.nz) is_shared_upper = (lookup_l1_patch(gi, gj, gk, l1_map) >= 0);

            if (!is_ghost && !is_shared_upper) continue;

            int owner = lookup_l1_patch(gi, gj, gk, l1_map);
            if (owner < 0)
            {
                int ci = gi - (axis == 0 ? 1 : 0);
                int cj = gj - (axis == 1 ? 1 : 0);
                int ck = gk - (axis == 2 ? 1 : 0);
                if (gi > 0 || axis != 0)
                    owner = lookup_l1_patch(ci, cj, ck, l1_map);
            }

            RealT val = RealT(0.0);
            if (owner >= 0 && owner != pat)
            {
                const Box3i qb = boxes[owner];
                int qli = gi - qb.lo[0], qlj = gj - qb.lo[1], qlk = gk - qb.lo[2];
                int qnxa = qb.nx + (axis == 0 ? 1 : 0);
                int qnya = qb.ny + (axis == 1 ? 1 : 0);
                int qnza = qb.nz + (axis == 2 ? 1 : 0);
                if (qli >= 0 && qli < qnxa && qlj >= 0 && qlj < qnya && qlk >= 0 && qlk < qnza)
                {
                    val = u_l1_a[qb.off_s_g1[axis] + sidx_g1(axis, qli, qlj, qlk, qb.nx, qb.ny)];
                    u_l1_a[b.off_s_g1[axis] + sidx_g1(axis, li, lj, lk, b.nx, b.ny)] = val;
                    continue;
                }
            }
            if (owner == pat) continue;

            int g_coord[3] = {gi, gj, gk};
            int c_idx[3];
            RealT off_c[3];
            for (int d = 0; d < 3; ++d)
            {
                c_idx[d] = g_coord[d] >> 1;
                if (d == axis) off_c[d] = (g_coord[d] & 1) ? RealT(0.5) : RealT(0.0);
                else off_c[d] = (g_coord[d] & 1) ? RealT(0.25) : RealT(-0.25);
            }
            int i0 = c_idx[0], j0 = c_idx[1], k0 = c_idx[2];
            if (i0 < 0 || i0 > c_nx0 || j0 < 0 || j0 > c_ny0 || k0 < 0 || k0 > c_nz0) continue;
            RealT vc = u_l0_a[sidx_g1(axis, i0, j0, k0, c_nx0, c_ny0)];
            if (axis == 0 && (gi & 1))
                val = RealT(0.5) * (vc + u_l0_a[sidx_g1(0, i0 + 1, j0, k0, c_nx0, c_ny0)]);
            else if (axis == 1 && (gj & 1))
                val = RealT(0.5) * (vc + u_l0_a[sidx_g1(1, i0, j0 + 1, k0, c_nx0, c_ny0)]);
            else if (axis == 2 && (gk & 1))
                val = RealT(0.5) * (vc + u_l0_a[sidx_g1(2, i0, j0, k0 + 1, c_nx0, c_ny0)]);
            else
            {
                RealT dx0_g = RealT(0.5) * (u_l0_a[sidx_g1(axis, i0 + 1, j0, k0, c_nx0, c_ny0)] -
                                            u_l0_a[sidx_g1(axis, i0 - 1, j0, k0, c_nx0, c_ny0)]);
                RealT dy0_g = RealT(0.5) * (u_l0_a[sidx_g1(axis, i0, j0 + 1, k0, c_nx0, c_ny0)] -
                                            u_l0_a[sidx_g1(axis, i0, j0 - 1, k0, c_nx0, c_ny0)]);
                RealT dz0_g = RealT(0.5) * (u_l0_a[sidx_g1(axis, i0, j0, k0 + 1, c_nx0, c_ny0)] -
                                            u_l0_a[sidx_g1(axis, i0, j0, k0 - 1, c_nx0, c_ny0)]);
                if (axis == 0) dx0_g = RealT(0.0);
                if (axis == 1) dy0_g = RealT(0.0);
                if (axis == 2) dz0_g = RealT(0.0);
                val = vc + off_c[0] * dx0_g + off_c[1] * dy0_g + off_c[2] * dz0_g;
            }
            u_l1_a[b.off_s_g1[axis] + sidx_g1(axis, li, lj, lk, b.nx, b.ny)] = val;
        }
    }
}


// Apply constant-mu (-0.5*c_mu0 or +0.5*c_mu0) staggered vc_laplace + C_sc*u + grad(p), matching VCStaggeredStokesOperator::apply and viscous_rhs!
// Note: when compute_grad_p == true, computes A_U = C_sc*u + d_ec*Lap(u) + Grad(p), and enforces copyDataAtDirichletBoundaries (u1 on y_lo/y_up, u2 on z_lo/z_up).
// When compute_grad_p == false (viscous_rhs), computes u_rhs = d_ec*Lap(u) + c_sc*u (with c_sc = rho/dt, d_ec = +0.5*c_mu0).
__global__ void apply_momentum_stokes_kernel(
    const double* u0, const double* u1, const double* u2,
    const double* p_g1,
    const double* c0, const double* c1, const double* c2,
    double d_ec, double c_scale,
    double* out0, double* out1, double* out2,
    int nx, int ny, int nz, double h,
    bool compute_grad_p, bool dirichlet_copy_bc, bool zero_covered_l0,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    double inv_h = 1.0 / h;
    double inv_h2 = inv_h * inv_h;
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;

    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o_c1 = is_l1 ? b.off_c_g1 : 0;
        size_t o0_g0 = is_l1 ? b.off_s_g0[0] : 0;
        size_t o1_g0 = is_l1 ? b.off_s_g0[1] : 0;
        size_t o2_g0 = is_l1 ? b.off_s_g0[2] : 0;
        size_t o0_g1 = is_l1 ? b.off_s_g1[0] : 0;
        size_t o1_g1 = is_l1 ? b.off_s_g1[1] : 0;
        size_t o2_g1 = is_l1 ? b.off_s_g1[2] : 0;

        // Axis 0
        int nx0 = b.nx + 1, ny0 = b.ny, nz0 = b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx0 * ny0 * nz0; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx0, t = idx / nx0, j = t % ny0, k = t / ny0;
            if (is_l1 && i == b.nx && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            if (!is_l1 && zero_covered_l0)
            {
                bool cr = (i < c_nx0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (i > 0) && (l1_map[(i - 1) + c_nx0 * (j + c_ny0 * k)] >= 0);
                if (cr && cl) { out0[idx] = 0.0; continue; }
            }
            size_t s = o0_g1 + sidx_g1(0, i, j, k, b.nx, b.ny);
            double uc = u0[s];
            double t_xx = 2.0 * (u0[o0_g1 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - 2.0 * uc + u0[o0_g1 + sidx_g1(0, i - 1, j, k, b.nx, b.ny)]);
            double t_yy = u0[o0_g1 + sidx_g1(0, i, j + 1, k, b.nx, b.ny)] - 2.0 * uc + u0[o0_g1 + sidx_g1(0, i, j - 1, k, b.nx, b.ny)];
            double t_xy = (u1[o1_g1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[o1_g1 + sidx_g1(1, i - 1, j + 1, k, b.nx, b.ny)]) -
                          (u1[o1_g1 + sidx_g1(1, i, j,     k, b.nx, b.ny)] - u1[o1_g1 + sidx_g1(1, i - 1, j,     k, b.nx, b.ny)]);
            double t_zz = u0[o0_g1 + sidx_g1(0, i, j, k + 1, b.nx, b.ny)] - 2.0 * uc + u0[o0_g1 + sidx_g1(0, i, j, k - 1, b.nx, b.ny)];
            double t_xz = (u2[o2_g1 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[o2_g1 + sidx_g1(2, i - 1, j, k + 1, b.nx, b.ny)]) -
                          (u2[o2_g1 + sidx_g1(2, i, j, k,     b.nx, b.ny)] - u2[o2_g1 + sidx_g1(2, i - 1, j, k,     b.nx, b.ny)]);
            double lap = (t_xx + t_yy + t_xy + t_zz + t_xz) * inv_h2;
            double val = c_scale * c0[s] * uc + d_ec * lap;
            if (compute_grad_p)
                val += (p_g1[o_c1 + cidx_g1(i, j, k, b.nx, b.ny)] - p_g1[o_c1 + cidx_g1(i - 1, j, k, b.nx, b.ny)]) * inv_h;
            out0[o0_g0 + idx] = val;
        }
        // Axis 1
        int nx1 = b.nx, ny1 = b.ny + 1, nz1 = b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx1 * ny1 * nz1; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx1, t = idx / nx1, j = t % ny1, k = t / ny1;
            if (is_l1 && j == b.ny && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            if (!is_l1 && zero_covered_l0)
            {
                bool cr = (j < c_ny0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (j > 0) && (l1_map[i + c_nx0 * ((j - 1) + c_ny0 * k)] >= 0);
                if (cr && cl) { out1[idx] = 0.0; continue; }
            }
            size_t s = o1_g1 + sidx_g1(1, i, j, k, b.nx, b.ny);
            int gj = b.lo[1] + j;
            if (dirichlet_copy_bc && (gj == 0 || gj == gy_max)) { out1[o1_g0 + idx] = u1[s]; continue; }
            double uc = u1[s];
            double t_yy = 2.0 * (u1[o1_g1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - 2.0 * uc + u1[o1_g1 + sidx_g1(1, i, j - 1, k, b.nx, b.ny)]);
            double t_xx = u1[o1_g1 + sidx_g1(1, i + 1, j, k, b.nx, b.ny)] - 2.0 * uc + u1[o1_g1 + sidx_g1(1, i - 1, j, k, b.nx, b.ny)];
            double t_yx = (u0[o0_g1 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[o0_g1 + sidx_g1(0, i + 1, j - 1, k, b.nx, b.ny)]) -
                          (u0[o0_g1 + sidx_g1(0, i,     j, k, b.nx, b.ny)] - u0[o0_g1 + sidx_g1(0, i,     j - 1, k, b.nx, b.ny)]);
            double t_zz = u1[o1_g1 + sidx_g1(1, i, j, k + 1, b.nx, b.ny)] - 2.0 * uc + u1[o1_g1 + sidx_g1(1, i, j, k - 1, b.nx, b.ny)];
            double t_yz = (u2[o2_g1 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[o2_g1 + sidx_g1(2, i, j - 1, k + 1, b.nx, b.ny)]) -
                          (u2[o2_g1 + sidx_g1(2, i, j, k,     b.nx, b.ny)] - u2[o2_g1 + sidx_g1(2, i, j - 1, k,     b.nx, b.ny)]);
            double lap = (t_yy + t_xx + t_yx + t_zz + t_yz) * inv_h2;
            double val = c_scale * c1[s] * uc + d_ec * lap;
            if (compute_grad_p)
                val += (p_g1[o_c1 + cidx_g1(i, j, k, b.nx, b.ny)] - p_g1[o_c1 + cidx_g1(i, j - 1, k, b.nx, b.ny)]) * inv_h;
            out1[o1_g0 + idx] = val;
        }
        // Axis 2
        int nx2 = b.nx, ny2 = b.ny, nz2 = b.nz + 1;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx2 * ny2 * nz2; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx2, t = idx / nx2, j = t % ny2, k = t / ny2;
            if (is_l1 && k == b.nz && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            if (!is_l1 && zero_covered_l0)
            {
                bool cr = (k < c_nz0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (k > 0) && (l1_map[i + c_nx0 * (j + c_ny0 * (k - 1))] >= 0);
                if (cr && cl) { out2[idx] = 0.0; continue; }
            }
            size_t s = o2_g1 + sidx_g1(2, i, j, k, b.nx, b.ny);
            int gk = b.lo[2] + k;
            if (dirichlet_copy_bc && (gk == 0 || gk == gz_max)) { out2[o2_g0 + idx] = u2[s]; continue; }
            double uc = u2[s];
            double t_zz = 2.0 * (u2[o2_g1 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - 2.0 * uc + u2[o2_g1 + sidx_g1(2, i, j, k - 1, b.nx, b.ny)]);
            double t_yy = u2[o2_g1 + sidx_g1(2, i, j + 1, k, b.nx, b.ny)] - 2.0 * uc + u2[o2_g1 + sidx_g1(2, i, j - 1, k, b.nx, b.ny)];
            double t_zy = (u1[o1_g1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[o1_g1 + sidx_g1(1, i, j + 1, k - 1, b.nx, b.ny)]) -
                          (u1[o1_g1 + sidx_g1(1, i, j,     k, b.nx, b.ny)] - u1[o1_g1 + sidx_g1(1, i, j,     k - 1, b.nx, b.ny)]);
            double t_xx = u2[o2_g1 + sidx_g1(2, i + 1, j, k, b.nx, b.ny)] - 2.0 * uc + u2[o2_g1 + sidx_g1(2, i - 1, j, k, b.nx, b.ny)];
            double t_zx = (u0[o0_g1 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[o0_g1 + sidx_g1(0, i + 1, j, k - 1, b.nx, b.ny)]) -
                          (u0[o0_g1 + sidx_g1(0, i,     j, k, b.nx, b.ny)] - u0[o0_g1 + sidx_g1(0, i,     j, k - 1, b.nx, b.ny)]);
            double lap = (t_zz + t_yy + t_zy + t_xx + t_zx) * inv_h2;
            double val = c_scale * c2[s] * uc + d_ec * lap;
            if (compute_grad_p)
                val += (p_g1[o_c1 + cidx_g1(i, j, k, b.nx, b.ny)] - p_g1[o_c1 + cidx_g1(i, j, k - 1, b.nx, b.ny)]) * inv_h;
            out2[o2_g0 + idx] = val;
        }
    }
}

__global__ void compute_div_kernel(
    const double* u0, const double* u1, const double* u2,
    double* div_out, double scale, int nx, int ny, int nz, double h,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    double fac = scale / h;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o_c0 = is_l1 ? b.off_c_g0 : 0;
        size_t o0 = is_l1 ? b.off_s_g1[0] : 0;
        size_t o1 = is_l1 ? b.off_s_g1[1] : 0;
        size_t o2 = is_l1 ? b.off_s_g1[2] : 0;
        int tot = b.nx * b.ny * b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            if (!is_l1 && l1_map != nullptr && l1_map[idx] >= 0) { div_out[idx] = 0.0; continue; }
            int i = idx % b.nx, t = idx / b.nx, j = t % b.ny, k = t / b.ny;
            double dv = (u0[o0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[o0 + sidx_g1(0, i, j, k, b.nx, b.ny)]) +
                        (u1[o1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[o1 + sidx_g1(1, i, j, k, b.nx, b.ny)]) +
                        (u2[o2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[o2 + sidx_g1(2, i, j, k, b.nx, b.ny)]);
            div_out[o_c0 + idx] = fac * dv;
        }
    }
}

// ============================================================================
// Dual-Precision (FP32 / FP64) 7-Level Symmetric-Tensor Matrix-Free Multigrid
// ============================================================================
template <typename RealT>
__global__ void p_poisson_residual_kernel(
    const RealT* phi, const RealT* b, const RealT* d0, const RealT* d1, const RealT* d2,
    RealT* res, int nx, int ny, int nz, RealT dx0, RealT dx1, RealT dx2,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    RealT idx0 = RealT(1.0) / (dx0 * dx0);
    RealT idx1 = RealT(1.0) / (dx1 * dx1);
    RealT idx2 = RealT(1.0) / (dx2 * dx2);

    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b_box{};
        if (is_l1) b_box = boxes[pat];
        else { b_box.lo[0] = b_box.lo[1] = b_box.lo[2] = 0; b_box.nx = nx; b_box.ny = ny; b_box.nz = nz; }
        size_t off_c = is_l1 ? b_box.off_c_g1 : 0;
        size_t off0 = is_l1 ? b_box.off_s_g1[0] : 0;
        size_t off1 = is_l1 ? b_box.off_s_g1[1] : 0;
        size_t off2 = is_l1 ? b_box.off_s_g1[2] : 0;
        int gy_max = is_l1 ? 2 * c_ny0 : ny;
        int gz_max = is_l1 ? 2 * c_nz0 : nz;

        int tot = b_box.nx * b_box.ny * b_box.nz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int i = lin % b_box.nx, t = lin / b_box.nx, j = t % b_box.ny, k = t / b_box.ny;
            size_t c = off_c + cidx_g1(i, j, k, b_box.nx, b_box.ny);
            if (!is_l1 && l1_map != nullptr && l1_map[lin] >= 0) { res[c] = RealT(0.0); continue; }

            RealT pc = phi[c];
            RealT d_xp = d0[off0 + sidx_g1(0, i + 1, j, k, b_box.nx, b_box.ny)];
            RealT d_xm = d0[off0 + sidx_g1(0, i, j, k, b_box.nx, b_box.ny)];
            RealT d_yp = (b_box.lo[1] + j + 1 == gy_max) ? RealT(0.0) : d1[off1 + sidx_g1(1, i, j + 1, k, b_box.nx, b_box.ny)];
            RealT d_ym = (b_box.lo[1] + j == 0)          ? RealT(0.0) : d1[off1 + sidx_g1(1, i, j, k, b_box.nx, b_box.ny)];
            RealT d_zp = (b_box.lo[2] + k + 1 == gz_max) ? RealT(0.0) : d2[off2 + sidx_g1(2, i, j, k + 1, b_box.nx, b_box.ny)];
            RealT d_zm = (b_box.lo[2] + k == 0)          ? RealT(0.0) : d2[off2 + sidx_g1(2, i, j, k, b_box.nx, b_box.ny)];

            RealT div_d_grad =
                idx0 * (d_xp * (phi[off_c + cidx_g1(i + 1, j, k, b_box.nx, b_box.ny)] - pc) -
                        d_xm * (pc - phi[off_c + cidx_g1(i - 1, j, k, b_box.nx, b_box.ny)])) +
                idx1 * (d_yp * (phi[off_c + cidx_g1(i, j + 1, k, b_box.nx, b_box.ny)] - pc) -
                        d_ym * (pc - phi[off_c + cidx_g1(i, j - 1, k, b_box.nx, b_box.ny)])) +
                idx2 * (d_zp * (phi[off_c + cidx_g1(i, j, k + 1, b_box.nx, b_box.ny)] - pc) -
                        d_zm * (pc - phi[off_c + cidx_g1(i, j, k - 1, b_box.nx, b_box.ny)]));

            res[c] = b[c] - div_d_grad;
        }
    }
}

template <typename RealT>
__global__ void p_poisson_rbgs_kernel(
    RealT* phi, const RealT* b, const RealT* d0, const RealT* d1, const RealT* d2,
    int nx, int ny, int nz, RealT dx0, RealT dx1, RealT dx2, int color,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    RealT idx0 = RealT(1.0) / (dx0 * dx0);
    RealT idx1 = RealT(1.0) / (dx1 * dx1);
    RealT idx2 = RealT(1.0) / (dx2 * dx2);

    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b_box{};
        if (is_l1) b_box = boxes[pat];
        else { b_box.lo[0] = b_box.lo[1] = b_box.lo[2] = 0; b_box.nx = nx; b_box.ny = ny; b_box.nz = nz; }
        size_t off_c = is_l1 ? b_box.off_c_g1 : 0;
        size_t off0 = is_l1 ? b_box.off_s_g1[0] : 0;
        size_t off1 = is_l1 ? b_box.off_s_g1[1] : 0;
        size_t off2 = is_l1 ? b_box.off_s_g1[2] : 0;
        int gx_max = is_l1 ? 2 * c_nx0 : nx;
        int gy_max = is_l1 ? 2 * c_ny0 : ny;
        int gz_max = is_l1 ? 2 * c_nz0 : nz;

        int tot = b_box.nx * b_box.ny * b_box.nz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int i = lin % b_box.nx, t = lin / b_box.nx, j = t % b_box.ny, k = t / b_box.ny;
            if (((b_box.lo[0] + i + b_box.lo[1] + j + b_box.lo[2] + k) & 1) != color) continue;
            if (!is_l1 && l1_map != nullptr && l1_map[lin] >= 0) continue;

            size_t c = off_c + cidx_g1(i, j, k, b_box.nx, b_box.ny);
            int gi = b_box.lo[0] + i;
            int gj = b_box.lo[1] + j;
            int gk = b_box.lo[2] + k;

            RealT md_xp = -d0[off0 + sidx_g1(0, i + 1, j, k, b_box.nx, b_box.ny)] * idx0;
            RealT md_xm = -d0[off0 + sidx_g1(0, i, j, k, b_box.nx, b_box.ny)] * idx0;
            RealT md_yp = (gj + 1 == gy_max) ? RealT(0.0) : -d1[off1 + sidx_g1(1, i, j + 1, k, b_box.nx, b_box.ny)] * idx1;
            RealT md_ym = (gj == 0)          ? RealT(0.0) : -d1[off1 + sidx_g1(1, i, j, k, b_box.nx, b_box.ny)] * idx1;
            RealT md_zp = (gk + 1 == gz_max) ? RealT(0.0) : -d2[off2 + sidx_g1(2, i, j, k + 1, b_box.nx, b_box.ny)] * idx2;
            RealT md_zm = (gk == 0)          ? RealT(0.0) : -d2[off2 + sidx_g1(2, i, j, k, b_box.nx, b_box.ny)] * idx2;

            RealT diag = RealT(0.0), off_rhs = b[c];
            if (gi == 0) { diag += RealT(2.0) * md_xm; }
            else { diag += md_xm; off_rhs += md_xm * phi[off_c + cidx_g1(i - 1, j, k, b_box.nx, b_box.ny)]; }
            if (gi + 1 == gx_max) { diag += RealT(2.0) * md_xp; }
            else { diag += md_xp; off_rhs += md_xp * phi[off_c + cidx_g1(i + 1, j, k, b_box.nx, b_box.ny)]; }

            if (gj > 0) { diag += md_ym; off_rhs += md_ym * phi[off_c + cidx_g1(i, j - 1, k, b_box.nx, b_box.ny)]; }
            if (gj + 1 < gy_max) { diag += md_yp; off_rhs += md_yp * phi[off_c + cidx_g1(i, j + 1, k, b_box.nx, b_box.ny)]; }
            if (gk > 0) { diag += md_zm; off_rhs += md_zm * phi[off_c + cidx_g1(i, j, k - 1, b_box.nx, b_box.ny)]; }
            if (gk + 1 < gz_max) { diag += md_zp; off_rhs += md_zp * phi[off_c + cidx_g1(i, j, k + 1, b_box.nx, b_box.ny)]; }

            if (diag > RealT(0.0)) phi[c] = off_rhs / diag;
        }
    }
}

template <typename RealT>
__global__ void restrict_cell_kernel(
    const RealT* rf, RealT* bc, RealT* xc, int nx_c, int ny_c, int nz_c, int nx_f, int ny_f, int rx, int ry, int rz)
{
    int tot = nx_c * ny_c * nz_c;
    for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
    {
        int ic = lin % nx_c, t = lin / nx_c, jc = t % ny_c, kc = t / ny_c;
        RealT s = RealT(0.0);
        for (int dk = 0; dk < rz; ++dk)
            for (int dj = 0; dj < ry; ++dj)
                for (int di = 0; di < rx; ++di)
                    s += rf[cidx_g1(ic * rx + di, jc * ry + dj, kc * rz + dk, nx_f, ny_f)];
        int c = cidx_g1(ic, jc, kc, nx_c, ny_c);
        bc[c] = s / static_cast<RealT>(rx * ry * rz);
        xc[c] = RealT(0.0);
    }
}

template <typename RealT>
__global__ void prolong_cell_add_kernel(
    const RealT* xc, RealT* xf, int nx_f, int ny_f, int nz_f, int nx_c, int ny_c, int rx, int ry, int rz,
    bool is_l1, const Box3i* boxes, int n_l1)
{
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx_f; b.ny = ny_f; b.nz = nz_f; }
        size_t off = is_l1 ? b.off_c_g1 : 0;
        int tot = b.nx * b.ny * b.nz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int li = lin % b.nx, t = lin / b.nx, lj = t % b.ny, lk = t / b.ny;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            int ic = gi / rx, jc = gj / ry, kc = gk / rz;
            xf[off + cidx_g1(li, lj, lk, b.nx, b.ny)] += xc[cidx_g1(ic, jc, kc, nx_c, ny_c)];
        }
    }
}

template <typename RealT>
__global__ void p_poisson_c6_pcg_kernel(
    RealT* phi, const RealT* b, const RealT* d0, const RealT* d1, const RealT* d2,
    int nx, int ny, int nz, RealT dx0, RealT dx1, RealT dx2)
{
    constexpr int MAX_C6 = 256;
    __shared__ double s_x[MAX_C6];
    __shared__ double s_r[MAX_C6];
    __shared__ double s_z[MAX_C6];
    __shared__ double s_p[MAX_C6];
    __shared__ double s_ap[MAX_C6];
    __shared__ double s_invd[MAX_C6];
    __shared__ double s_red[MAX_C6];

    int tid = threadIdx.x;
    int n = nx * ny * nz;
    double idx0 = 1.0 / (static_cast<double>(dx0) * static_cast<double>(dx0));
    double idx1 = 1.0 / (static_cast<double>(dx1) * static_cast<double>(dx1));
    double idx2 = 1.0 / (static_cast<double>(dx2) * static_cast<double>(dx2));

    auto shared_dot = [&](double local_v) -> double {
        s_red[tid] = local_v;
        __syncthreads();
        for (int stride = 128; stride > 0; stride >>= 1)
        {
            if (tid < stride) s_red[tid] += s_red[tid + stride];
            __syncthreads();
        }
        return s_red[0];
    };

    int i = 0, j = 0, k = 0;
    double md_xm = 0.0, md_xp = 0.0, md_ym = 0.0, md_yp = 0.0, md_zm = 0.0, md_zp = 0.0, diag = 0.0;
    if (tid < n)
    {
        i = tid % nx;
        int t = tid / nx;
        j = t % ny;
        k = t / ny;
        md_xp = -static_cast<double>(d0[sidx_g1(0, i + 1, j, k, nx, ny)]) * idx0;
        md_xm = -static_cast<double>(d0[sidx_g1(0, i,     j, k, nx, ny)]) * idx0;
        md_yp = (j + 1 == ny) ? 0.0 : -static_cast<double>(d1[sidx_g1(1, i, j + 1, k, nx, ny)]) * idx1;
        md_ym = (j == 0)      ? 0.0 : -static_cast<double>(d1[sidx_g1(1, i, j,     k, nx, ny)]) * idx1;
        md_zp = (k + 1 == nz) ? 0.0 : -static_cast<double>(d2[sidx_g1(2, i, j, k + 1, nx, ny)]) * idx2;
        md_zm = (k == 0)      ? 0.0 : -static_cast<double>(d2[sidx_g1(2, i, j, k,     nx, ny)]) * idx2;
        diag = (i == 0 ? 2.0 * md_xm : md_xm) + (i + 1 == nx ? 2.0 * md_xp : md_xp) +
               md_ym + md_yp + md_zm + md_zp;
        s_invd[tid] = (diag > 0.0) ? (1.0 / diag) : 0.0;
        s_x[tid] = static_cast<double>(phi[cidx_g1(i, j, k, nx, ny)]);
    }
    __syncthreads();

    auto eval_A = [&](const double* v) -> double {
        double vc = v[tid];
        double out = diag * vc;
        if (i > 0)      out -= md_xm * v[tid - 1];
        if (i + 1 < nx) out -= md_xp * v[tid + 1];
        if (j > 0)      out -= md_ym * v[tid - nx];
        if (j + 1 < ny) out -= md_yp * v[tid + nx];
        if (k > 0)      out -= md_zm * v[tid - nx * ny];
        if (k + 1 < nz) out -= md_zp * v[tid + nx * ny];
        return out;
    };

    double bval = 0.0;
    if (tid < n)
    {
        bval = static_cast<double>(b[cidx_g1(i, j, k, nx, ny)]);
        double ax0 = eval_A(s_x);
        double r0 = bval - ax0;
        double z0 = s_invd[tid] * r0;
        s_r[tid] = r0;
        s_z[tid] = z0;
        s_p[tid] = z0;
    }
    __syncthreads();

    double b2 = shared_dot(tid < n ? bval * bval : 0.0);
    double rz = shared_dot(tid < n ? s_r[tid] * s_z[tid] : 0.0);
    double tol2 = 1e-24 * b2;

    for (int it = 0; it < 200; ++it)
    {
        if (rz <= 0.0 || !isfinite(rz)) break;
        if (tid < n) s_ap[tid] = eval_A(s_p);
        __syncthreads();

        double pAp = shared_dot(tid < n ? s_p[tid] * s_ap[tid] : 0.0);
        if (pAp <= 0.0 || !isfinite(pAp)) break;
        double alpha = rz / pAp;

        if (tid < n)
        {
            s_x[tid] += alpha * s_p[tid];
            s_r[tid] -= alpha * s_ap[tid];
        }
        __syncthreads();

        double r2 = shared_dot(tid < n ? s_r[tid] * s_r[tid] : 0.0);
        if (r2 <= tol2) break;

        if (tid < n) s_z[tid] = s_invd[tid] * s_r[tid];
        __syncthreads();

        double rz_new = shared_dot(tid < n ? s_r[tid] * s_z[tid] : 0.0);
        double beta = rz_new / rz;
        if (tid < n) s_p[tid] = s_z[tid] + beta * s_p[tid];
        rz = rz_new;
        __syncthreads();
    }

    if (tid < n) phi[cidx_g1(i, j, k, nx, ny)] = static_cast<RealT>(s_x[tid]);
}

template <typename RealT>
__global__ void vel_helmholtz_residual_kernel(
    const RealT* u0, const RealT* u1, const RealT* u2,
    const RealT* b0, const RealT* b1, const RealT* b2,
    const RealT* c0, const RealT* c1, const RealT* c2,
    RealT* r0, RealT* r1, RealT* r2,
    int nx, int ny, int nz, RealT dx0, RealT dx1, RealT dx2,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    RealT f0 = RealT(1.0) / dx0, f1 = RealT(1.0) / dx1, f2 = RealT(1.0) / dx2;
    RealT dvel = static_cast<RealT>(c_d_vel_const);
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;

    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t off0 = is_l1 ? b.off_s_g1[0] : 0;
        size_t off1 = is_l1 ? b.off_s_g1[1] : 0;
        size_t off2 = is_l1 ? b.off_s_g1[2] : 0;

        int nx0 = b.nx + 1, ny0 = b.ny, nz0 = b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx0 * ny0 * nz0; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx0, t = idx / nx0, j = t % ny0, k = t / ny0;
            size_t s = off0 + sidx_g1(0, i, j, k, b.nx, b.ny);
            if (!is_l1 && l1_map != nullptr)
            {
                bool cr = (i < c_nx0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (i > 0) && (l1_map[(i - 1) + c_nx0 * (j + c_ny0 * k)] >= 0);
                if (cr || cl) { r0[s] = RealT(0.0); continue; }
            }
            RealT uc = u0[s];
            RealT au = c0[s] * uc + dvel * (
                RealT(2.0) * f0 * f0 * (u0[off0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - RealT(2.0) * uc + u0[off0 + sidx_g1(0, i - 1, j, k, b.nx, b.ny)]) +
                f1 * f1 * (u0[off0 + sidx_g1(0, i, j + 1, k, b.nx, b.ny)] - RealT(2.0) * uc + u0[off0 + sidx_g1(0, i, j - 1, k, b.nx, b.ny)]) +
                f2 * f2 * (u0[off0 + sidx_g1(0, i, j, k + 1, b.nx, b.ny)] - RealT(2.0) * uc + u0[off0 + sidx_g1(0, i, j, k - 1, b.nx, b.ny)]) +
                f0 * f1 * ((u1[off1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[off1 + sidx_g1(1, i - 1, j + 1, k, b.nx, b.ny)]) -
                           (u1[off1 + sidx_g1(1, i, j, k, b.nx, b.ny)]     - u1[off1 + sidx_g1(1, i - 1, j, k, b.nx, b.ny)])) +
                f0 * f2 * ((u2[off2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[off2 + sidx_g1(2, i - 1, j, k + 1, b.nx, b.ny)]) -
                           (u2[off2 + sidx_g1(2, i, j, k, b.nx, b.ny)]     - u2[off2 + sidx_g1(2, i - 1, j, k, b.nx, b.ny)])));
            r0[s] = b0[s] - au;
        }
        int nx1 = b.nx, ny1 = b.ny + 1, nz1 = b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx1 * ny1 * nz1; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx1, t = idx / nx1, j = t % ny1, k = t / ny1;
            size_t s = off1 + sidx_g1(1, i, j, k, b.nx, b.ny);
            if (!is_l1 && l1_map != nullptr)
            {
                bool cr = (j < c_ny0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (j > 0) && (l1_map[i + c_nx0 * ((j - 1) + c_ny0 * k)] >= 0);
                if (cr || cl) { r1[s] = RealT(0.0); continue; }
            }
            int gj = b.lo[1] + j;
            if (gj == 0 || gj == gy_max) { r1[s] = b1[s] - u1[s]; continue; }
            RealT vc = u1[s];
            RealT av = c1[s] * vc + dvel * (
                f0 * f0 * (u1[off1 + sidx_g1(1, i + 1, j, k, b.nx, b.ny)] - RealT(2.0) * vc + u1[off1 + sidx_g1(1, i - 1, j, k, b.nx, b.ny)]) +
                RealT(2.0) * f1 * f1 * (u1[off1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - RealT(2.0) * vc + u1[off1 + sidx_g1(1, i, j - 1, k, b.nx, b.ny)]) +
                f2 * f2 * (u1[off1 + sidx_g1(1, i, j, k + 1, b.nx, b.ny)] - RealT(2.0) * vc + u1[off1 + sidx_g1(1, i, j, k - 1, b.nx, b.ny)]) +
                f1 * f0 * ((u0[off0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[off0 + sidx_g1(0, i + 1, j - 1, k, b.nx, b.ny)]) -
                           (u0[off0 + sidx_g1(0, i, j, k, b.nx, b.ny)]     - u0[off0 + sidx_g1(0, i, j - 1, k, b.nx, b.ny)])) +
                f1 * f2 * ((u2[off2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[off2 + sidx_g1(2, i, j - 1, k + 1, b.nx, b.ny)]) -
                           (u2[off2 + sidx_g1(2, i, j, k, b.nx, b.ny)]     - u2[off2 + sidx_g1(2, i, j - 1, k, b.nx, b.ny)])));
            r1[s] = b1[s] - av;
        }
        int nx2 = b.nx, ny2 = b.ny, nz2 = b.nz + 1;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx2 * ny2 * nz2; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx2, t = idx / nx2, j = t % ny2, k = t / ny2;
            size_t s = off2 + sidx_g1(2, i, j, k, b.nx, b.ny);
            if (!is_l1 && l1_map != nullptr)
            {
                bool cr = (k < c_nz0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (k > 0) && (l1_map[i + c_nx0 * (j + c_ny0 * (k - 1))] >= 0);
                if (cr || cl) { r2[s] = RealT(0.0); continue; }
            }
            int gk = b.lo[2] + k;
            if (gk == 0 || gk == gz_max) { r2[s] = b2[s] - u2[s]; continue; }
            RealT wc = u2[s];
            RealT aw = c2[s] * wc + dvel * (
                f0 * f0 * (u2[off2 + sidx_g1(2, i + 1, j, k, b.nx, b.ny)] - RealT(2.0) * wc + u2[off2 + sidx_g1(2, i - 1, j, k, b.nx, b.ny)]) +
                f1 * f1 * (u2[off2 + sidx_g1(2, i, j + 1, k, b.nx, b.ny)] - RealT(2.0) * wc + u2[off2 + sidx_g1(2, i, j - 1, k, b.nx, b.ny)]) +
                RealT(2.0) * f2 * f2 * (u2[off2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - RealT(2.0) * wc + u2[off2 + sidx_g1(2, i, j, k - 1, b.nx, b.ny)]) +
                f2 * f0 * ((u0[off0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[off0 + sidx_g1(0, i + 1, j, k - 1, b.nx, b.ny)]) -
                           (u0[off0 + sidx_g1(0, i, j, k, b.nx, b.ny)]     - u0[off0 + sidx_g1(0, i, j, k - 1, b.nx, b.ny)])) +
                f2 * f1 * ((u1[off1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[off1 + sidx_g1(1, i, j + 1, k - 1, b.nx, b.ny)]) -
                           (u1[off1 + sidx_g1(1, i, j, k, b.nx, b.ny)]     - u1[off1 + sidx_g1(1, i, j, k - 1, b.nx, b.ny)])));
            r2[s] = b2[s] - aw;
        }
    }
}

template <typename RealT>
__global__ void vel_helmholtz_rbgs_axis_kernel(
    RealT* u0, RealT* u1, RealT* u2,
    const RealT* b0, const RealT* b1, const RealT* b2,
    const RealT* c0, const RealT* c1, const RealT* c2,
    int axis, int nx, int ny, int nz, RealT dx0, RealT dx1, RealT dx2, int color,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    RealT f0 = RealT(1.0) / dx0, f1 = RealT(1.0) / dx1, f2 = RealT(1.0) / dx2;
    RealT md = static_cast<RealT>(-c_d_vel_const);
    int gx_max = is_l1 ? 2 * c_nx0 : nx;
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;

    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t off0 = is_l1 ? b.off_s_g1[0] : 0;
        size_t off1 = is_l1 ? b.off_s_g1[1] : 0;
        size_t off2 = is_l1 ? b.off_s_g1[2] : 0;

        int nxa = b.nx + (axis == 0 ? 1 : 0);
        int nya = b.ny + (axis == 1 ? 1 : 0);
        int nza = b.nz + (axis == 2 ? 1 : 0);
        int tot = nxa * nya * nza;

        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
            int gi = b.lo[0] + i, gj = b.lo[1] + j, gk = b.lo[2] + k;
            if (((gi + gj + gk) & 1) != color) continue;

            if (is_l1)
            {
                if (axis == 0 && i == b.nx && lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
                if (axis == 1 && j == b.ny && lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
                if (axis == 2 && k == b.nz && lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
            }
            else if (l1_map != nullptr)
            {
                if (axis == 0)
                {
                    if ((i < c_nx0 && l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0) || (i > 0 && l1_map[(i - 1) + c_nx0 * (j + c_ny0 * k)] >= 0)) continue;
                }
                else if (axis == 1)
                {
                    if ((j < c_ny0 && l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0) || (j > 0 && l1_map[i + c_nx0 * ((j - 1) + c_ny0 * k)] >= 0)) continue;
                }
                else
                {
                    if ((k < c_nz0 && l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0) || (k > 0 && l1_map[i + c_nx0 * (j + c_ny0 * (k - 1))] >= 0)) continue;
                }
            }

            if (axis == 0)
            {
                size_t s = off0 + sidx_g1(0, i, j, k, b.nx, b.ny);
                RealT uc = u0[s];
                RealT au = c0[s] * uc - md * (
                    RealT(2.0) * f0 * f0 * (u0[off0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - RealT(2.0) * uc + u0[off0 + sidx_g1(0, i - 1, j, k, b.nx, b.ny)]) +
                    f1 * f1 * (u0[off0 + sidx_g1(0, i, j + 1, k, b.nx, b.ny)] - RealT(2.0) * uc + u0[off0 + sidx_g1(0, i, j - 1, k, b.nx, b.ny)]) +
                    f2 * f2 * (u0[off0 + sidx_g1(0, i, j, k + 1, b.nx, b.ny)] - RealT(2.0) * uc + u0[off0 + sidx_g1(0, i, j, k - 1, b.nx, b.ny)]) +
                    f0 * f1 * ((u1[off1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[off1 + sidx_g1(1, i - 1, j + 1, k, b.nx, b.ny)]) -
                               (u1[off1 + sidx_g1(1, i, j, k, b.nx, b.ny)]     - u1[off1 + sidx_g1(1, i - 1, j, k, b.nx, b.ny)])) +
                    f0 * f2 * ((u2[off2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[off2 + sidx_g1(2, i - 1, j, k + 1, b.nx, b.ny)]) -
                               (u2[off2 + sidx_g1(2, i, j, k, b.nx, b.ny)]     - u2[off2 + sidx_g1(2, i - 1, j, k, b.nx, b.ny)])));
                RealT diag = c0[s] + md * (RealT(4.0) * f0 * f0 + RealT(2.0) * f1 * f1 + RealT(2.0) * f2 * f2);
                if (gj == 0) diag += md * f1 * f1;
                if (gj == gy_max - 1) diag += md * f1 * f1;
                if (gk == 0) diag -= md * f2 * f2;
                if (gk == gz_max - 1) diag += md * f2 * f2;
                u0[s] = uc + (b0[s] - au) / diag;
            }
            else if (axis == 1)
            {
                size_t s = off1 + sidx_g1(1, i, j, k, b.nx, b.ny);
                if (gj == 0 || gj == gy_max) { u1[s] = b1[s]; continue; }
                RealT vc = u1[s];
                RealT av = c1[s] * vc - md * (
                    f0 * f0 * (u1[off1 + sidx_g1(1, i + 1, j, k, b.nx, b.ny)] - RealT(2.0) * vc + u1[off1 + sidx_g1(1, i - 1, j, k, b.nx, b.ny)]) +
                    RealT(2.0) * f1 * f1 * (u1[off1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - RealT(2.0) * vc + u1[off1 + sidx_g1(1, i, j - 1, k, b.nx, b.ny)]) +
                    f2 * f2 * (u1[off1 + sidx_g1(1, i, j, k + 1, b.nx, b.ny)] - RealT(2.0) * vc + u1[off1 + sidx_g1(1, i, j, k - 1, b.nx, b.ny)]) +
                    f1 * f0 * ((u0[off0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[off0 + sidx_g1(0, i + 1, j - 1, k, b.nx, b.ny)]) -
                               (u0[off0 + sidx_g1(0, i, j, k, b.nx, b.ny)]     - u0[off0 + sidx_g1(0, i, j - 1, k, b.nx, b.ny)])) +
                    f1 * f2 * ((u2[off2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - u2[off2 + sidx_g1(2, i, j - 1, k + 1, b.nx, b.ny)]) -
                               (u2[off2 + sidx_g1(2, i, j, k, b.nx, b.ny)]     - u2[off2 + sidx_g1(2, i, j - 1, k, b.nx, b.ny)])));
                RealT diag = c1[s] + md * (RealT(2.0) * f0 * f0 + RealT(4.0) * f1 * f1 + RealT(2.0) * f2 * f2);
                if (gi == 0) diag += md * f0 * f0;
                if (gi == gx_max - 1) diag += md * f0 * f0;
                if (gk == 0) diag -= md * f2 * f2;
                if (gk == gz_max - 1) diag += md * f2 * f2;
                u1[s] = vc + (b1[s] - av) / diag;
            }
            else
            {
                size_t s = off2 + sidx_g1(2, i, j, k, b.nx, b.ny);
                if (gk == 0 || gk == gz_max) { u2[s] = b2[s]; continue; }
                RealT wc = u2[s];
                RealT aw = c2[s] * wc - md * (
                    f0 * f0 * (u2[off2 + sidx_g1(2, i + 1, j, k, b.nx, b.ny)] - RealT(2.0) * wc + u2[off2 + sidx_g1(2, i - 1, j, k, b.nx, b.ny)]) +
                    f1 * f1 * (u2[off2 + sidx_g1(2, i, j + 1, k, b.nx, b.ny)] - RealT(2.0) * wc + u2[off2 + sidx_g1(2, i, j - 1, k, b.nx, b.ny)]) +
                    RealT(2.0) * f2 * f2 * (u2[off2 + sidx_g1(2, i, j, k + 1, b.nx, b.ny)] - RealT(2.0) * wc + u2[off2 + sidx_g1(2, i, j, k - 1, b.nx, b.ny)]) +
                    f2 * f0 * ((u0[off0 + sidx_g1(0, i + 1, j, k, b.nx, b.ny)] - u0[off0 + sidx_g1(0, i + 1, j, k - 1, b.nx, b.ny)]) -
                               (u0[off0 + sidx_g1(0, i, j, k, b.nx, b.ny)]     - u0[off0 + sidx_g1(0, i, j, k - 1, b.nx, b.ny)])) +
                    f2 * f1 * ((u1[off1 + sidx_g1(1, i, j + 1, k, b.nx, b.ny)] - u1[off1 + sidx_g1(1, i, j + 1, k - 1, b.nx, b.ny)]) -
                               (u1[off1 + sidx_g1(1, i, j, k, b.nx, b.ny)]     - u1[off1 + sidx_g1(1, i, j, k - 1, b.nx, b.ny)])));
                RealT diag = c2[s] + md * (RealT(2.0) * f0 * f0 + RealT(2.0) * f1 * f1 + RealT(4.0) * f2 * f2);
                if (gi == 0) diag += md * f0 * f0;
                if (gi == gx_max - 1) diag += md * f0 * f0;
                if (gj == 0) diag += md * f1 * f1;
                if (gj == gy_max - 1) diag += md * f1 * f1;
                u2[s] = wc + (b2[s] - aw) / diag;
            }
        }
    }
}

template <typename RealT>
__global__ void restrict_side_kernel(
    const RealT* rf0, const RealT* rf1, const RealT* rf2,
    RealT* bc0, RealT* bc1, RealT* bc2,
    RealT* xc0, RealT* xc1, RealT* xc2,
    int nx_c, int ny_c, int nz_c, int nx_f, int ny_f, int rx, int ry, int rz)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        int cnx = nx_c + (axis == 0 ? 1 : 0);
        int cny = ny_c + (axis == 1 ? 1 : 0);
        int cnz = nz_c + (axis == 2 ? 1 : 0);
        int tot = cnx * cny * cnz;
        const RealT* rf = (axis == 0) ? rf0 : (axis == 1 ? rf1 : rf2);
        RealT* bc = (axis == 0) ? bc0 : (axis == 1 ? bc1 : bc2);
        RealT* xc = (axis == 0) ? xc0 : (axis == 1 ? xc1 : xc2);
        int r_t1 = (axis == 0) ? ry : rx;
        int r_t2 = (axis == 2) ? ry : rz;

        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int ic = lin % cnx, t = lin / cnx, jc = t % cny, kc = t / cny;
            RealT s = RealT(0.0);
            for (int dt2 = 0; dt2 < r_t2; ++dt2)
                for (int dt1 = 0; dt1 < r_t1; ++dt1)
                {
                    int fi = (axis == 0) ? (ic * rx) : (ic * rx + dt1);
                    int fj = (axis == 1) ? (jc * ry) : (jc * ry + (axis == 0 ? dt1 : dt2));
                    int fk = (axis == 2) ? (kc * rz) : (kc * rz + dt2);
                    s += rf[sidx_g1(axis, fi, fj, fk, nx_f, ny_f)];
                }
            int sc = sidx_g1(axis, ic, jc, kc, nx_c, ny_c);
            bc[sc] = s / static_cast<RealT>(r_t1 * r_t2);
            xc[sc] = RealT(0.0);
        }
    }
}

template <typename RealT>
__global__ void prolong_side_add_kernel(
    const RealT* xc0, const RealT* xc1, const RealT* xc2,
    RealT* xf0, RealT* xf1, RealT* xf2,
    int nx_f, int ny_f, int nz_f, int nx_c, int ny_c, int nz_c, int rx, int ry, int rz,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx_f; b.ny = ny_f; b.nz = nz_f; }
        for (int axis = 0; axis < 3; ++axis)
        {
            const RealT* xc = (axis == 0) ? xc0 : (axis == 1 ? xc1 : xc2);
            RealT* xf = (axis == 0) ? xf0 : (axis == 1 ? xf1 : xf2);
            size_t off = is_l1 ? b.off_s_g1[axis] : 0;
            int nxa = b.nx + (axis == 0 ? 1 : 0);
            int nya = b.ny + (axis == 1 ? 1 : 0);
            int nza = b.nz + (axis == 2 ? 1 : 0);
            int tot = nxa * nya * nza;
            for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
            {
                int li = lin % nxa, t = lin / nxa, lj = t % nya, lk = t / nya;
                int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
                if (is_l1 && l1_map != nullptr)
                {
                    if (axis == 0 && li == b.nx && lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
                    if (axis == 1 && lj == b.ny && lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
                    if (axis == 2 && lk == b.nz && lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
                }
                RealT add_val = RealT(0.0);
                if (axis == 0)
                {
                    int ic = gi / rx, jc = gj / ry, kc = gk / rz;
                    int rem = gi - ic * rx;
                    if (rem == 0 || ic >= nx_c) add_val = xc[sidx_g1(0, ic, jc, kc, nx_c, ny_c)];
                    else
                    {
                        RealT w = static_cast<RealT>(rem) / static_cast<RealT>(rx);
                        add_val = (RealT(1.0) - w) * xc[sidx_g1(0, ic, jc, kc, nx_c, ny_c)] +
                                  w * xc[sidx_g1(0, ic + 1, jc, kc, nx_c, ny_c)];
                    }
                }
                else if (axis == 1)
                {
                    int ic = gi / rx, jc = gj / ry, kc = gk / rz;
                    int rem = gj - jc * ry;
                    if (rem == 0 || jc >= ny_c) add_val = xc[sidx_g1(1, ic, jc, kc, nx_c, ny_c)];
                    else
                    {
                        RealT w = static_cast<RealT>(rem) / static_cast<RealT>(ry);
                        add_val = (RealT(1.0) - w) * xc[sidx_g1(1, ic, jc, kc, nx_c, ny_c)] +
                                  w * xc[sidx_g1(1, ic, jc + 1, kc, nx_c, ny_c)];
                    }
                }
                else
                {
                    int ic = gi / rx, jc = gj / ry, kc = gk / rz;
                    int rem = gk - kc * rz;
                    if (rem == 0 || kc >= nz_c) add_val = xc[sidx_g1(2, ic, jc, kc, nx_c, ny_c)];
                    else
                    {
                        RealT w = static_cast<RealT>(rem) / static_cast<RealT>(rz);
                        add_val = (RealT(1.0) - w) * xc[sidx_g1(2, ic, jc, kc, nx_c, ny_c)] +
                                  w * xc[sidx_g1(2, ic, jc, kc + 1, nx_c, ny_c)];
                    }
                }
                xf[off + sidx_g1(axis, li, lj, lk, b.nx, b.ny)] += add_val;
            }
        }
    }
}

template <typename RealT>
__global__ void coarsen_coef_level_kernel(
    const double* fc, const double* fd, double* cc, double* cd, RealT* ccf, RealT* cdf,
    int axis, int nxc, int nyc, int nzc, int nxf, int nyf, int rx, int ry, int rz, double d_val)
{
    int cnxa = nxc + (axis == 0), cnya = nyc + (axis == 1), cnza = nzc + (axis == 2);
    int tot = cnxa * cnya * cnza;
    int rt1 = (axis == 0) ? ry : rx, rt2 = (axis == 2) ? ry : rz;
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
    {
        int ic = idx % cnxa, t = idx / cnxa, jc = t % cnya, kc = t / cnya;
        double cs = 0.0, ds = 0.0;
        for (int d2 = 0; d2 < rt2; ++d2)
            for (int d1 = 0; d1 < rt1; ++d1)
            {
                int fi = (axis == 0) ? (ic * rx) : (ic * rx + d1);
                int fj = (axis == 1) ? (jc * ry) : (jc * ry + (axis == 0 ? d1 : d2));
                int fk = (axis == 2) ? (kc * rz) : (kc * rz + d2);
                int fsi = sidx_g1(axis, fi, fj, fk, nxf, nyf);
                cs += fc[fsi];
                ds += fd[fsi];
            }
        double cv = cs / (rt1 * rt2);
        int rn = (axis == 0) ? rx : (axis == 1 ? ry : rz);
        double dv = static_cast<double>(rn) * (ds / (rt1 * rt2));
        int si = sidx_g1(axis, ic, jc, kc, nxc, nyc);
        cc[si] = cv;
        cd[si] = dv;
        ccf[si] = static_cast<RealT>(cv);
        cdf[si] = static_cast<RealT>(dv);
    }
    (void)d_val;
}

template <typename RealT>
__global__ void cast_fine_coefs_fp32_kernel(
    const double* c0, const double* c1, const double* c2,
    const double* d0, const double* d1, const double* d2,
    RealT* c0f, RealT* c1f, RealT* c2f,
    RealT* d0f, RealT* d1f, RealT* d2f, size_t n)
{
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += static_cast<size_t>(gridDim.x) * blockDim.x)
    {
        c0f[i] = static_cast<RealT>(c0[i]); c1f[i] = static_cast<RealT>(c1[i]); c2f[i] = static_cast<RealT>(c2[i]);
        d0f[i] = static_cast<RealT>(d0[i]); d1f[i] = static_cast<RealT>(d1[i]); d2f[i] = static_cast<RealT>(d2[i]);
    }
}

template <typename RealT>
__global__ void scale_copy_f64_to_real_kernel(const double* src, RealT* dst, double scale, size_t n)
{
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += static_cast<size_t>(gridDim.x) * blockDim.x)
        dst[i] = static_cast<RealT>(scale * src[i]);
}

template <typename RealT>
__global__ void scale_copy_real_to_f64_kernel(const RealT* src, double* dst, double scale, size_t n)
{
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n; i += static_cast<size_t>(gridDim.x) * blockDim.x)
        dst[i] = scale * static_cast<double>(src[i]);
}

__global__ void projection_build_fphi_kernel(
    const double* fp_g0, const double* div_u, double* bphi_g1,
    int nx, int ny, int nz, bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o0 = is_l1 ? b.off_c_g0 : 0;
        size_t o1 = is_l1 ? b.off_c_g1 : 0;
        int tot = b.nx * b.ny * b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            if (!is_l1 && l1_map != nullptr && l1_map[idx] >= 0)
            { bphi_g1[o1 + cidx_g1(idx % nx, (idx / nx) % ny, idx / (nx * ny), nx, ny)] = 0.0; continue; }
            int i = idx % b.nx, t = idx / b.nx, j = t % b.ny, k = t / b.ny;
            bphi_g1[o1 + cidx_g1(i, j, k, b.nx, b.ny)] = -(div_u[o0 + idx] + fp_g0[o0 + idx]);
        }
    }
}

__global__ void projection_update_u_p_kernel(
    const double* us0, const double* us1, const double* us2,
    const double* phi_g1, const double* bphi_g1,
    const double* dp0, const double* dp1, const double* dp2,
    double* out_u0, double* out_u1, double* out_u2, double* out_p,
    double dt, double h, int nx, int ny, int nz,
    bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    double inv_dt = 1.0 / dt, inv_h = 1.0 / h;
    int gy_max = is_l1 ? 2 * c_ny0 : ny, gz_max = is_l1 ? 2 * c_nz0 : nz;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t oc0 = is_l1 ? b.off_c_g0 : 0, oc1 = is_l1 ? b.off_c_g1 : 0;
        size_t os0_0 = is_l1 ? b.off_s_g0[0] : 0, os1_0 = is_l1 ? b.off_s_g0[1] : 0, os2_0 = is_l1 ? b.off_s_g0[2] : 0;
        size_t os0_1 = is_l1 ? b.off_s_g1[0] : 0, os1_1 = is_l1 ? b.off_s_g1[1] : 0, os2_1 = is_l1 ? b.off_s_g1[2] : 0;
        int totc = b.nx * b.ny * b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < totc; idx += gridDim.x * blockDim.x)
        {
            if (!is_l1 && l1_map[idx] >= 0) { out_p[idx] = 0.0; continue; }
            int i = idx % b.nx, t = idx / b.nx, j = t % b.ny, k = t / b.ny;
            size_t cg = oc1 + cidx_g1(i, j, k, b.nx, b.ny);
            out_p[oc0 + idx] = inv_dt * phi_g1[cg] + bphi_g1[cg];
        }
        int nx0 = b.nx + 1, ny0 = b.ny, nz0 = b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx0 * ny0 * nz0; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx0, t = idx / nx0, j = t % ny0, k = t / ny0;
            if (is_l1 && i == b.nx && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            if (!is_l1)
            {
                bool cr = (i < c_nx0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (i > 0) && (l1_map[(i - 1) + c_nx0 * (j + c_ny0 * k)] >= 0);
                if (cr && cl) { out_u0[idx] = 0.0; continue; }
            }
            size_t s = os0_1 + sidx_g1(0, i, j, k, b.nx, b.ny);
            double gp = (phi_g1[oc1 + cidx_g1(i, j, k, b.nx, b.ny)] - phi_g1[oc1 + cidx_g1(i - 1, j, k, b.nx, b.ny)]) * inv_h;
            out_u0[os0_0 + idx] = us0[s] + dp0[s] * gp;
        }
        int nx1 = b.nx, ny1 = b.ny + 1, nz1 = b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx1 * ny1 * nz1; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx1, t = idx / nx1, j = t % ny1, k = t / ny1;
            if (is_l1 && j == b.ny && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            if (!is_l1)
            {
                bool cr = (j < c_ny0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (j > 0) && (l1_map[i + c_nx0 * ((j - 1) + c_ny0 * k)] >= 0);
                if (cr && cl) { out_u1[idx] = 0.0; continue; }
            }
            if (b.lo[1] + j == 0 || b.lo[1] + j == gy_max) { out_u1[os1_0 + idx] = 0.0; continue; }
            size_t s = os1_1 + sidx_g1(1, i, j, k, b.nx, b.ny);
            double gp = (phi_g1[oc1 + cidx_g1(i, j, k, b.nx, b.ny)] - phi_g1[oc1 + cidx_g1(i, j - 1, k, b.nx, b.ny)]) * inv_h;
            out_u1[os1_0 + idx] = us1[s] + dp1[s] * gp;
        }
        int nx2 = b.nx, ny2 = b.ny, nz2 = b.nz + 1;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nx2 * ny2 * nz2; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nx2, t = idx / nx2, j = t % ny2, k = t / ny2;
            if (is_l1 && k == b.nz && lookup_l1_patch(b.lo[0] + i, b.lo[1] + j, b.lo[2] + k, l1_map) >= 0) continue;
            if (!is_l1)
            {
                bool cr = (k < c_nz0) && (l1_map[i + c_nx0 * (j + c_ny0 * k)] >= 0);
                bool cl = (k > 0) && (l1_map[i + c_nx0 * (j + c_ny0 * (k - 1))] >= 0);
                if (cr && cl) { out_u2[idx] = 0.0; continue; }
            }
            if (b.lo[2] + k == 0 || b.lo[2] + k == gz_max) { out_u2[os2_0 + idx] = 0.0; continue; }
            size_t s = os2_1 + sidx_g1(2, i, j, k, b.nx, b.ny);
            double gp = (phi_g1[oc1 + cidx_g1(i, j, k, b.nx, b.ny)] - phi_g1[oc1 + cidx_g1(i, j, k - 1, b.nx, b.ny)]) * inv_h;
            out_u2[os2_0 + idx] = us2[s] + dp2[s] * gp;
        }
    }
}

// ============================================================================
// OtherInAdv Fluid Operations Kernels
// ============================================================================

__global__ void coef_update_g1_kernel(
    const double* rho_g0, const double* lb_g0, double* c_g1, double* dp_g1,
    double dt, int axis, int nx, int ny, int nz, bool is_l1, const Box3i* boxes, int n_l1)
{
    double inv_dt = 1.0 / dt;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o0 = is_l1 ? b.off_s_g0[axis] : 0;
        size_t o1 = is_l1 ? b.off_s_g1[axis] : 0;
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        int tot = nxa * nya * nza;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, t = idx / nxa, j = t % nya, k = t / nya;
            double r = rho_g0[o0 + idx], lb = lb_g0[o0 + idx];
            size_t sg1 = o1 + sidx_g1(axis, i, j, k, b.nx, b.ny);
            c_g1[sg1] = inv_dt * r + lb;
            dp_g1[sg1] = -1.0 / (r + dt * lb);
        }
    }
}

// Fill g=1 physical ghost values on C_sc and P_D (mirror along boundaries)
__global__ void mirror_g1_side_l0_kernel(double* f, int axis, int nx, int ny, int nz)
{
    int nxa = nx + (axis == 0), nya = ny + (axis == 1), nza = nz + (axis == 2);
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nya * nza; idx += gridDim.x * blockDim.x)
    {
        int j = idx % nya, k = idx / nya;
        f[sidx_g1(axis, -1, j, k, nx, ny)] = f[sidx_g1(axis, 0, j, k, nx, ny)];
        f[sidx_g1(axis, nxa, j, k, nx, ny)] = f[sidx_g1(axis, nxa - 1, j, k, nx, ny)];
    }
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nza; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nxa, k = idx / nxa;
        f[sidx_g1(axis, i, -1, k, nx, ny)] = f[sidx_g1(axis, i, 0, k, nx, ny)];
        f[sidx_g1(axis, i, nya, k, nx, ny)] = f[sidx_g1(axis, i, nya - 1, k, nx, ny)];
    }
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nya; idx += gridDim.x * blockDim.x)
    {
        int i = idx % nxa, j = idx / nxa;
        f[sidx_g1(axis, i, j, -1, nx, ny)] = f[sidx_g1(axis, i, j, 0, nx, ny)];
        f[sidx_g1(axis, i, j, nza, nx, ny)] = f[sidx_g1(axis, i, j, nza - 1, nx, ny)];
    }
}

__global__ void mirror_g1_side_l1_kernel(double* f, const Box3i* boxes, int n_l1, int axis)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t o1 = b.off_s_g1[axis];
        int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nya * nza; idx += gridDim.x * blockDim.x)
        {
            int j = idx % nya, k = idx / nya;
            f[o1 + sidx_g1(axis, -1, j, k, b.nx, b.ny)] = f[o1 + sidx_g1(axis, 0, j, k, b.nx, b.ny)];
            f[o1 + sidx_g1(axis, nxa, j, k, b.nx, b.ny)] = f[o1 + sidx_g1(axis, nxa - 1, j, k, b.nx, b.ny)];
        }
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nza; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, k = idx / nxa;
            f[o1 + sidx_g1(axis, i, -1, k, b.nx, b.ny)] = f[o1 + sidx_g1(axis, i, 0, k, b.nx, b.ny)];
            f[o1 + sidx_g1(axis, i, nya, k, b.nx, b.ny)] = f[o1 + sidx_g1(axis, i, nya - 1, k, b.nx, b.ny)];
        }
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < nxa * nya; idx += gridDim.x * blockDim.x)
        {
            int i = idx % nxa, j = idx / nxa;
            f[o1 + sidx_g1(axis, i, j, -1, b.nx, b.ny)] = f[o1 + sidx_g1(axis, i, j, 0, b.nx, b.ny)];
            f[o1 + sidx_g1(axis, i, j, nza, b.nx, b.ny)] = f[o1 + sidx_g1(axis, i, j, nza - 1, b.nx, b.ny)];
        }
    }
}

__global__ void compute_stoccurl_kernel(
    const double* u0, const double* u1, const double* u2,
    double* w0, double* w1, double* w2,
    int nx, int ny, int nz, double h, bool is_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    double fac = 0.25 / h;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o_c0 = is_l1 ? b.off_c_g0 : 0;
        size_t o0 = is_l1 ? b.off_s_g1[0] : 0;
        size_t o1 = is_l1 ? b.off_s_g1[1] : 0;
        size_t o2 = is_l1 ? b.off_s_g1[2] : 0;
        int tot = b.nx * b.ny * b.nz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            if (!is_l1 && l1_map != nullptr && l1_map[idx] >= 0) { w0[idx] = w1[idx] = w2[idx] = 0.0; continue; }
            int i = idx % b.nx, t = idx / b.nx, j = t % b.ny, k = t / b.ny;
            double du1_dz = fac * (u1[o1 + sidx_g1(1, i, j, k + 1, b.nx, b.ny)] + u1[o1 + sidx_g1(1, i, j + 1, k + 1, b.nx, b.ny)] -
                                   u1[o1 + sidx_g1(1, i, j, k - 1, b.nx, b.ny)] - u1[o1 + sidx_g1(1, i, j + 1, k - 1, b.nx, b.ny)]);
            double du2_dy = fac * (u2[o2 + sidx_g1(2, i, j + 1, k, b.nx, b.ny)] + u2[o2 + sidx_g1(2, i, j + 1, k + 1, b.nx, b.ny)] -
                                   u2[o2 + sidx_g1(2, i, j - 1, k, b.nx, b.ny)] - u2[o2 + sidx_g1(2, i, j - 1, k + 1, b.nx, b.ny)]);
            w0[o_c0 + idx] = du2_dy - du1_dz;

            double du0_dz = fac * (u0[o0 + sidx_g1(0, i, j, k + 1, b.nx, b.ny)] + u0[o0 + sidx_g1(0, i + 1, j, k + 1, b.nx, b.ny)] -
                                   u0[o0 + sidx_g1(0, i, j, k - 1, b.nx, b.ny)] - u0[o0 + sidx_g1(0, i + 1, j, k - 1, b.nx, b.ny)]);
            double du2_dx = fac * (u2[o2 + sidx_g1(2, i + 1, j, k, b.nx, b.ny)] + u2[o2 + sidx_g1(2, i + 1, j, k + 1, b.nx, b.ny)] -
                                   u2[o2 + sidx_g1(2, i - 1, j, k, b.nx, b.ny)] - u2[o2 + sidx_g1(2, i - 1, j, k + 1, b.nx, b.ny)]);
            w1[o_c0 + idx] = du0_dz - du2_dx;

            double du0_dy = fac * (u0[o0 + sidx_g1(0, i, j + 1, k, b.nx, b.ny)] + u0[o0 + sidx_g1(0, i + 1, j + 1, k, b.nx, b.ny)] -
                                   u0[o0 + sidx_g1(0, i, j - 1, k, b.nx, b.ny)] - u0[o0 + sidx_g1(0, i + 1, j - 1, k, b.nx, b.ny)]);
            double du1_dx = fac * (u1[o1 + sidx_g1(1, i + 1, j, k, b.nx, b.ny)] + u1[o1 + sidx_g1(1, i + 1, j + 1, k, b.nx, b.ny)] -
                                   u1[o1 + sidx_g1(1, i - 1, j, k, b.nx, b.ny)] - u1[o1 + sidx_g1(1, i - 1, j + 1, k, b.nx, b.ny)]);
            w2[o_c0 + idx] = du1_dx - du0_dy;
        }
    }
}

// Regrid repack kernel (verify L1 overlap unchanged during SSPRK3->FGMRES->Regrid test)
__global__ void regrid_repack_cell_l1_kernel(
    const double* old_l1, const int16_t* old_l1_map, const Box3i* old_boxes,
    double* new_l1, const Box3i* new_boxes, int n_new)
{
    for (int pat = blockIdx.y; pat < n_new; pat += gridDim.y)
    {
        Box3i nb = new_boxes[pat];
        int tot = nb.nx * nb.ny * nb.nz;
        for (int lin = blockIdx.x * blockDim.x + threadIdx.x; lin < tot; lin += gridDim.x * blockDim.x)
        {
            int li = lin % nb.nx, t = lin / nb.nx, lj = t % nb.ny, lk = t / nb.ny;
            int gi = nb.lo[0] + li, gj = nb.lo[1] + lj, gk = nb.lo[2] + lk;
            int q = lookup_l1_patch(gi, gj, gk, old_l1_map);
            double v = 0.0;
            if (q >= 0)
            {
                const Box3i& ob = old_boxes[q];
                v = old_l1[ob.off_c_g0 + static_cast<size_t>(gi - ob.lo[0]) +
                           static_cast<size_t>(ob.nx) * ((gj - ob.lo[1]) + static_cast<size_t>(ob.ny) * (gk - ob.lo[2]))];
            }
            new_l1[nb.off_c_g0 + lin] = v;
        }
    }
}

// Native modifyRhsForBcs on GPU:
//  - enforces homogeneous Dirichlet velocity rows on y_lo/y_up (u1) and z_lo/z_up (u2) -> rhs = 0
//  - accounts for x_lo & x_up TRACTION pressure BC shift: y := y - A(u=0, p_bc_inhom)
//    On x_lo (i=0): p[-1] has +2*gamma_inhom, where gamma_inhom = mu*du0_cur_dx + p_in_half
//      Grad P contribution = (p[0] - p[-1])/h -> shift in y = +(2*gamma_inhom)/h!
//    On x_up (i=c_nx0): p[c_nx0] has +2*gamma_inhom, where gamma_inhom = mu*du0_cur_dx
//      Grad P contribution = (p[c_nx0] - p[c_nx0-1])/h -> shift in y = -(2*gamma_inhom)/h!
__global__ void modify_rhs_bcs_gpu_kernel(
    double* rhs0, double* rhs1, double* rhs2,
    const double* u0_cur_g0,
    double p_in_half, double mu, double h, int nx, int ny, int nz,
    bool is_l1, const Box3i* boxes, int n_l1)
{
    int gy_max = is_l1 ? 2 * c_ny0 : ny, gz_max = is_l1 ? 2 * c_nz0 : nz;
    if (!is_l1)
    {
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < ny * nz; idx += gridDim.x * blockDim.x)
        {
            int j = idx % ny, k = idx / ny;
            double uc0 = u0_cur_g0[sidx_g0(0, 0, j, k, nx, ny)];
            double uc1 = u0_cur_g0[sidx_g0(0, 1, j, k, nx, ny)];
            double uc2 = u0_cur_g0[sidx_g0(0, 2, j, k, nx, ny)];
            double g_lo = mu * (2.0 * uc1 - 1.5 * uc0 - 0.5 * uc2) / h + p_in_half;
            rhs0[sidx_g0(0, 0, j, k, nx, ny)] += 2.0 * g_lo / h;

            double ur0 = u0_cur_g0[sidx_g0(0, nx,     j, k, nx, ny)];
            double ur1 = u0_cur_g0[sidx_g0(0, nx - 1, j, k, nx, ny)];
            double ur2 = u0_cur_g0[sidx_g0(0, nx - 2, j, k, nx, ny)];
            double g_up = mu * (-(2.0 * ur1 - 1.5 * ur0 - 0.5 * ur2) / h);
            rhs0[sidx_g0(0, nx, j, k, nx, ny)] -= 2.0 * g_up / h;
        }
    }
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o1 = is_l1 ? b.off_s_g0[1] : 0, o2 = is_l1 ? b.off_s_g0[2] : 0;
        if (b.lo[1] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; rhs1[o1 + sidx_g0(1, i, 0, k, b.nx, b.ny)] = 0.0; }
        if (b.lo[1] + b.ny == gy_max)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; rhs1[o1 + sidx_g0(1, i, b.ny, k, b.nx, b.ny)] = 0.0; }
        if (b.lo[2] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; rhs2[o2 + sidx_g0(2, i, j, 0, b.nx, b.ny)] = 0.0; }
        if (b.lo[2] + b.nz == gz_max)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; rhs2[o2 + sidx_g0(2, i, j, b.nz, b.nx, b.ny)] = 0.0; }
    }
}


template <typename RealT>
__global__ void fill_cf_cell_const_l1_kernel(
    const RealT* phi_l0, RealT* phi_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map)
{
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const int sx = b.nx + 2, sy = b.ny + 2, sz = b.nz + 2;
        const int tot = sx * sy * sz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % sx) - 1;
            int t = idx / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;
            if (li >= 0 && li < b.nx && lj >= 0 && lj < b.ny && lk >= 0 && lk < b.nz) continue;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            if (gi < 0 || gi >= 2 * c_nx0 || gj < 0 || gj >= 2 * c_ny0 || gk < 0 || gk >= 2 * c_nz0) continue;
            if (lookup_l1_patch(gi, gj, gk, l1_map) >= 0) continue;
            phi_l1[b.off_c_g1 + cidx_g1(li, lj, lk, b.nx, b.ny)] = phi_l0[cidx_g1(gi >> 1, gj >> 1, gk >> 1, c_nx0, c_ny0)];
        }
    }
}

template <typename RealT>
__global__ void apply_phi_phys_bc_l1_kernel(RealT* phi_l1, const Box3i* boxes, int n_l1)
{
    int d_hi0 = 2 * c_nx0 - 1, d_hi1 = 2 * c_ny0 - 1, d_hi2 = 2 * c_nz0 - 1;
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        Box3i b = boxes[pat];
        size_t o1 = b.off_c_g1;
        if (b.lo[0] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
            { int j = idx % b.ny, k = idx / b.ny; phi_l1[o1 + cidx_g1(-1, j, k, b.nx, b.ny)] = -phi_l1[o1 + cidx_g1(0, j, k, b.nx, b.ny)]; }
        if (b.hi[0] == d_hi0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.ny * b.nz; idx += gridDim.x * blockDim.x)
            { int j = idx % b.ny, k = idx / b.ny; phi_l1[o1 + cidx_g1(b.nx, j, k, b.nx, b.ny)] = -phi_l1[o1 + cidx_g1(b.nx - 1, j, k, b.nx, b.ny)]; }
        if (b.lo[1] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; phi_l1[o1 + cidx_g1(i, -1, k, b.nx, b.ny)] = phi_l1[o1 + cidx_g1(i, 0, k, b.nx, b.ny)]; }
        if (b.hi[1] == d_hi1)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; phi_l1[o1 + cidx_g1(i, b.ny, k, b.nx, b.ny)] = phi_l1[o1 + cidx_g1(i, b.ny - 1, k, b.nx, b.ny)]; }
        if (b.lo[2] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; phi_l1[o1 + cidx_g1(i, j, -1, b.nx, b.ny)] = phi_l1[o1 + cidx_g1(i, j, 0, b.nx, b.ny)]; }
        if (b.hi[2] == d_hi2)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; phi_l1[o1 + cidx_g1(i, j, b.nz, b.nx, b.ny)] = phi_l1[o1 + cidx_g1(i, j, b.nz - 1, b.nx, b.ny)]; }
    }
}

template <typename RealT>
__global__ void fill_cf_side_const_l1_kernel(
    const RealT* u_l0, RealT* u_l1, const Box3i* boxes, int n_l1, const int16_t* l1_map, int axis)
{
    int max_gi = 2 * c_nx0 - (axis == 0 ? 0 : 1);
    int max_gj = 2 * c_ny0 - (axis == 1 ? 0 : 1);
    int max_gk = 2 * c_nz0 - (axis == 2 ? 0 : 1);
    for (int pat = blockIdx.y; pat < n_l1; pat += gridDim.y)
    {
        const Box3i b = boxes[pat];
        const int nxa = b.nx + (axis == 0), nya = b.ny + (axis == 1), nza = b.nz + (axis == 2);
        const int sx = nxa + 2, sy = nya + 2, sz = nza + 2;
        const int tot = sx * sy * sz;
        for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < tot; idx += gridDim.x * blockDim.x)
        {
            int li = (idx % sx) - 1;
            int t = idx / sx;
            int lj = (t % sy) - 1;
            int lk = (t / sy) - 1;
            if (li >= 0 && li < nxa && lj >= 0 && lj < nya && lk >= 0 && lk < nza) continue;
            int gi = b.lo[0] + li, gj = b.lo[1] + lj, gk = b.lo[2] + lk;
            if (gi < 0 || gi > max_gi || gj < 0 || gj > max_gj || gk < 0 || gk > max_gk) continue;
            int q, ql, qj, qk;
            if (lookup_l1_side_owner(axis, gi, gj, gk, l1_map, boxes, q, ql, qj, qk)) continue;
            int ic = gi >> 1, jc = gj >> 1, kc = gk >> 1;
            RealT val;
            if (axis == 0 && (gi & 1)) val = RealT(0.5) * (u_l0[sidx_g1(0, ic, jc, kc, c_nx0, c_ny0)] + u_l0[sidx_g1(0, ic + 1, jc, kc, c_nx0, c_ny0)]);
            else if (axis == 1 && (gj & 1)) val = RealT(0.5) * (u_l0[sidx_g1(1, ic, jc, kc, c_nx0, c_ny0)] + u_l0[sidx_g1(1, ic, jc + 1, kc, c_nx0, c_ny0)]);
            else if (axis == 2 && (gk & 1)) val = RealT(0.5) * (u_l0[sidx_g1(2, ic, jc, kc, c_nx0, c_ny0)] + u_l0[sidx_g1(2, ic, jc, kc + 1, c_nx0, c_ny0)]);
            else val = u_l0[sidx_g1(axis, ic, jc, kc, c_nx0, c_ny0)];
            u_l1[b.off_s_g1[axis] + sidx_g1(axis, li, lj, lk, b.nx, b.ny)] = val;
        }
    }
}


// Zero homogeneous Dirichlet side velocities on walls/symmetry planes (g=0)
__global__ void impose_dirichlet_zero_g0_kernel(
    double* u1, double* u2,
    int nx, int ny, int nz,
    bool is_l1, const Box3i* boxes, int n_l1)
{
    int gy_max = is_l1 ? 2 * c_ny0 : ny;
    int gz_max = is_l1 ? 2 * c_nz0 : nz;
    for (int pat = blockIdx.y; pat < (is_l1 ? n_l1 : 1); pat += gridDim.y)
    {
        Box3i b{};
        if (is_l1) b = boxes[pat];
        else { b.lo[0] = b.lo[1] = b.lo[2] = 0; b.nx = nx; b.ny = ny; b.nz = nz; }
        size_t o1 = is_l1 ? b.off_s_g0[1] : 0;
        size_t o2 = is_l1 ? b.off_s_g0[2] : 0;
        if (b.lo[1] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; u1[o1 + sidx_g0(1, i, 0, k, b.nx, b.ny)] = 0.0; }
        if (b.lo[1] + b.ny == gy_max)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.nz; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, k = idx / b.nx; u1[o1 + sidx_g0(1, i, b.ny, k, b.nx, b.ny)] = 0.0; }
        if (b.lo[2] == 0)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; u2[o2 + sidx_g0(2, i, j, 0, b.nx, b.ny)] = 0.0; }
        if (b.lo[2] + b.nz == gz_max)
            for (int idx = blockIdx.x * blockDim.x + threadIdx.x; idx < b.nx * b.ny; idx += gridDim.x * blockDim.x)
            { int i = idx % b.nx, j = idx / b.nx; u2[o2 + sidx_g0(2, i, j, b.nz, b.nx, b.ny)] = 0.0; }
    }
}

__global__ void check_finite_kernel(const double* arr, size_t n, int* nonfinite_flag)
{
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < n;
         i += static_cast<size_t>(gridDim.x) * blockDim.x)
    {
        if (!isfinite(arr[i]))
        {
            atomicOr(nonfinite_flag, 1);
            return;
        }
    }
}

} // namespace c2_gpu

#endif
