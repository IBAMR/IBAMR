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

#include "ibamr_c2_gpu.h"
#include "c2_gpu_kernels.cuh"

#include <dlfcn.h>
#include <nvml.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace c2_gpu
{

static constexpr int kMaxCoarseLevels = 10;
static constexpr int kMaxL1Patches = 512;

struct CoarseLevelMeta
{
    int nx = 0, ny = 0, nz = 0;
    int ratio[3] = {1, 1, 1};
    double dx[3] = {0.0, 0.0, 0.0};
    size_t n_c_g0 = 0, n_c_g1 = 0;
    size_t n_s_g0[3] = {}, n_s_g1[3] = {};
    size_t off_c1_c_g1 = 0;
    size_t off_c1_s_g1[3] = {};
    size_t off_slot_c_g1 = 0;
    size_t off_slot_s_g1[3] = {};
};

struct NvmlApi
{
    bool tried = false;
    bool ok = false;
    void* lib = nullptr;
    using FnInit = nvmlReturn_t (*)();
    using FnGetHandle = nvmlReturn_t (*)(unsigned int, nvmlDevice_t*);
    using FnGetProcsV3 = nvmlReturn_t (*)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_t*);
    using FnGetProcsV2 = nvmlReturn_t (*)(nvmlDevice_t, unsigned int*, nvmlProcessInfo_v2_t*);
    FnInit init_fn = nullptr;
    FnGetHandle handle_fn = nullptr;
    FnGetProcsV3 procs_v3_fn = nullptr;
    FnGetProcsV2 procs_v2_fn = nullptr;
    nvmlDevice_t dev = nullptr;
    unsigned int target_pid = 0;
    std::vector<unsigned int> pre_pids;
};

struct C2Context
{
    pid_t pid = 0;
    C2GridConfig cfg{};
    int nx0 = 0, ny0 = 0, nz0 = 0;
    size_t n_l0_cells = 0;
    double dx0 = 0.0, dx1 = 0.0;
    double rho0 = 1.0, mu0 = 1.0, d_vel_const = -0.5;
    double domain_vol = 0.0;

    size_t max_slot_g0 = 0;
    size_t max_slot_g1 = 0;
    int max_l1_patches = kMaxL1Patches;

    // Device state (Stage 1)
    //  0..2: U_cur (g=0)   3..5: U_old (g=0)   6..8: U_new (g=0)
    //  9..11: C_sc (g=1)  12..14: P_D (g=1)   15..17: rho (g=0)
    // 18..20: RHS_u (g=0) 21: P_cur (g=0)     22: P_new (g=0)   23: RHS_p (g=0)
    double* d_state[24] = {};

    // Krylov + SSPRK3 alias arena (Stage 2)
    double* d_krylov_v[6][4] = {};
    double* d_krylov_z[5][4] = {};

    // Dual-prec FAC MG + coarse levels (Stage 3)
    double* d_mg_fp64[10] = {};
    float* d_mg_c_fp32[3] = {};
    float* d_mg_d_fp32[3] = {};
    double* d_cl_fp64_vec[3][3] = {};
    double* d_cl_fp64_c[3] = {};
    double* d_cl_fp64_d[3] = {};
    float* d_cl_fp32_c[3] = {};
    float* d_cl_fp32_d[3] = {};
    CoarseLevelMeta cl_meta[kMaxCoarseLevels] = {};
    int num_coarse_levels = 0;
    size_t cl_total_g1 = 0;
    size_t cl_c1_total_g1 = 0;

    int16_t* d_l1_map = nullptr;
    Box3i* d_l1_boxes = nullptr;
    double* d_wgt_c = nullptr;
    double* d_wgt_s[3] = {};
    double* d_red_scratch = nullptr;
    double* h_red_pinned = nullptr;

    // Stage 4: Persistent Brinkman penalty L_b[3] + pinned H2D/D2H staging + finiteness flag
    double* d_L_b[3] = {};
    double* h_pcie_pinned = nullptr;
    int* d_nonfinite = nullptr;

    // Dynamic L1 snapshot layout
    int n_l1 = 0;
    size_t l0_c_g0 = 0, l0_c_g1 = 0;
    size_t l1_c_g0 = 0, l1_c_g1 = 0;
    size_t l0_s_g0[3] = {}, l0_s_g1[3] = {}, l0_s_g3[3] = {};
    size_t l1_s_g0[3] = {}, l1_s_g1[3] = {}, l1_s_g2[3] = {}, l1_s_g3[3] = {};
    std::vector<Box3i> h_l1_boxes;
    std::vector<int16_t> h_l1_map;

    bool has_valid_grid = false;
    bool has_brinkman = false;
    int sticky_fp64_active = 0;
    C2Telemetry tel{};
    NvmlApi nvml{};

    size_t tot_c_g0() const { return l0_c_g0 + l1_c_g0; }
    size_t tot_c_g1() const { return l0_c_g1 + l1_c_g1; }
    size_t tot_s_g0(int a) const { return l0_s_g0[a] + l1_s_g0[a]; }
    size_t tot_s_g1(int a) const { return l0_s_g1[a] + l1_s_g1[a]; }
};

static void set_fb_reason(C2Context* ctx, const char* msg)
{
    if (!ctx) return;
    std::snprintf(ctx->tel.last_fb_reason, sizeof(ctx->tel.last_fb_reason), "%s", msg ? msg : "unknown");
}

static bool cuda_ok(C2Context* ctx, cudaError_t err, const char* what)
{
    if (err == cudaSuccess) return true;
    if (ctx)
    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "cuda_%s:%s", what, cudaGetErrorString(err));
        set_fb_reason(ctx, buf);
    }
    return false;
}

#define C2_CUDA_CALL(ctx, call, what) \
    do {                              \
        if (!cuda_ok((ctx), (call), (what))) return C2_ERR_CUDA; \
    } while (0)

#define C2_CUDA_LAUNCH_CHK(ctx, what) \
    do {                              \
        cudaError_t lerr__ = cudaGetLastError(); \
        if (!cuda_ok((ctx), lerr__, (what))) return C2_ERR_CUDA; \
    } while (0)

static bool ensure_nvml_init(C2Context* ctx)
{
    if (!ctx) return false;
    NvmlApi& api = ctx->nvml;
    if (!api.tried)
    {
        api.tried = true;
        api.lib = dlopen("libnvidia-ml.so.1", RTLD_LAZY);
        if (api.lib)
        {
            api.init_fn = reinterpret_cast<NvmlApi::FnInit>(dlsym(api.lib, "nvmlInit_v2"));
            api.handle_fn = reinterpret_cast<NvmlApi::FnGetHandle>(dlsym(api.lib, "nvmlDeviceGetHandleByIndex_v2"));
            api.procs_v3_fn = reinterpret_cast<NvmlApi::FnGetProcsV3>(
                dlsym(api.lib, "nvmlDeviceGetComputeRunningProcesses_v3"));
            api.procs_v2_fn = reinterpret_cast<NvmlApi::FnGetProcsV2>(
                dlsym(api.lib, "nvmlDeviceGetComputeRunningProcesses_v2"));
            if (api.init_fn && api.handle_fn && (api.procs_v3_fn || api.procs_v2_fn))
            {
                if (api.init_fn() == NVML_SUCCESS && api.handle_fn(0, &api.dev) == NVML_SUCCESS)
                    api.ok = true;
            }
        }
    }
    return api.ok;
}

static bool get_nvml_procs(C2Context* ctx, std::vector<nvmlProcessInfo_t>& out)
{
    if (!ensure_nvml_init(ctx)) return false;
    NvmlApi& api = ctx->nvml;
    out.resize(64);
    unsigned int count = static_cast<unsigned int>(out.size());
    nvmlReturn_t rc = api.procs_v3_fn ? api.procs_v3_fn(api.dev, &count, out.data())
                                      : api.procs_v2_fn(api.dev, &count, out.data());
    if (rc == NVML_ERROR_INSUFFICIENT_SIZE && count > 0)
    {
        out.resize(count);
        rc = api.procs_v3_fn ? api.procs_v3_fn(api.dev, &count, out.data())
                             : api.procs_v2_fn(api.dev, &count, out.data());
    }
    if (rc != NVML_SUCCESS) { out.clear(); return false; }
    out.resize(count);
    return true;
}

static void snapshot_pre_cuda_nvml_pids(C2Context* ctx)
{
    std::vector<nvmlProcessInfo_t> procs;
    if (!get_nvml_procs(ctx, procs)) return;
    ctx->nvml.pre_pids.clear();
    for (const auto& p : procs) ctx->nvml.pre_pids.push_back(p.pid);
}

static unsigned int s_cached_nvml_host_pid = 0;

static long query_nvml_proc_mib_dl(C2Context* ctx)
{
    std::vector<nvmlProcessInfo_t> procs;
    if (!get_nvml_procs(ctx, procs)) return -1;
    if (ctx->nvml.target_pid == 0 && s_cached_nvml_host_pid != 0)
        ctx->nvml.target_pid = s_cached_nvml_host_pid;

    // 1. Direct container/host PID match
    for (const auto& p : procs)
    {
        if ((static_cast<pid_t>(p.pid) == ctx->pid || (ctx->nvml.target_pid != 0 && p.pid == ctx->nvml.target_pid)) &&
            p.usedGpuMemory != static_cast<unsigned long long>(NVML_VALUE_NOT_AVAILABLE))
        {
            ctx->nvml.target_pid = p.pid;
            s_cached_nvml_host_pid = p.pid;
            return static_cast<long>(p.usedGpuMemory / (1024ULL * 1024ULL));
        }
    }
    // 2. Inside unshared Docker PID namespace: match newly created host NVML compute PID absent before CUDA init
    if (ctx->nvml.target_pid == 0)
    {
        for (const auto& p : procs)
        {
            if (std::find(ctx->nvml.pre_pids.begin(), ctx->nvml.pre_pids.end(), p.pid) == ctx->nvml.pre_pids.end() &&
                p.usedGpuMemory != static_cast<unsigned long long>(NVML_VALUE_NOT_AVAILABLE))
            {
                ctx->nvml.target_pid = p.pid;
                s_cached_nvml_host_pid = p.pid;
                return static_cast<long>(p.usedGpuMemory / (1024ULL * 1024ULL));
            }
        }
    }
    return -1;
}

static void sample_vram(C2Context* ctx, int stage_idx)
{
    size_t f = 0, t = 0;
    if (cudaMemGetInfo(&f, &t) == cudaSuccess)
    {
        double free_mib = static_cast<double>(f) / (1024.0 * 1024.0);
        if (ctx->tel.min_dev_free_mib <= 0.0 || free_mib < ctx->tel.min_dev_free_mib)
            ctx->tel.min_dev_free_mib = free_mib;
    }
    long mib = query_nvml_proc_mib_dl(ctx);
    ctx->tel.last_nvml_mib = mib;
    if (stage_idx >= 0 && stage_idx < 5) ctx->tel.stage_nvml_mib[stage_idx] = mib;
    if (mib > ctx->tel.peak_nvml_mib) ctx->tel.peak_nvml_mib = mib;
}

static void free_ctx_buffers(C2Context* ctx)
{
    if (!ctx) return;
    for (int i = 0; i < 24; ++i) { if (ctx->d_state[i]) cudaFree(ctx->d_state[i]); ctx->d_state[i] = nullptr; }
    for (int j = 0; j < 6; ++j)
        for (int c = 0; c < 4; ++c) { if (ctx->d_krylov_v[j][c]) cudaFree(ctx->d_krylov_v[j][c]); ctx->d_krylov_v[j][c] = nullptr; }
    for (int j = 0; j < 5; ++j)
        for (int c = 0; c < 4; ++c) { if (ctx->d_krylov_z[j][c]) cudaFree(ctx->d_krylov_z[j][c]); ctx->d_krylov_z[j][c] = nullptr; }
    for (int i = 0; i < 10; ++i) { if (ctx->d_mg_fp64[i]) cudaFree(ctx->d_mg_fp64[i]); ctx->d_mg_fp64[i] = nullptr; }
    for (int a = 0; a < 3; ++a)
    {
        if (ctx->d_mg_c_fp32[a]) cudaFree(ctx->d_mg_c_fp32[a]);
        if (ctx->d_mg_d_fp32[a]) cudaFree(ctx->d_mg_d_fp32[a]);
        ctx->d_mg_c_fp32[a] = nullptr; ctx->d_mg_d_fp32[a] = nullptr;
        for (int k = 0; k < 3; ++k) { if (ctx->d_cl_fp64_vec[k][a]) cudaFree(ctx->d_cl_fp64_vec[k][a]); ctx->d_cl_fp64_vec[k][a] = nullptr; }
        if (ctx->d_cl_fp64_c[a]) cudaFree(ctx->d_cl_fp64_c[a]);
        if (ctx->d_cl_fp64_d[a]) cudaFree(ctx->d_cl_fp64_d[a]);
        if (ctx->d_cl_fp32_c[a]) cudaFree(ctx->d_cl_fp32_c[a]);
        if (ctx->d_cl_fp32_d[a]) cudaFree(ctx->d_cl_fp32_d[a]);
        ctx->d_cl_fp64_c[a] = nullptr; ctx->d_cl_fp64_d[a] = nullptr;
        ctx->d_cl_fp32_c[a] = nullptr; ctx->d_cl_fp32_d[a] = nullptr;
        if (ctx->d_wgt_s[a]) cudaFree(ctx->d_wgt_s[a]);
        ctx->d_wgt_s[a] = nullptr;
        if (ctx->d_L_b[a]) cudaFree(ctx->d_L_b[a]);
        ctx->d_L_b[a] = nullptr;
    }
    if (ctx->d_l1_map) { cudaFree(ctx->d_l1_map); ctx->d_l1_map = nullptr; }
    if (ctx->d_l1_boxes) { cudaFree(ctx->d_l1_boxes); ctx->d_l1_boxes = nullptr; }
    if (ctx->d_wgt_c) { cudaFree(ctx->d_wgt_c); ctx->d_wgt_c = nullptr; }
    if (ctx->d_red_scratch) { cudaFree(ctx->d_red_scratch); ctx->d_red_scratch = nullptr; }
    if (ctx->h_red_pinned) { cudaFreeHost(ctx->h_red_pinned); ctx->h_red_pinned = nullptr; }
    if (ctx->h_pcie_pinned) { cudaFreeHost(ctx->h_pcie_pinned); ctx->h_pcie_pinned = nullptr; }
    if (ctx->d_nonfinite) { cudaFree(ctx->d_nonfinite); ctx->d_nonfinite = nullptr; }
    if (ctx->nvml.lib) { dlclose(ctx->nvml.lib); ctx->nvml.lib = nullptr; }
}

static bool build_coarse_mg_chain(C2Context* ctx)
{
    int dims[kMaxCoarseLevels][3] = {};
    int ratios[kMaxCoarseLevels][3] = {};
    double dxs[kMaxCoarseLevels][3] = {};

    dims[0][0] = ctx->nx0; dims[0][1] = ctx->ny0; dims[0][2] = ctx->nz0;
    dxs[0][0] = ctx->dx0;  dxs[0][1] = ctx->dx0;  dxs[0][2] = ctx->dx0;

    int lv = 0;
    while (lv + 1 < kMaxCoarseLevels)
    {
        size_t ncells = static_cast<size_t>(dims[lv][0]) * dims[lv][1] * dims[lv][2];
        if (ncells <= 256 || (dims[lv][0] <= 8 && dims[lv][1] <= 8 && dims[lv][2] <= 8))
            break;

        int r[3] = {1, 1, 1};
        bool any2 = false;
        for (int d = 0; d < 3; ++d)
        {
            if (dims[lv][d] % 2 == 0 && dims[lv][d] / 2 >= 4)
            {
                r[d] = 2;
                any2 = true;
            }
        }
        if (any2)
        {
            double max_next_dx = 0.0;
            double min_curr_dx = std::min({dxs[lv][0], dxs[lv][1], dxs[lv][2]});
            for (int d = 0; d < 3; ++d)
                max_next_dx = std::max(max_next_dx, dxs[lv][d] * r[d]);
            if (max_next_dx / min_curr_dx > 4.1)
            {
                for (int d = 0; d < 3; ++d)
                {
                    if (r[d] == 1 && dims[lv][d] > 7 && (max_next_dx / dxs[lv][d] >= 3.9))
                    {
                        for (int f : {3, 5, 7})
                        {
                            if (dims[lv][d] % f == 0 && dims[lv][d] / f >= 4)
                            {
                                r[d] = f;
                                break;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            bool progressed = false;
            for (int f : {3, 5, 7})
            {
                for (int d = 0; d < 3; ++d)
                {
                    if (dims[lv][d] % f == 0 && dims[lv][d] / f >= 4)
                    {
                        r[d] = f;
                        progressed = true;
                    }
                }
                if (progressed) break;
            }
            if (!progressed) break;
        }

        for (int d = 0; d < 3; ++d)
        {
            ratios[lv][d] = r[d];
            dims[lv + 1][d] = dims[lv][d] / r[d];
            dxs[lv + 1][d] = dxs[lv][d] * r[d];
        }
        lv++;
    }
    ratios[lv][0] = ratios[lv][1] = ratios[lv][2] = 1;
    ctx->num_coarse_levels = lv + 1;
    if (ctx->num_coarse_levels < 2 || ctx->num_coarse_levels > kMaxCoarseLevels) return false;

    size_t off_c = 0, off_s[3] = {0, 0, 0};
    size_t c1_off_c = 0, c1_off_s[3] = {0, 0, 0};
    for (int i = 0; i < ctx->num_coarse_levels; ++i)
    {
        CoarseLevelMeta& m = ctx->cl_meta[i];
        m.nx = dims[i][0]; m.ny = dims[i][1]; m.nz = dims[i][2];
        for (int a = 0; a < 3; ++a)
        {
            m.ratio[a] = ratios[i][a];
            m.dx[a] = dxs[i][a];
        }
        m.n_c_g0 = static_cast<size_t>(m.nx) * m.ny * m.nz;
        m.n_c_g1 = static_cast<size_t>(m.nx + 2) * (m.ny + 2) * (m.nz + 2);
        m.off_slot_c_g1 = off_c;
        off_c += m.n_c_g1;
        m.off_c1_c_g1 = (i >= 1) ? c1_off_c : 0;
        if (i >= 1) c1_off_c += m.n_c_g1;

        for (int a = 0; a < 3; ++a)
        {
            int nx1 = m.nx + 2 + (a == 0 ? 1 : 0);
            int ny1 = m.ny + 2 + (a == 1 ? 1 : 0);
            int nz1 = m.nz + 2 + (a == 2 ? 1 : 0);
            m.n_s_g0[a] = static_cast<size_t>(m.nx + (a == 0 ? 1 : 0)) * (m.ny + (a == 1 ? 1 : 0)) * (m.nz + (a == 2 ? 1 : 0));
            m.n_s_g1[a] = static_cast<size_t>(nx1) * ny1 * nz1;
            m.off_slot_s_g1[a] = off_s[a];
            off_s[a] += m.n_s_g1[a];
            m.off_c1_s_g1[a] = (i >= 1) ? c1_off_s[a] : 0;
            if (i >= 1) c1_off_s[a] += m.n_s_g1[a];
        }
    }
    ctx->cl_total_g1 = std::max(off_c, std::max({off_s[0], off_s[1], off_s[2]}));
    ctx->cl_c1_total_g1 = std::max(c1_off_c, std::max({c1_off_s[0], c1_off_s[1], c1_off_s[2]}));
    return true;
}

static int gpu_wgt_norm4(C2Context* ctx, double* const v[4], double* out_norm)
{
    const int nblk = 256;
    wgt_dot_ab_kernel<<<nblk, 256>>>(
        v[0], v[1], v[2], v[3], v[0], v[1], v[2], v[3],
        ctx->d_wgt_s[0], ctx->d_wgt_s[1], ctx->d_wgt_s[2], ctx->d_wgt_c,
        ctx->tot_s_g0(0), ctx->tot_s_g0(1), ctx->tot_s_g0(2), ctx->tot_c_g0(),
        ctx->d_red_scratch);
    C2_CUDA_LAUNCH_CHK(ctx, "wgt_dot_ab_norm");
    C2_CUDA_CALL(ctx, cudaMemcpy(ctx->h_red_pinned, ctx->d_red_scratch, nblk * sizeof(double), cudaMemcpyDeviceToHost), "norm4_d2h");
    double s = 0.0;
    for (int i = 0; i < nblk; ++i) s += ctx->h_red_pinned[i];
    if (!std::isfinite(s)) return C2_ERR_NONFINITE;
    *out_norm = std::sqrt(std::max(0.0, s));
    return C2_OK;
}

static int gpu_wgt_dot4(C2Context* ctx, double* const a[4], double* const b[4], double* out_dot)
{
    const int nblk = 256;
    wgt_dot_ab_kernel<<<nblk, 256>>>(
        a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3],
        ctx->d_wgt_s[0], ctx->d_wgt_s[1], ctx->d_wgt_s[2], ctx->d_wgt_c,
        ctx->tot_s_g0(0), ctx->tot_s_g0(1), ctx->tot_s_g0(2), ctx->tot_c_g0(),
        ctx->d_red_scratch);
    C2_CUDA_LAUNCH_CHK(ctx, "wgt_dot_ab");
    C2_CUDA_CALL(ctx, cudaMemcpy(ctx->h_red_pinned, ctx->d_red_scratch, nblk * sizeof(double), cudaMemcpyDeviceToHost), "dot4_d2h");
    double s = 0.0;
    for (int i = 0; i < nblk; ++i) s += ctx->h_red_pinned[i];
    if (!std::isfinite(s)) return C2_ERR_NONFINITE;
    *out_dot = s;
    return C2_OK;
}

template <typename RealT>
static int sync_composite_side(C2Context* ctx, RealT* l0_a[3], RealT* l1_a[3], bool sync_l1_to_l0)
{
    const dim3 l1g(16, std::max(1, ctx->n_l1));
    if (sync_l1_to_l0 && ctx->n_l1 > 0)
    {
        for (int a = 0; a < 3; ++a)
            restrict_l1_to_l0_side_kernel<<<l1g, 256>>>(l1_a[a], l0_a[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
    }
    enforce_homog_phys_bcs_side_kernel<<<dim3(64, 1), 256>>>(
        l0_a[0], l0_a[1], l0_a[2], ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1);
    if (ctx->n_l1 > 0)
    {
        for (int a = 0; a < 3; ++a)
            fill_l1_side_ghosts_kernel<<<l1g, 256>>>(l1_a[a], l0_a[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        enforce_homog_phys_bcs_side_kernel<<<l1g, 256>>>(
            l1_a[0], l1_a[1], l1_a[2], ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1);
    }
    C2_CUDA_LAUNCH_CHK(ctx, "sync_composite_side");
    return C2_OK;
}

template <typename RealT>
static int sync_composite_cell(C2Context* ctx, RealT* p_l0, RealT* p_l1, bool sync_l1_to_l0)
{
    const dim3 l1g(16, std::max(1, ctx->n_l1));
    if (sync_l1_to_l0 && ctx->n_l1 > 0)
        restrict_l1_to_l0_cell_kernel<<<l1g, 256>>>(p_l1, p_l0, ctx->d_l1_boxes, ctx->n_l1);
    enforce_homog_phys_bcs_cell_kernel<<<dim3(64, 1), 256>>>(p_l0, ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1);
    if (ctx->n_l1 > 0)
    {
        fill_l1_cell_ghosts_kernel<<<l1g, 256>>>(p_l1, p_l0, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        enforce_homog_phys_bcs_cell_kernel<<<l1g, 256>>>(
            p_l1, ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1);
    }
    C2_CUDA_LAUNCH_CHK(ctx, "sync_composite_cell");
    return C2_OK;
}

static int fill_ghost_side_ssprk_g3(
    C2Context* ctx,
    double* d_g0_vec[3],
    double* d_l0_g3[3],
    double* d_l1_g3[3],
    bool is_rho)
{
    const int g = 3;
    const dim3 l1g(16, std::max(1, ctx->n_l1));
    for (int a = 0; a < 3; ++a)
    {
        pack_g0_to_g_side_l0_kernel<<<64, 256>>>(d_g0_vec[a], d_l0_g3[a], a, ctx->nx0, ctx->ny0, ctx->nz0, g);
        if (ctx->n_l1 > 0)
            pack_g0_to_g_side_l1_kernel<<<l1g, 256>>>(d_g0_vec[a] + ctx->l0_s_g0[a], d_l1_g3[a], ctx->d_l1_boxes, ctx->n_l1, a, g);
    }
    for (int a = 0; a < 3; ++a)
    {
        if (ctx->n_l1 > 0)
        {
            sync_l1_intralevel_side_kernel<<<l1g, 256>>>(d_l1_g3[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, g);
            restrict_l1_to_l0_side_kernel<<<l1g, 256>>>(d_l1_g3[a], d_l0_g3[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, g);
        }
        cart_extrap_phys_side_kernel<<<64, 256>>>(
            d_l0_g3[a], a, ctx->nx0, ctx->ny0, ctx->nz0, 0, ctx->nx0 - 1, 0, ctx->ny0 - 1, 0, ctx->nz0 - 1,
            0, 0, 0, ctx->nx0, ctx->ny0, ctx->nz0, g, 0, 0);
        if (ctx->n_l1 > 0)
        {
            refine_l0_to_l1_cons_linear_side_g3_kernel<<<l1g, 256>>>(
                d_l0_g3[a], d_l1_g3[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            cart_extrap_phys_side_l1_kernel<<<l1g, 256>>>(d_l1_g3[a], a, ctx->d_l1_boxes, ctx->n_l1, g, 0);
        }
        apply_channel_robin_side_phys_kernel<<<32, 256>>>(
            d_l0_g3[a], a, 0, 0, 0, ctx->nx0, ctx->ny0, ctx->nz0, ctx->nx0 - 1, ctx->ny0 - 1, ctx->nz0 - 1, g, is_rho, 0);
        if (ctx->n_l1 > 0)
            apply_channel_robin_side_phys_l1_kernel<<<l1g, 256>>>(d_l1_g3[a], ctx->d_l1_boxes, ctx->n_l1, a, g, is_rho);
        if (!is_rho)
        {
            apply_channel_robin_side_codim23_kernel<<<32, 256>>>(
                d_l0_g3[a], a, ctx->nx0, ctx->ny0, ctx->nz0, 0, 0, 0, ctx->nx0 - 1, ctx->ny0 - 1, ctx->nz0 - 1, g, 0);
            if (ctx->n_l1 > 0)
                apply_channel_robin_side_codim23_l1_kernel<<<l1g, 256>>>(d_l1_g3[a], ctx->d_l1_boxes, ctx->n_l1, a, g);
        }
        if (ctx->n_l1 > 0)
            sync_l1_intralevel_side_kernel<<<l1g, 256>>>(d_l1_g3[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, g);
    }
    if (!is_rho)
    {
        enforce_div_free_phys_bdry_l0_kernel<<<32, 256>>>(d_l0_g3[0], d_l0_g3[1], d_l0_g3[2], ctx->nx0, ctx->ny0, ctx->nz0, g, true);
        if (ctx->n_l1 > 0)
        {
            enforce_div_free_phys_bdry_l1_kernel<<<l1g, 256>>>(d_l1_g3[0], d_l1_g3[1], d_l1_g3[2], ctx->d_l1_boxes, ctx->n_l1, g);
            enforce_div_free_cf_l1_kernel<<<l1g, 256>>>(d_l1_g3[0], d_l1_g3[1], d_l1_g3[2], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, g);
        }
    }
    C2_CUDA_LAUNCH_CHK(ctx, "fill_ghost_side_ssprk_g3");
    return C2_OK;
}

static int fill_native_stokes_ghosts(C2Context* ctx, double* const x_g0[4], bool enforce_phys_bc)
{
    const dim3 l1g(16, std::max(1, ctx->n_l1));
    double* u_l0[3] = {ctx->d_mg_fp64[0], ctx->d_mg_fp64[1], ctx->d_mg_fp64[2]};
    double* u_l1[3] = {u_l0[0] + ctx->l0_s_g1[0], u_l0[1] + ctx->l0_s_g1[1], u_l0[2] + ctx->l0_s_g1[2]};
    double* w_l1[3] = {ctx->d_mg_fp64[6], ctx->d_mg_fp64[7], ctx->d_mg_fp64[8]};
    double* p_l0 = ctx->d_mg_fp64[3];
    double* p_l1 = p_l0 + ctx->l0_c_g1;

    for (int a = 0; a < 3; ++a)
    {
        pack_g0_to_g_side_l0_kernel<<<64, 256>>>(x_g0[a], u_l0[a], a, ctx->nx0, ctx->ny0, ctx->nz0, 1);
        if (ctx->n_l1 > 0)
        {
            pack_g0_to_g_side_l1_kernel<<<l1g, 256>>>(x_g0[a] + ctx->l0_s_g0[a], u_l1[a], ctx->d_l1_boxes, ctx->n_l1, a, 1);
            sync_l1_intralevel_side_kernel<<<l1g, 256>>>(u_l1[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
            restrict_l1_to_l0_side_kernel<<<l1g, 256>>>(u_l1[a], u_l0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
        }
    }
    pack_g0_to_g1_cell_l0_kernel<<<64, 256>>>(x_g0[3], p_l0, ctx->nx0, ctx->ny0, ctx->nz0);
    if (ctx->n_l1 > 0)
    {
        pack_g0_to_g1_cell_l1_kernel<<<l1g, 256>>>(x_g0[3] + ctx->l0_c_g0, p_l1, ctx->d_l1_boxes, ctx->n_l1);
        sync_l1_intralevel_cell_kernel<<<l1g, 256>>>(p_l1, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        restrict_l1_to_l0_cell_kernel<<<l1g, 256>>>(p_l1, p_l0, ctx->d_l1_boxes, ctx->n_l1);
    }

    for (int a = 0; a < 3; ++a)
        cart_extrap_phys_side_kernel<<<64, 256>>>(
            u_l0[a], a, ctx->nx0, ctx->ny0, ctx->nz0, 0, ctx->nx0 - 1, 0, ctx->ny0 - 1, 0, ctx->nz0 - 1,
            0, 0, 0, ctx->nx0, ctx->ny0, ctx->nz0, 1, 1, 0);
    cart_extrap_phys_cell_l0_kernel<<<64, 256>>>(p_l0, ctx->nx0, ctx->ny0, ctx->nz0);

    if (ctx->n_l1 > 0)
    {
        quad_cf_tangential_cell_l1_kernel<<<l1g, 256>>>(p_l0, p_l1, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        quad_cf_normal_cell_l1_kernel<<<l1g, 256>>>(p_l1, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        for (int a = 0; a < 3; ++a)
            C2_CUDA_CALL(ctx, cudaMemcpy(w_l1[a], u_l1[a], ctx->l1_s_g1[a] * sizeof(double), cudaMemcpyDeviceToDevice), "stokes_cf_cp1");
        quad_cf_tangential_side_l1_kernel<<<l1g, 256>>>(
            u_l0[0], u_l0[1], u_l0[2], w_l1[0], w_l1[1], w_l1[2], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        for (int a = 0; a < 3; ++a)
            C2_CUDA_CALL(ctx, cudaMemcpy(u_l1[a], w_l1[a], ctx->l1_s_g1[a] * sizeof(double), cudaMemcpyDeviceToDevice), "stokes_cf_cp2");
        quad_cf_normal_side_l1_kernel<<<l1g, 256>>>(
            w_l1[0], w_l1[1], w_l1[2], u_l1[0], u_l1[1], u_l1[2], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
    }

    if (enforce_phys_bc)
    {
        for (int a = 0; a < 3; ++a)
        {
            if (ctx->n_l1 > 0)
                cart_extrap_phys_side_l1_kernel<<<l1g, 256>>>(u_l1[a], a, ctx->d_l1_boxes, ctx->n_l1, 1, 1);
            apply_channel_robin_side_phys_kernel<<<32, 256>>>(
                u_l0[a], a, 0, 0, 0, ctx->nx0, ctx->ny0, ctx->nz0, ctx->nx0 - 1, ctx->ny0 - 1, ctx->nz0 - 1, 1, false, 0);
            if (ctx->n_l1 > 0)
                apply_channel_robin_side_phys_l1_kernel<<<l1g, 256>>>(u_l1[a], ctx->d_l1_boxes, ctx->n_l1, a, 1, false);
            apply_channel_robin_side_codim23_kernel<<<32, 256>>>(
                u_l0[a], a, ctx->nx0, ctx->ny0, ctx->nz0, 0, 0, 0, ctx->nx0 - 1, ctx->ny0 - 1, ctx->nz0 - 1, 1, 0);
            if (ctx->n_l1 > 0)
            {
                apply_channel_robin_side_codim23_l1_kernel<<<l1g, 256>>>(u_l1[a], ctx->d_l1_boxes, ctx->n_l1, a, 1);
                sync_l1_intralevel_side_kernel<<<l1g, 256>>>(u_l1[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
            }
        }
        enforce_div_free_phys_bdry_l0_kernel<<<32, 256>>>(u_l0[0], u_l0[1], u_l0[2], ctx->nx0, ctx->ny0, ctx->nz0, 1, false);
        apply_pressure_phys_bc_l0_kernel<<<32, 256>>>(
            p_l0, u_l0[0], ctx->d_state[0], ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, ctx->mu0, 0.0, true);
        if (ctx->n_l1 > 0)
            apply_pressure_neumann_l1_kernel<<<l1g, 256>>>(p_l1, ctx->d_l1_boxes, ctx->n_l1);
    }
    C2_CUDA_LAUNCH_CHK(ctx, "fill_native_stokes_ghosts");
    return C2_OK;
}

static int apply_stokes_fp64(C2Context* ctx, double* const x_g0[4], double* const y_g0[4])
{
    const dim3 l1g(16, std::max(1, ctx->n_l1));
    int rc = fill_native_stokes_ghosts(ctx, x_g0, true);
    if (rc != C2_OK) return rc;

    double* u_l0[3] = {ctx->d_mg_fp64[0], ctx->d_mg_fp64[1], ctx->d_mg_fp64[2]};
    double* u_l1[3] = {u_l0[0] + ctx->l0_s_g1[0], u_l0[1] + ctx->l0_s_g1[1], u_l0[2] + ctx->l0_s_g1[2]};
    double* p_l0 = ctx->d_mg_fp64[3];
    double* p_l1 = p_l0 + ctx->l0_c_g1;
    double* c_l0[3] = {ctx->d_state[9], ctx->d_state[10], ctx->d_state[11]};
    double* c_l1[3] = {c_l0[0] + ctx->l0_s_g1[0], c_l0[1] + ctx->l0_s_g1[1], c_l0[2] + ctx->l0_s_g1[2]};
    double* y_l1[3] = {y_g0[0] + ctx->l0_s_g0[0], y_g0[1] + ctx->l0_s_g0[1], y_g0[2] + ctx->l0_s_g0[2]};

    apply_momentum_stokes_kernel<<<dim3(64, 1), 256>>>(
        u_l0[0], u_l0[1], u_l0[2], p_l0, c_l0[0], c_l0[1], c_l0[2],
        ctx->d_vel_const, 1.0, y_g0[0], y_g0[1], y_g0[2],
        ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, true, false, true, false, nullptr, 1, ctx->d_l1_map);
    if (ctx->n_l1 > 0)
    {
        apply_momentum_stokes_kernel<<<l1g, 256>>>(
            u_l1[0], u_l1[1], u_l1[2], p_l1, c_l1[0], c_l1[1], c_l1[2],
            ctx->d_vel_const, 1.0, y_l1[0], y_l1[1], y_l1[2],
            ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, true, false, false, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        for (int a = 0; a < 3; ++a)
        {
            sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(y_l1[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            outerside_restrict_cf_l1_to_l0_side_kernel<<<l1g, 256>>>(y_l1[a], y_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }
    copy_dirichlet_bc_g0_kernel<<<dim3(32, 1), 256>>>(
        x_g0[1], x_g0[2], y_g0[1], y_g0[2], ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1);
    if (ctx->n_l1 > 0)
    {
        copy_dirichlet_bc_g0_kernel<<<l1g, 256>>>(
            x_g0[1] + ctx->l0_s_g0[1], x_g0[2] + ctx->l0_s_g0[2], y_l1[1], y_l1[2],
            ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1);
    }

    compute_div_kernel<<<dim3(64, 1), 256>>>(
        u_l0[0], u_l0[1], u_l0[2], y_g0[3], -1.0, ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, false, nullptr, 1, ctx->d_l1_map);
    if (ctx->n_l1 > 0)
    {
        compute_div_kernel<<<l1g, 256>>>(
            u_l1[0], u_l1[1], u_l1[2], y_g0[3] + ctx->l0_c_g0, -1.0,
            ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
    }
    C2_CUDA_LAUNCH_CHK(ctx, "apply_stokes_fp64");
    return C2_OK;
}

template <typename RealT>
static int run_pressure_coarse_wcycle(
    C2Context* ctx, int lv, const RealT* const d_l0[3],
    RealT* cl_x[kMaxCoarseLevels], RealT* cl_b[kMaxCoarseLevels],
    RealT* cl_r[kMaxCoarseLevels], RealT* cl_d[kMaxCoarseLevels][3])
{
    const CoarseLevelMeta& m = ctx->cl_meta[lv];
    const RealT* lv_d[3] = {lv == 0 ? d_l0[0] : cl_d[lv][0], lv == 0 ? d_l0[1] : cl_d[lv][1], lv == 0 ? d_l0[2] : cl_d[lv][2]};
    if (lv == ctx->num_coarse_levels - 1)
    {
        p_poisson_c6_pcg_kernel<<<1, 256>>>(
            cl_x[lv], cl_b[lv], lv_d[0], lv_d[1], lv_d[2],
            m.nx, m.ny, m.nz, RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]));
        C2_CUDA_LAUNCH_CHK(ctx, "p_poisson_coarsest");
        return C2_OK;
    }

    int ns = 2;
    for (int sw = 0; sw < ns; ++sw)
        for (int col = 0; col < 2; ++col)
            p_poisson_rbgs_kernel<<<dim3(32, 1), 256>>>(
                cl_x[lv], cl_b[lv], lv_d[0], lv_d[1], lv_d[2], m.nx, m.ny, m.nz,
                RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]), col, false, nullptr, 1, nullptr);

    enforce_homog_phys_bcs_cell_kernel<<<dim3(32, 1), 256>>>(cl_x[lv], m.nx, m.ny, m.nz, false, nullptr, 1);
    p_poisson_residual_kernel<<<dim3(32, 1), 256>>>(
        cl_x[lv], cl_b[lv], lv_d[0], lv_d[1], lv_d[2], cl_r[lv], m.nx, m.ny, m.nz,
        RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]), false, nullptr, 1, nullptr);

    const CoarseLevelMeta& mc = ctx->cl_meta[lv + 1];
    restrict_cell_kernel<<<16, 256>>>(
        cl_r[lv], cl_b[lv + 1], cl_x[lv + 1], mc.nx, mc.ny, mc.nz, m.nx, m.ny,
        m.ratio[0], m.ratio[1], m.ratio[2]);

    int gamma = (lv <= 2) ? 2 : 1;
    for (int g = 0; g < gamma; ++g)
    {
        int rc = run_pressure_coarse_wcycle<RealT>(ctx, lv + 1, d_l0, cl_x, cl_b, cl_r, cl_d);
        if (rc != C2_OK) return rc;
    }

    prolong_cell_add_kernel<<<dim3(32, 1), 256>>>(
        cl_x[lv + 1], cl_x[lv], m.nx, m.ny, m.nz, mc.nx, mc.ny, m.ratio[0], m.ratio[1], m.ratio[2], false, nullptr, 1);

    for (int sw = 0; sw < ns; ++sw)
        for (int col = 1; col >= 0; --col)
            p_poisson_rbgs_kernel<<<dim3(32, 1), 256>>>(
                cl_x[lv], cl_b[lv], lv_d[0], lv_d[1], lv_d[2], m.nx, m.ny, m.nz,
                RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]), col, false, nullptr, 1, nullptr);

    C2_CUDA_LAUNCH_CHK(ctx, "p_wcycle");
    return C2_OK;
}

template <typename RealT>
static int run_pressure_fac_mg(
    C2Context* ctx, RealT* phi_l0, RealT* phi_l1, const RealT* b_l0, const RealT* b_l1, RealT* r_l0, RealT* r_l1,
    const RealT* const d_l0[3], const RealT* const d_l1[3], RealT* cl_x[kMaxCoarseLevels], RealT* cl_b[kMaxCoarseLevels],
    RealT* cl_r[kMaxCoarseLevels], RealT* cl_d[kMaxCoarseLevels][3], int n_vcycles)
{
    dim3 gr_l1(16, std::max(1, ctx->n_l1));
    for (int cyc = 0; cyc < n_vcycles; ++cyc)
    {
        if (ctx->n_l1 > 0)
        {
            for (int sw = 0; sw < 2; ++sw)
            {
                int rc = sync_composite_cell(ctx, phi_l0, phi_l1, true);
                if (rc != C2_OK) return rc;
                for (int col = 0; col < 2; ++col)
                {
                    p_poisson_rbgs_kernel<<<dim3(64, 1), 256>>>(
                        phi_l0, b_l0, d_l0[0], d_l0[1], d_l0[2], ctx->nx0, ctx->ny0, ctx->nz0,
                        RealT(ctx->dx0), RealT(ctx->dx0), RealT(ctx->dx0), col, false, nullptr, 1, ctx->d_l1_map);
                    p_poisson_rbgs_kernel<<<gr_l1, 256>>>(
                        phi_l1, b_l1, d_l1[0], d_l1[1], d_l1[2], ctx->nx0, ctx->ny0, ctx->nz0,
                        RealT(ctx->dx1), RealT(ctx->dx1), RealT(ctx->dx1), col, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
                    rc = sync_composite_cell(ctx, phi_l0, phi_l1, false);
                    if (rc != C2_OK) return rc;
                }
            }
            int rc = sync_composite_cell(ctx, phi_l0, phi_l1, true);
            if (rc != C2_OK) return rc;
            p_poisson_residual_kernel<<<gr_l1, 256>>>(
                phi_l1, b_l1, d_l1[0], d_l1[1], d_l1[2], r_l1, ctx->nx0, ctx->ny0, ctx->nz0,
                RealT(ctx->dx1), RealT(ctx->dx1), RealT(ctx->dx1), true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
            p_poisson_residual_kernel<<<dim3(64, 1), 256>>>(
                phi_l0, b_l0, d_l0[0], d_l0[1], d_l0[2], r_l0, ctx->nx0, ctx->ny0, ctx->nz0,
                RealT(ctx->dx0), RealT(ctx->dx0), RealT(ctx->dx0), false, nullptr, 1, ctx->d_l1_map);
            restrict_l1_to_l0_cell_kernel<<<gr_l1, 256>>>(r_l1, r_l0, ctx->d_l1_boxes, ctx->n_l1);
        }
        else
        {
            p_poisson_residual_kernel<<<dim3(64, 1), 256>>>(
                phi_l0, b_l0, d_l0[0], d_l0[1], d_l0[2], r_l0, ctx->nx0, ctx->ny0, ctx->nz0,
                RealT(ctx->dx0), RealT(ctx->dx0), RealT(ctx->dx0), false, nullptr, 1, ctx->d_l1_map);
        }

        for (int lv = 0; lv < ctx->num_coarse_levels; ++lv)
            C2_CUDA_CALL(ctx, cudaMemset(cl_x[lv], 0, ctx->cl_meta[lv].n_c_g1 * sizeof(RealT)), "p_fac_memset");
        cl_b[0] = r_l0;
        int rc = run_pressure_coarse_wcycle<RealT>(ctx, 0, d_l0, cl_x, cl_b, cl_r, cl_d);
        if (rc != C2_OK) return rc;

        prolong_cell_add_kernel<<<dim3(32, 1), 256>>>(
            cl_x[0], phi_l0, ctx->nx0, ctx->ny0, ctx->nz0, ctx->nx0, ctx->ny0, 1, 1, 1, false, nullptr, 1);
        if (ctx->n_l1 > 0)
        {
            prolong_cell_add_kernel<<<gr_l1, 256>>>(
                cl_x[0], phi_l1, 2 * ctx->nx0, 2 * ctx->ny0, 2 * ctx->nz0, ctx->nx0, ctx->ny0,
                2, 2, 2, true, ctx->d_l1_boxes, ctx->n_l1);
            for (int sw = 0; sw < 2; ++sw)
            {
                rc = sync_composite_cell(ctx, phi_l0, phi_l1, true);
                if (rc != C2_OK) return rc;
                for (int col = 1; col >= 0; --col)
                {
                    p_poisson_rbgs_kernel<<<dim3(64, 1), 256>>>(
                        phi_l0, b_l0, d_l0[0], d_l0[1], d_l0[2], ctx->nx0, ctx->ny0, ctx->nz0,
                        RealT(ctx->dx0), RealT(ctx->dx0), RealT(ctx->dx0), col, false, nullptr, 1, ctx->d_l1_map);
                    p_poisson_rbgs_kernel<<<gr_l1, 256>>>(
                        phi_l1, b_l1, d_l1[0], d_l1[1], d_l1[2], ctx->nx0, ctx->ny0, ctx->nz0,
                        RealT(ctx->dx1), RealT(ctx->dx1), RealT(ctx->dx1), col, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
                    rc = sync_composite_cell(ctx, phi_l0, phi_l1, false);
                    if (rc != C2_OK) return rc;
                }
            }
        }
        rc = sync_composite_cell(ctx, phi_l0, phi_l1, ctx->n_l1 > 0);
        if (rc != C2_OK) return rc;
    }
    return C2_OK;
}

template <typename RealT>
static int run_vel_fac_mg(
    C2Context* ctx, RealT* u_l0[3], RealT* u_l1[3], const RealT* const b_l0[3], const RealT* const b_l1[3], RealT* r_l0[3], RealT* r_l1[3],
    const RealT* const c_l0[3], const RealT* const c_l1[3], RealT* cl_x[kMaxCoarseLevels][3], RealT* cl_b[kMaxCoarseLevels][3],
    RealT* cl_r[kMaxCoarseLevels][3], RealT* cl_c[kMaxCoarseLevels][3])
{
    dim3 gr_l1(16, std::max(1, ctx->n_l1));
    if (ctx->n_l1 > 0)
    {
        for (int sw = 0; sw < 2; ++sw)
        {
            int rc = sync_composite_side(ctx, u_l0, u_l1, true);
            if (rc != C2_OK) return rc;
            for (int col = 0; col < 2; ++col)
            {
                for (int ax = 0; ax < 3; ++ax)
                    vel_helmholtz_rbgs_axis_kernel<<<gr_l1, 256>>>(
                        u_l1[0], u_l1[1], u_l1[2], b_l1[0], b_l1[1], b_l1[2], c_l1[0], c_l1[1], c_l1[2],
                        ax, ctx->nx0, ctx->ny0, ctx->nz0, RealT(ctx->dx1), RealT(ctx->dx1), RealT(ctx->dx1),
                        col, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
                rc = sync_composite_side(ctx, u_l0, u_l1, false);
                if (rc != C2_OK) return rc;
            }
        }
        int rc = sync_composite_side(ctx, u_l0, u_l1, true);
        if (rc != C2_OK) return rc;
        vel_helmholtz_residual_kernel<<<gr_l1, 256>>>(
            u_l1[0], u_l1[1], u_l1[2], b_l1[0], b_l1[1], b_l1[2], c_l1[0], c_l1[1], c_l1[2], r_l1[0], r_l1[1], r_l1[2],
            ctx->nx0, ctx->ny0, ctx->nz0, RealT(ctx->dx1), RealT(ctx->dx1), RealT(ctx->dx1),
            true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        vel_helmholtz_residual_kernel<<<dim3(64, 1), 256>>>(
            u_l0[0], u_l0[1], u_l0[2], b_l0[0], b_l0[1], b_l0[2], c_l0[0], c_l0[1], c_l0[2], r_l0[0], r_l0[1], r_l0[2],
            ctx->nx0, ctx->ny0, ctx->nz0, RealT(ctx->dx0), RealT(ctx->dx0), RealT(ctx->dx0),
            false, nullptr, 1, ctx->d_l1_map);
        for (int ax = 0; ax < 3; ++ax)
            restrict_l1_to_l0_side_kernel<<<gr_l1, 256>>>(r_l1[ax], r_l0[ax], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, ax, 1);
    }
    else
    {
        vel_helmholtz_residual_kernel<<<dim3(64, 1), 256>>>(
            u_l0[0], u_l0[1], u_l0[2], b_l0[0], b_l0[1], b_l0[2], c_l0[0], c_l0[1], c_l0[2], r_l0[0], r_l0[1], r_l0[2],
            ctx->nx0, ctx->ny0, ctx->nz0, RealT(ctx->dx0), RealT(ctx->dx0), RealT(ctx->dx0),
            false, nullptr, 1, ctx->d_l1_map);
    }

    for (int lv = 0; lv < ctx->num_coarse_levels; ++lv)
        for (int ax = 0; ax < 3; ++ax)
            C2_CUDA_CALL(ctx, cudaMemset(cl_x[lv][ax], 0, ctx->cl_meta[lv].n_s_g1[ax] * sizeof(RealT)), "v_fac_memset");
    for (int ax = 0; ax < 3; ++ax)
        cl_b[0][ax] = r_l0[ax];
    const int ncl = ctx->num_coarse_levels;
    for (int lv = 0; lv < ncl - 1; ++lv)
    {
        const CoarseLevelMeta& m = ctx->cl_meta[lv];
        const RealT* lv_c[3] = {lv == 0 ? c_l0[0] : cl_c[lv][0], lv == 0 ? c_l0[1] : cl_c[lv][1], lv == 0 ? c_l0[2] : cl_c[lv][2]};
        int ns = (lv >= 3) ? 4 : 2;
        for (int sw = 0; sw < ns; ++sw)
            for (int col = 0; col < 2; ++col)
            {
                enforce_homog_phys_bcs_side_kernel<<<dim3(32, 1), 256>>>(cl_x[lv][0], cl_x[lv][1], cl_x[lv][2], m.nx, m.ny, m.nz, false, nullptr, 1);
                for (int ax = 0; ax < 3; ++ax)
                    vel_helmholtz_rbgs_axis_kernel<<<dim3(32, 1), 256>>>(
                        cl_x[lv][0], cl_x[lv][1], cl_x[lv][2], cl_b[lv][0], cl_b[lv][1], cl_b[lv][2], lv_c[0], lv_c[1], lv_c[2],
                        ax, m.nx, m.ny, m.nz, RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]), col, false, nullptr, 1, nullptr);
            }
        enforce_homog_phys_bcs_side_kernel<<<dim3(32, 1), 256>>>(cl_x[lv][0], cl_x[lv][1], cl_x[lv][2], m.nx, m.ny, m.nz, false, nullptr, 1);
        vel_helmholtz_residual_kernel<<<dim3(32, 1), 256>>>(
            cl_x[lv][0], cl_x[lv][1], cl_x[lv][2], cl_b[lv][0], cl_b[lv][1], cl_b[lv][2], lv_c[0], lv_c[1], lv_c[2],
            cl_r[lv][0], cl_r[lv][1], cl_r[lv][2], m.nx, m.ny, m.nz, RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]), false, nullptr, 1, nullptr);

        const CoarseLevelMeta& mc = ctx->cl_meta[lv + 1];
        restrict_side_kernel<<<16, 256>>>(
            cl_r[lv][0], cl_r[lv][1], cl_r[lv][2],
            cl_b[lv + 1][0], cl_b[lv + 1][1], cl_b[lv + 1][2],
            cl_x[lv + 1][0], cl_x[lv + 1][1], cl_x[lv + 1][2],
            mc.nx, mc.ny, mc.nz, m.nx, m.ny, m.ratio[0], m.ratio[1], m.ratio[2]);
    }

    {
        const CoarseLevelMeta& mc = ctx->cl_meta[ncl - 1];
        for (int sw = 0; sw < 24; ++sw)
            for (int col = 0; col < 2; ++col)
            {
                enforce_homog_phys_bcs_side_kernel<<<1, 256>>>(
                    cl_x[ncl - 1][0], cl_x[ncl - 1][1], cl_x[ncl - 1][2],
                    mc.nx, mc.ny, mc.nz, false, nullptr, 1);
                for (int ax = 0; ax < 3; ++ax)
                    vel_helmholtz_rbgs_axis_kernel<<<1, 256>>>(
                        cl_x[ncl - 1][0], cl_x[ncl - 1][1], cl_x[ncl - 1][2],
                        cl_b[ncl - 1][0], cl_b[ncl - 1][1], cl_b[ncl - 1][2],
                        cl_c[ncl - 1][0], cl_c[ncl - 1][1], cl_c[ncl - 1][2],
                        ax, mc.nx, mc.ny, mc.nz, RealT(mc.dx[0]), RealT(mc.dx[1]), RealT(mc.dx[2]), col, false, nullptr, 1, nullptr);
            }
    }

    for (int lv = ncl - 2; lv >= 0; --lv)
    {
        const CoarseLevelMeta& m = ctx->cl_meta[lv];
        const CoarseLevelMeta& mc = ctx->cl_meta[lv + 1];
        prolong_side_add_kernel<<<dim3(32, 1), 256>>>(
            cl_x[lv + 1][0], cl_x[lv + 1][1], cl_x[lv + 1][2],
            cl_x[lv][0], cl_x[lv][1], cl_x[lv][2],
            m.nx, m.ny, m.nz, mc.nx, mc.ny, mc.nz, m.ratio[0], m.ratio[1], m.ratio[2], false, nullptr, 1, nullptr);
        const RealT* lv_c[3] = {lv == 0 ? c_l0[0] : cl_c[lv][0], lv == 0 ? c_l0[1] : cl_c[lv][1], lv == 0 ? c_l0[2] : cl_c[lv][2]};
        int ns = (lv >= 3) ? 4 : 2;
        for (int sw = 0; sw < ns; ++sw)
            for (int col = 0; col < 2; ++col)
            {
                enforce_homog_phys_bcs_side_kernel<<<dim3(32, 1), 256>>>(cl_x[lv][0], cl_x[lv][1], cl_x[lv][2], m.nx, m.ny, m.nz, false, nullptr, 1);
                for (int ax = 0; ax < 3; ++ax)
                    vel_helmholtz_rbgs_axis_kernel<<<dim3(32, 1), 256>>>(
                        cl_x[lv][0], cl_x[lv][1], cl_x[lv][2], cl_b[lv][0], cl_b[lv][1], cl_b[lv][2], lv_c[0], lv_c[1], lv_c[2],
                        ax, m.nx, m.ny, m.nz, RealT(m.dx[0]), RealT(m.dx[1]), RealT(m.dx[2]), col, false, nullptr, 1, nullptr);
            }
    }

    prolong_side_add_kernel<<<dim3(32, 1), 256>>>(
        cl_x[0][0], cl_x[0][1], cl_x[0][2], u_l0[0], u_l0[1], u_l0[2],
        ctx->nx0, ctx->ny0, ctx->nz0, ctx->nx0, ctx->ny0, ctx->nz0, 1, 1, 1, false, nullptr, 1, nullptr);
    if (ctx->n_l1 > 0)
    {
        prolong_side_add_kernel<<<gr_l1, 256>>>(
            cl_x[0][0], cl_x[0][1], cl_x[0][2], u_l1[0], u_l1[1], u_l1[2],
            2 * ctx->nx0, 2 * ctx->ny0, 2 * ctx->nz0, ctx->nx0, ctx->ny0, ctx->nz0,
            2, 2, 2, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        for (int sw = 0; sw < 2; ++sw)
        {
            int rc = sync_composite_side(ctx, u_l0, u_l1, true);
            if (rc != C2_OK) return rc;
            for (int col = 0; col < 2; ++col)
            {
                for (int ax = 0; ax < 3; ++ax)
                    vel_helmholtz_rbgs_axis_kernel<<<gr_l1, 256>>>(
                        u_l1[0], u_l1[1], u_l1[2], b_l1[0], b_l1[1], b_l1[2], c_l1[0], c_l1[1], c_l1[2],
                        ax, ctx->nx0, ctx->ny0, ctx->nz0, RealT(ctx->dx1), RealT(ctx->dx1), RealT(ctx->dx1),
                        col, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
                rc = sync_composite_side(ctx, u_l0, u_l1, false);
                if (rc != C2_OK) return rc;
            }
        }
    }
    return sync_composite_side(ctx, u_l0, u_l1, ctx->n_l1 > 0);
}

static int apply_projection_preconditioner(
    C2Context* ctx, double* const rhs_g0[4], double* const z_g0[4], double dt, bool use_fp32)
{
    const dim3 l1g(16, std::max(1, ctx->n_l1));
    bool use_fp64_mg = !use_fp32;

    double bu_norm = 0.0;
    {
        const int nblk = 256;
        for (int a = 0; a < 3; ++a)
        {
            wgt_dot1_kernel<<<nblk, 256>>>(rhs_g0[a], ctx->d_wgt_s[a], ctx->tot_s_g0(a), ctx->d_red_scratch);
            C2_CUDA_LAUNCH_CHK(ctx, "bu_norm");
            C2_CUDA_CALL(ctx, cudaMemcpy(ctx->h_red_pinned, ctx->d_red_scratch, nblk * sizeof(double), cudaMemcpyDeviceToHost), "bu_norm_d2h");
            for (int i = 0; i < nblk; ++i) bu_norm += ctx->h_red_pinned[i];
        }
        if (!std::isfinite(bu_norm)) return C2_ERR_NONFINITE;
        bu_norm = std::sqrt(std::max(0.0, bu_norm));
    }
    if (bu_norm == 0.0) bu_norm = 1.0;

    if (use_fp64_mg)
    {
        double *x0[3], *x1[3], *b0[3], *b1[3], *r0[3], *r1[3];
        const double *c0[3], *c1[3];
        for (int a = 0; a < 3; ++a)
        {
            C2_CUDA_CALL(ctx, cudaMemset(ctx->d_mg_fp64[a], 0, ctx->tot_s_g1(a) * sizeof(double)), "pc_x_zero");
            C2_CUDA_CALL(ctx, cudaMemset(ctx->d_mg_fp64[3 + a], 0, ctx->tot_s_g1(a) * sizeof(double)), "pc_b_zero");
            x0[a] = ctx->d_mg_fp64[a]; x1[a] = ctx->d_mg_fp64[a] + ctx->l0_s_g1[a];
            b0[a] = ctx->d_mg_fp64[3 + a]; b1[a] = ctx->d_mg_fp64[3 + a] + ctx->l0_s_g1[a];
            r0[a] = ctx->d_mg_fp64[6 + a]; r1[a] = ctx->d_mg_fp64[6 + a] + ctx->l0_s_g1[a];
            c0[a] = ctx->d_state[9 + a]; c1[a] = ctx->d_state[9 + a] + ctx->l0_s_g1[a];

            pack_g0_to_g_side_l0_kernel<<<64, 256>>>(rhs_g0[a], b0[a], a, ctx->nx0, ctx->ny0, ctx->nz0, 1);
            if (ctx->n_l1 > 0)
                pack_g0_to_g_side_l1_kernel<<<l1g, 256>>>(rhs_g0[a] + ctx->l0_s_g0[a], b1[a], ctx->d_l1_boxes, ctx->n_l1, a, 1);
        }
        double* cl_x[kMaxCoarseLevels][3] = {};
        double* cl_b[kMaxCoarseLevels][3] = {};
        double* cl_r[kMaxCoarseLevels][3] = {};
        double* cl_c[kMaxCoarseLevels][3] = {};
        for (int lv = 0; lv < ctx->num_coarse_levels; ++lv)
            for (int a = 0; a < 3; ++a)
            {
                cl_x[lv][a] = ctx->d_cl_fp64_vec[0][a] + ctx->cl_meta[lv].off_slot_s_g1[a];
                cl_b[lv][a] = (lv == 0) ? r0[a] : (ctx->d_cl_fp64_vec[1][a] + ctx->cl_meta[lv].off_c1_s_g1[a]);
                cl_r[lv][a] = ctx->d_cl_fp64_vec[2][a] + ctx->cl_meta[lv].off_slot_s_g1[a];
                cl_c[lv][a] = (lv == 0) ? ctx->d_state[9 + a] : (ctx->d_cl_fp64_c[a] + ctx->cl_meta[lv].off_c1_s_g1[a]);
            }
        int rc = run_vel_fac_mg<double>(ctx, x0, x1, b0, b1, r0, r1, c0, c1, cl_x, cl_b, cl_r, cl_c);
        if (rc != C2_OK) return rc;
    }
    else
    {
        float* x_ptrs[3] = {
            reinterpret_cast<float*>(ctx->d_mg_fp64[0]),
            reinterpret_cast<float*>(ctx->d_mg_fp64[0]) + ctx->max_slot_g1,
            reinterpret_cast<float*>(ctx->d_mg_fp64[1])};
        float* b_ptrs[3] = {
            reinterpret_cast<float*>(ctx->d_mg_fp64[1]) + ctx->max_slot_g1,
            reinterpret_cast<float*>(ctx->d_mg_fp64[2]),
            reinterpret_cast<float*>(ctx->d_mg_fp64[2]) + ctx->max_slot_g1};
        float* r_ptrs[3] = {
            reinterpret_cast<float*>(ctx->d_mg_fp64[3]),
            reinterpret_cast<float*>(ctx->d_mg_fp64[3]) + ctx->max_slot_g1,
            reinterpret_cast<float*>(ctx->d_mg_fp64[4])};

        float* x0[3]; float* x1[3]; float* b0[3]; float* b1[3]; float* r0[3]; float* r1[3];
        const float* c0[3]; const float* c1[3];
        for (int a = 0; a < 3; ++a)
        {
            float* xp = x_ptrs[a]; float* bp = b_ptrs[a]; float* rp = r_ptrs[a];
            C2_CUDA_CALL(ctx, cudaMemset(xp, 0, ctx->tot_s_g1(a) * sizeof(float)), "pc_x32_zero");
            x0[a] = xp; x1[a] = xp + ctx->l0_s_g1[a];
            b0[a] = bp; b1[a] = bp + ctx->l0_s_g1[a];
            r0[a] = rp; r1[a] = rp + ctx->l0_s_g1[a];
            c0[a] = ctx->d_mg_c_fp32[a];
            c1[a] = ctx->d_mg_c_fp32[a] + ctx->l0_s_g1[a];

            double* tmp_d = ctx->d_mg_fp64[6 + a];
            C2_CUDA_CALL(ctx, cudaMemset(tmp_d, 0, ctx->tot_s_g1(a) * sizeof(double)), "pc_tmp_zero");
            pack_g0_to_g_side_l0_kernel<<<64, 256>>>(rhs_g0[a], tmp_d, a, ctx->nx0, ctx->ny0, ctx->nz0, 1);
            if (ctx->n_l1 > 0)
                pack_g0_to_g_side_l1_kernel<<<l1g, 256>>>(rhs_g0[a] + ctx->l0_s_g0[a], tmp_d + ctx->l0_s_g1[a], ctx->d_l1_boxes, ctx->n_l1, a, 1);
            scale_copy_f64_to_real_kernel<<<128, 256>>>(tmp_d, bp, 1.0 / bu_norm, ctx->tot_s_g1(a));
        }
        float* cl_x[kMaxCoarseLevels][3] = {};
        float* cl_b[kMaxCoarseLevels][3] = {};
        float* cl_r[kMaxCoarseLevels][3] = {};
        float* cl_c[kMaxCoarseLevels][3] = {};
        for (int lv = 0; lv < ctx->num_coarse_levels; ++lv)
            for (int a = 0; a < 3; ++a)
            {
                cl_x[lv][a] = reinterpret_cast<float*>(ctx->d_cl_fp64_vec[0][a]) + ctx->cl_meta[lv].off_slot_s_g1[a];
                cl_b[lv][a] = (lv == 0) ? r0[a] : (reinterpret_cast<float*>(ctx->d_cl_fp64_vec[1][a]) + ctx->cl_meta[lv].off_c1_s_g1[a]);
                cl_r[lv][a] = reinterpret_cast<float*>(ctx->d_cl_fp64_vec[2][a]) + ctx->cl_meta[lv].off_slot_s_g1[a];
                cl_c[lv][a] = (lv == 0) ? ctx->d_mg_c_fp32[a] : (ctx->d_cl_fp32_c[a] + ctx->cl_meta[lv].off_c1_s_g1[a]);
            }
        int rc = run_vel_fac_mg<float>(ctx, x0, x1, b0, b1, r0, r1, c0, c1, cl_x, cl_b, cl_r, cl_c);
        if (rc != C2_OK) return rc;
        for (int a = 0; a < 3; ++a)
            scale_copy_real_to_f64_kernel<<<128, 256>>>(x0[a], ctx->d_mg_fp64[6 + a], bu_norm, ctx->tot_s_g1(a));
        for (int a = 0; a < 3; ++a)
            C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_mg_fp64[a], ctx->d_mg_fp64[6 + a], ctx->tot_s_g1(a) * sizeof(double), cudaMemcpyDeviceToDevice), "v32_to_64");
    }

    double* ustar_l0[3] = {ctx->d_mg_fp64[0], ctx->d_mg_fp64[1], ctx->d_mg_fp64[2]};
    double* ustar_l1[3] = {ctx->d_mg_fp64[0] + ctx->l0_s_g1[0], ctx->d_mg_fp64[1] + ctx->l0_s_g1[1], ctx->d_mg_fp64[2] + ctx->l0_s_g1[2]};
    int rc = sync_composite_side(ctx, ustar_l0, ustar_l1, ctx->n_l1 > 0);
    if (rc != C2_OK) return rc;

    double* div_ustar_g0 = ctx->d_mg_fp64[6];
    compute_div_kernel<<<dim3(64, 1), 256>>>(
        ustar_l0[0], ustar_l0[1], ustar_l0[2], div_ustar_g0, 1.0,
        ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, false, nullptr, 1, ctx->d_l1_map);
    if (ctx->n_l1 > 0)
        compute_div_kernel<<<l1g, 256>>>(
            ustar_l1[0], ustar_l1[1], ustar_l1[2], div_ustar_g0 + ctx->l0_c_g0, 1.0,
            ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);

    double* bphi_l0 = ctx->d_mg_fp64[5];
    double* bphi_l1 = ctx->d_mg_fp64[5] + ctx->l0_c_g1;
    C2_CUDA_CALL(ctx, cudaMemset(bphi_l0, 0, ctx->tot_c_g1() * sizeof(double)), "bphi_zero");
    projection_build_fphi_kernel<<<dim3(64, 1), 256>>>(
        rhs_g0[3], div_ustar_g0, bphi_l0, ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1, ctx->d_l1_map);
    if (ctx->n_l1 > 0)
        projection_build_fphi_kernel<<<l1g, 256>>>(
            rhs_g0[3] + ctx->l0_c_g0, div_ustar_g0 + ctx->l0_c_g0, bphi_l1,
            ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);

    linear_sum_side_kernel<<<64, 256>>>(div_ustar_g0, 1.0, div_ustar_g0, 1.0, rhs_g0[3], ctx->tot_c_g0());
    double bp_norm = 0.0;
    {
        const int nblk = 256;
        wgt_dot1_kernel<<<nblk, 256>>>(div_ustar_g0, ctx->d_wgt_c, ctx->tot_c_g0(), ctx->d_red_scratch);
        C2_CUDA_LAUNCH_CHK(ctx, "bp_norm");
        C2_CUDA_CALL(ctx, cudaMemcpy(ctx->h_red_pinned, ctx->d_red_scratch, nblk * sizeof(double), cudaMemcpyDeviceToHost), "bp_norm_d2h");
        for (int i = 0; i < nblk; ++i) bp_norm += ctx->h_red_pinned[i];
        if (!std::isfinite(bp_norm)) return C2_ERR_NONFINITE;
        bp_norm = std::sqrt(std::max(0.0, bp_norm));
    }
    if (bp_norm == 0.0) bp_norm = 1.0;

    if (use_fp64_mg)
    {
        double* phi_l0 = ctx->d_mg_fp64[9];
        double* phi_l1 = ctx->d_mg_fp64[9] + ctx->l0_c_g1;
        C2_CUDA_CALL(ctx, cudaMemset(phi_l0, 0, ctx->tot_c_g1() * sizeof(double)), "phi64_zero");
        double* r_l0 = ctx->d_mg_fp64[4];
        double* r_l1 = ctx->d_mg_fp64[4] + ctx->l0_c_g1;
        const double* d_l0[3] = {ctx->d_state[12], ctx->d_state[13], ctx->d_state[14]};
        const double* d_l1[3] = {ctx->d_state[12] + ctx->l0_s_g1[0], ctx->d_state[13] + ctx->l0_s_g1[1], ctx->d_state[14] + ctx->l0_s_g1[2]};

        double* cl_x[kMaxCoarseLevels] = {};
        double* cl_b[kMaxCoarseLevels] = {};
        double* cl_r[kMaxCoarseLevels] = {};
        double* cl_d[kMaxCoarseLevels][3] = {};
        for (int lv = 0; lv < ctx->num_coarse_levels; ++lv)
        {
            cl_x[lv] = ctx->d_cl_fp64_vec[0][0] + ctx->cl_meta[lv].off_slot_c_g1;
            cl_b[lv] = (lv == 0) ? r_l0 : (ctx->d_cl_fp64_vec[1][0] + ctx->cl_meta[lv].off_c1_c_g1);
            cl_r[lv] = ctx->d_cl_fp64_vec[2][0] + ctx->cl_meta[lv].off_slot_c_g1;
            for (int a = 0; a < 3; ++a)
                cl_d[lv][a] = (lv == 0) ? ctx->d_state[12 + a] : (ctx->d_cl_fp64_d[a] + ctx->cl_meta[lv].off_c1_s_g1[a]);
        }
        rc = run_pressure_fac_mg<double>(ctx, phi_l0, phi_l1, bphi_l0, bphi_l1, r_l0, r_l1, d_l0, d_l1, cl_x, cl_b, cl_r, cl_d, 2);
        if (rc != C2_OK) return rc;
    }
    else
    {
        float* phi_f = reinterpret_cast<float*>(ctx->d_mg_fp64[3]);
        float* bp_f = reinterpret_cast<float*>(ctx->d_mg_fp64[3]) + ctx->max_slot_g1;
        float* r_f = reinterpret_cast<float*>(ctx->d_mg_fp64[4]);
        C2_CUDA_CALL(ctx, cudaMemset(phi_f, 0, ctx->tot_c_g1() * sizeof(float)), "phi32_zero");
        scale_copy_f64_to_real_kernel<<<128, 256>>>(bphi_l0, bp_f, 1.0 / bp_norm, ctx->tot_c_g1());

        float* phi_l0 = phi_f; float* phi_l1 = phi_f + ctx->l0_c_g1;
        float* b_l0_f = bp_f;  float* b_l1_f = bp_f + ctx->l0_c_g1;
        float* r_l0_f = r_f;   float* r_l1_f = r_f + ctx->l0_c_g1;
        const float* d_l0[3] = {ctx->d_mg_d_fp32[0], ctx->d_mg_d_fp32[1], ctx->d_mg_d_fp32[2]};
        const float* d_l1[3] = {ctx->d_mg_d_fp32[0] + ctx->l0_s_g1[0], ctx->d_mg_d_fp32[1] + ctx->l0_s_g1[1], ctx->d_mg_d_fp32[2] + ctx->l0_s_g1[2]};

        float* cl_x[kMaxCoarseLevels] = {};
        float* cl_b[kMaxCoarseLevels] = {};
        float* cl_r[kMaxCoarseLevels] = {};
        float* cl_d[kMaxCoarseLevels][3] = {};
        for (int lv = 0; lv < ctx->num_coarse_levels; ++lv)
        {
            cl_x[lv] = reinterpret_cast<float*>(ctx->d_cl_fp64_vec[0][0]) + ctx->cl_meta[lv].off_slot_c_g1;
            cl_b[lv] = (lv == 0) ? r_l0_f : (reinterpret_cast<float*>(ctx->d_cl_fp64_vec[1][0]) + ctx->cl_meta[lv].off_c1_c_g1);
            cl_r[lv] = reinterpret_cast<float*>(ctx->d_cl_fp64_vec[2][0]) + ctx->cl_meta[lv].off_slot_c_g1;
            for (int a = 0; a < 3; ++a)
                cl_d[lv][a] = (lv == 0) ? ctx->d_mg_d_fp32[a] : (ctx->d_cl_fp32_d[a] + ctx->cl_meta[lv].off_c1_s_g1[a]);
        }
        rc = run_pressure_fac_mg<float>(ctx, phi_l0, phi_l1, b_l0_f, b_l1_f, r_l0_f, r_l1_f, d_l0, d_l1, cl_x, cl_b, cl_r, cl_d, 2);
        if (rc != C2_OK) return rc;
        scale_copy_real_to_f64_kernel<<<128, 256>>>(phi_f, ctx->d_mg_fp64[9], bp_norm, ctx->tot_c_g1());
    }

    double* phi_l0 = ctx->d_mg_fp64[9];
    double* phi_l1 = ctx->d_mg_fp64[9] + ctx->l0_c_g1;
    rc = sync_composite_cell(ctx, phi_l0, phi_l1, ctx->n_l1 > 0);
    if (rc != C2_OK) return rc;
    const double* dp_l0[3] = {ctx->d_state[12], ctx->d_state[13], ctx->d_state[14]};
    const double* dp_l1[3] = {ctx->d_state[12] + ctx->l0_s_g1[0], ctx->d_state[13] + ctx->l0_s_g1[1], ctx->d_state[14] + ctx->l0_s_g1[2]};

    projection_update_u_p_kernel<<<dim3(64, 1), 256>>>(
        ustar_l0[0], ustar_l0[1], ustar_l0[2], phi_l0, bphi_l0, dp_l0[0], dp_l0[1], dp_l0[2],
        z_g0[0], z_g0[1], z_g0[2], z_g0[3], dt, ctx->dx0, ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1, ctx->d_l1_map);
    if (ctx->n_l1 > 0)
    {
        projection_update_u_p_kernel<<<l1g, 256>>>(
            ustar_l1[0], ustar_l1[1], ustar_l1[2], phi_l1, bphi_l1, dp_l1[0], dp_l1[1], dp_l1[2],
            z_g0[0] + ctx->l0_s_g0[0], z_g0[1] + ctx->l0_s_g0[1], z_g0[2] + ctx->l0_s_g0[2], z_g0[3] + ctx->l0_c_g0,
            dt, ctx->dx1, ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
        for (int a = 0; a < 3; ++a)
        {
            sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(z_g0[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            restrict_g0_l1_to_l0_side_kernel<<<l1g, 256>>>(z_g0[a] + ctx->l0_s_g0[a], z_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }
    C2_CUDA_LAUNCH_CHK(ctx, "proj_update");
    return C2_OK;
}

static int refresh_mg_coefs_from_state(C2Context* ctx)
{
    double* c_l0[3] = {ctx->d_state[9], ctx->d_state[10], ctx->d_state[11]};
    double* d_l0[3] = {ctx->d_state[12], ctx->d_state[13], ctx->d_state[14]};
    size_t max_s1 = std::max({ctx->tot_s_g1(0), ctx->tot_s_g1(1), ctx->tot_s_g1(2)});
    cast_fine_coefs_fp32_kernel<<<128, 256>>>(
        c_l0[0], c_l0[1], c_l0[2], d_l0[0], d_l0[1], d_l0[2],
        ctx->d_mg_c_fp32[0], ctx->d_mg_c_fp32[1], ctx->d_mg_c_fp32[2],
        ctx->d_mg_d_fp32[0], ctx->d_mg_d_fp32[1], ctx->d_mg_d_fp32[2],
        max_s1);
    for (int lv = 1; lv < ctx->num_coarse_levels; ++lv)
    {
        const CoarseLevelMeta& fm = ctx->cl_meta[lv - 1];
        const CoarseLevelMeta& cm = ctx->cl_meta[lv];
        for (int a = 0; a < 3; ++a)
        {
            const double* fc = (lv == 1) ? c_l0[a] : (ctx->d_cl_fp64_c[a] + fm.off_c1_s_g1[a]);
            const double* fd = (lv == 1) ? d_l0[a] : (ctx->d_cl_fp64_d[a] + fm.off_c1_s_g1[a]);
            double* cc = ctx->d_cl_fp64_c[a] + cm.off_c1_s_g1[a];
            double* cd = ctx->d_cl_fp64_d[a] + cm.off_c1_s_g1[a];
            float* ccf = ctx->d_cl_fp32_c[a] + cm.off_c1_s_g1[a];
            float* cdf = ctx->d_cl_fp32_d[a] + cm.off_c1_s_g1[a];
            coarsen_coef_level_kernel<<<32, 256>>>(
                fc, fd, cc, cd, ccf, cdf, a, cm.nx, cm.ny, cm.nz, fm.nx, fm.ny,
                fm.ratio[0], fm.ratio[1], fm.ratio[2], ctx->d_vel_const);
        }
    }
    C2_CUDA_LAUNCH_CHK(ctx, "refresh_mg_coefs");
    return C2_OK;
}

static int run_gpu_ssprk3_step(
    C2Context* ctx, double dt, double dt_prev, bool is_first_step,
    double* const u_old_g0[3], double* const u_cur_g0[3], double* const rho_cur_g0[3],
    double* const rho_new_g0[3], double* const n_full_g0[3], double* out_total_ms)
{
    double* v_comp[3] = {ctx->d_krylov_v[0][0], ctx->d_krylov_v[0][1], ctx->d_krylov_v[0][2]};
    double* v_scr_l0[3] = {ctx->d_krylov_v[1][0], ctx->d_krylov_v[1][1], ctx->d_krylov_v[1][2]};
    double* v_scr_l1[3] = {ctx->d_krylov_v[2][0], ctx->d_krylov_v[2][1], ctx->d_krylov_v[2][2]};
    double* r_pre_l0[3] = {ctx->d_krylov_v[3][0], ctx->d_krylov_v[3][1], ctx->d_krylov_v[3][2]};
    double* r_pre_l1[3] = {ctx->d_krylov_v[4][0], ctx->d_krylov_v[4][1], ctx->d_krylov_v[4][2]};
    double* r_old_work[3] = {ctx->d_krylov_z[2][0], ctx->d_krylov_z[2][1], ctx->d_krylov_z[2][2]};

    const dim3 l1g(16, std::max(1, ctx->n_l1));
    cudaEvent_t et0{}, et1{};
    C2_CUDA_CALL(ctx, cudaEventCreate(&et0), "ev_ss0");
    C2_CUDA_CALL(ctx, cudaEventCreate(&et1), "ev_ss1");
    C2_CUDA_CALL(ctx, cudaEventRecord(et0), "ev_ss_rec0");

    // Step 0
    for (int a = 0; a < 3; ++a)
    {
        C2_CUDA_CALL(ctx, cudaMemcpy(r_old_work[a], rho_cur_g0[a], ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyDeviceToDevice), "ssprk_rold");
        if (ctx->n_l1 > 0)
        {
            sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(r_old_work[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            restrict_g0_l1_to_l0_side_kernel<<<l1g, 256>>>(r_old_work[a] + ctx->l0_s_g0[a], r_old_work[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }
    int rc = fill_ghost_side_ssprk_g3(ctx, const_cast<double**>(rho_cur_g0), r_pre_l0, r_pre_l1, true);
    if (rc != C2_OK) { cudaEventDestroy(et0); cudaEventDestroy(et1); return rc; }
    rc = fill_ghost_side_ssprk_g3(ctx, const_cast<double**>(u_cur_g0), v_scr_l0, v_scr_l1, false);
    if (rc != C2_OK) { cudaEventDestroy(et0); cudaEventDestroy(et1); return rc; }

    for (int a = 0; a < 3; ++a)
    {
        ssprk3_stage_update_kernel<<<dim3(64, 1), 256>>>(
            a, dt, 0.5, 0.5, 1.0, false,
            r_old_work[a], r_pre_l0[a], v_scr_l0[0], v_scr_l0[1], v_scr_l0[2], u_cur_g0[a],
            rho_new_g0[a], nullptr, ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, false, nullptr, 1, 0, 0, 0, 0, 0);
        if (ctx->n_l1 > 0)
            ssprk3_stage_update_kernel<<<l1g, 256>>>(
                a, dt, 0.5, 0.5, 1.0, false,
                r_old_work[a] + ctx->l0_s_g0[a], r_pre_l1[a], v_scr_l1[0], v_scr_l1[1], v_scr_l1[2], u_cur_g0[a] + ctx->l0_s_g0[a],
                rho_new_g0[a] + ctx->l0_s_g0[a], nullptr, ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, true, ctx->d_l1_boxes, ctx->n_l1, 0, 0, 0, 0, 0);
    }

    const double omega = (!is_first_step && dt_prev > 0.0) ? (dt / dt_prev) : 1.0;

    // Step 1
    rc = fill_ghost_side_ssprk_g3(ctx, const_cast<double**>(rho_new_g0), r_pre_l0, r_pre_l1, true);
    if (rc != C2_OK) { cudaEventDestroy(et0); cudaEventDestroy(et1); return rc; }
    {
        double w0 = is_first_step ? 0.0 : -omega;
        double w1 = is_first_step ? 1.0 : (1.0 + omega);
        for (int a = 0; a < 3; ++a)
            linear_sum_side_kernel<<<128, 256>>>(v_comp[a], w0, u_old_g0[a], w1, u_cur_g0[a], ctx->tot_s_g0(a));
    }
    rc = fill_ghost_side_ssprk_g3(ctx, v_comp, v_scr_l0, v_scr_l1, false);
    if (rc != C2_OK) { cudaEventDestroy(et0); cudaEventDestroy(et1); return rc; }
    for (int a = 0; a < 3; ++a)
    {
        ssprk3_stage_update_kernel<<<dim3(64, 1), 256>>>(
            a, dt, 0.75, 0.25, 0.25, false,
            r_old_work[a], r_pre_l0[a], v_scr_l0[0], v_scr_l0[1], v_scr_l0[2], u_cur_g0[a],
            rho_new_g0[a], nullptr, ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, false, nullptr, 1, 0, 0, 0, 0, 0);
        if (ctx->n_l1 > 0)
            ssprk3_stage_update_kernel<<<l1g, 256>>>(
                a, dt, 0.75, 0.25, 0.25, false,
                r_old_work[a] + ctx->l0_s_g0[a], r_pre_l1[a], v_scr_l1[0], v_scr_l1[1], v_scr_l1[2], u_cur_g0[a] + ctx->l0_s_g0[a],
                rho_new_g0[a] + ctx->l0_s_g0[a], nullptr, ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, true, ctx->d_l1_boxes, ctx->n_l1, 0, 0, 0, 0, 0);
    }

    // Step 2
    rc = fill_ghost_side_ssprk_g3(ctx, const_cast<double**>(rho_new_g0), r_pre_l0, r_pre_l1, true);
    if (rc != C2_OK) { cudaEventDestroy(et0); cudaEventDestroy(et1); return rc; }
    {
        double w0 = is_first_step ? 0.0 : (-0.5 * omega);
        double w1 = is_first_step ? 1.0 : (1.0 + 0.5 * omega);
        for (int a = 0; a < 3; ++a)
            linear_sum_side_kernel<<<128, 256>>>(v_comp[a], w0, u_old_g0[a], w1, u_cur_g0[a], ctx->tot_s_g0(a));
    }
    rc = fill_ghost_side_ssprk_g3(ctx, v_comp, v_scr_l0, v_scr_l1, false);
    if (rc != C2_OK) { cudaEventDestroy(et0); cudaEventDestroy(et1); return rc; }
    for (int a = 0; a < 3; ++a)
    {
        ssprk3_stage_update_kernel<<<dim3(64, 1), 256>>>(
            a, dt, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, true,
            r_old_work[a], r_pre_l0[a], v_scr_l0[0], v_scr_l0[1], v_scr_l0[2], u_cur_g0[a],
            rho_new_g0[a], n_full_g0[a], ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, false, nullptr, 1, 0, 0, 0, 0, 0);
        if (ctx->n_l1 > 0)
            ssprk3_stage_update_kernel<<<l1g, 256>>>(
                a, dt, 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0, true,
                r_old_work[a] + ctx->l0_s_g0[a], r_pre_l1[a], v_scr_l1[0], v_scr_l1[1], v_scr_l1[2], u_cur_g0[a] + ctx->l0_s_g0[a],
                rho_new_g0[a] + ctx->l0_s_g0[a], n_full_g0[a] + ctx->l0_s_g0[a],
                ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, true, ctx->d_l1_boxes, ctx->n_l1, 0, 0, 0, 0, 0);
    }

    if (ctx->n_l1 > 0)
    {
        for (int a = 0; a < 3; ++a)
        {
            sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(rho_new_g0[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            restrict_g0_l1_to_l0_side_kernel<<<l1g, 256>>>(rho_new_g0[a] + ctx->l0_s_g0[a], rho_new_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(n_full_g0[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }

    C2_CUDA_CALL(ctx, cudaEventRecord(et1), "ev_ss_rec1");
    C2_CUDA_CALL(ctx, cudaEventSynchronize(et1), "ev_ss_sync");
    float tot = 0.0f;
    C2_CUDA_CALL(ctx, cudaEventElapsedTime(&tot, et0, et1), "ev_ss_el");
    if (out_total_ms) *out_total_ms = static_cast<double>(tot);
    cudaEventDestroy(et0); cudaEventDestroy(et1);
    return C2_OK;
}

struct SolveReport
{
    bool converged = false;
    int iters_total = 0;
    int fp32_iters = 0;
    int fp64_fb_iters = 0;
    bool fallback_triggered = false;
    std::string fallback_reason = "none";
    double final_true_rel_res = 0.0;
    double final_true_abs_res = 0.0;
    double rhs_norm = 0.0;
    double wall_ms = 0.0;
};

static int solve_fgmres_m5(
    C2Context* ctx, double* const rhs_g0[4], double* const x_g0[4],
    double dt, double rel_tol, double abs_tol, int max_iterations,
    bool zero_init_guess, bool start_in_fp32, bool allow_dynamic_fp64_fallback,
    SolveReport& rep)
{
    rep = SolveReport{};
    const int M = 5;
    if (zero_init_guess)
    {
        for (int c = 0; c < 4; ++c)
        {
            size_t n = (c < 3) ? ctx->tot_s_g0(c) : ctx->tot_c_g0();
            C2_CUDA_CALL(ctx, cudaMemset(x_g0[c], 0, n * sizeof(double)), "solve_zero");
        }
    }

    double b_norm = 0.0;
    int rc = gpu_wgt_norm4(ctx, rhs_g0, &b_norm);
    if (rc != C2_OK) return rc;
    rep.rhs_norm = b_norm;
    if (b_norm == 0.0 || b_norm <= abs_tol)
    {
        rep.converged = true;
        rep.final_true_rel_res = 0.0;
        rep.final_true_abs_res = b_norm;
        return C2_OK;
    }

    cudaEvent_t ev0{}, ev1{};
    C2_CUDA_CALL(ctx, cudaEventCreate(&ev0), "ev_st0");
    C2_CUDA_CALL(ctx, cudaEventCreate(&ev1), "ev_st1");
    C2_CUDA_CALL(ctx, cudaEventRecord(ev0), "ev_st_rec0");

    bool using_fp64 = !start_in_fp32;
    int low_prog = 0;
    int fp32_cycles = 0;

    while (rep.iters_total < max_iterations)
    {
        double* r0[4] = {ctx->d_krylov_v[0][0], ctx->d_krylov_v[0][1], ctx->d_krylov_v[0][2], ctx->d_krylov_v[0][3]};
        rc = apply_stokes_fp64(ctx, x_g0, r0);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
        for (int c = 0; c < 4; ++c)
        {
            size_t n = (c < 3) ? ctx->tot_s_g0(c) : ctx->tot_c_g0();
            linear_sum_side_kernel<<<128, 256>>>(r0[c], 1.0, rhs_g0[c], -1.0, r0[c], n);
        }
        double beta = 0.0;
        rc = gpu_wgt_norm4(ctx, r0, &beta);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
        rep.final_true_abs_res = beta;
        rep.final_true_rel_res = beta / b_norm;
        if (rep.final_true_rel_res <= rel_tol || beta <= abs_tol) { rep.converged = true; break; }

        double inv_beta = 1.0 / beta;
        for (int c = 0; c < 4; ++c)
        {
            size_t n = (c < 3) ? ctx->tot_s_g0(c) : ctx->tot_c_g0();
            linear_sum_side_kernel<<<128, 256>>>(r0[c], inv_beta, r0[c], 0.0, r0[c], n);
        }

        double H[M + 1][M] = {};
        double g[M + 1] = {}; g[0] = beta;
        double cs[M] = {}, sn[M] = {};
        double cycle_start_rel = rep.final_true_rel_res;
        int inner_done = 0;
        bool abort_to_fp64 = false;

        int max_j = std::min(M, max_iterations - rep.iters_total);
        for (int j = 0; j < max_j; ++j)
        {
            double* vj[4] = {ctx->d_krylov_v[j][0], ctx->d_krylov_v[j][1], ctx->d_krylov_v[j][2], ctx->d_krylov_v[j][3]};
            double* zj[4] = {ctx->d_krylov_z[j][0], ctx->d_krylov_z[j][1], ctx->d_krylov_z[j][2], ctx->d_krylov_z[j][3]};
            double* vj1[4] = {ctx->d_krylov_v[j + 1][0], ctx->d_krylov_v[j + 1][1], ctx->d_krylov_v[j + 1][2], ctx->d_krylov_v[j + 1][3]};

            rc = apply_projection_preconditioner(ctx, vj, zj, dt, !using_fp64);
            if (rc == C2_ERR_NONFINITE && !using_fp64 && allow_dynamic_fp64_fallback)
            {
                abort_to_fp64 = true;
                break;
            }
            if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }

            rc = apply_stokes_fp64(ctx, zj, vj1);
            if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }

            rep.iters_total++;
            if (!using_fp64) rep.fp32_iters++; else rep.fp64_fb_iters++;

            for (int i = 0; i <= j; ++i)
            {
                double* vi[4] = {ctx->d_krylov_v[i][0], ctx->d_krylov_v[i][1], ctx->d_krylov_v[i][2], ctx->d_krylov_v[i][3]};
                double h_ij = 0.0;
                int drc = gpu_wgt_dot4(ctx, vj1, vi, &h_ij);
                if (drc != C2_OK || !std::isfinite(h_ij)) { abort_to_fp64 = true; break; }
                H[i][j] += h_ij;
                axpy4_kernel<<<128, 256>>>(-h_ij, vi[0], vi[1], vi[2], vi[3], vj1[0], vj1[1], vj1[2], vj1[3],
                    ctx->tot_s_g0(0), ctx->tot_s_g0(1), ctx->tot_s_g0(2), ctx->tot_c_g0());
            }
            if (abort_to_fp64) break;
            double h_next = 0.0;
            int nrc = gpu_wgt_norm4(ctx, vj1, &h_next);
            if (nrc != C2_OK || !std::isfinite(h_next)) { abort_to_fp64 = true; break; }
            H[j + 1][j] = h_next;
            if (h_next > 1.0e-30)
            {
                double inv_h = 1.0 / h_next;
                for (int c = 0; c < 4; ++c)
                {
                    size_t n = (c < 3) ? ctx->tot_s_g0(c) : ctx->tot_c_g0();
                    linear_sum_side_kernel<<<128, 256>>>(vj1[c], inv_h, vj1[c], 0.0, vj1[c], n);
                }
            }
            for (int i = 0; i < j; ++i)
            {
                double t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                H[i][j] = t;
            }
            double den = std::hypot(H[j][j], H[j + 1][j]);
            cs[j] = (den > 0.0) ? H[j][j] / den : 1.0;
            sn[j] = (den > 0.0) ? H[j + 1][j] / den : 0.0;
            H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
            H[j + 1][j] = 0.0;
            double g_j = g[j];
            g[j] = cs[j] * g_j;
            g[j + 1] = -sn[j] * g_j;
            inner_done = j + 1;

            double est_abs = std::fabs(g[j + 1]);
            double est_rel = est_abs / b_norm;
            if (est_rel <= 0.9 * rel_tol || est_abs <= 0.9 * abs_tol) break;
        }
        if (abort_to_fp64)
        {
            if (!using_fp64 && allow_dynamic_fp64_fallback)
            {
                rep.fallback_triggered = true;
                rep.fallback_reason = "nan_inf_in_fp32";
                using_fp64 = true;
                continue;
            }
            cudaEventDestroy(ev0); cudaEventDestroy(ev1);
            return C2_ERR_NONFINITE;
        }
        if (inner_done > 0)
        {
            double y[M] = {};
            for (int i = inner_done - 1; i >= 0; --i)
            {
                double s = g[i];
                for (int k = i + 1; k < inner_done; ++k) s -= H[i][k] * y[k];
                y[i] = (std::fabs(H[i][i]) > 1.0e-30) ? (s / H[i][i]) : 0.0;
            }
            for (int i = 0; i < inner_done; ++i)
                axpy4_kernel<<<128, 256>>>(y[i], ctx->d_krylov_z[i][0], ctx->d_krylov_z[i][1], ctx->d_krylov_z[i][2], ctx->d_krylov_z[i][3],
                    x_g0[0], x_g0[1], x_g0[2], x_g0[3], ctx->tot_s_g0(0), ctx->tot_s_g0(1), ctx->tot_s_g0(2), ctx->tot_c_g0());
        }
        rc = apply_stokes_fp64(ctx, x_g0, r0);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
        for (int c = 0; c < 4; ++c)
        {
            size_t n = (c < 3) ? ctx->tot_s_g0(c) : ctx->tot_c_g0();
            linear_sum_side_kernel<<<128, 256>>>(r0[c], 1.0, rhs_g0[c], -1.0, r0[c], n);
        }
        double cur_abs = 0.0;
        rc = gpu_wgt_norm4(ctx, r0, &cur_abs);
        if (rc != C2_OK)
        {
            if (!using_fp64 && allow_dynamic_fp64_fallback)
            {
                rep.fallback_triggered = true; rep.fallback_reason = "nan_true_res"; using_fp64 = true; continue;
            }
            cudaEventDestroy(ev0); cudaEventDestroy(ev1);
            return C2_ERR_NONFINITE;
        }
        double cur_rel = cur_abs / b_norm;
        rep.final_true_abs_res = cur_abs;
        rep.final_true_rel_res = cur_rel;
        if (cur_rel <= rel_tol || cur_abs <= abs_tol) { rep.converged = true; break; }
        if (!using_fp64 && allow_dynamic_fp64_fallback)
        {
            fp32_cycles++;
            if (cur_rel > 0.90 * cycle_start_rel) low_prog++; else low_prog = 0;
            if (!std::isfinite(cur_rel)) { rep.fallback_triggered = true; rep.fallback_reason = "nan_true_res"; using_fp64 = true; }
            else if (low_prog >= 2) { rep.fallback_triggered = true; rep.fallback_reason = "stall_ge_0p9"; using_fp64 = true; }
            else if (fp32_cycles >= 4) { rep.fallback_triggered = true; rep.fallback_reason = "fp32_restarts_ge_4"; using_fp64 = true; }
        }
    }

    C2_CUDA_CALL(ctx, cudaEventRecord(ev1), "ev_st_rec1");
    C2_CUDA_CALL(ctx, cudaEventSynchronize(ev1), "ev_st_sync");
    float ms = 0.0f;
    C2_CUDA_CALL(ctx, cudaEventElapsedTime(&ms, ev0, ev1), "ev_st_el");
    rep.wall_ms = static_cast<double>(ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
    return C2_OK;
}

static int upload_side_array(C2Context* ctx, double* dst, const double* l0, const double* l1, int axis)
{
    if (!l0) return C2_ERR_INVALID;
    C2_CUDA_CALL(ctx, cudaMemcpy(dst, l0, ctx->l0_s_g0[axis] * sizeof(double), cudaMemcpyHostToDevice), "up_s_l0");
    if (ctx->n_l1 > 0 && ctx->l1_s_g0[axis] > 0)
    {
        if (!l1) return C2_ERR_INVALID;
        C2_CUDA_CALL(ctx, cudaMemcpy(dst + ctx->l0_s_g0[axis], l1, ctx->l1_s_g0[axis] * sizeof(double), cudaMemcpyHostToDevice), "up_s_l1");
    }
    return C2_OK;
}

static int upload_cell_array(C2Context* ctx, double* dst, const double* l0, const double* l1)
{
    if (!l0) return C2_ERR_INVALID;
    C2_CUDA_CALL(ctx, cudaMemcpy(dst, l0, ctx->l0_c_g0 * sizeof(double), cudaMemcpyHostToDevice), "up_c_l0");
    if (ctx->n_l1 > 0 && ctx->l1_c_g0 > 0)
    {
        if (!l1) return C2_ERR_INVALID;
        C2_CUDA_CALL(ctx, cudaMemcpy(dst + ctx->l0_c_g0, l1, ctx->l1_c_g0 * sizeof(double), cudaMemcpyHostToDevice), "up_c_l1");
    }
    return C2_OK;
}

static int download_side_array(C2Context* ctx, const double* src, double* l0, double* l1, int axis)
{
    if (l0)
        C2_CUDA_CALL(ctx, cudaMemcpy(l0, src, ctx->l0_s_g0[axis] * sizeof(double), cudaMemcpyDeviceToHost), "dn_s_l0");
    if (l1 && ctx->n_l1 > 0 && ctx->l1_s_g0[axis] > 0)
        C2_CUDA_CALL(ctx, cudaMemcpy(l1, src + ctx->l0_s_g0[axis], ctx->l1_s_g0[axis] * sizeof(double), cudaMemcpyDeviceToHost), "dn_s_l1");
    return C2_OK;
}

static int download_cell_array(C2Context* ctx, const double* src, double* l0, double* l1)
{
    if (l0)
        C2_CUDA_CALL(ctx, cudaMemcpy(l0, src, ctx->l0_c_g0 * sizeof(double), cudaMemcpyDeviceToHost), "dn_c_l0");
    if (l1 && ctx->n_l1 > 0 && ctx->l1_c_g0 > 0)
        C2_CUDA_CALL(ctx, cudaMemcpy(l1, src + ctx->l0_c_g0, ctx->l1_c_g0 * sizeof(double), cudaMemcpyDeviceToHost), "dn_c_l1");
    return C2_OK;
}

static int verify_device_finite(C2Context* ctx, const double* arr, size_t n)
{
    if (n == 0) return C2_OK;
    C2_CUDA_CALL(ctx, cudaMemset(ctx->d_nonfinite, 0, sizeof(int)), "nf_zero");
    check_finite_kernel<<<128, 256>>>(arr, n, ctx->d_nonfinite);
    C2_CUDA_LAUNCH_CHK(ctx, "check_finite");
    int flag = 0;
    C2_CUDA_CALL(ctx, cudaMemcpy(&flag, ctx->d_nonfinite, sizeof(int), cudaMemcpyDeviceToHost), "nf_d2h");
    if (flag != 0)
    {
        set_fb_reason(ctx, "nonfinite_gpu_output");
        return C2_ERR_NONFINITE;
    }
    return C2_OK;
}

} // namespace c2_gpu

extern "C" {

void* c2_gpu_create(const C2GridConfig* cfg, int* out_err)
{
    if (out_err) *out_err = C2_OK;
    if (!cfg || cfg->nx0 < 8 || cfg->ny0 < 8 || cfg->nz0 < 8 ||
        cfg->dx0 <= 0.0 || cfg->mu0 <= 0.0 || cfg->rho0 <= 0.0 || cfg->max_iterations < 1)
    {
        if (out_err) *out_err = C2_ERR_INVALID;
        return nullptr;
    }

    auto ctx = std::unique_ptr<c2_gpu::C2Context>(new (std::nothrow) c2_gpu::C2Context());
    if (!ctx)
    {
        if (out_err) *out_err = C2_ERR_CAPACITY;
        return nullptr;
    }
    ctx->pid = getpid();
    ctx->cfg = *cfg;
    ctx->nx0 = cfg->nx0; ctx->ny0 = cfg->ny0; ctx->nz0 = cfg->nz0;
    ctx->n_l0_cells = static_cast<size_t>(cfg->nx0) * cfg->ny0 * cfg->nz0;
    ctx->dx0 = cfg->dx0; ctx->dx1 = 0.5 * cfg->dx0;
    ctx->rho0 = cfg->rho0; ctx->mu0 = cfg->mu0; ctx->d_vel_const = -0.5 * cfg->mu0;
    ctx->domain_vol = (cfg->nx0 * cfg->dx0) * (cfg->ny0 * cfg->dx0) * (cfg->nz0 * cfg->dx0);
    c2_gpu::set_fb_reason(ctx.get(), "none");

    ctx->l0_c_g0 = ctx->n_l0_cells;
    ctx->l0_c_g1 = static_cast<size_t>(cfg->nx0 + 2) * (cfg->ny0 + 2) * (cfg->nz0 + 2);
    size_t max_l0_s1 = 0, max_l0_s3 = 0;
    for (int a = 0; a < 3; ++a)
    {
        ctx->l0_s_g0[a] = static_cast<size_t>(cfg->nx0 + (a == 0)) * (cfg->ny0 + (a == 1)) * (cfg->nz0 + (a == 2));
        ctx->l0_s_g1[a] = static_cast<size_t>(cfg->nx0 + 2 + (a == 0)) * (cfg->ny0 + 2 + (a == 1)) * (cfg->nz0 + 2 + (a == 2));
        ctx->l0_s_g3[a] = static_cast<size_t>(cfg->nx0 + 6 + (a == 0)) * (cfg->ny0 + 6 + (a == 1)) * (cfg->nz0 + 6 + (a == 2));
        max_l0_s1 = std::max(max_l0_s1, ctx->l0_s_g1[a]);
        max_l0_s3 = std::max(max_l0_s3, ctx->l0_s_g3[a]);
    }
    if (ctx->n_l0_cells >= 2000000ULL)
    {
        ctx->max_slot_g0 = std::max<size_t>(4480000ULL, max_l0_s3 + 65536ULL);
        ctx->max_slot_g1 = std::max<size_t>(4820000ULL, max_l0_s1 + 65536ULL);
    }
    else
    {
        ctx->max_slot_g0 = std::max<size_t>(4ULL * ctx->n_l0_cells + max_l0_s3, 800000ULL);
        ctx->max_slot_g1 = std::max<size_t>(ctx->max_slot_g0 + 540000ULL, max_l0_s1 + 262144ULL);
    }
    ctx->tel.max_l1_patches = ctx->max_l1_patches;
    ctx->tel.max_slot_g0 = static_cast<unsigned long long>(ctx->max_slot_g0);
    ctx->tel.max_slot_g1 = static_cast<unsigned long long>(ctx->max_slot_g1);

    if (!c2_gpu::build_coarse_mg_chain(ctx.get()))
    {
        if (out_err) *out_err = C2_ERR_INVALID;
        return nullptr;
    }

    c2_gpu::snapshot_pre_cuda_nvml_pids(ctx.get());

    if (cudaSetDevice(0) != cudaSuccess)
    {
        if (out_err) *out_err = C2_ERR_CUDA;
        return nullptr;
    }

    auto fail_cuda = [&]() -> void* {
        c2_gpu::free_ctx_buffers(ctx.get());
        if (out_err) *out_err = C2_ERR_CUDA;
        return nullptr;
    };

    if (cudaMemcpyToSymbol(c2_gpu::c_nx0, &ctx->nx0, sizeof(int)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_ny0, &ctx->ny0, sizeof(int)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_nz0, &ctx->nz0, sizeof(int)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_n_l0_cells, &ctx->n_l0_cells, sizeof(size_t)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_dx0, &ctx->dx0, sizeof(double)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_dx1, &ctx->dx1, sizeof(double)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_rho0, &ctx->rho0, sizeof(double)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_mu0, &ctx->mu0, sizeof(double)) != cudaSuccess ||
        cudaMemcpyToSymbol(c2_gpu::c_d_vel_const, &ctx->d_vel_const, sizeof(double)) != cudaSuccess)
        return fail_cuda();

    // Stage 0
    int* d_dummy = nullptr;
    if (cudaMalloc(&d_dummy, 256) != cudaSuccess) return fail_cuda();
    cudaMemset(d_dummy, 0, 256);
    c2_gpu::sample_vram(ctx.get(), 0);
    cudaFree(d_dummy);

    // Stage 1: 24 Persistent FP64 State arrays
    for (int i = 0; i < 24; ++i)
    {
        if (cudaMalloc(&ctx->d_state[i], ctx->max_slot_g1 * sizeof(double)) != cudaSuccess ||
            cudaMemset(ctx->d_state[i], 0, ctx->max_slot_g1 * sizeof(double)) != cudaSuccess)
            return fail_cuda();
    }
    c2_gpu::sample_vram(ctx.get(), 1);

    // Stage 2: 44 Krylov + SSPRK3 arena slots
    for (int j = 0; j < 6; ++j)
        for (int c = 0; c < 4; ++c)
            if (cudaMalloc(&ctx->d_krylov_v[j][c], ctx->max_slot_g0 * sizeof(double)) != cudaSuccess ||
                cudaMemset(ctx->d_krylov_v[j][c], 0, ctx->max_slot_g0 * sizeof(double)) != cudaSuccess)
                return fail_cuda();
    for (int j = 0; j < 5; ++j)
        for (int c = 0; c < 4; ++c)
            if (cudaMalloc(&ctx->d_krylov_z[j][c], ctx->max_slot_g0 * sizeof(double)) != cudaSuccess ||
                cudaMemset(ctx->d_krylov_z[j][c], 0, ctx->max_slot_g0 * sizeof(double)) != cudaSuccess)
                return fail_cuda();
    c2_gpu::sample_vram(ctx.get(), 2);

    // Stage 3: Dual-Prec MG + Coarse + Maps + Weights
    for (int i = 0; i < 10; ++i)
        if (cudaMalloc(&ctx->d_mg_fp64[i], ctx->max_slot_g1 * sizeof(double)) != cudaSuccess ||
            cudaMemset(ctx->d_mg_fp64[i], 0, ctx->max_slot_g1 * sizeof(double)) != cudaSuccess)
            return fail_cuda();
    for (int a = 0; a < 3; ++a)
    {
        if (cudaMalloc(&ctx->d_mg_c_fp32[a], ctx->max_slot_g1 * sizeof(float)) != cudaSuccess ||
            cudaMalloc(&ctx->d_mg_d_fp32[a], ctx->max_slot_g1 * sizeof(float)) != cudaSuccess ||
            cudaMemset(ctx->d_mg_c_fp32[a], 0, ctx->max_slot_g1 * sizeof(float)) != cudaSuccess ||
            cudaMemset(ctx->d_mg_d_fp32[a], 0, ctx->max_slot_g1 * sizeof(float)) != cudaSuccess)
            return fail_cuda();
        for (int k = 0; k < 3; ++k)
        {
            size_t cl_sz = (k == 1) ? (ctx->cl_c1_total_g1 + 64) : (ctx->cl_total_g1 + 64);
            if (cudaMalloc(&ctx->d_cl_fp64_vec[k][a], cl_sz * sizeof(double)) != cudaSuccess ||
                cudaMemset(ctx->d_cl_fp64_vec[k][a], 0, cl_sz * sizeof(double)) != cudaSuccess)
                return fail_cuda();
        }
        size_t c1_sz = ctx->cl_c1_total_g1 + 64;
        if (cudaMalloc(&ctx->d_cl_fp64_c[a], c1_sz * sizeof(double)) != cudaSuccess ||
            cudaMalloc(&ctx->d_cl_fp64_d[a], c1_sz * sizeof(double)) != cudaSuccess ||
            cudaMalloc(&ctx->d_cl_fp32_c[a], c1_sz * sizeof(float)) != cudaSuccess ||
            cudaMalloc(&ctx->d_cl_fp32_d[a], c1_sz * sizeof(float)) != cudaSuccess ||
            cudaMemset(ctx->d_cl_fp64_c[a], 0, c1_sz * sizeof(double)) != cudaSuccess ||
            cudaMemset(ctx->d_cl_fp64_d[a], 0, c1_sz * sizeof(double)) != cudaSuccess ||
            cudaMemset(ctx->d_cl_fp32_c[a], 0, c1_sz * sizeof(float)) != cudaSuccess ||
            cudaMemset(ctx->d_cl_fp32_d[a], 0, c1_sz * sizeof(float)) != cudaSuccess)
            return fail_cuda();
    }
    if (cudaMalloc(&ctx->d_l1_map, ctx->n_l0_cells * sizeof(int16_t)) != cudaSuccess ||
        cudaMalloc(&ctx->d_l1_boxes, ctx->max_l1_patches * sizeof(c2_gpu::Box3i)) != cudaSuccess ||
        cudaMalloc(&ctx->d_wgt_c, ctx->max_slot_g0 * sizeof(double)) != cudaSuccess)
        return fail_cuda();
    for (int a = 0; a < 3; ++a)
        if (cudaMalloc(&ctx->d_wgt_s[a], ctx->max_slot_g0 * sizeof(double)) != cudaSuccess)
            return fail_cuda();
    if (cudaMalloc(&ctx->d_red_scratch, 256 * sizeof(double)) != cudaSuccess ||
        cudaMallocHost(&ctx->h_red_pinned, 512 * sizeof(double)) != cudaSuccess)
        return fail_cuda();
    c2_gpu::sample_vram(ctx.get(), 3);

    // Stage 4: Persistent Brinkman penalty d_L_b[3] + pinned H2D/D2H staging + finiteness flag (replaces d_pcie_crop + d_pcie_full)
    for (int a = 0; a < 3; ++a)
        if (cudaMalloc(&ctx->d_L_b[a], ctx->max_slot_g0 * sizeof(double)) != cudaSuccess ||
            cudaMemset(ctx->d_L_b[a], 0, ctx->max_slot_g0 * sizeof(double)) != cudaSuccess)
            return fail_cuda();
    if (cudaMallocHost(&ctx->h_pcie_pinned, ctx->max_slot_g0 * sizeof(double)) != cudaSuccess ||
        cudaMalloc(&ctx->d_nonfinite, sizeof(int)) != cudaSuccess)
        return fail_cuda();
    c2_gpu::sample_vram(ctx.get(), 4);

    ctx->h_l1_map.assign(ctx->n_l0_cells, -1);
    return ctx.release();
}

void c2_gpu_destroy(void* handle)
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx) return;
    c2_gpu::free_ctx_buffers(ctx);
    delete ctx;
}

int c2_gpu_regrid(void* handle, int n_l1, const C2PatchBox* l1_boxes)
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx || n_l1 < 0 || (n_l1 > 0 && !l1_boxes)) return C2_ERR_INVALID;
    ctx->tel.regrid_count += 1;
    ctx->tel.last_n_l1 = n_l1;
    ctx->tel.peak_n_l1 = std::max(ctx->tel.peak_n_l1, n_l1);
    ctx->tel.max_l1_patches = ctx->max_l1_patches;
    ctx->tel.max_slot_g0 = static_cast<unsigned long long>(ctx->max_slot_g0);
    ctx->tel.max_slot_g1 = static_cast<unsigned long long>(ctx->max_slot_g1);
    if (n_l1 > ctx->max_l1_patches)
    {
        c2_gpu::set_fb_reason(ctx, "capacity_n_l1_patches");
        return C2_ERR_CAPACITY;
    }

    cudaEvent_t ev0{}, ev1{};
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);

    ctx->has_valid_grid = false;
    ctx->n_l1 = n_l1;
    ctx->h_l1_boxes.resize(n_l1);
    ctx->h_l1_map.assign(ctx->n_l0_cells, -1);

    size_t cg0 = 0, cg1 = 0;
    size_t sg0[3] = {0, 0, 0}, sg1[3] = {0, 0, 0}, sg2[3] = {0, 0, 0}, sg3[3] = {0, 0, 0};
    for (int p = 0; p < n_l1; ++p)
    {
        const C2PatchBox& raw = l1_boxes[p];
        if (raw.lo[0] < 0 || raw.lo[1] < 0 || raw.lo[2] < 0 ||
            raw.hi[0] < raw.lo[0] || raw.hi[1] < raw.lo[1] || raw.hi[2] < raw.lo[2] ||
            raw.hi[0] >= 2 * ctx->nx0 || raw.hi[1] >= 2 * ctx->ny0 || raw.hi[2] >= 2 * ctx->nz0 ||
            (raw.lo[0] & 1) || (raw.lo[1] & 1) || (raw.lo[2] & 1) ||
            ((raw.hi[0] + 1) & 1) || ((raw.hi[1] + 1) & 1) || ((raw.hi[2] + 1) & 1))
        {
            c2_gpu::set_fb_reason(ctx, "invalid_l1_box_bounds");
            cudaEventDestroy(ev0); cudaEventDestroy(ev1);
            return C2_ERR_INVALID;
        }
        c2_gpu::Box3i b{};
        b.lo[0] = raw.lo[0]; b.lo[1] = raw.lo[1]; b.lo[2] = raw.lo[2];
        b.hi[0] = raw.hi[0]; b.hi[1] = raw.hi[1]; b.hi[2] = raw.hi[2];
        b.nx = b.hi[0] - b.lo[0] + 1;
        b.ny = b.hi[1] - b.lo[1] + 1;
        b.nz = b.hi[2] - b.lo[2] + 1;
        b.off_c_g0 = cg0;
        b.off_c_g1 = cg1;
        cg0 += static_cast<size_t>(b.nx) * b.ny * b.nz;
        cg1 += static_cast<size_t>(b.nx + 2) * (b.ny + 2) * (b.nz + 2);
        for (int a = 0; a < 3; ++a)
        {
            const int sx0 = b.nx + (a == 0 ? 1 : 0);
            const int sy0 = b.ny + (a == 1 ? 1 : 0);
            const int sz0 = b.nz + (a == 2 ? 1 : 0);
            b.off_s_g0[a] = sg0[a];
            b.off_s_g1[a] = sg1[a];
            b.off_s_g2[a] = sg2[a];
            b.off_s_g3[a] = sg3[a];
            sg0[a] += static_cast<size_t>(sx0) * sy0 * sz0;
            sg1[a] += static_cast<size_t>(sx0 + 2) * (sy0 + 2) * (sz0 + 2);
            sg2[a] += static_cast<size_t>(sx0 + 4) * (sy0 + 4) * (sz0 + 4);
            sg3[a] += static_cast<size_t>(sx0 + 6) * (sy0 + 6) * (sz0 + 6);
        }
        ctx->h_l1_boxes[p] = b;
        for (int k0 = b.lo[2] >> 1; k0 <= (b.hi[2] >> 1); ++k0)
            for (int j0 = b.lo[1] >> 1; j0 <= (b.hi[1] >> 1); ++j0)
                for (int i0 = b.lo[0] >> 1; i0 <= (b.hi[0] >> 1); ++i0)
                    ctx->h_l1_map[i0 + ctx->nx0 * (j0 + ctx->ny0 * k0)] = static_cast<int16_t>(p);
    }

    ctx->l1_c_g0 = cg0;
    ctx->l1_c_g1 = cg1;
    for (int a = 0; a < 3; ++a)
    {
        ctx->l1_s_g0[a] = sg0[a];
        ctx->l1_s_g1[a] = sg1[a];
        ctx->l1_s_g2[a] = sg2[a];
        ctx->l1_s_g3[a] = sg3[a];
    }

    const unsigned long long tot_cg0 = static_cast<unsigned long long>(ctx->tot_c_g0());
    const unsigned long long tot_cg1 = static_cast<unsigned long long>(ctx->tot_c_g1());
    ctx->tel.last_tot_c_g0 = tot_cg0;
    ctx->tel.last_tot_c_g1 = tot_cg1;
    ctx->tel.peak_tot_c_g0 = std::max(ctx->tel.peak_tot_c_g0, tot_cg0);
    ctx->tel.peak_tot_c_g1 = std::max(ctx->tel.peak_tot_c_g1, tot_cg1);
    unsigned long long max_sg0 = 0ULL, max_sg1 = 0ULL, max_sg3 = 0ULL;
    for (int a = 0; a < 3; ++a)
    {
        const unsigned long long sg0_a = static_cast<unsigned long long>(ctx->tot_s_g0(a));
        const unsigned long long sg1_a = static_cast<unsigned long long>(ctx->tot_s_g1(a));
        const unsigned long long l0_sg3_a = static_cast<unsigned long long>(ctx->l0_s_g3[a]);
        const unsigned long long l1_sg3_a = static_cast<unsigned long long>(ctx->l1_s_g3[a]);
        ctx->tel.last_tot_s_g0[a] = sg0_a;
        ctx->tel.last_tot_s_g1[a] = sg1_a;
        ctx->tel.last_l0_s_g3[a] = l0_sg3_a;
        ctx->tel.last_l1_s_g3[a] = l1_sg3_a;
        ctx->tel.peak_tot_s_g0[a] = std::max(ctx->tel.peak_tot_s_g0[a], sg0_a);
        ctx->tel.peak_tot_s_g1[a] = std::max(ctx->tel.peak_tot_s_g1[a], sg1_a);
        ctx->tel.peak_l0_s_g3[a] = std::max(ctx->tel.peak_l0_s_g3[a], l0_sg3_a);
        ctx->tel.peak_l1_s_g3[a] = std::max(ctx->tel.peak_l1_s_g3[a], l1_sg3_a);
        max_sg0 = std::max(max_sg0, sg0_a);
        max_sg1 = std::max(max_sg1, sg1_a);
        max_sg3 = std::max(max_sg3, std::max(l0_sg3_a, l1_sg3_a));
    }
    ctx->tel.last_max_s_g0 = max_sg0;
    ctx->tel.last_max_s_g1 = max_sg1;
    ctx->tel.last_max_s_g3 = max_sg3;
    ctx->tel.peak_max_s_g0 = std::max(ctx->tel.peak_max_s_g0, max_sg0);
    ctx->tel.peak_max_s_g1 = std::max(ctx->tel.peak_max_s_g1, max_sg1);
    ctx->tel.peak_max_s_g3 = std::max(ctx->tel.peak_max_s_g3, max_sg3);

    bool fit = (ctx->tot_c_g0() <= ctx->max_slot_g0) && (ctx->tot_c_g1() <= ctx->max_slot_g1);
    for (int a = 0; a < 3; ++a)
    {
        fit = fit && (ctx->tot_s_g0(a) <= ctx->max_slot_g0) && (ctx->tot_s_g1(a) <= ctx->max_slot_g1) &&
              (ctx->l0_s_g3[a] <= ctx->max_slot_g0) && (ctx->l1_s_g3[a] <= ctx->max_slot_g0);
    }
    if (!fit)
    {
        c2_gpu::set_fb_reason(ctx, "capacity_slot_overflow");
        cudaEventDestroy(ev0); cudaEventDestroy(ev1);
        return C2_ERR_CAPACITY;
    }

    // Build exact IBAMR composite integration weights wc, ws[3]
    std::vector<double> wc(ctx->tot_c_g0(), 0.0), ws[3];
    for (int a = 0; a < 3; ++a) ws[a].assign(ctx->tot_s_g0(a), 0.0);
    const double v0 = ctx->dx0 * ctx->dx0 * ctx->dx0;
    const double v1 = ctx->dx1 * ctx->dx1 * ctx->dx1;
    const int nx0 = ctx->nx0, ny0 = ctx->ny0, nz0 = ctx->nz0;
    for (size_t i = 0; i < ctx->n_l0_cells; ++i)
        wc[i] = (ctx->h_l1_map[i] >= 0) ? 0.0 : v0;
    for (size_t i = 0; i < ctx->l1_c_g0; ++i)
        wc[ctx->l0_c_g0 + i] = v1;

    auto cov0 = [&](int i0, int j0, int k0) -> bool {
        return ctx->h_l1_map[i0 + nx0 * (j0 + ny0 * k0)] >= 0;
    };
    for (int k = 0; k < nz0; ++k)
        for (int j = 0; j < ny0; ++j)
            for (int i = 0; i <= nx0; ++i)
            {
                bool cr = (i < nx0) && cov0(i, j, k);
                bool cl = (i > 0) && cov0(i - 1, j, k);
                size_t off = c2_gpu::sidx_g0(0, i, j, k, nx0, ny0);
                ws[0][off] = (cr && cl) ? 0.0 : ((cr || cl) ? 0.5 * v0 : ((i == 0 || i == nx0 ? 0.5 : 1.0) * v0));
            }
    for (int k = 0; k < nz0; ++k)
        for (int j = 0; j <= ny0; ++j)
            for (int i = 0; i < nx0; ++i)
            {
                bool cr = (j < ny0) && cov0(i, j, k);
                bool cl = (j > 0) && cov0(i, j - 1, k);
                size_t off = c2_gpu::sidx_g0(1, i, j, k, nx0, ny0);
                ws[1][off] = (cr && cl) ? 0.0 : ((cr || cl) ? 0.5 * v0 : ((j == 0 || j == ny0 ? 0.5 : 1.0) * v0));
            }
    for (int k = 0; k <= nz0; ++k)
        for (int j = 0; j < ny0; ++j)
            for (int i = 0; i < nx0; ++i)
            {
                bool cr = (k < nz0) && cov0(i, j, k);
                bool cl = (k > 0) && cov0(i, j, k - 1);
                size_t off = c2_gpu::sidx_g0(2, i, j, k, nx0, ny0);
                ws[2][off] = (cr && cl) ? 0.0 : ((cr || cl) ? 0.5 * v0 : ((k == 0 || k == nz0 ? 0.5 : 1.0) * v0));
            }

    auto lookup_h = [&](int gi, int gj, int gk) -> int {
        if (gi < 0 || gi >= 2 * nx0 || gj < 0 || gj >= 2 * ny0 || gk < 0 || gk >= 2 * nz0) return -1;
        return ctx->h_l1_map[(gi >> 1) + nx0 * ((gj >> 1) + ny0 * (gk >> 1))];
    };
    for (int p = 0; p < n_l1; ++p)
    {
        const c2_gpu::Box3i& b = ctx->h_l1_boxes[p];
        for (int k = 0; k < b.nz; ++k)
            for (int j = 0; j < b.ny; ++j)
                for (int i = 0; i <= b.nx; ++i)
                {
                    int gi = b.lo[0] + i, gj = b.lo[1] + j, gk = b.lo[2] + k;
                    if (i == b.nx && gi < 2 * nx0 && lookup_h(gi, gj, gk) >= 0) continue;
                    bool hr = lookup_h(gi, gj, gk) >= 0, hl = lookup_h(gi - 1, gj, gk) >= 0;
                    ws[0][ctx->l0_s_g0[0] + b.off_s_g0[0] + c2_gpu::sidx_g0(0, i, j, k, b.nx, b.ny)] = (hr && hl) ? v1 : 0.5 * v1;
                }
        for (int k = 0; k < b.nz; ++k)
            for (int j = 0; j <= b.ny; ++j)
                for (int i = 0; i < b.nx; ++i)
                {
                    int gi = b.lo[0] + i, gj = b.lo[1] + j, gk = b.lo[2] + k;
                    if (j == b.ny && gj < 2 * ny0 && lookup_h(gi, gj, gk) >= 0) continue;
                    bool hr = lookup_h(gi, gj, gk) >= 0, hl = lookup_h(gi, gj - 1, gk) >= 0;
                    ws[1][ctx->l0_s_g0[1] + b.off_s_g0[1] + c2_gpu::sidx_g0(1, i, j, k, b.nx, b.ny)] = (hr && hl) ? v1 : 0.5 * v1;
                }
        for (int k = 0; k <= b.nz; ++k)
            for (int j = 0; j < b.ny; ++j)
                for (int i = 0; i < b.nx; ++i)
                {
                    int gi = b.lo[0] + i, gj = b.lo[1] + j, gk = b.lo[2] + k;
                    if (k == b.nz && gk < 2 * nz0 && lookup_h(gi, gj, gk) >= 0) continue;
                    bool hr = lookup_h(gi, gj, gk) >= 0, hl = lookup_h(gi, gj, gk - 1) >= 0;
                    ws[2][ctx->l0_s_g0[2] + b.off_s_g0[2] + c2_gpu::sidx_g0(2, i, j, k, b.nx, b.ny)] = (hr && hl) ? v1 : 0.5 * v1;
                }
    }

    if (cudaMemcpy(ctx->d_l1_map, ctx->h_l1_map.data(), ctx->n_l0_cells * sizeof(int16_t), cudaMemcpyHostToDevice) != cudaSuccess ||
        (n_l1 > 0 && cudaMemcpy(ctx->d_l1_boxes, ctx->h_l1_boxes.data(), n_l1 * sizeof(c2_gpu::Box3i), cudaMemcpyHostToDevice) != cudaSuccess) ||
        cudaMemcpy(ctx->d_wgt_c, wc.data(), ctx->tot_c_g0() * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
    {
        c2_gpu::set_fb_reason(ctx, "cuda_regrid_upload");
        cudaEventDestroy(ev0); cudaEventDestroy(ev1);
        return C2_ERR_CUDA;
    }
    for (int a = 0; a < 3; ++a)
    {
        if (cudaMemcpy(ctx->d_wgt_s[a], ws[a].data(), ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyHostToDevice) != cudaSuccess)
        {
            c2_gpu::set_fb_reason(ctx, "cuda_regrid_wgt_s");
            cudaEventDestroy(ev0); cudaEventDestroy(ev1);
            return C2_ERR_CUDA;
        }
    }

    ctx->has_valid_grid = true;
    ctx->has_brinkman = false;
    cudaEventRecord(ev1); cudaEventSynchronize(ev1);
    float ms = 0.0f; cudaEventElapsedTime(&ms, ev0, ev1);
    ctx->tel.ms_regrid = static_cast<double>(ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
    c2_gpu::sample_vram(ctx, -1);
    return C2_OK;
}

int c2_gpu_upload_state(
    void* handle,
    const double* const u_cur_l0[3],
    const double* const u_cur_l1[3],
    const double* const u_old_l0[3],
    const double* const u_old_l1[3],
    const double* p_cur_l0,
    const double* p_cur_l1,
    const double* const rho_cur_l0[3],
    const double* const rho_cur_l1[3])
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx || !ctx->has_valid_grid || !u_cur_l0 || !p_cur_l0 || !rho_cur_l0) return C2_ERR_INVALID;

    cudaEvent_t ev0{}, ev1{};
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);

    for (int a = 0; a < 3; ++a)
    {
        int rc = c2_gpu::upload_side_array(ctx, ctx->d_state[a], u_cur_l0[a], u_cur_l1 ? u_cur_l1[a] : nullptr, a);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
        if (u_old_l0 && u_old_l0[a])
        {
            rc = c2_gpu::upload_side_array(ctx, ctx->d_state[3 + a], u_old_l0[a], u_old_l1 ? u_old_l1[a] : nullptr, a);
            if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
        }
        else
        {
            C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[3 + a], ctx->d_state[a], ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyDeviceToDevice), "up_old_dup");
        }
        rc = c2_gpu::upload_side_array(ctx, ctx->d_state[15 + a], rho_cur_l0[a], rho_cur_l1 ? rho_cur_l1[a] : nullptr, a);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
    }
    int rc = c2_gpu::upload_cell_array(ctx, ctx->d_state[21], p_cur_l0, p_cur_l1);
    if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }

    cudaEventRecord(ev1); cudaEventSynchronize(ev1);
    float ms = 0.0f; cudaEventElapsedTime(&ms, ev0, ev1);
    ctx->tel.ms_h2d += static_cast<double>(ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
    return C2_OK;
}

int c2_gpu_upload_brinkman(
    void* handle,
    const double* const L_l0[3],
    const double* const L_l1[3])
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx || !ctx->has_valid_grid || !L_l0) return C2_ERR_INVALID;

    cudaEvent_t ev0{}, ev1{};
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);

    const dim3 l1g(16, std::max(1, ctx->n_l1));
    for (int a = 0; a < 3; ++a)
    {
        int rc = c2_gpu::upload_side_array(ctx, ctx->d_L_b[a], L_l0[a], L_l1 ? L_l1[a] : nullptr, a);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
        if (ctx->n_l1 > 0)
        {
            c2_gpu::sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(ctx->d_L_b[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            c2_gpu::restrict_g0_l1_to_l0_side_kernel<<<l1g, 256>>>(ctx->d_L_b[a] + ctx->l0_s_g0[a], ctx->d_L_b[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }
    C2_CUDA_LAUNCH_CHK(ctx, "upload_brinkman_sync");
    ctx->has_brinkman = true;

    cudaEventRecord(ev1); cudaEventSynchronize(ev1);
    float ms = 0.0f; cudaEventElapsedTime(&ms, ev0, ev1);
    ctx->tel.ms_h2d += static_cast<double>(ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
    return C2_OK;
}

int c2_gpu_upload_rho(
    void* handle,
    const double* const rho_l0[3],
    const double* const rho_l1[3])
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx || !ctx->has_valid_grid || !rho_l0) return C2_ERR_INVALID;

    cudaEvent_t ev0{}, ev1{};
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);

    for (int a = 0; a < 3; ++a)
    {
        int rc = c2_gpu::upload_side_array(ctx, ctx->d_state[15 + a], rho_l0[a], rho_l1 ? rho_l1[a] : nullptr, a);
        if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
    }

    cudaEventRecord(ev1); cudaEventSynchronize(ev1);
    float ms = 0.0f; cudaEventElapsedTime(&ms, ev0, ev1);
    ctx->tel.ms_h2d += static_cast<double>(ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
    return C2_OK;
}

int c2_gpu_preprocess(
    void* handle,
    double dt,
    double dt_prev,
    int is_first_step,
    double constant_rho)
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx || !ctx->has_valid_grid || dt <= 0.0) return C2_ERR_INVALID;

    cudaEvent_t ev0{}, ev1{};
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);

    const dim3 l1g(16, std::max(1, ctx->n_l1));
    double* rho_cur_p[3] = {ctx->d_state[15], ctx->d_state[16], ctx->d_state[17]};
    if (constant_rho > 0.0)
    {
        for (int a = 0; a < 3; ++a)
            c2_gpu::fill_const_kernel<<<128, 256>>>(rho_cur_p[a], constant_rho, ctx->tot_s_g0(a));
    }

    // 1. Viscous RHS: urhs_p = +0.5 * mu0 * laplace(U_cur) + (1/dt) * rho_cur * U_cur
    double* g4[4] = {ctx->d_state[0], ctx->d_state[1], ctx->d_state[2], ctx->d_state[21]};
    int rc = c2_gpu::fill_native_stokes_ghosts(ctx, g4, true);
    if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }

    double* ul0[3] = {ctx->d_mg_fp64[0], ctx->d_mg_fp64[1], ctx->d_mg_fp64[2]};
    double* ul1[3] = {ul0[0] + ctx->l0_s_g1[0], ul0[1] + ctx->l0_s_g1[1], ul0[2] + ctx->l0_s_g1[2]};
    double* urhs_p[3] = {ctx->d_state[18], ctx->d_state[19], ctx->d_state[20]};
    double* urhs_l1[3] = {urhs_p[0] + ctx->l0_s_g0[0], urhs_p[1] + ctx->l0_s_g0[1], urhs_p[2] + ctx->l0_s_g0[2]};

    c2_gpu::apply_momentum_stokes_kernel<<<dim3(64, 1), 256>>>(
        ul0[0], ul0[1], ul0[2], ctx->d_mg_fp64[3], ctx->d_state[9], ctx->d_state[10], ctx->d_state[11],
        0.5 * ctx->mu0, 0.0, urhs_p[0], urhs_p[1], urhs_p[2],
        ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx0, false, false, false, false, nullptr, 1, ctx->d_l1_map);
    if (ctx->n_l1 > 0)
    {
        c2_gpu::apply_momentum_stokes_kernel<<<l1g, 256>>>(
            ul1[0], ul1[1], ul1[2], ctx->d_mg_fp64[3] + ctx->l0_c_g1,
            ctx->d_state[9] + ctx->l0_s_g1[0], ctx->d_state[10] + ctx->l0_s_g1[1], ctx->d_state[11] + ctx->l0_s_g1[2],
            0.5 * ctx->mu0, 0.0, urhs_l1[0], urhs_l1[1], urhs_l1[2],
            ctx->nx0, ctx->ny0, ctx->nz0, ctx->dx1, false, false, false, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
    }
    for (int a = 0; a < 3; ++a)
    {
        if (ctx->n_l1 > 0)
        {
            c2_gpu::sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(urhs_l1[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            c2_gpu::outerside_restrict_cf_l1_to_l0_side_kernel<<<l1g, 256>>>(urhs_l1[a], urhs_p[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
        c2_gpu::add_scaled_prod_side_g1_to_g0_kernel<<<dim3(64, 1), 256>>>(
            urhs_p[a], rho_cur_p[a], ul0[a], 1.0 / dt, a, ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1, ctx->d_l1_map);
        if (ctx->n_l1 > 0)
        {
            c2_gpu::add_scaled_prod_side_g1_to_g0_kernel<<<l1g, 256>>>(
                urhs_l1[a], rho_cur_p[a] + ctx->l0_s_g0[a], ul1[a], 1.0 / dt, a,
                ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map);
            c2_gpu::sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(urhs_l1[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }
    C2_CUDA_LAUNCH_CHK(ctx, "viscous_rhs");

    // Set initial guess (U_new = U_cur, P_new = P_cur)
    for (int a = 0; a < 3; ++a)
        C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[6 + a], ctx->d_state[a], ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyDeviceToDevice), "init_unew");
    C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[22], ctx->d_state[21], ctx->tot_c_g0() * sizeof(double), cudaMemcpyDeviceToDevice), "init_pnew");

    // 2. SSPRK3 density + convective derivative N_full (into d_krylov_z[0][0..2]) and rho_new (into d_state[15..17])
    double* u_cur_p[3] = {ctx->d_state[0], ctx->d_state[1], ctx->d_state[2]};
    double* u_old_p[3] = {ctx->d_state[3], ctx->d_state[4], ctx->d_state[5]};
    double* rho_tmp_p[3] = {ctx->d_krylov_v[5][0], ctx->d_krylov_v[5][1], ctx->d_krylov_v[5][2]};
    double* n_full_p[3] = {ctx->d_krylov_z[0][0], ctx->d_krylov_z[0][1], ctx->d_krylov_z[0][2]};
    double ssprk_ms = 0.0;
    rc = c2_gpu::run_gpu_ssprk3_step(
        ctx, dt, dt_prev, (is_first_step != 0), u_old_p, u_cur_p, rho_cur_p, rho_tmp_p, n_full_p, &ssprk_ms);
    if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }

    for (int a = 0; a < 3; ++a)
        C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[15 + a], rho_tmp_p[a], ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyDeviceToDevice), "persist_rhonew");

    cudaEventRecord(ev1); cudaEventSynchronize(ev1);
    float tot_ms = 0.0f; cudaEventElapsedTime(&tot_ms, ev0, ev1);
    ctx->tel.ms_ssprk = ssprk_ms;
    ctx->tel.ms_preprocess = static_cast<double>(tot_ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);
    c2_gpu::sample_vram(ctx, -1);
    return C2_OK;
}

int c2_gpu_integrate_solve(
    void* handle,
    double dt,
    double p_in_new,
    const double* const f_l0[3],
    const double* const f_l1[3],
    int has_f,
    double* const u_new_l0[3],
    double* const u_new_l1[3],
    double* p_new_l0,
    double* p_new_l1,
    double* const rho_new_l0[3],
    double* const rho_new_l1[3],
    int download_rho,
    int* out_converged,
    int* out_iters,
    double* out_rel_res)
{
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx || !ctx->has_valid_grid || !ctx->has_brinkman || dt <= 0.0 || !u_new_l0 || !p_new_l0)
        return C2_ERR_INVALID;

    if (out_converged) *out_converged = 0;
    if (out_iters) *out_iters = 0;
    if (out_rel_res) *out_rel_res = 0.0;

    cudaEvent_t ev0{}, ev1{};
    cudaEventCreate(&ev0); cudaEventCreate(&ev1);
    cudaEventRecord(ev0);

    const dim3 l1g(16, std::max(1, ctx->n_l1));
    double* rho_new_p[3] = {ctx->d_state[15], ctx->d_state[16], ctx->d_state[17]};
    double* urhs_p[3] = {ctx->d_state[18], ctx->d_state[19], ctx->d_state[20]};
    double* n_full_p[3] = {ctx->d_krylov_z[0][0], ctx->d_krylov_z[0][1], ctx->d_krylov_z[0][2]};

    // 1. Sync rho_new, compute C_sc & P_D from (rho_new, d_L_b), and coarsen MG coefs
    for (int a = 0; a < 3; ++a)
    {
        if (ctx->n_l1 > 0)
        {
            c2_gpu::sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(rho_new_p[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            c2_gpu::restrict_g0_l1_to_l0_side_kernel<<<l1g, 256>>>(rho_new_p[a] + ctx->l0_s_g0[a], rho_new_p[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
        c2_gpu::coef_update_g1_kernel<<<dim3(64, 1), 256>>>(
            rho_new_p[a], ctx->d_L_b[a], ctx->d_state[9 + a], ctx->d_state[12 + a],
            dt, a, ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1);
        if (ctx->n_l1 > 0)
        {
            c2_gpu::coef_update_g1_kernel<<<l1g, 256>>>(
                rho_new_p[a] + ctx->l0_s_g0[a], ctx->d_L_b[a] + ctx->l0_s_g0[a],
                ctx->d_state[9 + a] + ctx->l0_s_g1[a], ctx->d_state[12 + a] + ctx->l0_s_g1[a],
                dt, a, ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1);
            double* c_l1 = ctx->d_state[9 + a] + ctx->l0_s_g1[a];
            double* d_l1 = ctx->d_state[12 + a] + ctx->l0_s_g1[a];
            c2_gpu::sync_l1_intralevel_side_kernel<<<l1g, 256>>>(c_l1, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
            c2_gpu::sync_l1_intralevel_side_kernel<<<l1g, 256>>>(d_l1, ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
            c2_gpu::restrict_l1_to_l0_side_kernel<<<l1g, 256>>>(c_l1, ctx->d_state[9 + a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
            c2_gpu::restrict_l1_to_l0_side_kernel<<<l1g, 256>>>(d_l1, ctx->d_state[12 + a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a, 1);
            c2_gpu::mirror_g1_side_l1_kernel<<<l1g, 256>>>(c_l1, ctx->d_l1_boxes, ctx->n_l1, a);
            c2_gpu::mirror_g1_side_l1_kernel<<<l1g, 256>>>(d_l1, ctx->d_l1_boxes, ctx->n_l1, a);
        }
        c2_gpu::mirror_g1_side_l0_kernel<<<32, 256>>>(ctx->d_state[9 + a], a, ctx->nx0, ctx->ny0, ctx->nz0);
        c2_gpu::mirror_g1_side_l0_kernel<<<32, 256>>>(ctx->d_state[12 + a], a, ctx->nx0, ctx->ny0, ctx->nz0);
    }
    int rc = c2_gpu::refresh_mg_coefs_from_state(ctx);
    if (rc != C2_OK) { cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }

    // 2. Build RHS: urhs_p = viscous_rhs - N_full + F (+ sync + modifyRhsForBcs)
    for (int a = 0; a < 3; ++a)
        c2_gpu::linear_sum_side_kernel<<<128, 256>>>(urhs_p[a], 1.0, urhs_p[a], -1.0, n_full_p[a], ctx->tot_s_g0(a));

    if (has_f != 0 && f_l0 != nullptr)
    {
        cudaEvent_t fh0{}, fh1{};
        cudaEventCreate(&fh0); cudaEventCreate(&fh1); cudaEventRecord(fh0);
        for (int a = 0; a < 3; ++a)
        {
            if (f_l0[a])
            {
                rc = c2_gpu::upload_side_array(ctx, ctx->d_krylov_z[1][a], f_l0[a], f_l1 ? f_l1[a] : nullptr, a);
                if (rc != C2_OK) { cudaEventDestroy(fh0); cudaEventDestroy(fh1); cudaEventDestroy(ev0); cudaEventDestroy(ev1); return rc; }
                c2_gpu::linear_sum_side_kernel<<<128, 256>>>(urhs_p[a], 1.0, urhs_p[a], 1.0, ctx->d_krylov_z[1][a], ctx->tot_s_g0(a));
            }
        }
        cudaEventRecord(fh1); cudaEventSynchronize(fh1);
        float h2d_ms = 0.0f; cudaEventElapsedTime(&h2d_ms, fh0, fh1);
        ctx->tel.ms_h2d += static_cast<double>(h2d_ms);
        cudaEventDestroy(fh0); cudaEventDestroy(fh1);
    }

    if (ctx->n_l1 > 0)
    {
        for (int a = 0; a < 3; ++a)
        {
            c2_gpu::sync_l1_shared_faces_g0_kernel<<<l1g, 256>>>(urhs_p[a] + ctx->l0_s_g0[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
            c2_gpu::restrict_g0_l1_to_l0_side_kernel<<<l1g, 256>>>(urhs_p[a] + ctx->l0_s_g0[a], urhs_p[a], ctx->d_l1_boxes, ctx->n_l1, ctx->d_l1_map, a);
        }
    }
    c2_gpu::modify_rhs_bcs_gpu_kernel<<<dim3(64, 1), 256>>>(
        urhs_p[0], urhs_p[1], urhs_p[2], ctx->d_state[0],
        p_in_new, ctx->mu0, ctx->dx0, ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1);
    if (ctx->n_l1 > 0)
        c2_gpu::modify_rhs_bcs_gpu_kernel<<<l1g, 256>>>(
            urhs_p[0] + ctx->l0_s_g0[0], urhs_p[1] + ctx->l0_s_g0[1], urhs_p[2] + ctx->l0_s_g0[2],
            ctx->d_state[0] + ctx->l0_s_g0[0], p_in_new, ctx->mu0, ctx->dx1,
            ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1);
    C2_CUDA_CALL(ctx, cudaMemset(ctx->d_state[23], 0, ctx->tot_c_g0() * sizeof(double)), "rhsp_zero");

    cudaEventRecord(ev1); cudaEventSynchronize(ev1);
    float other_ms = 0.0f; cudaEventElapsedTime(&other_ms, ev0, ev1);
    ctx->tel.ms_other = static_cast<double>(other_ms);
    cudaEventDestroy(ev0); cudaEventDestroy(ev1);

    // 3. FGMRES(5) solve with warm start from (U_cur, P_cur) in d_state[6..8, 22]
    double* rhs_solve[4] = {urhs_p[0], urhs_p[1], urhs_p[2], ctx->d_state[23]};
    double* sol_gpu[4] = {ctx->d_state[6], ctx->d_state[7], ctx->d_state[8], ctx->d_state[22]};
    c2_gpu::impose_dirichlet_zero_g0_kernel<<<dim3(32, 1), 256>>>(
        sol_gpu[1], sol_gpu[2], ctx->nx0, ctx->ny0, ctx->nz0, false, nullptr, 1);
    if (ctx->n_l1 > 0)
        c2_gpu::impose_dirichlet_zero_g0_kernel<<<l1g, 256>>>(
            sol_gpu[1] + ctx->l0_s_g0[1], sol_gpu[2] + ctx->l0_s_g0[2],
            ctx->nx0, ctx->ny0, ctx->nz0, true, ctx->d_l1_boxes, ctx->n_l1);

    bool use_fp32 = (ctx->cfg.enable_fp32_mg != 0) && (ctx->sticky_fp64_active == 0);
    c2_gpu::SolveReport rep{};
    rc = c2_gpu::solve_fgmres_m5(
        ctx, rhs_solve, sol_gpu, dt,
        ctx->cfg.rel_tol, ctx->cfg.abs_tol, ctx->cfg.max_iterations,
        false, use_fp32, true, rep);

    ctx->tel.last_iters = rep.iters_total;
    ctx->tel.last_fp32_iters = rep.fp32_iters;
    ctx->tel.last_fp64_fb_iters = rep.fp64_fb_iters;
    ctx->tel.last_fallback = rep.fallback_triggered ? 1 : 0;
    ctx->tel.last_true_rel_res = rep.final_true_rel_res;
    ctx->tel.last_true_abs_res = rep.final_true_abs_res;
    ctx->tel.last_rhs_norm = rep.rhs_norm;
    ctx->tel.ms_stokes = rep.wall_ms;
    c2_gpu::set_fb_reason(ctx, rep.fallback_reason.c_str());
    if (ctx->cfg.sticky_fp64_on_fb && rep.fallback_triggered)
        ctx->sticky_fp64_active = 1;
    ctx->tel.sticky_fp64_active = ctx->sticky_fp64_active;

    if (rc != C2_OK)
    {
        c2_gpu::sample_vram(ctx, -1);
        return rc;
    }

    // Finiteness check across all outputs
    for (int a = 0; a < 3; ++a)
    {
        rc = c2_gpu::verify_device_finite(ctx, sol_gpu[a], ctx->tot_s_g0(a));
        if (rc != C2_OK) return rc;
        rc = c2_gpu::verify_device_finite(ctx, rho_new_p[a], ctx->tot_s_g0(a));
        if (rc != C2_OK) return rc;
    }
    rc = c2_gpu::verify_device_finite(ctx, sol_gpu[3], ctx->tot_c_g0());
    if (rc != C2_OK) return rc;

    if (out_iters) *out_iters = rep.iters_total;
    if (out_rel_res) *out_rel_res = rep.final_true_rel_res;
    if (out_converged) *out_converged = rep.converged ? 1 : 0;

    if (rep.converged)
    {
        cudaEvent_t de0{}, de1{};
        cudaEventCreate(&de0); cudaEventCreate(&de1);
        cudaEventRecord(de0);
        for (int a = 0; a < 3; ++a)
        {
            rc = c2_gpu::download_side_array(ctx, sol_gpu[a], u_new_l0[a], u_new_l1 ? u_new_l1[a] : nullptr, a);
            if (rc != C2_OK) { cudaEventDestroy(de0); cudaEventDestroy(de1); return rc; }
            if (download_rho && rho_new_l0 && rho_new_l0[a])
            {
                rc = c2_gpu::download_side_array(ctx, rho_new_p[a], rho_new_l0[a], rho_new_l1 ? rho_new_l1[a] : nullptr, a);
                if (rc != C2_OK) { cudaEventDestroy(de0); cudaEventDestroy(de1); return rc; }
            }
        }
        rc = c2_gpu::download_cell_array(ctx, sol_gpu[3], p_new_l0, p_new_l1);
        if (rc != C2_OK) { cudaEventDestroy(de0); cudaEventDestroy(de1); return rc; }

        // Advance persistent GPU state (U_old <- U_cur <- U_new, P_cur <- P_new)
        for (int a = 0; a < 3; ++a)
        {
            C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[3 + a], ctx->d_state[a], ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyDeviceToDevice), "adv_uold");
            C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[a], sol_gpu[a], ctx->tot_s_g0(a) * sizeof(double), cudaMemcpyDeviceToDevice), "adv_ucur");
        }
        C2_CUDA_CALL(ctx, cudaMemcpy(ctx->d_state[21], sol_gpu[3], ctx->tot_c_g0() * sizeof(double), cudaMemcpyDeviceToDevice), "adv_pcur");

        cudaEventRecord(de1); cudaEventSynchronize(de1);
        float d2h_ms = 0.0f; cudaEventElapsedTime(&d2h_ms, de0, de1);
        ctx->tel.ms_d2h += static_cast<double>(d2h_ms);
        cudaEventDestroy(de0); cudaEventDestroy(de1);
    }
    else
    {
        c2_gpu::set_fb_reason(ctx, "stokes_not_converged");
    }

    c2_gpu::sample_vram(ctx, -1);
    return C2_OK;
}

void c2_gpu_get_telemetry(void* handle, C2Telemetry* out_tel)
{
    if (!out_tel) return;
    std::memset(out_tel, 0, sizeof(*out_tel));
    auto* ctx = static_cast<c2_gpu::C2Context*>(handle);
    if (!ctx) return;
    *out_tel = ctx->tel;
}

} // extern "C"
