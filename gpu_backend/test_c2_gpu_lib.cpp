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

#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace
{

static size_t sidx0(int a, int i, int j, int k, int nx, int ny)
{
    int sx = nx + (a == 0);
    int sy = ny + (a == 1);
    return static_cast<size_t>(i) + static_cast<size_t>(sx) * (j + static_cast<size_t>(sy) * k);
}

static double eval_brinkman_L_b(double x, double y, double penalty, double dt, double h_min)
{
    const double pad = 10.0 * 112.0;
    const double y_lo_ext = 0.0 - pad, y_hi_ext = 20.0;
    const double cx = 0.5 * (50.0 + 70.0), cy = 0.5 * (y_lo_ext + y_hi_ext);
    const double hx = 0.5 * (70.0 - 50.0), hy = 0.5 * (y_hi_ext - y_lo_ext);
    const double r_c = 1.0, hx_core = hx - r_c, hy_core = hy - r_c;
    double dx = std::fabs(x - cx) - hx_core, dy = std::fabs(y - cy) - hy_core;
    double ox = std::max(dx, 0.0), oy = std::max(dy, 0.0);
    double phi = std::hypot(ox, oy) + std::min(std::max(dx, dy), 0.0) - r_c;
    if (phi <= 1.0e-12 * h_min)
        return penalty * (1.0 / dt + 1.0 / (h_min * h_min));
    return 0.0;
}

struct Api
{
    void* lib = nullptr;
    decltype(&c2_gpu_create) create = nullptr;
    decltype(&c2_gpu_destroy) destroy = nullptr;
    decltype(&c2_gpu_regrid) regrid = nullptr;
    decltype(&c2_gpu_upload_state) upload_state = nullptr;
    decltype(&c2_gpu_upload_brinkman) upload_brinkman = nullptr;
    decltype(&c2_gpu_preprocess) preprocess = nullptr;
    decltype(&c2_gpu_integrate_solve) integrate_solve = nullptr;
    decltype(&c2_gpu_get_telemetry) get_telemetry = nullptr;
};

static Api load_api(const char* so_path)
{
    Api a;
    a.lib = dlopen(so_path, RTLD_NOW);
    if (!a.lib) { std::fprintf(stderr, "dlopen(%s) failed: %s\n", so_path, dlerror()); std::exit(EXIT_FAILURE); }
    a.create = reinterpret_cast<decltype(&c2_gpu_create)>(dlsym(a.lib, "c2_gpu_create"));
    a.destroy = reinterpret_cast<decltype(&c2_gpu_destroy)>(dlsym(a.lib, "c2_gpu_destroy"));
    a.regrid = reinterpret_cast<decltype(&c2_gpu_regrid)>(dlsym(a.lib, "c2_gpu_regrid"));
    a.upload_state = reinterpret_cast<decltype(&c2_gpu_upload_state)>(dlsym(a.lib, "c2_gpu_upload_state"));
    a.upload_brinkman = reinterpret_cast<decltype(&c2_gpu_upload_brinkman)>(dlsym(a.lib, "c2_gpu_upload_brinkman"));
    a.preprocess = reinterpret_cast<decltype(&c2_gpu_preprocess)>(dlsym(a.lib, "c2_gpu_preprocess"));
    a.integrate_solve = reinterpret_cast<decltype(&c2_gpu_integrate_solve)>(dlsym(a.lib, "c2_gpu_integrate_solve"));
    a.get_telemetry = reinterpret_cast<decltype(&c2_gpu_get_telemetry)>(dlsym(a.lib, "c2_gpu_get_telemetry"));
    if (!a.create || !a.destroy || !a.regrid || !a.upload_state || !a.upload_brinkman || !a.preprocess || !a.integrate_solve || !a.get_telemetry)
    {
        std::fprintf(stderr, "Missing symbols in %s\n", so_path);
        std::exit(EXIT_FAILURE);
    }
    return a;
}

static bool test_synthetic_channel_amr(const Api& api, int n_l1_case)
{
    const int nx0 = 224, ny0 = 50, nz0 = 25;
    const double dx0 = 0.5, dt = 0.00625;
    C2GridConfig cfg{};
    cfg.nx0 = nx0; cfg.ny0 = ny0; cfg.nz0 = nz0;
    cfg.dx0 = dx0; cfg.mu0 = 1.0; cfg.rho0 = 1.0;
    cfg.enable_fp32_mg = 1; cfg.sticky_fp64_on_fb = 1;
    cfg.max_iterations = 100; cfg.rel_tol = 1.0e-5; cfg.abs_tol = 1.0e-50;
    int err = 0;
    void* h = api.create(&cfg, &err);
    if (!h || err != C2_OK) return false;

    // Verify oversized L1 patch count cleanly returns C2_ERR_CAPACITY
    std::vector<C2PatchBox> too_many(600);
    if (api.regrid(h, 600, too_many.data()) != C2_ERR_CAPACITY)
    {
        std::fprintf(stderr, "[ChannelAMR] expected C2_ERR_CAPACITY for 600 patches\n");
        api.destroy(h); return false;
    }

    std::vector<C2PatchBox> boxes;
    if (n_l1_case == 4)
    {
        boxes = {
            {{184, 72, 0},  {243, 85, 49}},
            {{244, 72, 0},  {303, 85, 49}},
            {{184, 86, 0},  {243, 99, 49}},
            {{244, 86, 0},  {303, 99, 49}},
        };
    }
    if (api.regrid(h, static_cast<int>(boxes.size()), boxes.empty() ? nullptr : boxes.data()) != C2_OK)
    { api.destroy(h); return false; }

    size_t n_c0 = static_cast<size_t>(nx0) * ny0 * nz0;
    size_t l1_c0 = 0, l1_s0[3] = {0, 0, 0};
    for (const auto& b : boxes)
    {
        int bx = b.hi[0] - b.lo[0] + 1, by = b.hi[1] - b.lo[1] + 1, bz = b.hi[2] - b.lo[2] + 1;
        l1_c0 += static_cast<size_t>(bx) * by * bz;
        for (int a = 0; a < 3; ++a)
            l1_s0[a] += static_cast<size_t>(bx + (a == 0)) * (by + (a == 1)) * (bz + (a == 2));
    }

    std::array<std::vector<double>, 3> u0, u1, rho0, rho1, L0, L1, f0, f1, un0, un1, rn0, rn1;
    std::vector<double> p0(n_c0, 0.0), p1(l1_c0, 0.0), pn0(n_c0, 0.0), pn1(l1_c0, 0.0);
    for (int a = 0; a < 3; ++a)
    {
        int nxa = nx0 + (a == 0), nya = ny0 + (a == 1), nza = nz0 + (a == 2);
        size_t ns0 = static_cast<size_t>(nxa) * nya * nza;
        u0[a].assign(ns0, 0.0); rho0[a].assign(ns0, 1.0); L0[a].assign(ns0, 0.0);
        f0[a].assign(ns0, (a == 0) ? 0.01 : 0.0);
        u1[a].assign(l1_s0[a], 0.0); rho1[a].assign(l1_s0[a], 1.0); L1[a].assign(l1_s0[a], 0.0);
        f1[a].assign(l1_s0[a], (a == 0) ? 0.01 : 0.0);
        un0[a].assign(ns0, 0.0); rn0[a].assign(ns0, 0.0);
        un1[a].assign(l1_s0[a], 0.0); rn1[a].assign(l1_s0[a], 0.0);
        for (int k = 0; k < nza; ++k)
            for (int j = 0; j < nya; ++j)
                for (int i = 0; i < nxa; ++i)
                    L0[a][sidx0(a, i, j, k, nx0, ny0)] = eval_brinkman_L_b(dx0 * (i + (a == 0 ? 0.0 : 0.5)), dx0 * (j + (a == 1 ? 0.0 : 0.5)), 1.0e4, dt, dx0);
    }
    const double* up0[3] = {u0[0].data(), u0[1].data(), u0[2].data()};
    const double* up1[3] = {u1[0].data(), u1[1].data(), u1[2].data()};
    const double* rp0[3] = {rho0[0].data(), rho0[1].data(), rho0[2].data()};
    const double* rp1[3] = {rho1[0].data(), rho1[1].data(), rho1[2].data()};
    const double* Lp0[3] = {L0[0].data(), L0[1].data(), L0[2].data()};
    const double* Lp1[3] = {L1[0].data(), L1[1].data(), L1[2].data()};
    const double* fp0[3] = {f0[0].data(), f0[1].data(), f0[2].data()};
    const double* fp1[3] = {f1[0].data(), f1[1].data(), f1[2].data()};
    double* unp0[3] = {un0[0].data(), un0[1].data(), un0[2].data()};
    double* unp1[3] = {un1[0].data(), un1[1].data(), un1[2].data()};
    double* rnp0[3] = {rn0[0].data(), rn0[1].data(), rn0[2].data()};
    double* rnp1[3] = {rn1[0].data(), rn1[1].data(), rn1[2].data()};

    if (api.upload_state(h, up0, up1, nullptr, nullptr, p0.data(), p1.data(), rp0, rp1) != C2_OK ||
        api.upload_brinkman(h, Lp0, Lp1) != C2_OK ||
        api.preprocess(h, dt, -1.0, 1, 1.0) != C2_OK)
    { api.destroy(h); return false; }

    int c1 = 0, it1 = 0, c2 = 0, it2 = 0;
    double res1 = 0.0, res2 = 0.0;
    if (api.integrate_solve(h, dt, 0.5, fp0, fp1, 1, unp0, unp1, pn0.data(), pn1.data(), rnp0, rnp1, 1, &c1, &it1, &res1) != C2_OK)
    { api.destroy(h); return false; }

    if (api.preprocess(h, dt, dt, 0, 1.0) != C2_OK ||
        api.integrate_solve(h, dt, 1.5, fp0, fp1, 1, unp0, unp1, pn0.data(), pn1.data(), rnp0, rnp1, 0, &c2, &it2, &res2) != C2_OK)
    { api.destroy(h); return false; }

    double sym_u2 = 0.0;
    for (int j = 0; j < ny0; ++j)
        for (int i = 0; i < nx0; ++i)
            sym_u2 = std::max(sym_u2, std::fabs(un0[2][sidx0(2, i, j, 0, nx0, ny0)]));
    for (size_t p = 0; p < boxes.size(); ++p)
    {
        if (boxes[p].lo[2] == 0)
        {
            int bx = boxes[p].hi[0] - boxes[p].lo[0] + 1, by = boxes[p].hi[1] - boxes[p].lo[1] + 1;
            size_t off = 0;
            for (size_t q = 0; q < p; ++q)
                off += static_cast<size_t>(boxes[q].hi[0] - boxes[q].lo[0] + 1) * (boxes[q].hi[1] - boxes[q].lo[1] + 1) * (boxes[q].hi[2] - boxes[q].lo[2] + 2);
            for (int j = 0; j < by; ++j)
                for (int i = 0; i < bx; ++i)
                    sym_u2 = std::max(sym_u2, std::fabs(un1[2][off + sidx0(2, i, j, 0, bx, by)]));
        }
    }

    C2Telemetry tel{};
    api.get_telemetry(h, &tel);
    api.destroy(h);

    const bool tel_ok = (tel.regrid_count == 2) && (tel.last_n_l1 == n_l1_case) && (tel.peak_n_l1 == 600) &&
                        (tel.max_l1_patches == 512) && (tel.last_tot_c_g0 > 0) && (tel.last_tot_c_g1 >= tel.last_tot_c_g0) &&
                        (tel.last_max_s_g0 == tel.peak_max_s_g0) && (tel.last_max_s_g0 <= tel.max_slot_g0) &&
                        (tel.last_max_s_g1 == tel.peak_max_s_g1) && (tel.last_max_s_g1 <= tel.max_slot_g1) &&
                        (tel.last_max_s_g3 == tel.peak_max_s_g3) && (tel.last_max_s_g3 <= tel.max_slot_g0);

    std::printf("[ChannelAMR n_l1=%d] step1: conv=%d it=%d rel=%.3e | step2: conv=%d it=%d rel=%.3e (fp32_it=%d fp64_fb_it=%d fb=%d) | sym_u2=%.3e | NVML peak=%ld MiB | regrid=%d max_s_g0=%llu/%llu\n",
                n_l1_case, c1, it1, res1, c2, it2, res2, tel.last_fp32_iters, tel.last_fp64_fb_iters, tel.last_fallback,
                sym_u2, tel.peak_nvml_mib, tel.regrid_count, tel.peak_max_s_g0, tel.max_slot_g0);
    return c1 && c2 && tel_ok && (res1 <= 1.0e-5) && (res2 <= 1.0e-5) && (sym_u2 <= 1.0e-12) &&
           (tel.peak_nvml_mib > 0) && (tel.peak_nvml_mib <= 3584);
}

static bool test_capacity_telemetry(const Api& api)
{
    C2GridConfig cfg{224, 50, 25, 0.5, 1.0, 1.0, 1, 0, 100, 1.0e-5, 1.0e-50};
    int err = 0;
    void* h = api.create(&cfg, &err);
    if (!h || err != C2_OK) return false;

    std::vector<C2PatchBox> too_many(513, C2PatchBox{{0, 0, 0}, {1, 1, 1}});
    int rc_patches = api.regrid(h, static_cast<int>(too_many.size()), too_many.data());
    C2Telemetry tel_patches{};
    api.get_telemetry(h, &tel_patches);

    // Trigger capacity_slot_overflow on 224x50x25 channel geometry (max_slot_g0=1,532,160) with full-domain L1 (2,520,000 tot_c_g0)
    C2PatchBox big_box{{0, 0, 0}, {447, 99, 49}};
    int rc_slots = api.regrid(h, 1, &big_box);
    C2Telemetry tel_slots{};
    api.get_telemetry(h, &tel_slots);
    api.destroy(h);

    const bool ok = (rc_patches == C2_ERR_CAPACITY) && (tel_patches.regrid_count == 1) &&
                    (tel_patches.last_n_l1 == 513) && (tel_patches.peak_n_l1 == 513) &&
                    (std::string(tel_patches.last_fb_reason) == "capacity_n_l1_patches") &&
                    (rc_slots == C2_ERR_CAPACITY) && (tel_slots.regrid_count == 2) &&
                    (tel_slots.last_n_l1 == 1) && (tel_slots.peak_n_l1 == 513) &&
                    (tel_slots.last_tot_c_g0 > tel_slots.max_slot_g0) &&
                    (tel_slots.last_max_s_g0 > tel_slots.max_slot_g0) &&
                    (tel_slots.peak_max_s_g0 == tel_slots.last_max_s_g0) &&
                    (tel_slots.last_max_s_g1 == tel_slots.peak_max_s_g1) &&
                    (std::string(tel_slots.last_fb_reason) == "capacity_slot_overflow");
    std::printf("[CapacityTelemetry] rc_patches=%d n_l1=%d | rc_slots=%d tot_c_g0=%llu/%llu max_s_g0=%llu max_s_g1=%llu reason=%s\n",
                rc_patches, tel_patches.last_n_l1, rc_slots, tel_slots.last_tot_c_g0, tel_slots.max_slot_g0,
                tel_slots.last_max_s_g0, tel_slots.last_max_s_g1, tel_slots.last_fb_reason);
    return ok;
}

} // namespace

int main()
{
    Api api = load_api("gpu_backend/libibamr_c2_gpu.so");
    bool ok = true;
    ok = test_synthetic_channel_amr(api, 0) && ok;
    ok = test_synthetic_channel_amr(api, 4) && ok;
    ok = test_capacity_telemetry(api) && ok;
    std::printf("%s\n", ok ? "GPU_BACKEND_LIB_TEST_OK" : "GPU_BACKEND_LIB_TEST_FAIL");
    return ok ? 0 : 1;
}
