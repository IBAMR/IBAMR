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

#ifndef IBAMR_C2_GPU_H
#define IBAMR_C2_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

enum C2GpuStatus {
    C2_OK = 0,
    C2_ERR_INVALID = 1,
    C2_ERR_CAPACITY = 2,
    C2_ERR_CUDA = 3,
    C2_ERR_NONFINITE = 4
};

typedef struct {
    int nx0;
    int ny0;
    int nz0;
    double dx0;
    double mu0;
    double rho0;
    int enable_fp32_mg;
    int sticky_fp64_on_fb;
    int max_iterations;
    double rel_tol;
    double abs_tol;
} C2GridConfig;

typedef struct {
    int lo[3];
    int hi[3];
} C2PatchBox;

typedef struct {
    long stage_nvml_mib[5];
    long peak_nvml_mib;
    long last_nvml_mib;
    double min_dev_free_mib;
    int last_iters;
    int last_fp32_iters;
    int last_fp64_fb_iters;
    int last_fallback;
    double last_true_rel_res;
    char last_fb_reason[64];
    double ms_regrid;
    double ms_preprocess;
    double ms_ssprk;
    double ms_other;
    double ms_stokes;
    double ms_h2d;
    double ms_d2h;
    int sticky_fp64_active;
    double last_true_abs_res;
    double last_rhs_norm;
    int regrid_count;
    int last_n_l1;
    int max_l1_patches;
    int peak_n_l1;
    unsigned long long max_slot_g0;
    unsigned long long max_slot_g1;
    unsigned long long last_tot_c_g0;
    unsigned long long last_tot_c_g1;
    unsigned long long last_tot_s_g0[3];
    unsigned long long last_tot_s_g1[3];
    unsigned long long last_l0_s_g3[3];
    unsigned long long last_l1_s_g3[3];
    unsigned long long peak_tot_c_g0;
    unsigned long long peak_tot_c_g1;
    unsigned long long peak_tot_s_g0[3];
    unsigned long long peak_tot_s_g1[3];
    unsigned long long peak_l0_s_g3[3];
    unsigned long long peak_l1_s_g3[3];
    unsigned long long last_max_s_g0;
    unsigned long long last_max_s_g1;
    unsigned long long last_max_s_g3;
    unsigned long long peak_max_s_g0;
    unsigned long long peak_max_s_g1;
    unsigned long long peak_max_s_g3;
} C2Telemetry;

void* c2_gpu_create(const C2GridConfig* cfg, int* out_err);
void c2_gpu_destroy(void* handle);

int c2_gpu_regrid(void* handle, int n_l1, const C2PatchBox* l1_boxes);

int c2_gpu_upload_state(
    void* handle,
    const double* const u_cur_l0[3],
    const double* const u_cur_l1[3],
    const double* const u_old_l0[3],
    const double* const u_old_l1[3],
    const double* p_cur_l0,
    const double* p_cur_l1,
    const double* const rho_cur_l0[3],
    const double* const rho_cur_l1[3]);

int c2_gpu_upload_brinkman(
    void* handle,
    const double* const L_l0[3],
    const double* const L_l1[3]);

int c2_gpu_upload_rho(
    void* handle,
    const double* const rho_l0[3],
    const double* const rho_l1[3]);

int c2_gpu_preprocess(
    void* handle,
    double dt,
    double dt_prev,
    int is_first_step,
    double constant_rho);

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
    double* out_rel_res);

void c2_gpu_get_telemetry(void* handle, C2Telemetry* out_tel);

#ifdef __cplusplus
}
#endif

#endif
