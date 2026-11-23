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

#ifndef included_IBAMR_INSVCStaggeredConservativeGPUBridge
#define included_IBAMR_INSVCStaggeredConservativeGPUBridge

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#if (NDIM == 3)

#include "ibamr_c2_gpu.h"

#include <ibamr/StaggeredStokesSolver.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <Variable.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/MathUtilities.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <mpi.h>
#include <petscksp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <dlfcn.h>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace IBAMR
{
namespace detail
{
inline void
c2_dladdr_anchor()
{
}

inline double*
nonnull_data(std::vector<double>& v)
{
    static double dummy = 0.0;
    return v.empty() ? &dummy : v.data();
}

inline const double*
nonnull_cdata(const std::vector<double>& v)
{
    static const double dummy = 0.0;
    return v.empty() ? &dummy : v.data();
}

inline int*
nonnull_idata(std::vector<int>& v)
{
    static int dummy = 0;
    return v.empty() ? &dummy : v.data();
}
} // namespace detail

/*!
 * \brief Rank-0 dynamic loader and MPI gather/scatter bridge for the opt-in CUDA fluid solver backend.
 */
class INSVCStaggeredConservativeGPUBridge
{
public:
    struct Config
    {
        bool enabled = false;
        std::string backend = "candidate2";
        std::string library_path;
        bool enable_fp32_mg = true;
        bool sticky_fp64_on_fb = false;
        bool fallback_to_cpu = true;
        int force_cpu_fallback_step = -1;
        bool enable_logging = true;
    };

    INSVCStaggeredConservativeGPUBridge(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                        std::string integrator_name)
        : d_integrator_name(std::move(integrator_name))
    {
        if (input_db)
        {
            d_cfg.enabled = input_db->getBoolWithDefault("use_gpu_fluid_backend", false);
            if (input_db->keyExists("gpu_fluid_backend_db"))
            {
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db = input_db->getDatabase("gpu_fluid_backend_db");
                if (db)
                {
                    d_cfg.backend = db->getStringWithDefault("backend", d_cfg.backend);
                    d_cfg.library_path = db->getStringWithDefault("library_path", d_cfg.library_path);
                    d_cfg.enable_fp32_mg = db->getBoolWithDefault("enable_fp32_mg", d_cfg.enable_fp32_mg);
                    d_cfg.sticky_fp64_on_fb = db->getBoolWithDefault("sticky_fp64_on_fb", d_cfg.sticky_fp64_on_fb);
                    d_cfg.fallback_to_cpu = db->getBoolWithDefault("fallback_to_cpu", d_cfg.fallback_to_cpu);
                    d_cfg.force_cpu_fallback_step =
                        db->getIntegerWithDefault("force_cpu_fallback_step", d_cfg.force_cpu_fallback_step);
                    d_cfg.enable_logging = db->getBoolWithDefault("enable_logging", d_cfg.enable_logging);
                }
            }
        }
    }

    ~INSVCStaggeredConservativeGPUBridge()
    {
        if (IBTK::IBTK_MPI::getRank() == 0)
        {
            destroyHandleRank0();
            if (d_dl_handle)
            {
                dlclose(d_dl_handle);
                d_dl_handle = nullptr;
            }
        }
    }

    const Config& config() const
    {
        return d_cfg;
    }

    bool fallbackToCpu() const
    {
        return d_cfg.fallback_to_cpu;
    }

    int forceCpuFallbackStep() const
    {
        return d_cfg.force_cpu_fallback_step;
    }

    void invalidateBrinkmanCache()
    {
        d_cached_L_valid = false;
    }

    static bool extractEffectiveSolverSettings(
        const SAMRAI::tbox::Pointer<StaggeredStokesSolver>& stokes_solver,
        const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& stokes_solver_db,
        double& rel_tol,
        double& abs_tol,
        int& max_iters,
        std::string& ksp_type_out)
    {
        rel_tol = 1.0e-5;
        abs_tol = 1.0e-50;
        max_iters = 100;
        ksp_type_out = "fgmres";

        if (stokes_solver)
        {
            rel_tol = stokes_solver->getRelativeTolerance();
            abs_tol = stokes_solver->getAbsoluteTolerance();
            max_iters = stokes_solver->getMaxIterations();
        }
        else if (stokes_solver_db)
        {
            rel_tol = stokes_solver_db->getDoubleWithDefault("rel_residual_tol", rel_tol);
            abs_tol = stokes_solver_db->getDoubleWithDefault("abs_residual_tol", abs_tol);
            max_iters = stokes_solver_db->getIntegerWithDefault("max_iterations", max_iters);
        }

        if (stokes_solver_db && stokes_solver_db->keyExists("ksp_type"))
        {
            ksp_type_out = stokes_solver_db->getString("ksp_type");
        }

        std::string pfx = "stokes_";
        if (stokes_solver_db && stokes_solver_db->keyExists("options_prefix"))
        {
            pfx = stokes_solver_db->getString("options_prefix");
        }
        if (stokes_solver)
        {
            if (auto* petsc_solver = dynamic_cast<const IBTK::PETScKrylovLinearSolver*>(stokes_solver.getPointer()))
            {
                if (petsc_solver->getPETScKSP() != nullptr)
                {
                    const char* ksp_pfx = nullptr;
                    KSPGetOptionsPrefix(petsc_solver->getPETScKSP(), &ksp_pfx);
                    if (ksp_pfx) pfx = ksp_pfx;
                }
            }
        }

        PetscReal opt_rtol = 0.0, opt_atol = 0.0;
        PetscInt opt_maxit = 0;
        PetscBool flg = PETSC_FALSE;
        PetscOptionsGetReal(nullptr, pfx.c_str(), "-ksp_rtol", &opt_rtol, &flg);
        if (flg) rel_tol = static_cast<double>(opt_rtol);
        flg = PETSC_FALSE;
        PetscOptionsGetReal(nullptr, pfx.c_str(), "-ksp_atol", &opt_atol, &flg);
        if (flg) abs_tol = static_cast<double>(opt_atol);
        flg = PETSC_FALSE;
        PetscOptionsGetInt(nullptr, pfx.c_str(), "-ksp_max_it", &opt_maxit, &flg);
        if (flg) max_iters = static_cast<int>(opt_maxit);

        char ksp_buf[64] = { 0 };
        flg = PETSC_FALSE;
        PetscOptionsGetString(nullptr, pfx.c_str(), "-ksp_type", ksp_buf, sizeof(ksp_buf), &flg);
        if (flg && ksp_buf[0] != '\0')
        {
            ksp_type_out = ksp_buf;
        }
        return (ksp_type_out == "fgmres" || ksp_type_out == "gmres");
    }

    static bool inspectPressureDrivenHalfZBcs(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& u_var,
        const int u_scratch_idx,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        const double eval_time,
        double& p_in_new,
        std::string& reason)
    {
        if (!hierarchy || bc_coefs.size() != NDIM || !bc_coefs[0] || !bc_coefs[1] || !bc_coefs[2])
        {
            reason = "missing_velocity_bc_coefs";
            return false;
        }

        int local_ok = 1;
        std::string local_reason;
        int local_saw_inlet = 0;
        double local_p_in_min = std::numeric_limits<double>::infinity();
        double local_p_in_max = -std::numeric_limits<double>::infinity();
        const double tol = 1.0e-12;

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l0 = hierarchy->getPatchLevel(0);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(l0); p && local_ok; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = l0->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM>> codim1_boxes =
                IBTK::PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            for (int n = 0; n < codim1_boxes.size() && local_ok; ++n)
            {
                const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box = codim1_boxes[n];
                const unsigned int loc = bdry_box.getLocationIndex();
                if (loc >= 2 * NDIM) continue;
                const SAMRAI::hier::Box<NDIM> bc_fill_box =
                    pgeom->getBoundaryFillBox(bdry_box, patch_box, SAMRAI::hier::IntVector<NDIM>(1));
                const SAMRAI::hier::BoundaryBox<NDIM> trimmed_bdry_box(
                    bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), loc);
                const SAMRAI::hier::Box<NDIM> bc_coef_box =
                    IBTK::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                if (bc_coef_box.empty()) continue;

                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>> acoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM, double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>> bcoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM, double>(bc_coef_box, 1);
                SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>> gcoef_data =
                    new SAMRAI::pdat::ArrayData<NDIM, double>(bc_coef_box, 1);

                for (unsigned int comp = 0; comp < NDIM && local_ok; ++comp)
                {
                    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc = bc_coefs[comp];
                    auto* ext_bc = dynamic_cast<IBTK::ExtendedRobinBcCoefStrategy*>(bc);
                    if (ext_bc)
                    {
                        ext_bc->setTargetPatchDataIndex(u_scratch_idx);
                        ext_bc->setHomogeneousBc(false);
                    }
                    acoef_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                    bcoef_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                    gcoef_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                    bc->setBcCoefs(acoef_data, bcoef_data, gcoef_data, u_var, *patch, trimmed_bdry_box, eval_time);
                    if (ext_bc) ext_bc->clearTargetPatchDataIndex();

                    for (SAMRAI::hier::Box<NDIM>::Iterator it(bc_coef_box); it; it++)
                    {
                        const SAMRAI::hier::Index<NDIM>& idx = it();
                        const double a = (*acoef_data)(idx, 0);
                        const double b = (*bcoef_data)(idx, 0);
                        const double g = (*gcoef_data)(idx, 0);
                        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(g))
                        {
                            local_ok = 0;
                            local_reason = "nonfinite_bc_coef";
                            break;
                        }

                        if (loc == 0 && comp == 0)
                        {
                            // Inlet normal traction: sigma.n = -p_in -> g = -p_in
                            if (std::abs(a) > tol || std::abs(b - 1.0) > tol)
                            {
                                local_ok = 0;
                                local_reason = "bc_inlet_not_traction";
                                break;
                            }
                            const double p_val = -g;
                            local_saw_inlet = 1;
                            local_p_in_min = std::min(local_p_in_min, p_val);
                            local_p_in_max = std::max(local_p_in_max, p_val);
                        }
                        else if ((loc == 1 && comp == 0) || (loc == 4 && (comp == 0 || comp == 1)))
                        {
                            // Outlet zero normal traction OR Half-Z symmetry zero tangential shear
                            if (std::abs(a) > tol || std::abs(b - 1.0) > tol || std::abs(g) > tol)
                            {
                                local_ok = 0;
                                local_reason = (loc == 4) ? "bc_zlo_not_halfz_symmetry" : "bc_outlet_not_zero_traction";
                                break;
                            }
                        }
                        else
                        {
                            // Homogeneous no-slip / impermeable / zero-tangential velocity Dirichlet
                            if (std::abs(a - 1.0) > tol || std::abs(b) > tol || std::abs(g) > tol)
                            {
                                local_ok = 0;
                                std::ostringstream oss;
                                oss << "bc_face_" << loc << "_comp_" << comp << "_not_zero_dirichlet";
                                local_reason = oss.str();
                                break;
                            }
                        }
                    }
                }
            }
        }

        const int all_ok = IBTK::IBTK_MPI::minReduction(local_ok);
        const int saw_inlet = IBTK::IBTK_MPI::maxReduction(local_saw_inlet);
        const double global_p_min = IBTK::IBTK_MPI::minReduction(local_p_in_min);
        const double global_p_max = IBTK::IBTK_MPI::maxReduction(local_p_in_max);

        if (all_ok == 0 || saw_inlet == 0 || !std::isfinite(global_p_min) || !std::isfinite(global_p_max))
        {
            reason = local_reason.empty() ? "unsupported_bc_configuration" : local_reason;
            p_in_new = 0.0;
            return false;
        }
        const double p_scale = std::max(1.0, std::abs(global_p_max));
        if ((global_p_max - global_p_min) > 1.0e-10 * p_scale)
        {
            reason = "nonuniform_inlet_pressure_bc";
            p_in_new = 0.0;
            return false;
        }
        p_in_new = 0.5 * (global_p_min + global_p_max);
        return true;
    }

    bool runPreprocess(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int u_cur_idx,
        const int u_old_idx,
        const int p_cur_idx,
        const int rho_cur_idx,
        const double mu0,
        const double fallback_rho0,
        const double dt,
        const double dt_prev,
        const int step_num,
        const bool is_first_step,
        const bool hierarchy_dirty,
        const double rel_tol,
        const double abs_tol,
        const int max_iters,
        bool& out_uniform_rho,
        double& out_rho0,
        std::string& out_reason)
    {
        out_reason.clear();
        const int rank = IBTK::IBTK_MPI::getRank();

        if (hierarchy_dirty || !d_layout.valid)
        {
            buildLayout(hierarchy);
            d_cached_L_valid = false;
        }

        // Check if current density is uniform across the entire hierarchy
        double local_rho_min = std::numeric_limits<double>::infinity();
        double local_rho_max = -std::numeric_limits<double>::infinity();
        int local_rho_finite = 1;
        const int finest_ln = hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& pbox = patch->getBox();
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> rho_data = patch->getPatchData(rho_cur_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    const SAMRAI::hier::Box<NDIM> sbox = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(pbox, axis);
                    const SAMRAI::pdat::ArrayData<NDIM, double>& arr = rho_data->getArrayData(axis);
                    const SAMRAI::hier::Box<NDIM>& abox = arr.getBox();
                    const double* ptr = arr.getPointer(0);
                    const int gnx = abox.numberCells(0);
                    const int gny = abox.numberCells(1);
                    const int ox = sbox.lower(0) - abox.lower(0);
                    const int oy = sbox.lower(1) - abox.lower(1);
                    const int oz = sbox.lower(2) - abox.lower(2);
                    const int sx = sbox.numberCells(0);
                    const int sy = sbox.numberCells(1);
                    const int sz = sbox.numberCells(2);
                    for (int k = 0; k < sz; ++k)
                    {
                        for (int j = 0; j < sy; ++j)
                        {
                            const size_t base = static_cast<size_t>(ox) +
                                                static_cast<size_t>(gnx) * (static_cast<size_t>(oy + j) +
                                                                            static_cast<size_t>(gny) * static_cast<size_t>(oz + k));
                            for (int i = 0; i < sx; ++i)
                            {
                                const double val = ptr[base + static_cast<size_t>(i)];
                                if (!std::isfinite(val)) local_rho_finite = 0;
                                local_rho_min = std::min(local_rho_min, val);
                                local_rho_max = std::max(local_rho_max, val);
                            }
                        }
                    }
                }
            }
        }
        const int all_rho_finite = IBTK::IBTK_MPI::minReduction(local_rho_finite);
        const double global_rho_min = IBTK::IBTK_MPI::minReduction(local_rho_min);
        const double global_rho_max = IBTK::IBTK_MPI::maxReduction(local_rho_max);
        if (all_rho_finite == 0 || !std::isfinite(global_rho_min) || !std::isfinite(global_rho_max))
        {
            out_reason = "nonfinite_rho_current";
            return false;
        }
        const double rho_scale = std::max(1.0, std::abs(global_rho_max));
        out_uniform_rho = (global_rho_max > 0.0) && ((global_rho_max - global_rho_min) <= 1.0e-12 * rho_scale);
        out_rho0 = out_uniform_rho ? 0.5 * (global_rho_min + global_rho_max) : fallback_rho0;

        C2GridConfig desired{};
        desired.nx0 = d_layout.nx0;
        desired.ny0 = d_layout.ny0;
        desired.nz0 = d_layout.nz0;
        desired.dx0 = d_layout.dx0;
        desired.mu0 = mu0;
        desired.rho0 = out_rho0;
        desired.enable_fp32_mg = d_cfg.enable_fp32_mg ? 1 : 0;
        desired.sticky_fp64_on_fb = d_cfg.sticky_fp64_on_fb ? 1 : 0;
        desired.max_iterations = max_iters;
        desired.rel_tol = rel_tol;
        desired.abs_tol = abs_tol;

        int upload_mode = 0;
        if (rank == 0)
        {
            if (hierarchy_dirty || !d_has_active_gpu_state || !d_gpu_handle || !d_has_uploaded_grid ||
                gridConfigChanged(desired))
            {
                upload_mode = 1;
            }
            else if (!out_uniform_rho)
            {
                upload_mode = 2;
            }
        }
        upload_mode = IBTK::IBTK_MPI::bcast(upload_mode, 0);
        const bool need_state_upload = (upload_mode == 1);

        // Collect all required field arrays on Rank 0 BEFORE Rank-0 CUDA calls
        if (need_state_upload)
        {
            for (int a = 0; a < 3; ++a)
            {
                gatherSideLevelToGlobalL0(hierarchy, u_cur_idx, a, d_h_u_cur_l0[a]);
                gatherSideLevelToL1(hierarchy, u_cur_idx, a, d_h_u_cur_l1[a]);
                if (!is_first_step && u_old_idx >= 0)
                {
                    gatherSideLevelToGlobalL0(hierarchy, u_old_idx, a, d_h_u_old_l0[a]);
                    gatherSideLevelToL1(hierarchy, u_old_idx, a, d_h_u_old_l1[a]);
                }
                gatherSideLevelToGlobalL0(hierarchy, rho_cur_idx, a, d_h_rho_l0[a]);
                gatherSideLevelToL1(hierarchy, rho_cur_idx, a, d_h_rho_l1[a]);
            }
            gatherCellLevelToGlobalL0(hierarchy, p_cur_idx, d_h_p_cur_l0);
            gatherCellLevelToL1(hierarchy, p_cur_idx, d_h_p_cur_l1);
        }
        else if (upload_mode == 2)
        {
            for (int a = 0; a < 3; ++a)
            {
                gatherSideLevelToGlobalL0(hierarchy, rho_cur_idx, a, d_h_rho_l0[a]);
                gatherSideLevelToL1(hierarchy, rho_cur_idx, a, d_h_rho_l1[a]);
            }
        }

        int decision = 0;
        std::string rank0_reason;
        if (rank == 0)
        {
            int rc = C2_OK;
            if (!ensureLoadedRank0(rank0_reason))
            {
                rc = C2_ERR_INVALID;
            }
            if (rc == C2_OK && (!d_gpu_handle || gridConfigChanged(desired)))
            {
                destroyHandleRank0();
                int create_err = C2_OK;
                d_gpu_handle = d_fn.create(&desired, &create_err);
                if (!d_gpu_handle || create_err != C2_OK)
                {
                    rc = (create_err != C2_OK) ? create_err : C2_ERR_CUDA;
                    rank0_reason = (rc == C2_ERR_CAPACITY) ? "c2_gpu_create_capacity" : "c2_gpu_create_failed";
                }
                else
                {
                    d_active_cfg = desired;
                    d_has_active_cfg = true;
                    d_has_uploaded_grid = false;
                }
            }

            if (rc == C2_OK && (need_state_upload || !d_has_uploaded_grid))
            {
                rc = d_fn.regrid(
                    d_gpu_handle,
                    static_cast<int>(d_layout.l1_boxes.size()),
                    d_layout.l1_boxes.empty() ? nullptr : d_layout.l1_boxes.data());
                C2Telemetry regrid_tel{};
                d_fn.get_telemetry(d_gpu_handle, &regrid_tel);
                if (d_cfg.enable_logging)
                {
                    SAMRAI::tbox::plog << "[C2_GPU_REGRID] step=" << step_num
                                       << " regrid=" << regrid_tel.regrid_count
                                       << " rc=" << rc
                                       << " n_l1=" << regrid_tel.last_n_l1 << "/" << regrid_tel.max_l1_patches
                                       << " peak_n_l1=" << regrid_tel.peak_n_l1
                                       << " tot_c_g0=" << regrid_tel.last_tot_c_g0 << "/" << regrid_tel.max_slot_g0
                                       << " peak_c_g0=" << regrid_tel.peak_tot_c_g0
                                       << " tot_c_g1=" << regrid_tel.last_tot_c_g1 << "/" << regrid_tel.max_slot_g1
                                       << " peak_c_g1=" << regrid_tel.peak_tot_c_g1
                                       << " tot_s_g0=" << regrid_tel.last_tot_s_g0[0] << ","
                                       << regrid_tel.last_tot_s_g0[1] << "," << regrid_tel.last_tot_s_g0[2] << "/"
                                       << regrid_tel.max_slot_g0
                                       << " max_s_g0=" << regrid_tel.last_max_s_g0
                                       << " peak_s_g0=" << regrid_tel.peak_tot_s_g0[0] << ","
                                       << regrid_tel.peak_tot_s_g0[1] << "," << regrid_tel.peak_tot_s_g0[2]
                                       << " peak_max_s_g0=" << regrid_tel.peak_max_s_g0
                                       << " tot_s_g1=" << regrid_tel.last_tot_s_g1[0] << ","
                                       << regrid_tel.last_tot_s_g1[1] << "," << regrid_tel.last_tot_s_g1[2] << "/"
                                       << regrid_tel.max_slot_g1
                                       << " max_s_g1=" << regrid_tel.last_max_s_g1
                                       << " peak_s_g1=" << regrid_tel.peak_tot_s_g1[0] << ","
                                       << regrid_tel.peak_tot_s_g1[1] << "," << regrid_tel.peak_tot_s_g1[2]
                                       << " peak_max_s_g1=" << regrid_tel.peak_max_s_g1
                                       << " l0_s_g3=" << regrid_tel.last_l0_s_g3[0] << ","
                                       << regrid_tel.last_l0_s_g3[1] << "," << regrid_tel.last_l0_s_g3[2] << "/"
                                       << regrid_tel.max_slot_g0
                                       << " l1_s_g3=" << regrid_tel.last_l1_s_g3[0] << ","
                                       << regrid_tel.last_l1_s_g3[1] << "," << regrid_tel.last_l1_s_g3[2] << "/"
                                       << regrid_tel.max_slot_g0
                                       << " max_s_g3=" << regrid_tel.last_max_s_g3
                                       << " peak_max_s_g3=" << regrid_tel.peak_max_s_g3
                                       << " nvml_peak_mib=" << regrid_tel.peak_nvml_mib
                                       << "\n" << std::flush;
                }
                if (rc != C2_OK)
                {
                    rank0_reason = readRank0FailureReason("c2_gpu_regrid", rc);
                }
                else
                {
                    d_has_uploaded_grid = true;
                    const double* u_cur_l0_p[3] = { d_h_u_cur_l0[0].data(), d_h_u_cur_l0[1].data(), d_h_u_cur_l0[2].data() };
                    const double* u_cur_l1_p[3] = {
                        detail::nonnull_cdata(d_h_u_cur_l1[0]),
                        detail::nonnull_cdata(d_h_u_cur_l1[1]),
                        detail::nonnull_cdata(d_h_u_cur_l1[2])
                    };
                    const double* u_old_l0_p[3] = { d_h_u_old_l0[0].data(), d_h_u_old_l0[1].data(), d_h_u_old_l0[2].data() };
                    const double* u_old_l1_p[3] = {
                        detail::nonnull_cdata(d_h_u_old_l1[0]),
                        detail::nonnull_cdata(d_h_u_old_l1[1]),
                        detail::nonnull_cdata(d_h_u_old_l1[2])
                    };
                    const double* rho_l0_p[3] = { d_h_rho_l0[0].data(), d_h_rho_l0[1].data(), d_h_rho_l0[2].data() };
                    const double* rho_l1_p[3] = {
                        detail::nonnull_cdata(d_h_rho_l1[0]),
                        detail::nonnull_cdata(d_h_rho_l1[1]),
                        detail::nonnull_cdata(d_h_rho_l1[2])
                    };
                    const bool have_u_old = (!is_first_step && u_old_idx >= 0);
                    const bool has_l1 = !d_layout.l1_boxes.empty();
                    rc = d_fn.upload_state(
                        d_gpu_handle,
                        u_cur_l0_p,
                        has_l1 ? u_cur_l1_p : nullptr,
                        have_u_old ? u_old_l0_p : nullptr,
                        (have_u_old && has_l1) ? u_old_l1_p : nullptr,
                        d_h_p_cur_l0.data(),
                        has_l1 ? detail::nonnull_cdata(d_h_p_cur_l1) : nullptr,
                        rho_l0_p,
                        has_l1 ? rho_l1_p : nullptr);
                    if (rc != C2_OK)
                    {
                        rank0_reason = readRank0FailureReason("c2_gpu_upload_state", rc);
                    }
                }
            }
            else if (rc == C2_OK && !out_uniform_rho)
            {
                const double* rho_l0_p[3] = { d_h_rho_l0[0].data(), d_h_rho_l0[1].data(), d_h_rho_l0[2].data() };
                const double* rho_l1_p[3] = {
                    detail::nonnull_cdata(d_h_rho_l1[0]),
                    detail::nonnull_cdata(d_h_rho_l1[1]),
                    detail::nonnull_cdata(d_h_rho_l1[2])
                };
                const bool has_l1 = !d_layout.l1_boxes.empty();
                rc = d_fn.upload_rho(d_gpu_handle, rho_l0_p, has_l1 ? rho_l1_p : nullptr);
                if (rc != C2_OK)
                {
                    rank0_reason = readRank0FailureReason("c2_gpu_upload_rho", rc);
                }
            }

            if (rc == C2_OK)
            {
                rc = d_fn.preprocess(
                    d_gpu_handle,
                    dt,
                    dt_prev,
                    is_first_step ? 1 : 0,
                    out_uniform_rho ? out_rho0 : -1.0);
                if (rc != C2_OK)
                {
                    rank0_reason = readRank0FailureReason("c2_gpu_preprocess", rc);
                }
            }
            decision = (rc == C2_OK) ? 1 : 0;
        }

        decision = IBTK::IBTK_MPI::bcast(decision, 0);
        if (decision == 1)
        {
            d_active_cfg = desired;
            d_has_active_cfg = true;
            d_has_active_gpu_state = true;
            return true;
        }
        d_has_active_cfg = false;
        d_has_active_gpu_state = false;
        d_cached_L_valid = false;
        if (rank == 0)
        {
            out_reason = rank0_reason.empty() ? "preprocess_failed" : rank0_reason;
        }
        else
        {
            out_reason = "rank0_preprocess_failed";
        }
        return false;
    }

    bool runIntegrateSolve(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int step_num,
        const int L_idx,
        const int f_idx,
        const bool has_f,
        const double dt,
        const double p_in_new,
        const bool brinkman_dirty,
        const bool download_rho,
        const int u_new_idx,
        const int p_new_idx,
        const int rho_new_idx,
        const bool log_solver_iters,
        const bool log_residual_norm,
        std::string& out_reason)
    {
        out_reason.clear();
        const int rank = IBTK::IBTK_MPI::getRank();

        // Check whether Brinkman penalty array L changed on any rank
        int local_L_same = (!brinkman_dirty && d_cached_L_valid) ? 1 : 0;
        packLocalSidesForComparison(hierarchy, L_idx, d_local_L_scratch);
        if (local_L_same == 1 && d_local_L_scratch.size() == d_cached_local_L.size())
        {
            for (size_t idx = 0; idx < d_local_L_scratch.size(); ++idx)
            {
                if (d_local_L_scratch[idx] != d_cached_local_L[idx])
                {
                    local_L_same = 0;
                    break;
                }
            }
        }
        else
        {
            local_L_same = 0;
        }
        const int all_L_same = IBTK::IBTK_MPI::minReduction(local_L_same);

        if (all_L_same == 0)
        {
            for (int a = 0; a < 3; ++a)
            {
                gatherSideLevelToGlobalL0(hierarchy, L_idx, a, d_h_L_l0[a]);
                gatherSideLevelToL1(hierarchy, L_idx, a, d_h_L_l1[a]);
            }
        }

        if (has_f && f_idx >= 0)
        {
            for (int a = 0; a < 3; ++a)
            {
                gatherSideLevelToGlobalL0(hierarchy, f_idx, a, d_h_f_l0[a]);
                gatherSideLevelToL1(hierarchy, f_idx, a, d_h_f_l1[a]);
            }
        }

        int decision = 0;
        int out_converged = 0;
        int out_iters = 0;
        double out_rel_res = 0.0;
        C2Telemetry tel{};
        std::string rank0_reason;

        if (rank == 0)
        {
            int rc = C2_OK;
            if (!d_gpu_handle)
            {
                rc = C2_ERR_INVALID;
                rank0_reason = "missing_gpu_handle";
            }

            const bool has_l1 = !d_layout.l1_boxes.empty();
            if (rc == C2_OK && all_L_same == 0)
            {
                const double* L_l0_p[3] = { d_h_L_l0[0].data(), d_h_L_l0[1].data(), d_h_L_l0[2].data() };
                const double* L_l1_p[3] = {
                    detail::nonnull_cdata(d_h_L_l1[0]),
                    detail::nonnull_cdata(d_h_L_l1[1]),
                    detail::nonnull_cdata(d_h_L_l1[2])
                };
                rc = d_fn.upload_brinkman(d_gpu_handle, L_l0_p, has_l1 ? L_l1_p : nullptr);
                if (rc != C2_OK)
                {
                    rank0_reason = readRank0FailureReason("c2_gpu_upload_brinkman", rc);
                }
            }

            if (rc == C2_OK)
            {
                allocateScratchOutputsRank0();
                const double* f_l0_p[3] = { d_h_f_l0[0].data(), d_h_f_l0[1].data(), d_h_f_l0[2].data() };
                const double* f_l1_p[3] = {
                    detail::nonnull_cdata(d_h_f_l1[0]),
                    detail::nonnull_cdata(d_h_f_l1[1]),
                    detail::nonnull_cdata(d_h_f_l1[2])
                };
                double* u_new_l0_p[3] = { d_h_u_new_l0[0].data(), d_h_u_new_l0[1].data(), d_h_u_new_l0[2].data() };
                double* u_new_l1_p[3] = {
                    detail::nonnull_data(d_h_u_new_l1[0]),
                    detail::nonnull_data(d_h_u_new_l1[1]),
                    detail::nonnull_data(d_h_u_new_l1[2])
                };
                double* rho_new_l0_p[3] = {
                    d_h_rho_new_l0[0].data(), d_h_rho_new_l0[1].data(), d_h_rho_new_l0[2].data()
                };
                double* rho_new_l1_p[3] = {
                    detail::nonnull_data(d_h_rho_new_l1[0]),
                    detail::nonnull_data(d_h_rho_new_l1[1]),
                    detail::nonnull_data(d_h_rho_new_l1[2])
                };

                rc = d_fn.integrate_solve(
                    d_gpu_handle,
                    dt,
                    p_in_new,
                    has_f ? f_l0_p : nullptr,
                    (has_f && has_l1) ? f_l1_p : nullptr,
                    has_f ? 1 : 0,
                    u_new_l0_p,
                    has_l1 ? u_new_l1_p : nullptr,
                    d_h_p_new_l0.data(),
                    has_l1 ? detail::nonnull_data(d_h_p_new_l1) : nullptr,
                    download_rho ? rho_new_l0_p : nullptr,
                    (download_rho && has_l1) ? rho_new_l1_p : nullptr,
                    download_rho ? 1 : 0,
                    &out_converged,
                    &out_iters,
                    &out_rel_res);

                d_fn.get_telemetry(d_gpu_handle, &tel);
                if (rc != C2_OK)
                {
                    rank0_reason = readRank0FailureReason("c2_gpu_integrate_solve", rc);
                }
                else if (out_converged == 0)
                {
                    rank0_reason = (tel.last_fb_reason[0] != '\0')
                                       ? std::string(tel.last_fb_reason)
                                       : std::string("stokes_not_converged");
                }
            }
            decision = (rc == C2_OK && out_converged == 1) ? 1 : 0;
        }

        decision = IBTK::IBTK_MPI::bcast(decision, 0);
        if (decision == 1)
        {
            if (all_L_same == 0)
            {
                d_cached_local_L = d_local_L_scratch;
                d_cached_L_valid = true;
            }
            for (int a = 0; a < 3; ++a)
            {
                scatterGlobalL0ToSideLevel(hierarchy, d_h_u_new_l0[a], u_new_idx, a);
                scatterL1ToSideLevel(hierarchy, d_h_u_new_l1[a], u_new_idx, a);
                if (download_rho)
                {
                    scatterGlobalL0ToSideLevel(hierarchy, d_h_rho_new_l0[a], rho_new_idx, a);
                    scatterL1ToSideLevel(hierarchy, d_h_rho_new_l1[a], rho_new_idx, a);
                }
            }
            scatterGlobalL0ToCellLevel(hierarchy, d_h_p_new_l0, p_new_idx);
            scatterL1ToCellLevel(hierarchy, d_h_p_new_l1, p_new_idx);

            if (rank == 0)
            {
                if (log_solver_iters)
                {
                    SAMRAI::tbox::plog << d_integrator_name
                                       << "::integrateHierarchy(): stokes solve number of iterations = " << out_iters
                                       << "\n";
                }
                if (log_residual_norm)
                {
                    SAMRAI::tbox::plog << d_integrator_name
                                       << "::integrateHierarchy(): stokes solve residual norm        = "
                                       << tel.last_true_abs_res << "\n";
                }
                if (d_cfg.enable_logging)
                {
                    SAMRAI::tbox::plog << "[C2_GPU] step=" << step_num
                                       << " iters=" << tel.last_iters
                                       << " fp32_iters=" << tel.last_fp32_iters
                                       << " fp64_fb_iters=" << tel.last_fp64_fb_iters
                                       << " fp64_fb=" << tel.last_fallback
                                       << " rel_res=" << tel.last_true_rel_res
                                       << " abs_res=" << tel.last_true_abs_res
                                       << " rhs_norm=" << tel.last_rhs_norm
                                       << " nvml_peak_mib=" << tel.peak_nvml_mib
                                       << " ms_prep=" << tel.ms_preprocess
                                       << " ms_ssprk=" << tel.ms_ssprk
                                       << " ms_stokes=" << tel.ms_stokes
                                       << " ms_h2d=" << tel.ms_h2d
                                       << " ms_d2h=" << tel.ms_d2h
                                       << " fb_reason=" << (tel.last_fb_reason[0] ? tel.last_fb_reason : "none")
                                       << "\n" << std::flush;
                }
            }
            d_has_active_gpu_state = true;
            return true;
        }

        // Recoverable failure or non-convergence: leave SAMRAI new state untouched!
        d_has_active_gpu_state = false;
        d_cached_L_valid = false;
        d_has_uploaded_grid = false;
        if (rank == 0)
        {
            out_reason = rank0_reason.empty() ? "integrate_solve_failed" : rank0_reason;
        }
        else
        {
            out_reason = "rank0_integrate_solve_failed";
        }
        return false;
    }

private:
    struct C2FnTable
    {
        void* (*create)(const C2GridConfig*, int*) = nullptr;
        void (*destroy)(void*) = nullptr;
        int (*regrid)(void*, int, const C2PatchBox*) = nullptr;
        int (*upload_state)(
            void*, const double* const[3], const double* const[3], const double* const[3], const double* const[3],
            const double*, const double*, const double* const[3], const double* const[3]) = nullptr;
        int (*upload_brinkman)(void*, const double* const[3], const double* const[3]) = nullptr;
        int (*upload_rho)(void*, const double* const[3], const double* const[3]) = nullptr;
        int (*preprocess)(void*, double, double, int, double) = nullptr;
        int (*integrate_solve)(
            void*, double, double, const double* const[3], const double* const[3], int, double* const[3],
            double* const[3], double*, double*, double* const[3], double* const[3], int, int*, int*, double*) = nullptr;
        void (*get_telemetry)(void*, C2Telemetry*) = nullptr;
    };

    struct PatchBoxDesc
    {
        int lo[3] = { 0, 0, 0 };
        int hi[3] = { -1, -1, -1 };
        int nx = 0;
        int ny = 0;
        int nz = 0;
        size_t c_count = 0;
        size_t s_count[3] = { 0, 0, 0 };
    };

    struct HierarchyLayout
    {
        bool valid = false;
        int nx0 = 0;
        int ny0 = 0;
        int nz0 = 0;
        double dx0 = 0.0;
        size_t l0_c_tot = 0;
        size_t l0_s_tot[3] = { 0, 0, 0 };
        std::vector<PatchBoxDesc> l0_local_patches;
        std::vector<PatchBoxDesc> l0_global_patches;
        std::vector<int> l0_c_counts, l0_c_displs;
        std::vector<int> l0_s_counts[3], l0_s_displs[3];
        size_t l0_c_recv_tot = 0;
        size_t l0_s_recv_tot[3] = { 0, 0, 0 };
        size_t local_l0_c_tot = 0;
        size_t local_l0_s_tot[3] = { 0, 0, 0 };

        std::vector<PatchBoxDesc> l1_local_patches;
        std::vector<PatchBoxDesc> l1_global_patches;
        std::vector<C2PatchBox> l1_boxes;
        std::vector<int> l1_c_counts, l1_c_displs;
        std::vector<int> l1_s_counts[3], l1_s_displs[3];
        size_t l1_c_tot = 0;
        size_t l1_s_tot[3] = { 0, 0, 0 };
        size_t local_l1_c_tot = 0;
        size_t local_l1_s_tot[3] = { 0, 0, 0 };
    };

    void destroyHandleRank0()
    {
        if (d_gpu_handle && d_fn.destroy)
        {
            d_fn.destroy(d_gpu_handle);
        }
        d_gpu_handle = nullptr;
        d_has_active_cfg = false;
        d_has_uploaded_grid = false;
        d_has_active_gpu_state = false;
    }

    bool ensureLoadedRank0(std::string& err_msg)
    {
        if (d_dl_handle) return true;

        std::vector<std::string> candidates;
        if (!d_cfg.library_path.empty())
        {
            candidates.push_back(d_cfg.library_path);
        }
        else
        {
            Dl_info info{};
            if (dladdr(reinterpret_cast<const void*>(&detail::c2_dladdr_anchor), &info) != 0 && info.dli_fname)
            {
                std::string fname(info.dli_fname);
                const size_t slash = fname.find_last_of('/');
                if (slash != std::string::npos)
                {
                    candidates.push_back(fname.substr(0, slash) + "/libibamr_c2_gpu.so");
                }
            }
            candidates.push_back("libibamr_c2_gpu.so");
        }

        std::ostringstream errs;
        void* handle = nullptr;
        for (const std::string& cand : candidates)
        {
            dlerror();
            handle = dlopen(cand.c_str(), RTLD_NOW | RTLD_LOCAL);
            if (handle) break;
            const char* e = dlerror();
            errs << cand << ": " << (e ? e : "unknown dlopen error") << "; ";
        }
        if (!handle)
        {
            err_msg = "dlopen_failed: " + errs.str();
            return false;
        }

        auto load_sym = [&](auto& fn_ptr, const char* name) -> bool {
            dlerror();
            void* sym = dlsym(handle, name);
            const char* e = dlerror();
            if (!sym || e)
            {
                err_msg = std::string("dlsym_failed(") + name + "): " + (e ? e : "null");
                return false;
            }
            fn_ptr = reinterpret_cast<std::remove_reference_t<decltype(fn_ptr)>>(sym);
            return true;
        };

        C2FnTable tbl{};
        if (!load_sym(tbl.create, "c2_gpu_create") || !load_sym(tbl.destroy, "c2_gpu_destroy") ||
            !load_sym(tbl.regrid, "c2_gpu_regrid") || !load_sym(tbl.upload_state, "c2_gpu_upload_state") ||
            !load_sym(tbl.upload_brinkman, "c2_gpu_upload_brinkman") ||
            !load_sym(tbl.upload_rho, "c2_gpu_upload_rho") || !load_sym(tbl.preprocess, "c2_gpu_preprocess") ||
            !load_sym(tbl.integrate_solve, "c2_gpu_integrate_solve") ||
            !load_sym(tbl.get_telemetry, "c2_gpu_get_telemetry"))
        {
            dlclose(handle);
            return false;
        }
        d_dl_handle = handle;
        d_fn = tbl;
        return true;
    }

    std::string readRank0FailureReason(const char* stage, int rc) const
    {
        C2Telemetry tel{};
        if (d_gpu_handle && d_fn.get_telemetry)
        {
            d_fn.get_telemetry(d_gpu_handle, &tel);
        }
        if (tel.last_fb_reason[0] != '\0')
        {
            return std::string(stage) + ":" + tel.last_fb_reason;
        }
        std::ostringstream oss;
        oss << stage << ":rc=" << rc;
        return oss.str();
    }

    bool gridConfigChanged(const C2GridConfig& c) const
    {
        if (!d_has_active_cfg) return true;
        return d_active_cfg.nx0 != c.nx0 || d_active_cfg.ny0 != c.ny0 || d_active_cfg.nz0 != c.nz0 ||
               std::abs(d_active_cfg.dx0 - c.dx0) > 1.0e-13 * std::max(1.0, std::abs(c.dx0)) ||
               std::abs(d_active_cfg.mu0 - c.mu0) > 1.0e-13 * std::max(1.0, std::abs(c.mu0)) ||
               std::abs(d_active_cfg.rho0 - c.rho0) > 1.0e-13 * std::max(1.0, std::abs(c.rho0)) ||
               d_active_cfg.enable_fp32_mg != c.enable_fp32_mg ||
               d_active_cfg.sticky_fp64_on_fb != c.sticky_fp64_on_fb ||
               d_active_cfg.max_iterations != c.max_iterations ||
               d_active_cfg.rel_tol != c.rel_tol || d_active_cfg.abs_tol != c.abs_tol;
    }

    static PatchBoxDesc makeBoxDesc(const SAMRAI::hier::Box<NDIM>& box)
    {
        PatchBoxDesc d{};
        for (int a = 0; a < 3; ++a)
        {
            d.lo[a] = box.lower(a);
            d.hi[a] = box.upper(a);
        }
        d.nx = d.hi[0] - d.lo[0] + 1;
        d.ny = d.hi[1] - d.lo[1] + 1;
        d.nz = d.hi[2] - d.lo[2] + 1;
        d.c_count = static_cast<size_t>(d.nx) * d.ny * d.nz;
        for (int a = 0; a < 3; ++a)
        {
            const int sx = d.nx + (a == 0 ? 1 : 0);
            const int sy = d.ny + (a == 1 ? 1 : 0);
            const int sz = d.nz + (a == 2 ? 1 : 0);
            d.s_count[a] = static_cast<size_t>(sx) * sy * sz;
        }
        return d;
    }

    static void gatherLevelBoxesToRank0(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int ln,
        std::vector<PatchBoxDesc>& out_local,
        std::vector<PatchBoxDesc>& out_global,
        std::vector<int>& out_c_counts,
        std::vector<int>& out_c_displs,
        std::vector<int> out_s_counts[3],
        std::vector<int> out_s_displs[3],
        size_t& out_local_c_tot,
        size_t out_local_s_tot[3],
        size_t& out_global_c_tot,
        size_t out_global_s_tot[3])
    {
        MPI_Comm comm = IBTK::IBTK_MPI::getCommunicator();
        const int rank = IBTK::IBTK_MPI::getRank();
        const int nprocs = IBTK::IBTK_MPI::getNodes();

        out_local.clear();
        out_local_c_tot = 0;
        out_local_s_tot[0] = out_local_s_tot[1] = out_local_s_tot[2] = 0;

        if (ln <= hierarchy->getFinestLevelNumber())
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                PatchBoxDesc d = makeBoxDesc(level->getPatch(p())->getBox());
                out_local_c_tot += d.c_count;
                for (int a = 0; a < 3; ++a) out_local_s_tot[a] += d.s_count[a];
                out_local.push_back(d);
            }
        }

        std::vector<int> local_ints;
        local_ints.reserve(out_local.size() * 6);
        for (const auto& d : out_local)
        {
            local_ints.push_back(d.lo[0]);
            local_ints.push_back(d.lo[1]);
            local_ints.push_back(d.lo[2]);
            local_ints.push_back(d.hi[0]);
            local_ints.push_back(d.hi[1]);
            local_ints.push_back(d.hi[2]);
        }

        const int local_int_count = static_cast<int>(local_ints.size());
        std::vector<int> all_int_counts(rank == 0 ? nprocs : 0, 0);
        MPI_Gather(&local_int_count, 1, MPI_INT, all_int_counts.data(), 1, MPI_INT, 0, comm);

        std::vector<int> int_displs(rank == 0 ? nprocs : 0, 0);
        int tot_ints = 0;
        if (rank == 0)
        {
            for (int r = 0; r < nprocs; ++r)
            {
                int_displs[r] = tot_ints;
                tot_ints += all_int_counts[r];
            }
        }
        std::vector<int> all_ints(rank == 0 ? tot_ints : 0);
        MPI_Gatherv(
            detail::nonnull_idata(local_ints),
            local_int_count,
            MPI_INT,
            detail::nonnull_idata(all_ints),
            detail::nonnull_idata(all_int_counts),
            detail::nonnull_idata(int_displs),
            MPI_INT,
            0,
            comm);

        out_global.clear();
        out_global_c_tot = 0;
        out_global_s_tot[0] = out_global_s_tot[1] = out_global_s_tot[2] = 0;
        if (rank == 0)
        {
            out_c_counts.assign(nprocs, 0);
            out_c_displs.assign(nprocs, 0);
            for (int a = 0; a < 3; ++a)
            {
                out_s_counts[a].assign(nprocs, 0);
                out_s_displs[a].assign(nprocs, 0);
            }
            for (int r = 0; r < nprocs; ++r)
            {
                out_c_displs[r] = static_cast<int>(out_global_c_tot);
                for (int a = 0; a < 3; ++a) out_s_displs[a][r] = static_cast<int>(out_global_s_tot[a]);
                const int n_p = all_int_counts[r] / 6;
                const int off = int_displs[r];
                for (int k = 0; k < n_p; ++k)
                {
                    const int base = off + 6 * k;
                    PatchBoxDesc d{};
                    d.lo[0] = all_ints[base + 0];
                    d.lo[1] = all_ints[base + 1];
                    d.lo[2] = all_ints[base + 2];
                    d.hi[0] = all_ints[base + 3];
                    d.hi[1] = all_ints[base + 4];
                    d.hi[2] = all_ints[base + 5];
                    d.nx = d.hi[0] - d.lo[0] + 1;
                    d.ny = d.hi[1] - d.lo[1] + 1;
                    d.nz = d.hi[2] - d.lo[2] + 1;
                    d.c_count = static_cast<size_t>(d.nx) * d.ny * d.nz;
                    out_c_counts[r] += static_cast<int>(d.c_count);
                    out_global_c_tot += d.c_count;
                    for (int a = 0; a < 3; ++a)
                    {
                        const int sx = d.nx + (a == 0 ? 1 : 0);
                        const int sy = d.ny + (a == 1 ? 1 : 0);
                        const int sz = d.nz + (a == 2 ? 1 : 0);
                        d.s_count[a] = static_cast<size_t>(sx) * sy * sz;
                        out_s_counts[a][r] += static_cast<int>(d.s_count[a]);
                        out_global_s_tot[a] += d.s_count[a];
                    }
                    out_global.push_back(d);
                }
            }
        }
    }

    void buildLayout(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy)
    {
        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM>> grid_geom = hierarchy->getGridGeometry();
        const SAMRAI::hier::Box<NDIM>& domain_box = grid_geom->getPhysicalDomain()[0];
        d_layout.nx0 = domain_box.numberCells(0);
        d_layout.ny0 = domain_box.numberCells(1);
        d_layout.nz0 = domain_box.numberCells(2);
        d_layout.dx0 = grid_geom->getDx()[0];
        d_layout.l0_c_tot = static_cast<size_t>(d_layout.nx0) * d_layout.ny0 * d_layout.nz0;
        for (int a = 0; a < 3; ++a)
        {
            const size_t sx = static_cast<size_t>(d_layout.nx0 + (a == 0 ? 1 : 0));
            const size_t sy = static_cast<size_t>(d_layout.ny0 + (a == 1 ? 1 : 0));
            const size_t sz = static_cast<size_t>(d_layout.nz0 + (a == 2 ? 1 : 0));
            d_layout.l0_s_tot[a] = sx * sy * sz;
        }

        gatherLevelBoxesToRank0(
            hierarchy,
            0,
            d_layout.l0_local_patches,
            d_layout.l0_global_patches,
            d_layout.l0_c_counts,
            d_layout.l0_c_displs,
            d_layout.l0_s_counts,
            d_layout.l0_s_displs,
            d_layout.local_l0_c_tot,
            d_layout.local_l0_s_tot,
            d_layout.l0_c_recv_tot,
            d_layout.l0_s_recv_tot);

        gatherLevelBoxesToRank0(
            hierarchy,
            1,
            d_layout.l1_local_patches,
            d_layout.l1_global_patches,
            d_layout.l1_c_counts,
            d_layout.l1_c_displs,
            d_layout.l1_s_counts,
            d_layout.l1_s_displs,
            d_layout.local_l1_c_tot,
            d_layout.local_l1_s_tot,
            d_layout.l1_c_tot,
            d_layout.l1_s_tot);

        d_layout.l1_boxes.clear();
        if (IBTK::IBTK_MPI::getRank() == 0)
        {
            d_layout.l1_boxes.reserve(d_layout.l1_global_patches.size());
            for (const auto& p : d_layout.l1_global_patches)
            {
                C2PatchBox b{};
                b.lo[0] = p.lo[0];
                b.lo[1] = p.lo[1];
                b.lo[2] = p.lo[2];
                b.hi[0] = p.hi[0];
                b.hi[1] = p.hi[1];
                b.hi[2] = p.hi[2];
                d_layout.l1_boxes.push_back(b);
            }
            for (int a = 0; a < 3; ++a)
            {
                d_h_u_cur_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
                d_h_u_cur_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
                d_h_u_old_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
                d_h_u_old_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
                d_h_rho_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
                d_h_rho_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
                d_h_L_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
                d_h_L_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
                d_h_f_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
                d_h_f_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
            }
            d_h_p_cur_l0.resize(d_layout.l0_c_tot, 0.0);
            d_h_p_cur_l1.resize(d_layout.l1_c_tot, 0.0);
        }
        d_layout.valid = true;
    }

    void allocateScratchOutputsRank0()
    {
        for (int a = 0; a < 3; ++a)
        {
            d_h_u_new_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
            d_h_u_new_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
            d_h_rho_new_l0[a].resize(d_layout.l0_s_tot[a], 0.0);
            d_h_rho_new_l1[a].resize(d_layout.l1_s_tot[a], 0.0);
        }
        d_h_p_new_l0.resize(d_layout.l0_c_tot, 0.0);
        d_h_p_new_l1.resize(d_layout.l1_c_tot, 0.0);
    }

    static void packLocalSideData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& level,
        const int data_idx,
        const int axis,
        const size_t total_count,
        std::vector<double>& out_buf)
    {
        out_buf.resize(total_count);
        size_t pos = 0;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM> sbox = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> sdata = patch->getPatchData(data_idx);
            const SAMRAI::pdat::ArrayData<NDIM, double>& arr = sdata->getArrayData(axis);
            const SAMRAI::hier::Box<NDIM>& abox = arr.getBox();
            const double* ptr = arr.getPointer(0);
            const int gnx = abox.numberCells(0);
            const int gny = abox.numberCells(1);
            const int ox = sbox.lower(0) - abox.lower(0);
            const int oy = sbox.lower(1) - abox.lower(1);
            const int oz = sbox.lower(2) - abox.lower(2);
            const int sx = sbox.numberCells(0);
            const int sy = sbox.numberCells(1);
            const int sz = sbox.numberCells(2);
            for (int k = 0; k < sz; ++k)
            {
                for (int j = 0; j < sy; ++j)
                {
                    const size_t src_base = static_cast<size_t>(ox) +
                                            static_cast<size_t>(gnx) * (static_cast<size_t>(oy + j) +
                                                                        static_cast<size_t>(gny) * static_cast<size_t>(oz + k));
                    std::memcpy(out_buf.data() + pos, ptr + src_base, static_cast<size_t>(sx) * sizeof(double));
                    pos += static_cast<size_t>(sx);
                }
            }
        }
    }

    static void unpackLocalSideData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& level,
        const int data_idx,
        const int axis,
        const std::vector<double>& in_buf)
    {
        size_t pos = 0;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM> sbox = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> sdata = patch->getPatchData(data_idx);
            SAMRAI::pdat::ArrayData<NDIM, double>& arr = sdata->getArrayData(axis);
            const SAMRAI::hier::Box<NDIM>& abox = arr.getBox();
            double* ptr = arr.getPointer(0);
            const int gnx = abox.numberCells(0);
            const int gny = abox.numberCells(1);
            const int ox = sbox.lower(0) - abox.lower(0);
            const int oy = sbox.lower(1) - abox.lower(1);
            const int oz = sbox.lower(2) - abox.lower(2);
            const int sx = sbox.numberCells(0);
            const int sy = sbox.numberCells(1);
            const int sz = sbox.numberCells(2);
            for (int k = 0; k < sz; ++k)
            {
                for (int j = 0; j < sy; ++j)
                {
                    const size_t dst_base = static_cast<size_t>(ox) +
                                            static_cast<size_t>(gnx) * (static_cast<size_t>(oy + j) +
                                                                        static_cast<size_t>(gny) * static_cast<size_t>(oz + k));
                    std::memcpy(ptr + dst_base, in_buf.data() + pos, static_cast<size_t>(sx) * sizeof(double));
                    pos += static_cast<size_t>(sx);
                }
            }
        }
    }

    static void packLocalCellData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& level,
        const int data_idx,
        const size_t total_count,
        std::vector<double>& out_buf)
    {
        out_buf.resize(total_count);
        size_t pos = 0;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& cbox = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> cdata = patch->getPatchData(data_idx);
            const SAMRAI::pdat::ArrayData<NDIM, double>& arr = cdata->getArrayData();
            const SAMRAI::hier::Box<NDIM>& abox = arr.getBox();
            const double* ptr = arr.getPointer(0);
            const int gnx = abox.numberCells(0);
            const int gny = abox.numberCells(1);
            const int ox = cbox.lower(0) - abox.lower(0);
            const int oy = cbox.lower(1) - abox.lower(1);
            const int oz = cbox.lower(2) - abox.lower(2);
            const int nx = cbox.numberCells(0);
            const int ny = cbox.numberCells(1);
            const int nz = cbox.numberCells(2);
            for (int k = 0; k < nz; ++k)
            {
                for (int j = 0; j < ny; ++j)
                {
                    const size_t src_base = static_cast<size_t>(ox) +
                                            static_cast<size_t>(gnx) * (static_cast<size_t>(oy + j) +
                                                                        static_cast<size_t>(gny) * static_cast<size_t>(oz + k));
                    std::memcpy(out_buf.data() + pos, ptr + src_base, static_cast<size_t>(nx) * sizeof(double));
                    pos += static_cast<size_t>(nx);
                }
            }
        }
    }

    static void unpackLocalCellData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>& level,
        const int data_idx,
        const std::vector<double>& in_buf)
    {
        size_t pos = 0;
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& cbox = patch->getBox();
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double>> cdata = patch->getPatchData(data_idx);
            SAMRAI::pdat::ArrayData<NDIM, double>& arr = cdata->getArrayData();
            const SAMRAI::hier::Box<NDIM>& abox = arr.getBox();
            double* ptr = arr.getPointer(0);
            const int gnx = abox.numberCells(0);
            const int gny = abox.numberCells(1);
            const int ox = cbox.lower(0) - abox.lower(0);
            const int oy = cbox.lower(1) - abox.lower(1);
            const int oz = cbox.lower(2) - abox.lower(2);
            const int nx = cbox.numberCells(0);
            const int ny = cbox.numberCells(1);
            const int nz = cbox.numberCells(2);
            for (int k = 0; k < nz; ++k)
            {
                for (int j = 0; j < ny; ++j)
                {
                    const size_t dst_base = static_cast<size_t>(ox) +
                                            static_cast<size_t>(gnx) * (static_cast<size_t>(oy + j) +
                                                                        static_cast<size_t>(gny) * static_cast<size_t>(oz + k));
                    std::memcpy(ptr + dst_base, in_buf.data() + pos, static_cast<size_t>(nx) * sizeof(double));
                    pos += static_cast<size_t>(nx);
                }
            }
        }
    }

    static void packLocalSidesForComparison(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int data_idx,
        std::vector<double>& out_buf)
    {
        size_t total = 0;
        const int finest_ln = hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                const SAMRAI::hier::Box<NDIM>& pbox = level->getPatch(p())->getBox();
                for (int a = 0; a < 3; ++a)
                {
                    const SAMRAI::hier::Box<NDIM> sbox = SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(pbox, a);
                    total += static_cast<size_t>(sbox.size());
                }
            }
        }
        out_buf.resize(total);
        size_t pos = 0;
        std::vector<double> tmp;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (int a = 0; a < 3; ++a)
            {
                size_t level_a = 0;
                for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    level_a += static_cast<size_t>(
                        SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(level->getPatch(p())->getBox(), a).size());
                }
                packLocalSideData(level, data_idx, a, level_a, tmp);
                if (level_a > 0)
                {
                    std::memcpy(out_buf.data() + pos, tmp.data(), level_a * sizeof(double));
                    pos += level_a;
                }
            }
        }
    }

    void gatherSideLevelToGlobalL0(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int data_idx,
        const int axis,
        std::vector<double>& out_global_l0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l0 = hierarchy->getPatchLevel(0);
        packLocalSideData(l0, data_idx, axis, d_layout.local_l0_s_tot[axis], d_send_tmp);
        const int rank = IBTK::IBTK_MPI::getRank();
        if (rank == 0) d_recv_tmp.resize(d_layout.l0_s_recv_tot[axis]);
        MPI_Gatherv(
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l0_s_tot[axis]),
            MPI_DOUBLE,
            detail::nonnull_data(d_recv_tmp),
            detail::nonnull_idata(d_layout.l0_s_counts[axis]),
            detail::nonnull_idata(d_layout.l0_s_displs[axis]),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
        if (rank == 0)
        {
            out_global_l0.resize(d_layout.l0_s_tot[axis]);
            const int gsx = d_layout.nx0 + (axis == 0 ? 1 : 0);
            const int gsy = d_layout.ny0 + (axis == 1 ? 1 : 0);
            size_t pos = 0;
            for (const auto& d : d_layout.l0_global_patches)
            {
                const int psx = d.nx + (axis == 0 ? 1 : 0);
                const int psy = d.ny + (axis == 1 ? 1 : 0);
                const int psz = d.nz + (axis == 2 ? 1 : 0);
                for (int k = 0; k < psz; ++k)
                {
                    for (int j = 0; j < psy; ++j)
                    {
                        const size_t dst_off = static_cast<size_t>(d.lo[0]) +
                                               static_cast<size_t>(gsx) * (static_cast<size_t>(d.lo[1] + j) +
                                                                           static_cast<size_t>(gsy) *
                                                                               static_cast<size_t>(d.lo[2] + k));
                        std::memcpy(
                            out_global_l0.data() + dst_off, d_recv_tmp.data() + pos, static_cast<size_t>(psx) * sizeof(double));
                        pos += static_cast<size_t>(psx);
                    }
                }
            }
        }
    }

    void scatterGlobalL0ToSideLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const std::vector<double>& global_l0,
        const int data_idx,
        const int axis)
    {
        const int rank = IBTK::IBTK_MPI::getRank();
        if (rank == 0)
        {
            d_recv_tmp.resize(d_layout.l0_s_recv_tot[axis]);
            const int gsx = d_layout.nx0 + (axis == 0 ? 1 : 0);
            const int gsy = d_layout.ny0 + (axis == 1 ? 1 : 0);
            size_t pos = 0;
            for (const auto& d : d_layout.l0_global_patches)
            {
                const int psx = d.nx + (axis == 0 ? 1 : 0);
                const int psy = d.ny + (axis == 1 ? 1 : 0);
                const int psz = d.nz + (axis == 2 ? 1 : 0);
                for (int k = 0; k < psz; ++k)
                {
                    for (int j = 0; j < psy; ++j)
                    {
                        const size_t src_off = static_cast<size_t>(d.lo[0]) +
                                               static_cast<size_t>(gsx) * (static_cast<size_t>(d.lo[1] + j) +
                                                                           static_cast<size_t>(gsy) *
                                                                               static_cast<size_t>(d.lo[2] + k));
                        std::memcpy(
                            d_recv_tmp.data() + pos, global_l0.data() + src_off, static_cast<size_t>(psx) * sizeof(double));
                        pos += static_cast<size_t>(psx);
                    }
                }
            }
        }
        d_send_tmp.resize(d_layout.local_l0_s_tot[axis]);
        MPI_Scatterv(
            detail::nonnull_data(d_recv_tmp),
            detail::nonnull_idata(d_layout.l0_s_counts[axis]),
            detail::nonnull_idata(d_layout.l0_s_displs[axis]),
            MPI_DOUBLE,
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l0_s_tot[axis]),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l0 = hierarchy->getPatchLevel(0);
        unpackLocalSideData(l0, data_idx, axis, d_send_tmp);
    }

    void gatherCellLevelToGlobalL0(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int data_idx,
        std::vector<double>& out_global_l0)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l0 = hierarchy->getPatchLevel(0);
        packLocalCellData(l0, data_idx, d_layout.local_l0_c_tot, d_send_tmp);
        const int rank = IBTK::IBTK_MPI::getRank();
        if (rank == 0) d_recv_tmp.resize(d_layout.l0_c_recv_tot);
        MPI_Gatherv(
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l0_c_tot),
            MPI_DOUBLE,
            detail::nonnull_data(d_recv_tmp),
            detail::nonnull_idata(d_layout.l0_c_counts),
            detail::nonnull_idata(d_layout.l0_c_displs),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
        if (rank == 0)
        {
            out_global_l0.resize(d_layout.l0_c_tot);
            const int gnx = d_layout.nx0;
            const int gny = d_layout.ny0;
            size_t pos = 0;
            for (const auto& d : d_layout.l0_global_patches)
            {
                for (int k = 0; k < d.nz; ++k)
                {
                    for (int j = 0; j < d.ny; ++j)
                    {
                        const size_t dst_off = static_cast<size_t>(d.lo[0]) +
                                               static_cast<size_t>(gnx) * (static_cast<size_t>(d.lo[1] + j) +
                                                                           static_cast<size_t>(gny) *
                                                                               static_cast<size_t>(d.lo[2] + k));
                        std::memcpy(
                            out_global_l0.data() + dst_off, d_recv_tmp.data() + pos, static_cast<size_t>(d.nx) * sizeof(double));
                        pos += static_cast<size_t>(d.nx);
                    }
                }
            }
        }
    }

    void scatterGlobalL0ToCellLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const std::vector<double>& global_l0,
        const int data_idx)
    {
        const int rank = IBTK::IBTK_MPI::getRank();
        if (rank == 0)
        {
            d_recv_tmp.resize(d_layout.l0_c_recv_tot);
            const int gnx = d_layout.nx0;
            const int gny = d_layout.ny0;
            size_t pos = 0;
            for (const auto& d : d_layout.l0_global_patches)
            {
                for (int k = 0; k < d.nz; ++k)
                {
                    for (int j = 0; j < d.ny; ++j)
                    {
                        const size_t src_off = static_cast<size_t>(d.lo[0]) +
                                               static_cast<size_t>(gnx) * (static_cast<size_t>(d.lo[1] + j) +
                                                                           static_cast<size_t>(gny) *
                                                                               static_cast<size_t>(d.lo[2] + k));
                        std::memcpy(
                            d_recv_tmp.data() + pos, global_l0.data() + src_off, static_cast<size_t>(d.nx) * sizeof(double));
                        pos += static_cast<size_t>(d.nx);
                    }
                }
            }
        }
        d_send_tmp.resize(d_layout.local_l0_c_tot);
        MPI_Scatterv(
            detail::nonnull_data(d_recv_tmp),
            detail::nonnull_idata(d_layout.l0_c_counts),
            detail::nonnull_idata(d_layout.l0_c_displs),
            MPI_DOUBLE,
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l0_c_tot),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l0 = hierarchy->getPatchLevel(0);
        unpackLocalCellData(l0, data_idx, d_send_tmp);
    }

    void gatherSideLevelToL1(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int data_idx,
        const int axis,
        std::vector<double>& out_l1)
    {
        if (hierarchy->getFinestLevelNumber() < 1)
        {
            out_l1.clear();
            return;
        }
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l1 = hierarchy->getPatchLevel(1);
        packLocalSideData(l1, data_idx, axis, d_layout.local_l1_s_tot[axis], d_send_tmp);
        if (IBTK::IBTK_MPI::getRank() == 0) out_l1.resize(d_layout.l1_s_tot[axis]);
        MPI_Gatherv(
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l1_s_tot[axis]),
            MPI_DOUBLE,
            detail::nonnull_data(out_l1),
            detail::nonnull_idata(d_layout.l1_s_counts[axis]),
            detail::nonnull_idata(d_layout.l1_s_displs[axis]),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
    }

    void scatterL1ToSideLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const std::vector<double>& in_l1,
        const int data_idx,
        const int axis)
    {
        if (hierarchy->getFinestLevelNumber() < 1) return;
        d_send_tmp.resize(d_layout.local_l1_s_tot[axis]);
        std::vector<double>* send_ptr =
            (IBTK::IBTK_MPI::getRank() == 0) ? const_cast<std::vector<double>*>(&in_l1) : &d_recv_tmp;
        MPI_Scatterv(
            detail::nonnull_data(*send_ptr),
            detail::nonnull_idata(d_layout.l1_s_counts[axis]),
            detail::nonnull_idata(d_layout.l1_s_displs[axis]),
            MPI_DOUBLE,
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l1_s_tot[axis]),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l1 = hierarchy->getPatchLevel(1);
        unpackLocalSideData(l1, data_idx, axis, d_send_tmp);
    }

    void gatherCellLevelToL1(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const int data_idx,
        std::vector<double>& out_l1)
    {
        if (hierarchy->getFinestLevelNumber() < 1)
        {
            out_l1.clear();
            return;
        }
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l1 = hierarchy->getPatchLevel(1);
        packLocalCellData(l1, data_idx, d_layout.local_l1_c_tot, d_send_tmp);
        if (IBTK::IBTK_MPI::getRank() == 0) out_l1.resize(d_layout.l1_c_tot);
        MPI_Gatherv(
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l1_c_tot),
            MPI_DOUBLE,
            detail::nonnull_data(out_l1),
            detail::nonnull_idata(d_layout.l1_c_counts),
            detail::nonnull_idata(d_layout.l1_c_displs),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
    }

    void scatterL1ToCellLevel(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>>& hierarchy,
        const std::vector<double>& in_l1,
        const int data_idx)
    {
        if (hierarchy->getFinestLevelNumber() < 1) return;
        d_send_tmp.resize(d_layout.local_l1_c_tot);
        std::vector<double>* send_ptr =
            (IBTK::IBTK_MPI::getRank() == 0) ? const_cast<std::vector<double>*>(&in_l1) : &d_recv_tmp;
        MPI_Scatterv(
            detail::nonnull_data(*send_ptr),
            detail::nonnull_idata(d_layout.l1_c_counts),
            detail::nonnull_idata(d_layout.l1_c_displs),
            MPI_DOUBLE,
            detail::nonnull_data(d_send_tmp),
            static_cast<int>(d_layout.local_l1_c_tot),
            MPI_DOUBLE,
            0,
            IBTK::IBTK_MPI::getCommunicator());
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> l1 = hierarchy->getPatchLevel(1);
        unpackLocalCellData(l1, data_idx, d_send_tmp);
    }

    Config d_cfg;
    std::string d_integrator_name;
    void* d_dl_handle = nullptr;
    C2FnTable d_fn;
    void* d_gpu_handle = nullptr;
    C2GridConfig d_active_cfg{};
    bool d_has_active_cfg = false;
    bool d_has_uploaded_grid = false;
    bool d_has_active_gpu_state = false;

    HierarchyLayout d_layout;
    std::vector<double> d_send_tmp;
    std::vector<double> d_recv_tmp;

    std::vector<double> d_h_u_cur_l0[3], d_h_u_cur_l1[3];
    std::vector<double> d_h_u_old_l0[3], d_h_u_old_l1[3];
    std::vector<double> d_h_rho_l0[3], d_h_rho_l1[3];
    std::vector<double> d_h_L_l0[3], d_h_L_l1[3];
    std::vector<double> d_h_f_l0[3], d_h_f_l1[3];
    std::vector<double> d_h_p_cur_l0, d_h_p_cur_l1;

    std::vector<double> d_h_u_new_l0[3], d_h_u_new_l1[3];
    std::vector<double> d_h_rho_new_l0[3], d_h_rho_new_l1[3];
    std::vector<double> d_h_p_new_l0, d_h_p_new_l1;

    std::vector<double> d_local_L_scratch;
    std::vector<double> d_cached_local_L;
    bool d_cached_L_valid = false;
};

} // namespace IBAMR

#endif // (NDIM == 3)
#endif // included_IBAMR_INSVCStaggeredConservativeGPUBridge
