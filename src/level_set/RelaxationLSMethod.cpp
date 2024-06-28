// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/LSInitStrategy.h"
#include "ibamr/RelaxationLSMethod.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_MPI.h"

#include "BasePatchLevel.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

// FORTRAN ROUTINES
#if (NDIM == 2)
#define RELAXATION_LS_1ST_ORDER_FC IBAMR_FC_FUNC(relaxationls1storder2d, RELAXATIONLS1STORDER2D)
#define RELAXATION_LS_3RD_ORDER_ENO_FC IBAMR_FC_FUNC(relaxationls3rdordereno2d, RELAXATIONLS3RDORDERENO2D)
#define GODUNOV_HAMILTONIAN_3RD_ORDER_ENO_FC IBAMR_FC_FUNC(godunovhamiltonianeno2d, GODUNOVHAMILTONIANENO2D)
#define RELAXATION_LS_3RD_ORDER_WENO_FC IBAMR_FC_FUNC(relaxationls3rdorderweno2d, RELAXATIONLS3RDORDERWENO2D)
#define GODUNOV_HAMILTONIAN_3RD_ORDER_WENO_FC IBAMR_FC_FUNC(godunovhamiltonianweno2d, GODUNOVHAMILTONIANWENO2D)
#define RELAXATION_LS_5TH_ORDER_WENO_FC IBAMR_FC_FUNC(relaxationls5thorderweno2d, RELAXATIONLS5THORDERWENO2D)
#define GODUNOV_HAMILTONIAN_5TH_ORDER_WENO_FC                                                                          \
    IBAMR_FC_FUNC(godunovhamiltonian5thorderweno2d, GODUNOVHAMILTONIAN5THORDERWENO2D)
#define PROJECT_LS_MASS_CONSTRAINT_FC IBAMR_FC_FUNC(projectlsmassconstraint2d, PROJECTLSMASSCONSTRAINT2D)
#endif

#if (NDIM == 3)
#define RELAXATION_LS_1ST_ORDER_FC IBAMR_FC_FUNC(relaxationls1storder3d, RELAXATIONLS1STORDER3D)
#define RELAXATION_LS_3RD_ORDER_ENO_FC IBAMR_FC_FUNC(relaxationls3rdordereno3d, RELAXATIONLS3RDORDERENO3D)
#define GODUNOV_HAMILTONIAN_3RD_ORDER_ENO_FC IBAMR_FC_FUNC(godunovhamiltonianeno3d, GODUNOVHAMILTONIANENO3D)
#define RELAXATION_LS_3RD_ORDER_WENO_FC IBAMR_FC_FUNC(relaxationls3rdorderweno3d, RELAXATIONLS3RDORDERWENO3D)
#define GODUNOV_HAMILTONIAN_3RD_ORDER_WENO_FC IBAMR_FC_FUNC(godunovhamiltonianweno3d, GODUNOVHAMILTONIANWENO3D)
#define RELAXATION_LS_5TH_ORDER_WENO_FC IBAMR_FC_FUNC(relaxationls5thorderweno3d, RELAXATIONLS5THORDERWENO3D)
#define GODUNOV_HAMILTONIAN_5TH_ORDER_WENO_FC                                                                          \
    IBAMR_FC_FUNC(godunovhamiltonian5thorderweno3d, GODUNOVHAMILTONIAN5THORDERWENO3D)
#define PROJECT_LS_MASS_CONSTRAINT_FC IBAMR_FC_FUNC(projectlsmassconstraint3d, PROJECTLSMASSCONSTRAINT3D)
#endif

extern "C"
{
    void RELAXATION_LS_1ST_ORDER_FC(double* U,
                                    const int& U_gcw,
                                    const double* V,
                                    const int& V_gcw,
                                    const double* dt,
                                    const int& dt_gcw,
                                    const int& ilower0,
                                    const int& iupper0,
                                    const int& ilower1,
                                    const int& iupper1,
#if (NDIM == 3)
                                    const int& ilower2,
                                    const int& iupper2,
#endif
                                    const double* dx,
                                    const int& dir);

    void RELAXATION_LS_3RD_ORDER_ENO_FC(double* U,
                                        const int& U_gcw,
                                        const double* V,
                                        const int& V_gcw,
                                        const double* dt,
                                        const int& dt_gcw,
                                        const int& ilower0,
                                        const int& iupper0,
                                        const int& ilower1,
                                        const int& iupper1,
#if (NDIM == 3)
                                        const int& ilower2,
                                        const int& iupper2,
#endif
                                        const double* dx,
                                        const int& dir,
                                        const int& use_subcell,
                                        const int& use_sign_fix);

    void GODUNOV_HAMILTONIAN_3RD_ORDER_ENO_FC(double* U,
                                              const int& U_gcw,
                                              const double* V,
                                              const int& V_gcw,
                                              const int& ilower0,
                                              const int& iupper0,
                                              const int& ilower1,
                                              const int& iupper1,
#if (NDIM == 3)
                                              const int& ilower2,
                                              const int& iupper2,
#endif
                                              const double* dx,
                                              const int& use_subcell);

    void RELAXATION_LS_3RD_ORDER_WENO_FC(double* U,
                                         const int& U_gcw,
                                         const double* V,
                                         const int& V_gcw,
                                         const double* dt,
                                         const int& dt_gcw,
                                         const int& ilower0,
                                         const int& iupper0,
                                         const int& ilower1,
                                         const int& iupper1,
#if (NDIM == 3)
                                         const int& ilower2,
                                         const int& iupper2,
#endif
                                         const double* dx,
                                         const int& dir,
                                         const int& use_subcell,
                                         const int& use_sign_fix);

    void GODUNOV_HAMILTONIAN_3RD_ORDER_WENO_FC(double* U,
                                               const int& U_gcw,
                                               const double* V,
                                               const int& V_gcw,
                                               const int& ilower0,
                                               const int& iupper0,
                                               const int& ilower1,
                                               const int& iupper1,
#if (NDIM == 3)
                                               const int& ilower2,
                                               const int& iupper2,
#endif
                                               const double* dx,
                                               const int& use_subcell);

    void RELAXATION_LS_5TH_ORDER_WENO_FC(double* U,
                                         const int& U_gcw,
                                         const double* V,
                                         const int& V_gcw,
                                         const double* dt,
                                         const int& dt_gcw,
                                         const int& ilower0,
                                         const int& iupper0,
                                         const int& ilower1,
                                         const int& iupper1,
#if (NDIM == 3)
                                         const int& ilower2,
                                         const int& iupper2,
#endif
                                         const double* dx,
                                         const int& dir,
                                         const int& use_sign_fix);

    void GODUNOV_HAMILTONIAN_5TH_ORDER_WENO_FC(double* U,
                                               const int& U_gcw,
                                               const double* V,
                                               const int& V_gcw,
                                               const int& ilower0,
                                               const int& iupper0,
                                               const int& ilower1,
                                               const int& iupper1,
#if (NDIM == 3)
                                               const int& ilower2,
                                               const int& iupper2,
#endif
                                               const double* dx);

    void PROJECT_LS_MASS_CONSTRAINT_FC(double* U,
                                       const int& U_gcw,
                                       const double* C,
                                       const int& C_gcw,
                                       const double* V,
                                       const int& V_gcw,
                                       const double* H,
                                       const int& H_gcw,
                                       const int& ilower0,
                                       const int& iupper0,
                                       const int& ilower1,
                                       const int& iupper1,
#if (NDIM == 3)
                                       const int& ilower2,
                                       const int& iupper2,
#endif
                                       const double* dx);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RelaxationLSMethod::RelaxationLSMethod(std::string object_name, Pointer<Database> db, bool register_for_restart)
    : LSInitStrategy(std::move(object_name), register_for_restart)
{
    // Some default values.
    d_ls_order = THIRD_ORDER_ENO_LS;

    // Get any additional or overwrite base class options.
    if (d_registered_for_restart) getFromRestart();
    if (!db.isNull()) getFromInput(db);

    return;
} // RelaxationLSMethod

void
RelaxationLSMethod::initializeLSData(int D_idx,
                                     Pointer<HierarchyMathOps> hier_math_ops,
                                     int integrator_step,
                                     double time,
                                     bool initial_time)
{
    const bool initialize_ls =
        d_reinitialize_ls || initial_time || (d_reinit_interval && integrator_step % d_reinit_interval == 0);
    if (!initialize_ls) return;

    if (d_apply_mass_constraint && initial_time)
    {
        TBOX_WARNING(d_object_name << "::initializeLSData():\n"
                                   << " Mass constraint is automatically turned "
                                      "off for initial hierarchy time"
                                   << std::endl);
    }
    const bool constrain_ls_mass = (d_apply_mass_constraint && !initial_time);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > data_var;
    var_db->mapIndexToVariable(D_idx, data_var);
    Pointer<CellVariable<NDIM, double> > D_var = data_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!D_var.isNull());
#endif

    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Create a temporary variable to hold previous iteration values with
    // appropriate ghost cell width since it is not guaranteed that D_idx will
    // have proper ghost cell width.
    IntVector<NDIM> cell_ghosts;
    IntVector<NDIM> no_ghosts = 0;
    if (d_ls_order == FIRST_ORDER_LS)
    {
        TBOX_WARNING(d_object_name << "::initializeLSData():\n"
                                   << " First order relxation is known to cause significant "
                                      "interface volume loss \n"
                                   << " consider trying THIRD_ORDER_ENO or THIRD_ORDER_WENO." << std::endl);
        cell_ghosts = std::max(1, d_D_gcw);
    }
    else if (d_ls_order == THIRD_ORDER_ENO_LS || d_ls_order == THIRD_ORDER_WENO_LS)
    {
        cell_ghosts = std::max(2, d_D_gcw);
    }
    else if (d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        cell_ghosts = std::max(3, d_D_gcw);
    }
    else
    {
        TBOX_ERROR("RelaxationLSMethod does not support " << enum_to_string(d_ls_order) << std::endl);
    }
    const int D_scratch_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);
    const int D_iter_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::ITER"), cell_ghosts);
    const int D_init_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::INIT"), cell_ghosts);
    const int D_copy_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::COPY"), cell_ghosts);
    const int H_init_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::H_INIT"), cell_ghosts);
    const int H_scratch_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::H_SCRATCH"), no_ghosts);

    // Time step size variable
    const int dt_idx = var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::DT"), no_ghosts);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_scratch_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_iter_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_init_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(D_copy_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(H_init_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(H_scratch_idx, time);
        hierarchy->getPatchLevel(ln)->allocatePatchData(dt_idx, time);
    }

    // First, fill cells with some positive/negative values
    // away from the interface and actual distance value near the interface.
    for (unsigned k = 0; k < d_locate_interface_fcns.size(); ++k)
    {
        (*d_locate_interface_fcns[k])(D_scratch_idx, hier_math_ops, time, initial_time, d_locate_interface_fcns_ctx[k]);
    }

    // Set hierarchy objects.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent D_transaction(
        D_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, d_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> D_fill_op = new HierarchyGhostCellInterpolation();
    InterpolationTransactionComponent H_transcation(
        H_init_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, nullptr);
    Pointer<HierarchyGhostCellInterpolation> H_fill_op = new HierarchyGhostCellInterpolation();
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);

    // Carry out relaxation
    double diff_L2_norm = 1.0e12;
    int outer_iter = 0;
    const int cc_wgt_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Initialize operator states
    D_fill_op->initializeOperatorState(D_transaction, hierarchy);
    H_fill_op->initializeOperatorState(H_transcation, hierarchy);

    // Copy initial condition, including ghost cells
    D_fill_op->fillData(time);
    hier_cc_data_ops.copyData(D_init_idx, D_scratch_idx, /*interior_only*/ false);

    // Compute and store |grad phi_0| to apply the mass constraint
    if (constrain_ls_mass)
    {
        computeHamiltonian(hier_math_ops, H_init_idx, D_init_idx);
        H_fill_op->fillData(time);
    }

    while (diff_L2_norm > d_abs_tol && outer_iter < d_max_its)
    {
        // Save a copy of previous iterate before modifying the level set function.
        hier_cc_data_ops.copyData(D_iter_idx, D_scratch_idx);
        relax(D_fill_op, hier_math_ops, D_scratch_idx, D_init_idx, D_copy_idx, dt_idx, outer_iter, time);
        hier_cc_data_ops.linearSum(D_scratch_idx, d_alpha, D_scratch_idx, 1.0 - d_alpha, D_iter_idx);

        if (constrain_ls_mass)
        {
            D_fill_op->fillData(time);
            hier_cc_data_ops.copyData(D_copy_idx,
                                      D_scratch_idx,
                                      /*interior_only*/ false);
            applyMassConstraint(hier_math_ops, D_scratch_idx, D_copy_idx, D_init_idx, H_init_idx);
        }

        if (d_apply_vol_redistribution)
        {
            D_fill_op->fillData(time);
            computeHamiltonian(hier_math_ops, H_scratch_idx, D_scratch_idx);
            applyVolumeRedistribution(
                hier_math_ops, /*lambda*/ D_copy_idx, D_scratch_idx, D_init_idx, H_scratch_idx, dt_idx);
        }

        // Compute error, but copy previous iteration beforehand
        hier_cc_data_ops.copyData(D_copy_idx, D_iter_idx);
        hier_cc_data_ops.axmy(D_iter_idx, 1.0, D_iter_idx, D_scratch_idx);
        diff_L2_norm = hier_cc_data_ops.L2Norm(D_iter_idx, cc_wgt_idx);

        outer_iter += 1;

        if (d_enable_logging)
        {
            // Compute difference between |grad phi| and 1
            // Reuse Hamiltonian from the volume redistribution algorithm.
            if (!d_apply_vol_redistribution)
            {
                D_fill_op->fillData(time);
                computeHamiltonian(hier_math_ops, H_scratch_idx, D_scratch_idx);
            }
            hier_cc_data_ops.addScalar(H_scratch_idx, H_scratch_idx, -1.0);
            const double grad_norm = hier_cc_data_ops.L2Norm(H_scratch_idx, cc_wgt_idx);

            plog << d_object_name << "::initializeLSData(): After iteration # " << outer_iter << std::endl;
            plog << d_object_name << "::initializeLSData(): L2-norm between successive iterations = " << diff_L2_norm
                 << std::endl;
            plog << d_object_name << "::initializeLSData(): L2-Norm || |grad phi| - 1|| = " << grad_norm << std::endl;
        }

        if (diff_L2_norm <= d_abs_tol && d_enable_logging)
        {
            plog << d_object_name << "::initializeLSData(): Relaxation converged for entire domain" << std::endl;
        }
    }

    if (outer_iter >= d_max_its)
    {
        if (d_enable_logging)
        {
            plog << d_object_name << "::initializeLSData(): Reached maximum allowable outer iterations" << std::endl;
            plog << d_object_name << "::initializeLSData(): ||distance_new - distance_old||_2 = " << diff_L2_norm
                 << std::endl;
        }
    }

    // Copy signed distance into supplied patch data index
    hier_cc_data_ops.copyData(D_idx, D_scratch_idx);

    // Deallocate the temporary variable.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_scratch_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_iter_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_init_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(D_copy_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(H_init_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(H_scratch_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(dt_idx);
    }
    var_db->removePatchDataIndex(D_scratch_idx);
    var_db->removePatchDataIndex(D_iter_idx);
    var_db->removePatchDataIndex(D_init_idx);
    var_db->removePatchDataIndex(D_copy_idx);
    var_db->removePatchDataIndex(H_init_idx);
    var_db->removePatchDataIndex(H_scratch_idx);
    var_db->removePatchDataIndex(dt_idx);

    // Indicate that the LS has been initialized.
    d_reinitialize_ls = false;

    return;
} // initializeLSData

void
RelaxationLSMethod::setApplyMassConstraint(bool apply_mass_constraint)
{
    d_apply_mass_constraint = apply_mass_constraint;
    return;
} // setApplyMassConstraint

void
RelaxationLSMethod::setApplySubcellFix(bool apply_subcell_fix)
{
    d_apply_subcell_fix = apply_subcell_fix;
    return;
} // setApplySubcellFix

void
RelaxationLSMethod::setApplySignFix(bool apply_sign_fix)
{
    d_apply_sign_fix = apply_sign_fix;
    return;
} // setApplySignFix

void
RelaxationLSMethod::setLSGhostCellWidth(int D_gcw)
{
    d_D_gcw = D_gcw;
    return;
} // setLSGhostCellWidth

void
RelaxationLSMethod::setApplyVolumeRedistribution(bool apply_volume_redistribution)
{
    d_apply_vol_redistribution = apply_volume_redistribution;
    return;
} // setApplyVolumeRedistribution

/////////////////////////////// PRIVATE //////////////////////////////////////

void
RelaxationLSMethod::relax(Pointer<HierarchyGhostCellInterpolation> D_fill_op,
                          Pointer<HierarchyMathOps> hier_math_ops,
                          int dist_idx,
                          int dist_init_idx,
                          int dist_copy_idx,
                          int dt_idx,
                          const int iter,
                          const double time) const
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    if (d_ls_ts == GAUSS_SEIDEL_PSEUDO_TS || d_ls_ts == RK1_TS)
    {
        // Fill the ghost cells of dist_idx
        D_fill_op->fillData(time);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > dist_data = patch->getPatchData(dist_idx);
                const Pointer<CellData<NDIM, double> > dist_init_data = patch->getPatchData(dist_init_idx);
                Pointer<CellData<NDIM, double> > dt_data = patch->getPatchData(dt_idx);
                relax(dist_data, dist_init_data, dt_data, patch, iter);
            }
        }
    }
    else if (d_ls_ts == TVD_RK2_TS)
    {
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(hierarchy, coarsest_ln, finest_ln);

        for (int stage = 0; stage < 2; ++stage)
        {
            D_fill_op->fillData(time);
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > dist_data = patch->getPatchData(dist_idx);
                    const Pointer<CellData<NDIM, double> > dist_init_data = patch->getPatchData(dist_init_idx);
                    Pointer<CellData<NDIM, double> > dt_data = patch->getPatchData(dt_idx);
                    relax(dist_data, dist_init_data, dt_data, patch, iter);
                }
            }
            if (stage == 0) hier_cc_data_ops.copyData(dist_copy_idx, dist_idx);
        }
        hier_cc_data_ops.linearSum(dist_idx, 0.5, dist_idx, 0.5, dist_copy_idx);
    }
    else
    {
        TBOX_ERROR("Unknown time stepping scheme for the RelaxationLSMethod: " << enum_to_string(d_ls_ts) << std::endl);
    }
    return;

} // relax

void
RelaxationLSMethod::relax(Pointer<CellData<NDIM, double> > dist_data,
                          const Pointer<CellData<NDIM, double> > dist_init_data,
                          Pointer<CellData<NDIM, double> > dt_data,
                          const Pointer<Patch<NDIM> > patch,
                          const int iter) const
{
    double* const D = dist_data->getPointer(0);
    const double* const P = dist_init_data->getPointer(0);
    double* const dt = dt_data->getPointer(0);
    const int D_ghosts = (dist_data->getGhostCellWidth()).max();
    const int P_ghosts = (dist_init_data->getGhostCellWidth()).max();
    const int dt_ghosts = (dt_data->getGhostCellWidth()).max();

#if !defined(NDEBUG)
    TBOX_ASSERT(dist_data->getDepth() == 1);
    TBOX_ASSERT(dist_init_data->getDepth() == 1);
    if (d_ls_order == FIRST_ORDER_LS)
    {
        TBOX_ASSERT(D_ghosts >= 1);
        TBOX_ASSERT(P_ghosts >= 1);
    }
    if (d_ls_order == THIRD_ORDER_ENO_LS || d_ls_order == THIRD_ORDER_WENO_LS)
    {
        TBOX_ASSERT(D_ghosts >= 2);
        TBOX_ASSERT(P_ghosts >= 2);
    }
    if (d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        TBOX_ASSERT(D_ghosts >= 3);
    }
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Here we are solving Hamilton-Jacobi equation of the type
    //
    //         d (phi)/dt + sign(phi_0)(|grad phi| - 1) = 0
    //
    // till steady state. We use an explicit time stepping scheme to solve the
    // above equation, which reads as
    //
    //         phi_{ijk}^{n+1} =  phi_{ijk}^{n} - dt_{ijk} sign(phi_0)(|grad phi|^n - 1)
    //
    // Though the time stepping scheme is fictitious, but it can affect the
    // spatial order of accuracy. The above update (till steady state) can also be
    // viewed as Gauss-Seidel (GS) update at grid location (i,j,k). If we choose
    // to use the GS scheme, we use a roster-scan algorithm to permute (i,j,k)
    // indices at each iteration.

    // Get the direction of sweeping (alternates according to iteration number)
    // For RK integrators the sweep direction does not matter so we select the
    // the usual (i,j,k) direction at each iteration. For GS updates,
    // we permute (i,j,k) with each iteration.
    int dir = 0;
    if (d_ls_ts == GAUSS_SEIDEL_PSEUDO_TS)
    {
        const int num_dirs = (NDIM < 3) ? 4 : 8;
        dir = iter % num_dirs;
    }

    if (d_ls_order == FIRST_ORDER_LS)
    {
        RELAXATION_LS_1ST_ORDER_FC(D,
                                   D_ghosts,
                                   P,
                                   P_ghosts,
                                   dt,
                                   dt_ghosts,
                                   patch_box.lower(0),
                                   patch_box.upper(0),
                                   patch_box.lower(1),
                                   patch_box.upper(1),
#if (NDIM == 3)
                                   patch_box.lower(2),
                                   patch_box.upper(2),
#endif
                                   dx,
                                   dir);
    }
    else if (d_ls_order == THIRD_ORDER_ENO_LS)
    {
        RELAXATION_LS_3RD_ORDER_ENO_FC(D,
                                       D_ghosts,
                                       P,
                                       P_ghosts,
                                       dt,
                                       dt_ghosts,
                                       patch_box.lower(0),
                                       patch_box.upper(0),
                                       patch_box.lower(1),
                                       patch_box.upper(1),
#if (NDIM == 3)
                                       patch_box.lower(2),
                                       patch_box.upper(2),
#endif
                                       dx,
                                       dir,
                                       static_cast<int>(d_apply_subcell_fix),
                                       static_cast<int>(d_apply_sign_fix));
    }
    else if (d_ls_order == THIRD_ORDER_WENO_LS)
    {
        RELAXATION_LS_3RD_ORDER_WENO_FC(D,
                                        D_ghosts,
                                        P,
                                        P_ghosts,
                                        dt,
                                        dt_ghosts,
                                        patch_box.lower(0),
                                        patch_box.upper(0),
                                        patch_box.lower(1),
                                        patch_box.upper(1),
#if (NDIM == 3)
                                        patch_box.lower(2),
                                        patch_box.upper(2),
#endif
                                        dx,
                                        dir,
                                        static_cast<int>(d_apply_subcell_fix),
                                        static_cast<int>(d_apply_sign_fix));
    }
    else if (d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        RELAXATION_LS_5TH_ORDER_WENO_FC(D,
                                        D_ghosts,
                                        P,
                                        P_ghosts,
                                        dt,
                                        dt_ghosts,
                                        patch_box.lower(0),
                                        patch_box.upper(0),
                                        patch_box.lower(1),
                                        patch_box.upper(1),
#if (NDIM == 3)
                                        patch_box.lower(2),
                                        patch_box.upper(2),
#endif
                                        dx,
                                        dir,
                                        static_cast<int>(d_apply_sign_fix));
    }
    else
    {
        TBOX_ERROR("RelaxationLSMethod does not support " << enum_to_string(d_ls_order) << std::endl);
    }

    return;
} // relax

void
RelaxationLSMethod::computeHamiltonian(Pointer<HierarchyMathOps> hier_math_ops, int ham_idx, int dist_idx) const
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > ham_data = patch->getPatchData(ham_idx);
            const Pointer<CellData<NDIM, double> > dist_data = patch->getPatchData(dist_idx);
            computeHamiltonian(ham_data, dist_data, patch);
        }
    }
    return;

} // computeHamiltonian

void
RelaxationLSMethod::computeHamiltonian(Pointer<CellData<NDIM, double> > ham_data,
                                       const Pointer<CellData<NDIM, double> > dist_data,
                                       const Pointer<Patch<NDIM> > patch) const
{
    double* const H = ham_data->getPointer(0);
    const double* const P = dist_data->getPointer(0);
    const int H_ghosts = (ham_data->getGhostCellWidth()).max();
    const int P_ghosts = (dist_data->getGhostCellWidth()).max();

#if !defined(NDEBUG)
    TBOX_ASSERT(ham_data->getDepth() == 1);
    TBOX_ASSERT(dist_data->getDepth() == 1);
    if (d_ls_order == THIRD_ORDER_ENO_LS || d_ls_order == THIRD_ORDER_WENO_LS)
    {
        TBOX_ASSERT(P_ghosts >= 2);
    }
    else if (d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        TBOX_ASSERT(P_ghosts >= 3);
    }
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    if (d_ls_order == THIRD_ORDER_ENO_LS)
    {
        GODUNOV_HAMILTONIAN_3RD_ORDER_ENO_FC(H,
                                             H_ghosts,
                                             P,
                                             P_ghosts,
                                             patch_box.lower(0),
                                             patch_box.upper(0),
                                             patch_box.lower(1),
                                             patch_box.upper(1),
#if (NDIM == 3)
                                             patch_box.lower(2),
                                             patch_box.upper(2),
#endif
                                             dx,
                                             static_cast<int>(d_apply_subcell_fix));
    }
    else if (d_ls_order == THIRD_ORDER_WENO_LS)
    {
        GODUNOV_HAMILTONIAN_3RD_ORDER_WENO_FC(H,
                                              H_ghosts,
                                              P,
                                              P_ghosts,
                                              patch_box.lower(0),
                                              patch_box.upper(0),
                                              patch_box.lower(1),
                                              patch_box.upper(1),
#if (NDIM == 3)
                                              patch_box.lower(2),
                                              patch_box.upper(2),
#endif
                                              dx,
                                              static_cast<int>(d_apply_subcell_fix));
    }
    else if (d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        GODUNOV_HAMILTONIAN_5TH_ORDER_WENO_FC(H,
                                              H_ghosts,
                                              P,
                                              P_ghosts,
                                              patch_box.lower(0),
                                              patch_box.upper(0),
                                              patch_box.lower(1),
                                              patch_box.upper(1),
#if (NDIM == 3)
                                              patch_box.lower(2),
                                              patch_box.upper(2),
#endif
                                              dx);
    }
    else
    {
        TBOX_ERROR("RelaxationLSMethod does not support mass constraint for " << enum_to_string(d_ls_order)
                                                                              << std::endl);
    }

    return;
} // computeHamiltonian

void
RelaxationLSMethod::applyMassConstraint(Pointer<HierarchyMathOps> hier_math_ops,
                                        int dist_idx,
                                        int dist_copy_idx,
                                        int dist_init_idx,
                                        int ham_init_idx) const
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, double> > dist_data = patch->getPatchData(dist_idx);
            const Pointer<CellData<NDIM, double> > dist_copy_data = patch->getPatchData(dist_copy_idx);
            const Pointer<CellData<NDIM, double> > dist_init_data = patch->getPatchData(dist_init_idx);
            const Pointer<CellData<NDIM, double> > ham_init_data = patch->getPatchData(ham_init_idx);
            applyMassConstraint(dist_data, dist_copy_data, dist_init_data, ham_init_data, patch);
        }
    }
    return;

} // applyMassConstraint

void
RelaxationLSMethod::applyMassConstraint(Pointer<CellData<NDIM, double> > dist_data,
                                        const Pointer<CellData<NDIM, double> > dist_copy_data,
                                        const Pointer<CellData<NDIM, double> > dist_init_data,
                                        const Pointer<CellData<NDIM, double> > ham_init_data,
                                        const Pointer<Patch<NDIM> > patch) const
{
    double* const U = dist_data->getPointer(0);
    const double* const C = dist_copy_data->getPointer(0);
    const double* const V = dist_init_data->getPointer(0);
    const double* const H = ham_init_data->getPointer(0);
    const int U_ghosts = (dist_data->getGhostCellWidth()).max();
    const int C_ghosts = (dist_copy_data->getGhostCellWidth()).max();
    const int V_ghosts = (dist_init_data->getGhostCellWidth()).max();
    const int H_ghosts = (ham_init_data->getGhostCellWidth()).max();

#if !defined(NDEBUG)
    TBOX_ASSERT(dist_data->getDepth() == 1);
    TBOX_ASSERT(dist_copy_data->getDepth() == 1);
    TBOX_ASSERT(ham_init_data->getDepth() == 1);
    TBOX_ASSERT(dist_init_data->getDepth() == 1);
    if (d_ls_order == THIRD_ORDER_ENO_LS || d_ls_order == THIRD_ORDER_WENO_LS)
    {
        TBOX_ASSERT(U_ghosts >= 2);
        TBOX_ASSERT(C_ghosts >= 2);
        TBOX_ASSERT(H_ghosts >= 2);
    }
    if (d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        TBOX_ASSERT(U_ghosts >= 3);
        TBOX_ASSERT(C_ghosts >= 3);
        TBOX_ASSERT(H_ghosts >= 3);
    }
#endif

    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    if (d_ls_order == THIRD_ORDER_ENO_LS || d_ls_order == THIRD_ORDER_WENO_LS || d_ls_order == FIFTH_ORDER_WENO_LS)
    {
        PROJECT_LS_MASS_CONSTRAINT_FC(U,
                                      U_ghosts,
                                      C,
                                      C_ghosts,
                                      V,
                                      V_ghosts,
                                      H,
                                      H_ghosts,
                                      patch_box.lower(0),
                                      patch_box.upper(0),
                                      patch_box.lower(1),
                                      patch_box.upper(1),
#if (NDIM == 3)
                                      patch_box.lower(2),
                                      patch_box.upper(2),
#endif
                                      dx);
    }
    else
    {
        TBOX_ERROR("RelaxationLSMethod does not support mass constraint for " << enum_to_string(d_ls_order)
                                                                              << std::endl);
    }

    return;
} // applyMassConstraint

void
RelaxationLSMethod::applyVolumeRedistribution(Pointer<HierarchyMathOps> hier_math_ops,
                                              int lambda_idx,
                                              int dist_idx,
                                              int dist_init_idx,
                                              int ham_idx,
                                              int dt_idx) const
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    auto signof = [](const double x) { return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); };
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    double source_integral = 0.0;
    double weight_integral = 0.0;

    // First compute the Lagrange multiplier lambda
    // that enforces integral of H(phi) to be constant at all times
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > wgt_data = patch->getPatchData(wgt_cc_idx);
            const Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(dist_idx);
            const Pointer<CellData<NDIM, double> > phi_init_data = patch->getPatchData(dist_init_idx);
            const Pointer<CellData<NDIM, double> > ham_data = patch->getPatchData(ham_idx);
            Pointer<CellData<NDIM, double> > lambda_data = patch->getPatchData(lambda_idx);

            // Get grid spacing information
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d) cell_size *= patch_dx[d];
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double num_cells = 1.0;
            const double alpha = num_cells * cell_size;

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci);
                const double phi0 = (*phi_init_data)(ci);
                const double grad_phi = (*ham_data)(ci);
                const double dv = (*wgt_data)(ci);

                // smoothened delta and Heaviside functions
                const double h_phi = IBTK::smooth_heaviside(phi, alpha);
                const double h_prime = IBTK::smooth_delta(phi, alpha);

                (*lambda_data)(ci) = h_phi * (1.0 - h_phi);

                source_integral += h_prime * signof(phi0) * (grad_phi - 1.0) * dv;
                weight_integral += h_prime * h_phi * (1.0 - h_phi) * dv;
            }
        }
    }
    const double S = IBTK_MPI::sumReduction(source_integral);
    const double W = IBTK_MPI::sumReduction(weight_integral);

    // Correct the distance function using lambda
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(dist_idx);
            const Pointer<CellData<NDIM, double> > lambda_data = patch->getPatchData(lambda_idx);
            const Pointer<CellData<NDIM, double> > dt_data = patch->getPatchData(dt_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double dt = (*dt_data)(ci);

                if (dt <= 0.0) continue;

                (*phi_data)(ci) = (*phi_data)(ci) + (*dt_data)(ci) * (S / W) * (*lambda_data)(ci);
            }
        }
    }

    return;
} // applyVolumeRedistribution

void
RelaxationLSMethod::getFromInput(Pointer<Database> input_db)
{
    std::string ls_order = "THIRD_ORDER_ENO";
    ls_order = input_db->getStringWithDefault("order", ls_order);
    d_ls_order = string_to_enum<LevelSetOrder>(ls_order);

    std::string ls_ts = "GAUSS_SEIDEL";
    ls_ts = input_db->getStringWithDefault("time_stepping_scheme", ls_ts);
    d_ls_ts = string_to_enum<LevelSetTimeStepping>(ls_ts);

    d_max_its = input_db->getIntegerWithDefault("max_iterations", d_max_its);
    d_max_its = input_db->getIntegerWithDefault("max_its", d_max_its);

    d_abs_tol = input_db->getDoubleWithDefault("abs_tol", d_abs_tol);

    d_reinit_interval = input_db->getIntegerWithDefault("reinit_interval", d_reinit_interval);

    d_enable_logging = input_db->getBoolWithDefault("enable_logging", d_enable_logging);

    d_apply_mass_constraint = input_db->getBoolWithDefault("apply_mass_constraint", d_apply_mass_constraint);
    d_apply_subcell_fix = input_db->getBoolWithDefault("apply_subcell_fix", d_apply_subcell_fix);

    d_apply_sign_fix = input_db->getBoolWithDefault("apply_sign_fix", d_apply_sign_fix);
    d_apply_sign_fix = input_db->getBoolWithDefault("apply_sgn_fix", d_apply_sign_fix);

    d_apply_vol_redistribution =
        input_db->getBoolWithDefault("apply_volume_redistribution", d_apply_vol_redistribution);

    d_D_gcw = input_db->getIntegerWithDefault("ghost_cell_width", d_D_gcw);

    return;
} // getFromInput

void
RelaxationLSMethod::getFromRestart()
{
    // intentionally left-blank.
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
