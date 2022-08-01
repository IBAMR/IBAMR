// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/CarmanKozenyDragForce.h"
#include "ibamr/IBHydrodynamicSurfaceForceEvaluator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"

#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "Eigen/src/Geometry/Quaternion.h"
#include "Eigen/src/LU/FullPivLU.h"

#include <cmath>
#include <string>
#include <utility>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC //////////////////////////////////////
CarmanKozenyDragForce::CarmanKozenyDragForce(std::string object_name,
                                             Pointer<CellVariable<NDIM, double> > H_var,
                                             Pointer<CellVariable<NDIM, double> > lf_var,
                                             Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                             Pointer<INSVCStaggeredHierarchyIntegrator> fluid_solver,
                                             Pointer<Database> input_db,
                                             bool register_for_restart)
    : BrinkmanPenalizationStrategy(std::move(object_name), register_for_restart),
      d_adv_diff_solver(adv_diff_solver),
      d_fluid_solver(fluid_solver),
      d_H_var(H_var),
      d_lf_var(lf_var)
{
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // CarmanKozenyDragForce

void
CarmanKozenyDragForce::preprocessComputeBrinkmanPenalization(double /*current_time*/,
                                                             double /*new_time*/,
                                                             int /*num_cycles*/)
{
    //    d_trans_vel_new = d_trans_vel_current;
    //    d_rot_vel_new = d_rot_vel_current;
    //    d_center_of_mass_new = d_center_of_mass_current;
    //    d_quaternion_new = d_quaternion_current;

    return;
} // preprocessComputeBrinkmanPenalization

void
CarmanKozenyDragForce::computeBrinkmanVelocity(int u_idx, double time, int cycle_num)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(time, d_new_time));
#endif
    const double dt = d_new_time - d_current_time;

    // Get the interpolated density variable
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    // Ghost fill the heaviside and liquid fraction values.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getScratchContext());

    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getScratchContext());

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> phi_transaction_comps(2);
    phi_transaction_comps[0] = InterpolationTransactionComponent(H_scratch_idx,
                                                                 H_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_adv_diff_solver->getPhysicalBcCoefs(d_H_var));

    // Using levelset bc for liquid fraction.
    phi_transaction_comps[1] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                 lf_new_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_adv_diff_solver->getPhysicalBcCoefs(d_H_var));

    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    // Set the rigid body velocity in u_idx
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_scratch_idx);
        Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double H_lower = (*H_data)(s_i.toCell(0));
                const double H_upper = (*H_data)(s_i.toCell(1));
                const double H = 0.5 * (H_lower + H_upper);

                const double lf_lower = (*lf_data)(s_i.toCell(0));
                const double lf_upper = (*lf_data)(s_i.toCell(1));
                const double liquid_fraction = 0.5 * (lf_lower + lf_upper);

                const double alpha_s = H * (1.0 - liquid_fraction);
                const double penalty = (*rho_data)(s_i) / dt;
                const double solid_velocity = 0.0;

                (*u_data)(s_i) = solid_velocity * penalty * alpha_s * alpha_s / (std::pow(1.0 - alpha_s, 3.0) + d_ed);
            }
        }
    }

    return;
} // computeBrinkmanVelocity

void
CarmanKozenyDragForce::demarcateBrinkmanZone(int u_idx, double time, int /*cycle_num*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(MathUtilities<double>::equalEps(time, d_new_time));
#else
    NULL_USE(time);
#endif

    const double dt = d_new_time - d_current_time;

    // Get the interpolated density variable
    const int rho_ins_idx = d_fluid_solver->getLinearOperatorRhoPatchDataIndex();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    const int H_scratch_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getScratchContext());

    const int lf_new_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    const int lf_scratch_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getScratchContext());

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_adv_diff_solver->getPatchHierarchy();

    // Ghost fill the level set values.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> phi_transaction_comps(2);
    phi_transaction_comps[0] = InterpolationTransactionComponent(H_scratch_idx,
                                                                 H_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_adv_diff_solver->getPhysicalBcCoefs(d_H_var));

    // Using levelset bc for liquid fraction.
    phi_transaction_comps[1] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                 lf_new_idx,
                                                                 "CONSERVATIVE_LINEAR_REFINE",
                                                                 false,
                                                                 "CONSERVATIVE_COARSEN",
                                                                 "LINEAR",
                                                                 false,
                                                                 d_adv_diff_solver->getPhysicalBcCoefs(d_H_var));

    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_dx = patch_geom->getDx();
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_scratch_idx);
        Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_ins_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double H_lower = (*H_data)(s_i.toCell(0));
                const double H_upper = (*H_data)(s_i.toCell(1));
                const double H = 0.5 * (H_lower + H_upper);

                const double lf_lower = (*lf_data)(s_i.toCell(0));
                const double lf_upper = (*lf_data)(s_i.toCell(1));
                const double liquid_fraction = 0.5 * (lf_lower + lf_upper);

                const double alpha_s = H * (1.0 - liquid_fraction);
                const double penalty = (*rho_data)(s_i) / dt;
                (*u_data)(s_i) = penalty * alpha_s * alpha_s / (std::pow(1.0 - alpha_s, 3.0) + d_ed);
            }
        }
    }

    return;
} // demarcateBrinkmanZone

void
CarmanKozenyDragForce::postprocessComputeBrinkmanPenalization(double current_time, double new_time, int num_cycles)
{
    BrinkmanPenalizationStrategy::postprocessComputeBrinkmanPenalization(current_time, new_time, num_cycles);

    return;
} // postprocessComputeBrinkmanPenalization

void
CarmanKozenyDragForce::putToDatabase(Pointer<Database> db)
{
    return;
} // postprocessComputeBrinkmanVelocity

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CarmanKozenyDragForce::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
    if (input_db->keyExists("num_interface_cells"))
    {
        d_num_interface_cells = input_db->getDouble("num_interface_cells");
    }
    return;
} // getFromInput

void
CarmanKozenyDragForce::getFromRestart()
{
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
