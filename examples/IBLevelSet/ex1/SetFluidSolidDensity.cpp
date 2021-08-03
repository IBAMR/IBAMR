// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/HierarchyMathOps.h>

#include "SetFluidSolidDensity.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidSolidDensityCallbackFunction(int rho_idx,
                                         Pointer<Variable<NDIM> > rho_var,
                                         Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         const int cycle_num,
                                         const double time,
                                         const double current_time,
                                         const double new_time,
                                         void* ctx)
{
    // Set the density from the level set information
    static SetFluidSolidDensity* ptr_SetFluidSolidDensity = static_cast<SetFluidSolidDensity*>(ctx);
    ptr_SetFluidSolidDensity->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidSolidDensityCallBackFunction

// Various options to setting side-centered densities

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidSolidDensity::SetFluidSolidDensity(const std::string& object_name,
                                           Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                           Pointer<CellVariable<NDIM, double> > ls_solid_var,
                                           const double rho_fluid,
                                           const double rho_solid,
                                           const int /*ls_reinit_interval*/,
                                           const double num_solid_interface_cells)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_solid_var(ls_solid_var),
      d_rho_fluid(rho_fluid),
      d_rho_solid(rho_solid),
      d_num_solid_interface_cells(num_solid_interface_cells)
{
    // intentionally left blank
    return;
} // SetFluidGasSolidDensity

void
SetFluidSolidDensity::setDensityPatchData(int rho_idx,
                                          Pointer<Variable<NDIM> > rho_var,
                                          Pointer<HierarchyMathOps> hier_math_ops,
                                          const int /*cycle_num*/,
                                          const double time,
                                          const double current_time,
                                          const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int ls_solid_idx = -1, ls_gas_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        ls_solid_idx = var_db->mapVariableAndContextToIndex(d_ls_solid_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Set the density based on the cell centered level set
    Pointer<CellVariable<NDIM, double> > rho_cc_var = rho_var;

    if (rho_cc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                const double alpha = d_num_solid_interface_cells * patch_dx[0];

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
                const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_idx);
                Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                // Li et al, 2015
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi_s = (*ls_solid_data)(ci);
                    double Hphi_s;

                    if (phi_s < -alpha)
                    {
                        Hphi_s = 0.0;
                    }
                    else if (std::abs(phi_s) <= alpha)
                    {
                        Hphi_s = 0.5 + 0.5 * phi_s / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_s / alpha);
                    }
                    else
                    {
                        Hphi_s = 1.0;
                    }

                    // Next, set the density of the solid phase
                    (*rho_data)(ci) = d_rho_fluid * Hphi_s + (1.0 - Hphi_s) * d_rho_solid;
                }
            }
        }
    }

    // Setting side centered density directly
    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        // Note, this method requires ghost cells to be filled for the level set variable
        RobinBcCoefStrategy<NDIM>* ls_solid_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_ls_solid_var).front();
        IntVector<NDIM> cell_ghosts = 1;
        const int ls_solid_scratch_idx = var_db->registerVariableAndContext(
            d_ls_solid_var, var_db->getContext(d_object_name + "::SOLID::SCRATCH"), cell_ghosts);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_solid_scratch_idx, time);
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ls_transaction_comps(1);
        ls_transaction_comps[0] = InterpolationTransactionComponent(ls_solid_scratch_idx,
                                                                    ls_solid_idx,
                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                    false,
                                                                    "CONSERVATIVE_COARSEN",
                                                                    "LINEAR",
                                                                    false,
                                                                    ls_solid_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(ls_transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                // Compute the indicators for both level sets
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double phi_solid_lower = (*ls_solid_data)(si.toCell(0));
                        const double phi_solid_upper = (*ls_solid_data)(si.toCell(1));

                        // SETTING 3: Simple average of phi onto side centers and set rho_sc directly
                        double h_solid;
                        const double* const patch_dx = patch_geom->getDx();
                        const double alpha = d_num_solid_interface_cells * patch_dx[0];
                        const double phi_solid = 0.5 * (phi_solid_lower + phi_solid_upper);

                        if (phi_solid < -alpha)
                        {
                            h_solid = 0.0;
                        }
                        else if (std::abs(phi_solid) <= alpha)
                        {
                            h_solid =
                                0.5 + 0.5 * phi_solid / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_solid / alpha);
                        }
                        else
                        {
                            h_solid = 1.0;
                        }

                        (*rho_data)(si) = d_rho_fluid * h_solid + (1.0 - h_solid) * d_rho_solid;
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_solid_scratch_idx);
        }
        var_db->removePatchDataIndex(ls_solid_scratch_idx);
    }

    return;
} // setDensityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
