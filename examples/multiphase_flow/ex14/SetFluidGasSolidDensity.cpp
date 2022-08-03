// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// APPLICATION INCLUDES
#include <ibtk/FEDataManager.h>
#include <ibtk/HierarchyMathOps.h>

#include "SetFluidGasSolidDensity.h"

#include <libmesh/equation_systems.h>
#include <libmesh/petsc_vector.h>

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidGasSolidDensityCallbackFunction(int rho_idx,
                                            Pointer<hier::Variable<NDIM> > rho_var,
                                            Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                            const int cycle_num,
                                            const double time,
                                            const double current_time,
                                            const double new_time,
                                            void* ctx)
{
    // Set the density from the level set information
    static SetFluidGasSolidDensity* ptr_SetFluidGasSolidDensity = static_cast<SetFluidGasSolidDensity*>(ctx);
    ptr_SetFluidGasSolidDensity->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidGasSolidDensityCallBackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidGasSolidDensity::SetFluidGasSolidDensity(const std::string& object_name,
                                                 Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                 EulerianFEStructure* efes,
                                                 Pointer<CellVariable<NDIM, double> > ls_gas_var,
                                                 const double rho_fluid,
                                                 const double rho_gas,
                                                 const double rho_solid,
                                                 const double num_gas_interface_cells)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_efes(efes),
      d_ls_gas_var(ls_gas_var),
      d_rho_fluid(rho_fluid),
      d_rho_gas(rho_gas),
      d_rho_solid(rho_solid),
      d_num_gas_interface_cells(num_gas_interface_cells)
{
    // intentionally left blank
    return;
} // SetFluidGasSolidDensity

SetFluidGasSolidDensity::~SetFluidGasSolidDensity()
{
    // intentionally left blank
    return;

} //~SetFluidGasSolidDensity

void
SetFluidGasSolidDensity::setDensityPatchData(int rho_idx,
                                             Pointer<hier::Variable<NDIM> > rho_var,
                                             Pointer<HierarchyMathOps> hier_math_ops,
                                             const int /*cycle_num*/,
                                             const double time,
                                             const double current_time,
                                             const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int ls_gas_idx = -1;
    if (IBTK::rel_equal_eps(time, current_time))
    {
        ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(time, new_time))
    {
        ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Prolong the FE structure onto the Eulerian grid
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    const std::string& X_system_name = d_efes->d_fe_data_manager->COORDINATES_SYSTEM_NAME;
    EquationSystems* equation_systems = d_efes->d_fe_data_manager->getEquationSystems();
    auto& X_system = equation_systems->get_system(X_system_name);
    NumericVector<double>& X_current = *X_system.current_local_solution;
    std::unique_ptr<PetscVector<double> > X_IB = d_efes->d_fe_data_manager->buildIBGhostedVector(X_system_name);
    copy_and_synch(X_current, *X_IB);
    std::unique_ptr<NumericVector<double> > chi_num_vec = X_IB->zero_clone();
    PetscVector<double>& chi_petsc_vec = dynamic_cast<PetscVector<double>&>(*chi_num_vec);
    chi_petsc_vec = 1.0;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_efes->d_chi_idx)) level->allocatePatchData(d_efes->d_chi_idx);
    }
    d_efes->d_fe_data_manager->prolongData(d_efes->d_chi_idx,
                                           chi_petsc_vec,
                                           *X_IB,
                                           X_system_name,
                                           /*is_density*/ false,
                                           /*accumulate_on_grid*/ false,
                                           /*close_F*/ true,
                                           /*close_X*/ false);

    // Setting side centered density directly
    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    TBOX_ASSERT(!rho_sc_var.isNull());
    if (rho_sc_var)
    {
        // Note, this method requires ghost cells to be filled for the level set variable
        std::vector<RobinBcCoefStrategy<NDIM>*> ls_solid_bc_coef(NDIM, nullptr);
        RobinBcCoefStrategy<NDIM>* ls_gas_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_ls_gas_var).front();
        IntVector<NDIM> cell_ghosts = 1;
        const int ls_solid_scratch_idx = var_db->registerVariableAndContext(
            d_efes->d_chi_var, var_db->getContext(d_object_name + "::SOLID::SCRATCH"), cell_ghosts);
        const int ls_gas_scratch_idx = var_db->registerVariableAndContext(
            d_ls_gas_var, var_db->getContext(d_object_name + "::GAS::SCRATCH"), cell_ghosts);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_solid_scratch_idx, time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_gas_scratch_idx, time);
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ls_transaction_comps(2);
        ls_transaction_comps[0] = InterpolationTransactionComponent(ls_solid_scratch_idx,
                                                                    d_efes->d_chi_idx,
                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                    false,
                                                                    "CONSERVATIVE_COARSEN",
                                                                    "LINEAR",
                                                                    false,
                                                                    ls_solid_bc_coef);
        ls_transaction_comps[1] = InterpolationTransactionComponent(ls_gas_scratch_idx,
                                                                    ls_gas_idx,
                                                                    "CONSERVATIVE_LINEAR_REFINE",
                                                                    false,
                                                                    "CONSERVATIVE_COARSEN",
                                                                    "LINEAR",
                                                                    false,
                                                                    ls_gas_bc_coef);
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
                const Pointer<SideData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_scratch_idx);
                const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                // Compute the indicators for both level sets
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double phi_solid = (*ls_solid_data)(si);
                        const double phi_gas_lower = (*ls_gas_data)(si.toCell(0));
                        const double phi_gas_upper = (*ls_gas_data)(si.toCell(1));

                        // SETTING 3: Simple average of phi onto side centers and set rho_sc directly
                        double h_solid, h_gas;
                        const double* const patch_dx = patch_geom->getDx();
                        const double beta = d_num_gas_interface_cells * patch_dx[1];
                        const double phi_gas = 0.5 * (phi_gas_lower + phi_gas_upper);

                        if (phi_solid > 0.0)
                        {
                            h_solid = 0.0;
                        }
                        else
                        {
                            h_solid = 1.0;
                        }

                        if (phi_gas < -beta)
                        {
                            h_gas = 0.0;
                        }
                        else if (std::abs(phi_gas) <= beta)
                        {
                            h_gas = 0.5 + 0.5 * phi_gas / beta + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_gas / beta);
                        }
                        else
                        {
                            h_gas = 1.0;
                        }
                        const double rho_flow = (d_rho_fluid - d_rho_gas) * h_gas + d_rho_gas;
                        const double rho_full = (rho_flow - d_rho_solid) * h_solid + d_rho_solid;

                        (*rho_data)(si) = rho_full;
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_solid_scratch_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_gas_scratch_idx);
            // patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_efes->d_chi_idx);
        }
        var_db->removePatchDataIndex(ls_solid_scratch_idx);
        var_db->removePatchDataIndex(ls_gas_scratch_idx);
    }

    return;
} // setDensityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
