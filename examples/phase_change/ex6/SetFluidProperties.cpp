// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

///////////////////////////// INCLUDES ///////////////////////////////////
#include "SetFluidProperties.h"

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////
void
callSetLiquidSolidGasDensityCallbackFunction(int rho_idx,
                                             Pointer<Variable<NDIM> > rho_var,
                                             Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                             const int cycle_num,
                                             const double time,
                                             const double current_time,
                                             const double new_time,
                                             void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidGasSolidDensityCallBackFunction

void
callSetLiquidSolidGasSpecificHeatCallbackFunction(int Cp_idx,
                                                  Pointer<Variable<NDIM> > Cp_var,
                                                  Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                  const int cycle_num,
                                                  const double time,
                                                  const double current_time,
                                                  const double new_time,
                                                  void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setSpecificHeatPatchData(
        Cp_idx, Cp_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetSpecificHeatCallbackFunction

void
callSetLiquidSolidGasConductivityCallbackFunction(int D_idx,
                                                  Pointer<Variable<NDIM> > D_var,
                                                  Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                  const int cycle_num,
                                                  const double time,
                                                  const double current_time,
                                                  const double new_time,
                                                  void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setDiffusionCoefficientPatchData(
        D_idx, D_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetSpecificHeatCallbackFunction

void
callSetLiquidGasSolidViscosityCallbackFunction(int mu_idx,
                                               Pointer<Variable<NDIM> > mu_var,
                                               Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                               const int cycle_num,
                                               const double time,
                                               const double current_time,
                                               const double new_time,
                                               void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidGasSolidViscosityCallBackFunction
/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       const Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       const Pointer<CellVariable<NDIM, double> > lf_var,
                                       RobinBcCoefStrategy<NDIM>* lf_bc_coef,
                                       const Pointer<CellVariable<NDIM, double> > H_var,
                                       RobinBcCoefStrategy<NDIM>* H_bc_coef,
                                       const double rho_liquid,
                                       const double rho_solid,
                                       const double rho_gas,
                                       const double kappa_liquid,
                                       const double kappa_solid,
                                       const double kappa_gas,
                                       const double Cp_liquid,
                                       const double Cp_solid,
                                       const double Cp_gas,
                                       const double mu_liquid,
                                       const double mu_solid,
                                       const double mu_gas)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_lf_var(lf_var),
      d_lf_bc_coef(lf_bc_coef),
      d_H_var(H_var),
      d_H_bc_coef(H_bc_coef),
      d_rho_liquid(rho_liquid),
      d_rho_solid(rho_solid),
      d_rho_gas(rho_gas),
      d_kappa_liquid(kappa_liquid),
      d_kappa_solid(kappa_solid),
      d_kappa_gas(kappa_gas),
      d_Cp_liquid(Cp_liquid),
      d_Cp_solid(Cp_solid),
      d_Cp_gas(Cp_gas),
      d_mu_liquid(mu_liquid),
      d_mu_solid(mu_solid),
      d_mu_gas(mu_gas)
{
    // intentionally left blank
    return;
} // SetFluidProperties

SetFluidProperties::~SetFluidProperties()
{
    // intentionally left blank
    return;

} //~SetFluidProperties

void
SetFluidProperties::setDensityPatchData(int rho_idx,
                                        Pointer<Variable<NDIM> > rho_var,
                                        Pointer<HierarchyMathOps> hier_math_ops,
                                        const int /*cycle_num*/,
                                        const double time,
                                        const double current_time,
                                        const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    int lf_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int H_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
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

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
                const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
                const Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
                // Li et al, 2015
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double liquid_fraction = (*lf_data)(ci);
                    const double heaviside = (*H_data)(ci);

                    (*rho_data)(ci) = d_rho_gas + (d_rho_solid - d_rho_gas) * heaviside +
                                      (d_rho_liquid - d_rho_solid) * liquid_fraction * heaviside;
                }
            }
        }
    }

    // Setting side centered density directly
    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        RobinBcCoefStrategy<NDIM>* lf_bc_coef = d_lf_bc_coef;
        RobinBcCoefStrategy<NDIM>* H_bc_coef = d_H_bc_coef;
        IntVector<NDIM> cell_ghosts = 1;
        const int lf_scratch_idx =
            var_db->registerVariableAndContext(d_lf_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);
        const int H_scratch_idx =
            var_db->registerVariableAndContext(d_H_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(lf_scratch_idx, time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(H_scratch_idx, time);
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> interpolation_transaction(2);
        interpolation_transaction[0] = InterpolationTransactionComponent(lf_scratch_idx,
                                                                         lf_idx,
                                                                         "CONSERVATIVE_LINEAR_REFINE",
                                                                         false,
                                                                         "CONSERVATIVE_COARSEN",
                                                                         "LINEAR",
                                                                         false,
                                                                         lf_bc_coef);

        interpolation_transaction[1] = InterpolationTransactionComponent(H_scratch_idx,
                                                                         H_idx,
                                                                         "CONSERVATIVE_LINEAR_REFINE",
                                                                         false,
                                                                         "CONSERVATIVE_COARSEN",
                                                                         "LINEAR",
                                                                         false,
                                                                         H_bc_coef);

        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(interpolation_transaction, patch_hierarchy);
        hier_bdry_fill->fillData(time);

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_scratch_idx);
                const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double liquid_fraction_lower = (*lf_data)(si.toCell(0));
                        const double liquid_fraction_upper = (*lf_data)(si.toCell(1));
                        const double liquid_fraction = 0.5 * (liquid_fraction_lower + liquid_fraction_upper);

                        const double heaviside_lower = (*H_data)(si.toCell(0));
                        const double heaviside_upper = (*H_data)(si.toCell(1));
                        const double heaviside = 0.5 * (heaviside_lower + heaviside_upper);

                        (*rho_data)(si) = d_rho_gas + (d_rho_solid - d_rho_gas) * heaviside +
                                          (d_rho_liquid - d_rho_solid) * liquid_fraction * heaviside;
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(lf_scratch_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(H_scratch_idx);
        }
        var_db->removePatchDataIndex(lf_scratch_idx);
        var_db->removePatchDataIndex(H_scratch_idx);
    }

    return;
} // setDensityPatchData

void
SetFluidProperties::setDiffusionCoefficientPatchData(int D_idx,
                                                     Pointer<Variable<NDIM> > /*D_var*/,
                                                     SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                                     const int /*cycle_num*/,
                                                     const double time,
                                                     const double current_time,
                                                     const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int lf_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int H_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            const Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double liquid_fraction = (*lf_data)(ci);
                const double heaviside = (*H_data)(ci);
                (*D_data)(ci) = d_kappa_gas + (d_kappa_solid - d_kappa_gas) * heaviside +
                                (d_kappa_liquid - d_kappa_solid) * liquid_fraction * heaviside;
            }
        }
    }

    return;
} // setDiffusionCoefficientPatchData

void
SetFluidProperties::setSpecificHeatPatchData(int Cp_idx,
                                             Pointer<Variable<NDIM> > /*Cp_var*/,
                                             SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                             const int /*cycle_num*/,
                                             const double time,
                                             const double current_time,
                                             const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int lf_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int H_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > Cp_data = patch->getPatchData(Cp_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double liquid_fraction = (*lf_data)(ci);
                const double heaviside = (*H_data)(ci);
                (*Cp_data)(ci) = d_Cp_gas + (d_Cp_solid - d_Cp_gas) * heaviside +
                                 (d_Cp_liquid - d_Cp_solid) * liquid_fraction * heaviside;
            }
        }
    }

    return;
} // setSpecificHeatPatchData

void
SetFluidProperties::setViscosityPatchData(int mu_idx,
                                          Pointer<Variable<NDIM> > /*mu_var*/,
                                          SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                          const int /*cycle_num*/,
                                          const double time,
                                          const double current_time,
                                          const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int lf_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        lf_idx = var_db->mapVariableAndContextToIndex(d_lf_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    int H_idx = -1;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        H_idx = var_db->mapVariableAndContextToIndex(d_H_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
            Pointer<CellData<NDIM, double> > lf_data = patch->getPatchData(lf_idx);
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double liquid_fraction = (*lf_data)(ci);
                const double heaviside = (*H_data)(ci);
                (*mu_data)(ci) = d_mu_gas + (d_mu_solid - d_mu_gas) * heaviside +
                                 (d_mu_liquid - d_mu_solid) * liquid_fraction * heaviside;
            }
        }
    }

    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////
