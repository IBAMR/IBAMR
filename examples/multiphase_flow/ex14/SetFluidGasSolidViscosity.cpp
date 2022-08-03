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
#include <ibtk/HierarchyMathOps.h>

#include "SetFluidGasSolidViscosity.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidGasSolidViscosityCallbackFunction(int mu_idx,
                                              Pointer<hier::Variable<NDIM> > mu_var,
                                              Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                              const int cycle_num,
                                              const double time,
                                              const double current_time,
                                              const double new_time,
                                              void* ctx)
{
    // Set the density from the level set information
    static SetFluidGasSolidViscosity* ptr_SetFluidGasSolidViscosity = static_cast<SetFluidGasSolidViscosity*>(ctx);
    ptr_SetFluidGasSolidViscosity->setViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidGasSolidViscosityCallBackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidGasSolidViscosity::SetFluidGasSolidViscosity(const std::string& object_name,
                                                     Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                     Pointer<CellVariable<NDIM, double> > ls_gas_var,
                                                     const double mu_fluid,
                                                     const double mu_gas,
                                                     const double num_gas_interface_cells)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_gas_var(ls_gas_var),
      d_mu_fluid(mu_fluid),
      d_mu_gas(mu_gas),
      d_num_gas_interface_cells(num_gas_interface_cells)
{
    // intentionally left blank
    return;
} // SetFluidGasSolidViscosity

SetFluidGasSolidViscosity::~SetFluidGasSolidViscosity()
{
    // intentionally left blank
    return;

} //~SetFluidGasSolidViscosity

void
SetFluidGasSolidViscosity::setViscosityPatchData(int mu_idx,
                                                 Pointer<Variable<NDIM> > mu_var,
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

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    Pointer<CellVariable<NDIM, double> > mu_cc_var = mu_var;
    if (mu_cc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                const double beta = d_num_gas_interface_cells * patch_dx[1];

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_idx);
                Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

                // Calderer et al, 2014
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi_g = (*ls_gas_data)(ci);
                    double Hphi_g;

                    if (phi_g < -beta)
                    {
                        Hphi_g = 0.0;
                    }
                    else if (std::abs(phi_g) <= beta)
                    {
                        Hphi_g = 0.5 + 0.5 * phi_g / beta + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_g / beta);
                    }
                    else
                    {
                        Hphi_g = 1.0;
                    }

                    // Compute the viscosity of the "flowing" phases
                    (*mu_data)(ci) = (d_mu_fluid - d_mu_gas) * Hphi_g + d_mu_gas;
                }
            }
        }
    }
    else
    {
        // Erroring out if any other centered is used for mu
        TBOX_ERROR("This statement should not have been reached");
    }

    return;
} // setViscosityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
