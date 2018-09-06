// Filename: SetFluidGasSolidViscosity.cpp
// Created on Dec 14, 2017 by Nishant Nangia

// APPLICATION INCLUDES
#include "SetFluidGasSolidViscosity.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidGasSolidViscosityCallbackFunction(int mu_idx,
                                              Pointer<Variable<NDIM> > mu_var,
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
                                                     Pointer<CellVariable<NDIM, double> > ls_solid_var,
                                                     Pointer<CellVariable<NDIM, double> > ls_gas_var,
                                                     const double mu_fluid,
                                                     const double mu_gas,
                                                     const double mu_solid,
                                                     const int ls_reinit_interval,
                                                     const double num_solid_interface_cells,
                                                     const double num_gas_interface_cells,
                                                     const bool set_mu_solid)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_solid_var(ls_solid_var),
      d_ls_gas_var(ls_gas_var),
      d_mu_fluid(mu_fluid),
      d_mu_gas(mu_gas),
      d_mu_solid(mu_solid),
      d_ls_reinit_interval(ls_reinit_interval),
      d_num_solid_interface_cells(num_solid_interface_cells),
      d_num_gas_interface_cells(num_gas_interface_cells),
      d_set_mu_solid(set_mu_solid)
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
    int ls_solid_idx = -1;
    int ls_gas_idx = -1;

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

    if (MathUtilities<double>::equalEps(time, current_time))
    {
        ls_gas_idx = var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
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
                const double alpha = d_num_solid_interface_cells * patch_dx[0];
                const double beta = d_num_gas_interface_cells * patch_dx[1];

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
                const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_idx);
                Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

                // Calderer et al, 2014
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi_s = (*ls_solid_data)(ci);
                    const double phi_g = (*ls_gas_data)(ci);
                    double Hphi_s, Hphi_g;
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

                    // First, compute the viscosity of the "flowing" phases
                    const double mu_flow = (d_mu_fluid - d_mu_gas) * Hphi_g + d_mu_gas;

                    // Next, set the viscosity of the solid phase in the usual way
                    if (d_set_mu_solid)
                        (*mu_data)(ci) = (mu_flow - d_mu_solid) * Hphi_s + d_mu_solid;
                    else
                        (*mu_data)(ci) = mu_flow;
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
