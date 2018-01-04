// Filename: SetFluidProperties.cpp
// Created on Dec 17, 2017 by Nishant Nangia

// APPLICATION INCLUDES
#include "SetFluidProperties.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidDensityCallbackFunction(int rho_idx,
                                    Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int integrator_step,
                                    const double current_time,
                                    const bool initial_time,
                                    const bool regrid_time,
                                    void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setDensityPatchData(rho_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;

} // callSetFluidDensityCallbackFunction

void
callSetFluidViscosityCallbackFunction(int mu_idx,
                                      Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                      const int integrator_step,
                                      const double current_time,
                                      const bool initial_time,
                                      const bool regrid_time,
                                      void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setViscosityPatchData(mu_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;

} // callSetFluidViscosityCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       const std::string& ls_name,
                                       const double rho_outside,
                                       const double rho_inside,
                                       const double mu_outside,
                                       const double mu_inside,
                                       const int ls_reinit_interval,
                                       const double num_interface_cells)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_name(ls_name),
      d_rho_outside(rho_outside),
      d_rho_inside(rho_inside),
      d_mu_outside(mu_outside),
      d_mu_inside(mu_inside),
      d_ls_reinit_interval(ls_reinit_interval),
      d_num_interface_cells(num_interface_cells)
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
                                        SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                        const int integrator_step,
                                        const double /*current_time*/,
                                        const bool initial_time,
                                        const bool /*regrid_time*/)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double> > ls_var = d_adv_diff_solver->getLevelSetVariable(d_ls_name);
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(ls_var, d_adv_diff_solver->getCurrentContext());

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
            const double* const patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / (double)NDIM);

            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

            // Calderer et al, 2014
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double phi = (*ls_data)(ci);
                double h_phi;
                if (phi < -alpha)
                {
                    h_phi = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    h_phi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                }
                else
                {
                    h_phi = 1.0;
                }

                (*rho_data)(ci) = d_rho_inside + (d_rho_outside - d_rho_inside) * h_phi;
            }
        }
    }

    return;
} // setDensityPatchData

void
SetFluidProperties::setViscosityPatchData(int mu_idx,
                                          SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                          const int integrator_step,
                                          const double /*current_time*/,
                                          const bool initial_time,
                                          const bool /*regrid_time*/)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double> > ls_var = d_adv_diff_solver->getLevelSetVariable(d_ls_name);
    const int ls_current_idx = var_db->mapVariableAndContextToIndex(ls_var, d_adv_diff_solver->getCurrentContext());

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
            const double* const patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / (double)NDIM);

            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

            // Calderer et al, 2014
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double phi = (*ls_data)(ci);
                double h_phi;
                if (phi < -alpha)
                {
                    h_phi = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    h_phi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                }
                else
                {
                    h_phi = 1.0;
                }

                (*mu_data)(ci) = d_mu_inside + (d_mu_outside - d_mu_inside) * h_phi;
            }
        }
    }

    return;
} // setViscosityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
