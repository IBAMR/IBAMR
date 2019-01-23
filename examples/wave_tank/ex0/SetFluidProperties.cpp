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

} // callSetFluidDensityCallbackFunction

// Various options to setting side-centered densities
#define SMOOTH_SC_RHO 1
#define DESJARDINS_SC_RHO 0

void
callSetFluidViscosityCallbackFunction(int mu_idx,
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

} // callSetFluidViscosityCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                       Pointer<CellVariable<NDIM, double> > ls_var,
                                       const double rho_outside,
                                       const double rho_inside,
                                       const double mu_outside,
                                       const double mu_inside,
                                       const int ls_reinit_interval,
                                       const double num_interface_cells)
    : d_object_name(object_name),
      d_adv_diff_solver(adv_diff_solver),
      d_ls_var(ls_var),
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
                                        Pointer<Variable<NDIM> > rho_var,
                                        SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                        const int /*cycle_num*/,
                                        const double time,
                                        const double current_time,
                                        const double new_time)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int ls_idx;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        ls_idx = var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        ls_idx = var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getNewContext());
    }
    else
    {
        TBOX_ERROR("This statement should not be reached");
    }

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    // Normal way to set cell centered density
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
                double vol_cell = 1.0;
                for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
                double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / (double)NDIM);

                const Box<NDIM>& patch_box = patch->getBox();
                const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
                Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi = (*ls_data)(ci);

                    // Calderer et al, 2014
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
    }

    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        // Note, this method requires ghost cells to be filled for the level set variable
        RobinBcCoefStrategy<NDIM>* ls_bc_coef = d_adv_diff_solver->getPhysicalBcCoefs(d_ls_var).front();
        IntVector<NDIM> cell_ghosts = 1;
        const int ls_scratch_idx =
            var_db->registerVariableAndContext(d_ls_var, var_db->getContext(d_object_name + "::SCRATCH"), cell_ghosts);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(ls_scratch_idx, time);
        }
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent ls_transaction(ls_scratch_idx,
                                                         ls_idx,
                                                         "CONSERVATIVE_LINEAR_REFINE",
                                                         false,
                                                         "CONSERVATIVE_COARSEN",
                                                         "LINEAR",
                                                         false,
                                                         ls_bc_coef);
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(ls_transaction, patch_hierarchy);
        hier_bdry_fill->fillData(time);

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
                const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_scratch_idx);
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                    {
                        SideIndex<NDIM> si(it(), axis, 0);
                        const double phi_lower = (*ls_data)(si.toCell(0));
                        const double phi_upper = (*ls_data)(si.toCell(1));
                        double h;
#if (DESJARDINS_SC_RHO)
                        // Desjardins way to set side-centered density
                        if (phi_lower >= 0.0 && phi_upper >= 0.0)
                        {
                            h = 1.0;
                        }
                        else if (phi_lower < 0.0 && phi_upper < 0.0)
                        {
                            h = 0.0;
                        }
                        else
                        {
                            h = (std::max(phi_lower, 0.0) + std::max(phi_upper, 0.0)) /
                                (std::abs(phi_lower) + std::abs(phi_upper));
                        }
                        (*rho_data)(si) = d_rho_inside + (d_rho_outside - d_rho_inside) * h;
#endif
#if (SMOOTH_SC_RHO)
                        // Simple average of phi onto side centers and set rho_sc directly
                        const double phi = 0.5 * (phi_lower + phi_upper);

                        if (phi < -alpha)
                        {
                            h = 0.0;
                        }
                        else if (std::abs(phi) <= alpha)
                        {
                            h = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                        }
                        else
                        {
                            h = 1.0;
                        }

                        (*rho_data)(si) = (d_rho_outside - d_rho_inside) * h + d_rho_inside;
#endif
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_scratch_idx);
        }
        var_db->removePatchDataIndex(ls_scratch_idx);
    }

    return;
} // setDensityPatchData

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
    int ls_idx;
    if (MathUtilities<double>::equalEps(time, current_time))
    {
        ls_idx = var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
    }
    else if (MathUtilities<double>::equalEps(time, new_time))
    {
        ls_idx = var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getNewContext());
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
            const double* const patch_dx = patch_geom->getDx();
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / (double)NDIM);

            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);
            if (!ls_data) TBOX_ERROR("This statement should not be reached");
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double phi = (*ls_data)(ci);

                // Calderer et al, 2014
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
