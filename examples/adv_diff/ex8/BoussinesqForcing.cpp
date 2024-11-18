// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "BoussinesqForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

BoussinesqForcing::BoussinesqForcing(Pointer<SAMRAI::hier::Variable<NDIM> > T_var,
                                     Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator,
                                     Pointer<CellVariable<NDIM, double> > ls_inner_solid_var,
                                     Pointer<CellVariable<NDIM, double> > ls_outer_solid_var,
                                     Pointer<Database> input_db)
    : d_T_var(T_var), d_adv_diff_hier_integrator(adv_diff_hier_integrator)
{
    d_rayleigh_number = input_db->getDouble("rayleigh_number");
    d_prandtl_number = input_db->getDouble("prandtl_number");

    d_ls_solid_vars.push_back(ls_inner_solid_var);
    d_ls_solid_vars.push_back(ls_outer_solid_var);

    d_chi_var = new SideVariable<NDIM, double>("chi");

    return;
} // BoussinesqForcing

BoussinesqForcing::~BoussinesqForcing()
{
    // intentionally blank
    return;
} // ~BoussinesqForcing

bool
BoussinesqForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
BoussinesqForcing::setDataOnPatchHierarchy(const int data_idx,
                                           Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                           Pointer<PatchHierarchy<NDIM> > hierarchy,
                                           const double data_time,
                                           const bool initial_time,
                                           const int coarsest_ln_in,
                                           const int finest_ln_in)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_chi_idx = var_db->registerVariableAndContext(d_chi_var, var_db->getContext(d_object_name + "::chi"));

    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(d_chi_idx, data_time);
    }
    int T_scratch_idx = var_db->mapVariableAndContextToIndex(d_T_var, d_adv_diff_hier_integrator->getScratchContext());

    int T_new_idx = var_db->mapVariableAndContextToIndex(d_T_var, d_adv_diff_hier_integrator->getNewContext());
    HierarchyDataOpsManager<NDIM>* hier_data_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_cc_data_ops =
        hier_data_ops_manager->getOperationsDouble(d_T_var, hierarchy, /*get_unique*/ true);
    hier_cc_data_ops->copyData(T_scratch_idx, T_new_idx);

    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent ghost_fill_component(T_scratch_idx,
                                                           "CONSERVATIVE_LINEAR_REFINE",
                                                           false,
                                                           "CONSERVATIVE_COARSEN",
                                                           "LINEAR",
                                                           false,
                                                           d_adv_diff_hier_integrator->getPhysicalBcCoefs(d_T_var));
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_fill_component, hierarchy);

    // Filling the ghost cells for all the level set variables
    std::vector<InterpolationTransactionComponent> ls_transaction_comps(d_ls_solid_vars.size());
    int n = 0;
    for (auto ls_solid_var : d_ls_solid_vars)
    {
        int ls_solid_var_new_idx =
            var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_hier_integrator->getNewContext());
        int ls_solid_var_scratch_idx =
            var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_hier_integrator->getScratchContext());

        hier_cc_data_ops->copyData(ls_solid_var_scratch_idx, ls_solid_var_new_idx);

        ls_transaction_comps[n] =
            InterpolationTransactionComponent(ls_solid_var_scratch_idx,
                                              "CONSERVATIVE_LINEAR_REFINE",
                                              false,
                                              "CONSERVATIVE_COARSEN",
                                              "LINEAR",
                                              false,
                                              d_adv_diff_hier_integrator->getPhysicalBcCoefs(ls_solid_var));

        n = n + 1;
    }
    ghost_fill_op.initializeOperatorState(ls_transaction_comps, hierarchy);
    ghost_fill_op.fillData(data_time);

    for (int level_num = coarsest_ln; level_num <= finest_ln; ++level_num)
    {
        setDataOnPatchLevel(data_idx, var, hierarchy->getPatchLevel(level_num), data_time, initial_time);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(d_chi_idx);
    }
    return;
} // setDataOnPatchHierarchy

void
BoussinesqForcing::setDataOnPatch(const int data_idx,
                                  Pointer<SAMRAI::hier::Variable<NDIM> > /*var*/,
                                  Pointer<Patch<NDIM> > patch,
                                  const double /*data_time*/,
                                  const bool initial_time,
                                  Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    if (initial_time) return;
    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);

    Pointer<CellData<NDIM, double> > T_data =
        patch->getPatchData(d_T_var, d_adv_diff_hier_integrator->getScratchContext());
    Pointer<SideData<NDIM, double> > chi_data = patch->getPatchData(d_chi_idx);
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    const double* patch_dx = patch_geom->getDx();
    const double alpha = 2.0 * patch_dx[0];
    const Box<NDIM>& patch_box = patch->getBox();

    for (unsigned int d = 0; d < NDIM; d++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, d)); it; it++)
        {
            SideIndex<NDIM> si(it(), d, SideIndex<NDIM>::Lower);
            (*chi_data)(si) = 0.0;
            for (auto ls_solid_var : d_ls_solid_vars)
            {
                VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
                int ls_solid_scratch_idx =
                    var_db->mapVariableAndContextToIndex(ls_solid_var, d_adv_diff_hier_integrator->getScratchContext());
                Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_scratch_idx);
                const double phi_lower = (*ls_solid_data)(si.toCell(0));
                const double phi_upper = (*ls_solid_data)(si.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                double Hphi;
                if (phi < -alpha)
                {
                    Hphi = 0.0;
                }
                else if (std::abs(phi) <= alpha)
                {
                    Hphi = 0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha);
                }
                else
                {
                    Hphi = 1.0;
                }
                const double chi = 1.0 - Hphi;

                (*chi_data)(si) = (*chi_data)(si) + chi;
            }
        }
    }

    for (unsigned int d = 0; d < NDIM; d++)
    {
        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, d)); it; it++)
        {
            SideIndex<NDIM> si(it(), d, SideIndex<NDIM>::Lower);
            (*F_data)(si) = 0.0;
            if (d == 1)
            {
                (*F_data)(si) = (1.0 - (*chi_data)(si)) * d_rayleigh_number * d_prandtl_number * 0.5 *
                                ((*T_data)(si.toCell(0)) + (*T_data)(si.toCell(1)));
            }
        }
    }
    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
