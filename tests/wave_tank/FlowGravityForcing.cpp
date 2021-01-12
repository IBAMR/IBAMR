// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "FlowGravityForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FlowGravityForcing::FlowGravityForcing(const std::string& object_name,
                                       Pointer<Database> input_db,
                                       Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator,
                                       Pointer<CellVariable<NDIM, double> > ls_gas_var,
                                       std::vector<double> grav_const)
    : d_object_name(object_name),
      d_adv_diff_hierarchy_integrator(adv_diff_hierarchy_integrator),
      d_ls_gas_var(ls_gas_var),
      d_grav_const(grav_const)

{
    d_rho_neg = input_db->getDouble("rho_neg");
    d_rho_pos = input_db->getDouble("rho_pos");
    d_num_gas_interface_cells = input_db->getDouble("num_interface_cells");
    return;
} // FlowGravityForcing

FlowGravityForcing::~FlowGravityForcing()
{
    // intentionally blank
    return;
} // ~FlowGravityForcing

bool
FlowGravityForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
FlowGravityForcing::setDataOnPatchHierarchy(const int data_idx,
                                            Pointer<hier::Variable<NDIM> > /*var*/,
                                            Pointer<PatchHierarchy<NDIM> > hierarchy,
                                            const double data_time,
                                            const bool /*initial_time*/,
                                            const int coarsest_ln_in,
                                            const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);

    // Get level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int ls_gas_current_idx =
        var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_hierarchy_integrator->getCurrentContext());
    int ls_gas_new_idx =
        var_db->mapVariableAndContextToIndex(d_ls_gas_var, d_adv_diff_hierarchy_integrator->getNewContext());
    const bool ls_gas_new_is_allocated = d_adv_diff_hierarchy_integrator->isAllocatedPatchData(ls_gas_new_idx);
    int ls_gas_idx = ls_gas_new_is_allocated ? ls_gas_new_idx : ls_gas_current_idx;

    IntVector<NDIM> cell_ghosts = 1;
    int ls_gas_scratch_idx = var_db->registerVariableAndContext(
        d_ls_gas_var, var_db->getContext(d_object_name + "::LS_GAS_SCRATCH"), cell_ghosts);

#if !defined(NDEBUG)
    TBOX_ASSERT(ls_gas_idx >= 0);
    TBOX_ASSERT(ls_gas_scratch_idx >= 0);
#endif

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(ls_gas_scratch_idx, data_time);
    }
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ls_transaction_comp(1);
    ls_transaction_comp[0] =
        InterpolationTransactionComponent(ls_gas_scratch_idx,
                                          ls_gas_idx,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          false,
                                          "CONSERVATIVE_COARSEN",
                                          "LINEAR",
                                          false,
                                          d_adv_diff_hierarchy_integrator->getPhysicalBcCoefs(d_ls_gas_var));
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(ls_transaction_comp, hierarchy);
    hier_bdry_fill->fillData(data_time);

    // Set the gravity force. In this version, the gravity force is reconstructed from the flow density field.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
            const Pointer<CellData<NDIM, double> > ls_gas_data = patch->getPatchData(ls_gas_scratch_idx);

            double beta = 1.0;
            for (int d = 0; d < NDIM; ++d) beta *= patch_dx[d];
            beta = std::pow(beta, 1.0 / NDIM) * d_num_gas_interface_cells;

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                    // Reconstruct density
                    double phi_gas_lower = (*ls_gas_data)(s_i.toCell(0));
                    double phi_gas_upper = (*ls_gas_data)(s_i.toCell(1));
                    double h_gas_lower, h_gas_upper;
                    if (phi_gas_lower < -beta)
                    {
                        h_gas_lower = 0.0;
                    }
                    else if (std::abs(phi_gas_lower) <= beta)
                    {
                        h_gas_lower = 0.5 + 0.5 * phi_gas_lower / beta +
                                      1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_gas_lower / beta);
                    }
                    else
                    {
                        h_gas_lower = 1.0;
                    }
                    if (phi_gas_upper < -beta)
                    {
                        h_gas_upper = 0.0;
                    }
                    else if (std::abs(phi_gas_upper) <= beta)
                    {
                        h_gas_upper = 0.5 + 0.5 * phi_gas_upper / beta +
                                      1.0 / (2.0 * M_PI) * std::sin(M_PI * phi_gas_upper / beta);
                    }
                    else
                    {
                        h_gas_upper = 1.0;
                    }
                    const double rho_flow_lower = (d_rho_pos - d_rho_neg) * h_gas_lower + d_rho_neg;
                    const double rho_flow_upper = (d_rho_pos - d_rho_neg) * h_gas_upper + d_rho_neg;
                    (*f_data)(s_i) = d_grav_const[axis] * 2.0 * (rho_flow_lower * rho_flow_upper) /
                                     (rho_flow_lower + rho_flow_upper);
                }
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(ls_gas_scratch_idx);
    }
    var_db->removePatchDataIndex(ls_gas_scratch_idx);
    return;
} // setDataOnPatchHierarchy

void
FlowGravityForcing::setDataOnPatch(const int data_idx,
                                   Pointer<hier::Variable<NDIM> > /*var*/,
                                   Pointer<Patch<NDIM> > patch,
                                   const double /*data_time*/,
                                   const bool initial_time,
                                   Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    if (initial_time)
    {
        Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(data_idx);
        f_data->fillAll(0.0);
    }
    // Intentionally left blank

} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
