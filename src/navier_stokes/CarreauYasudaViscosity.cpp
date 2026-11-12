// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibamr/CarreauYasudaViscosity.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <HierarchySideDataOpsReal.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

CarreauYasudaViscosity::CarreauYasudaViscosity(std::string object_name,
                                               Pointer<Database> input_db,
                                               INSVCStaggeredHierarchyIntegrator& ins_integrator)
    : d_object_name(std::move(object_name)), d_ins_integrator(ins_integrator)
{
    if (d_ins_integrator.muIsConstant())
    {
        TBOX_ERROR(d_object_name << "::CarreauYasudaViscosity():\n"
                                 << "  unsupported mu_is_const = TRUE in " << d_ins_integrator.getName() << "\n"
                                 << "  set mu_is_const = FALSE or remove it\n");
    }
    if (d_ins_integrator.getViscosityVariable() || d_ins_integrator.getTransportedViscosityVariable())
    {
        TBOX_ERROR(d_object_name << "::CarreauYasudaViscosity():\n"
                                 << "  " << d_ins_integrator.getName() << " already has a viscosity variable\n");
    }

    const auto get_parameter = [&](const std::string& key, const bool positive)
    {
        const double value = input_db->getDouble(key);
        if (!(positive ? value > 0.0 : value >= 0.0))
        {
            TBOX_ERROR(d_object_name << "::CarreauYasudaViscosity():\n"
                                     << "  invalid " << key << " = " << value << "\n"
                                     << "  " << key << " must be " << (positive ? "positive" : "nonnegative") << "\n");
        }
        return value;
    };
    d_zero_shear_viscosity = get_parameter("zero_shear_viscosity", true);
    d_infinite_shear_viscosity = get_parameter("infinite_shear_viscosity", false);
    d_relaxation_time = get_parameter("relaxation_time", false);
    d_power_law_index = get_parameter("power_law_index", true);
    d_yasuda_exponent = get_parameter("yasuda_exponent", true);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> scratch_ctx = var_db->getContext(d_object_name + "::SCRATCH");
    d_U_scratch_var = new SideVariable<NDIM, double>(d_object_name + "::U");
    d_U_scratch_idx = var_db->registerVariableAndContext(d_U_scratch_var, scratch_ctx, IntVector<NDIM>(1));
    d_EE_var = new CellVariable<NDIM, double>(d_object_name + "::EE", NDIM * (NDIM + 1) / 2);
    d_EE_idx = var_db->registerVariableAndContext(d_EE_var, scratch_ctx, IntVector<NDIM>(0));

    d_mu_var = new CellVariable<NDIM, double>(d_object_name + "::mu");
    d_ins_integrator.registerViscosityVariable(d_mu_var);
    d_ins_integrator.registerResetFluidViscosityFcn(&reset_viscosity, this);
    return;
} // CarreauYasudaViscosity

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CarreauYasudaViscosity::computeViscosity(const int mu_idx, const int U_idx, const double time)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = d_ins_integrator.getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_U_scratch_idx, time);
        level->allocatePatchData(d_EE_idx, time);
    }

    // Fill the velocity ghost cells in the same way as the integrator. The operator is built on each call because
    // its schedules refer to the current patch levels, and a parent integrator may regrid the hierarchy without
    // notifying this object.
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = d_ins_integrator.getVelocityBoundaryConditions();
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent U_bc_component(d_U_scratch_idx,
                                                     /*refine_op_name*/ "NONE",
                                                     /*use_cf_bdry_interpolation*/ true,
                                                     /*coarsen_op_name*/ "CUBIC_COARSEN",
                                                     /*phys_bdry_extrap_type*/ "LINEAR",
                                                     /*consistent_type_2_bdry*/ false,
                                                     U_bc_coefs);
    HierarchyGhostCellInterpolation U_bdry_bc_fill_op;
    U_bdry_bc_fill_op.initializeOperatorState(U_bc_component, hierarchy);
    HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(hierarchy, coarsest_ln, finest_ln);
    hier_sc_data_ops.copyData(d_U_scratch_idx, U_idx);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(U_bc_coefs,
                                                              /*P_bc_coef*/ nullptr,
                                                              d_U_scratch_idx,
                                                              /*P_data_idx*/ -1,
                                                              /*homogeneous_bc*/ false);
    U_bdry_bc_fill_op.fillData(time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(U_bc_coefs, /*P_bc_coef*/ nullptr);
    d_ins_integrator.getHierarchyMathOps()->strain_rate(
        d_EE_idx, d_EE_var, d_U_scratch_idx, d_U_scratch_var, /*src_ghost_fill*/ nullptr, time);

    const double exponent = (d_power_law_index - 1.0) / d_yasuda_exponent;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);
            Pointer<CellData<NDIM, double>> EE_data = patch->getPatchData(d_EE_idx);
            for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                // E : E counts each off-diagonal component twice.
                double EE_norm_sq = 0.0;
                for (int d = 0; d < NDIM * (NDIM + 1) / 2; ++d)
                {
                    EE_norm_sq += (d < NDIM ? 1.0 : 2.0) * (*EE_data)(i, d) * (*EE_data)(i, d);
                }
                const double shear_rate = std::sqrt(2.0 * EE_norm_sq);
                (*mu_data)(i) =
                    d_infinite_shear_viscosity +
                    (d_zero_shear_viscosity - d_infinite_shear_viscosity) *
                        std::pow(1.0 + std::pow(d_relaxation_time * shear_rate, d_yasuda_exponent), exponent);
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_U_scratch_idx);
        level->deallocatePatchData(d_EE_idx);
    }
    return;
} // computeViscosity

void
CarreauYasudaViscosity::reset_viscosity(const int mu_idx,
                                        Pointer<Variable<NDIM>> /*mu_var*/,
                                        Pointer<HierarchyMathOps> /*hier_math_ops*/,
                                        const int /*cycle_num*/,
                                        const double time,
                                        const double /*current_time*/,
                                        const double /*new_time*/,
                                        void* ctx)
{
    auto* const viscosity = static_cast<CarreauYasudaViscosity*>(ctx);
    const INSVCStaggeredHierarchyIntegrator& integrator = viscosity->d_ins_integrator;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const bool is_current =
        mu_idx == var_db->mapVariableAndContextToIndex(viscosity->d_mu_var, integrator.getCurrentContext());
    const int U_idx = var_db->mapVariableAndContextToIndex(
        integrator.getVelocityVariable(), is_current ? integrator.getCurrentContext() : integrator.getNewContext());
    viscosity->computeViscosity(mu_idx, U_idx, time);
    return;
} // reset_viscosity

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
