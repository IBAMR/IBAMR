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

#include <ibamr/INSVCStaggeredCarreauYasudaHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <HierarchyDataOpsReal.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <VariableDatabase.h>

#include <cmath>
#include <string>
#include <utility>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredCarreauYasudaHierarchyIntegrator::INSVCStaggeredCarreauYasudaHierarchyIntegrator(
    std::string object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : INSVCStaggeredNonConservativeHierarchyIntegrator(std::move(object_name), input_db, register_for_restart)
{
    if (d_mu_is_const)
    {
        TBOX_ERROR(d_object_name << "::INSVCStaggeredCarreauYasudaHierarchyIntegrator():\n"
                                 << "  unsupported mu_is_const = TRUE\n"
                                 << "  set mu_is_const = FALSE or remove it\n");
    }

    const auto get_parameter = [&](const std::string& key, const bool positive)
    {
        const double value = input_db->getDouble(key);
        if (!(positive ? value > 0.0 : value >= 0.0))
        {
            TBOX_ERROR(d_object_name << "::INSVCStaggeredCarreauYasudaHierarchyIntegrator():\n"
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

    registerViscosityVariable(new CellVariable<NDIM, double>(d_object_name + "::mu"));
    registerResetFluidViscosityFcn(&reset_viscosity, this);
    return;
} // INSVCStaggeredCarreauYasudaHierarchyIntegrator

void
INSVCStaggeredCarreauYasudaHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM>> hierarchy,
    Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    if (d_integrator_is_initialized) return;

    INSVCStaggeredNonConservativeHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    d_strain_rate_var = new CellVariable<NDIM, double>(d_object_name + "::strain_rate", NDIM * (NDIM + 1) / 2);
    d_strain_rate_idx = VariableDatabase<NDIM>::getDatabase()->registerVariableAndContext(
        d_strain_rate_var, getScratchContext(), IntVector<NDIM>(0));
    return;
} // initializeHierarchyIntegrator

void
INSVCStaggeredCarreauYasudaHierarchyIntegrator::initializePatchHierarchy(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                                         Pointer<GriddingAlgorithm<NDIM>> gridding_alg)
{
    INSVCStaggeredNonConservativeHierarchyIntegrator::initializePatchHierarchy(hierarchy, gridding_alg);
    computeViscosity(d_mu_current_idx, d_U_current_idx, d_integrator_time);
    return;
} // initializePatchHierarchy

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSVCStaggeredCarreauYasudaHierarchyIntegrator::computeViscosity(const int mu_idx, const int U_idx, const double time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // The velocity scratch data are allocated only during a time step.
    const Pointer<PatchLevel<NDIM>> coarsest_level = d_hierarchy->getPatchLevel(coarsest_ln);
    const bool allocate_U_scratch = !coarsest_level->checkAllocated(d_U_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (allocate_U_scratch) level->allocatePatchData(d_U_scratch_idx, time);
        level->allocatePatchData(d_strain_rate_idx, time);
    }

    d_hier_sc_data_ops->copyData(d_U_scratch_idx, U_idx);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ nullptr,
                                                              d_U_scratch_idx,
                                                              /*P_data_idx*/ -1,
                                                              /*homogeneous_bc*/ false);
    d_U_bdry_bc_fill_op->fillData(time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs,
                                                              /*P_bc_coef*/ nullptr);
    d_hier_math_ops->strain_rate(d_strain_rate_idx, d_strain_rate_var, d_U_scratch_idx, d_U_var, d_no_fill_op, time);

    const double exponent = (d_power_law_index - 1.0) / d_yasuda_exponent;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> mu_data = patch->getPatchData(mu_idx);
            Pointer<CellData<NDIM, double>> EE_data = patch->getPatchData(d_strain_rate_idx);
            for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                // EE : EE counts each off-diagonal component twice.
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
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        if (allocate_U_scratch) level->deallocatePatchData(d_U_scratch_idx);
        level->deallocatePatchData(d_strain_rate_idx);
    }
    return;
} // computeViscosity

void
INSVCStaggeredCarreauYasudaHierarchyIntegrator::reset_viscosity(const int mu_idx,
                                                                Pointer<Variable<NDIM>> /*mu_var*/,
                                                                Pointer<HierarchyMathOps> /*hier_math_ops*/,
                                                                const int /*cycle_num*/,
                                                                const double time,
                                                                const double /*current_time*/,
                                                                const double /*new_time*/,
                                                                void* ctx)
{
    auto* const integrator = static_cast<INSVCStaggeredCarreauYasudaHierarchyIntegrator*>(ctx);
    const int U_idx = mu_idx == integrator->d_mu_current_idx ? integrator->d_U_current_idx : integrator->d_U_new_idx;
    integrator->computeViscosity(mu_idx, U_idx, time);
    return;
} // reset_viscosity

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
