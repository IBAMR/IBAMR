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

#ifndef included_PhaseChangeExamples_LiquidFractionForceMask
#define included_PhaseChangeExamples_LiquidFractionForceMask

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibtk/HierarchyMathOps.h>

#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>
#include <VariableContext.h>

namespace PhaseChangeExamples
{
// Fill the selected new liquid fraction into scratch before masking the interfacial force.
// The selected variable may be the integrator's extrapolated liquid fraction.
class LiquidFractionForceMask
{
public:
    LiquidFractionForceMask(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> liquid_fraction,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                            SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_integrator,
                            SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> ins_integrator,
                            double rho_liquid,
                            double rho_gas);
    static void mask_surface_tension_force(int F_idx,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           int integrator_step,
                                           double time,
                                           double current_time,
                                           double new_time,
                                           void* ctx);

private:
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_liquid_fraction;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_bc_coef;
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_integrator;
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_ins_integrator;
    double d_rho_liquid, d_rho_gas;
};

} // namespace PhaseChangeExamples

#endif
