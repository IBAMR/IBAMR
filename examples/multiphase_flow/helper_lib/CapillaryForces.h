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

#ifndef included_MultiphaseExamples_CapillaryForces
#define included_MultiphaseExamples_CapillaryForces

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibtk/HierarchyMathOps.h>

#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <Patch.h>
#include <RobinBcCoefStrategy.h>
#include <VariableContext.h>

namespace MultiphaseExamples
{
// Scale each patch-local side. A negative liquid-fraction index selects density-only masking.
void mask_force(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                int F_idx,
                int rho_idx,
                int lf_scratch_idx,
                double rho_liquid,
                double rho_gas);

class DensityForceMask
{
public:
    DensityForceMask(SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> integrator,
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
    SAMRAI::tbox::Pointer<IBAMR::INSVCStaggeredHierarchyIntegrator> d_integrator;
    double d_rho_liquid, d_rho_gas;
};

// The supplied context identifies the temperature data actually consumed by the force calculation.
class SurfaceTensionCoefficients
{
public:
    SurfaceTensionCoefficients(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> temperature,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context,
                               double sigma0,
                               double dsigma_dT0,
                               double T_ref);
    static void compute_surface_tension_coef_function(int F_idx,
                                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                      int integrator_step,
                                                      double time,
                                                      double current_time,
                                                      double new_time,
                                                      void* ctx);
    static void compute_marangoni_coef_function(int F_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                int integrator_step,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

private:
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_temperature;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    double d_sigma0, d_dsigma_dT0, d_T_ref;
};

} // namespace MultiphaseExamples

#endif
