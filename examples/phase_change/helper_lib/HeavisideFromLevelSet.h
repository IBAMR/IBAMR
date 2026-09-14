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

#ifndef included_PhaseChangeExamples_HeavisideFromLevelSet
#define included_PhaseChangeExamples_HeavisideFromLevelSet

#include <ibamr/config.h>

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
// Synchronize the Heaviside field with the integrator's current advected level set.
class HeavisideFromLevelSet
{
public:
    HeavisideFromLevelSet(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                          double num_interface_cells);
    static void synchronize_levelset_with_heaviside_fcn(int H_current_idx,
                                                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                        int integrator_step,
                                                        double time,
                                                        bool initial_time,
                                                        bool regrid_time,
                                                        void* ctx);

private:
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_ls_var;
    double d_num_interface_cells;
};

} // namespace PhaseChangeExamples

#endif
