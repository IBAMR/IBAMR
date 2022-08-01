// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/PhaseChangeDivUSourceFunction.h"
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

PhaseChangeDivUSourceFunction::PhaseChangeDivUSourceFunction(
    const std::string& object_name,
    const Pointer<PhaseChangeHierarchyIntegrator> pc_hier_integrator,
    const Pointer<HierarchyMathOps> hier_math_ops)
    : d_object_name(object_name), d_pc_hier_integrator(pc_hier_integrator), d_hier_math_ops(hier_math_ops)
{
    // intentionally blank
    return;
} // PhaseChangeDivUSourceFunction

bool
PhaseChangeDivUSourceFunction::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
PhaseChangeDivUSourceFunction::setDataOnPatch(const int data_idx,
                                              Pointer<Variable<NDIM> > var,
                                              Pointer<Patch<NDIM> > patch,
                                              const double data_time,
                                              const bool initial_time,
                                              Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    // set -Div U = -S where source term S is computed from PhaseChangeHierarchyIntegrator.
    const int S_idx = d_pc_hier_integrator->getDivergenceVelocitySourceTermIndex();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);

    // Subtract S_idx to the Div_U_F_idx.
    hier_cc_data_ops.axpy(data_idx, -1.0, S_idx, data_idx);

    return;
} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
