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

#ifndef included_example_LSLocateInterface
#define included_example_LSLocateInterface

#include <ibamr/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/CartGridFunction.h>

namespace IBTK
{
class HierarchyMathOps;
}

// Supply the initial interface to LSInitStrategy, and preserve the advected
// interface on subsequent resets. The initializer defines the example geometry.
class LSLocateInterface
{
public:
    /*! \brief Retain the integrator, advected variable, and initial geometry. */
    LSLocateInterface(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> integrator,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                      SAMRAI::tbox::Pointer<IBTK::CartGridFunction> initial_conditions);

    /*! \brief Initialize the geometry or copy the current advected field. */
    void setLevelSetPatchData(int data_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              double time,
                              bool initial_time);

private:
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_ls_var;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_initial_conditions;
};

void call_locate_interface(int data_idx,
                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                           double time,
                           bool initial_time,
                           void* ctx);

#endif
