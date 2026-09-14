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

#ifndef included_LSLocateInterface
#define included_LSLocateInterface

#include <ibamr/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/CartGridFunction.h>

namespace IBTK
{
class HierarchyMathOps;
}

void call_ls_locate_interface_callback(int D_idx,
                                       SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                       double time,
                                       bool initial_time,
                                       void* ctx);

// Supply the initial level set to LSInitStrategy, then preserve the integrator's
// advected level set on later interface resets.
class LSLocateInterface
{
public:
    LSLocateInterface(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                      SAMRAI::tbox::Pointer<IBTK::CartGridFunction> initial_condition);

    void setLevelSetPatchData(int D_idx,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              double time,
                              bool initial_time);

private:
    LSLocateInterface(const LSLocateInterface& from) = delete;
    LSLocateInterface& operator=(const LSLocateInterface& that) = delete;

    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_ls_var;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_initial_condition;
};

#endif
