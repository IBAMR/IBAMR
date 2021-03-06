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

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

// Initialize the neighborhood of a circular interface.
struct CircularInterface
{
    IBTK::Vector X0;
    double R;
};
void circular_interface_neighborhood(int D_idx,
                                     SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                     double time,
                                     bool initial_time,
                                     void* ctx);
