// Filename LSLocateInterface.h
// Created on Oct 10, 2017 by Nishant Nangia

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
