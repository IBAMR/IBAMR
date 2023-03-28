// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_TagInterfaceRefinementCells
#define included_TagInterfaceRefinementCells

///////////////////////////// INCLUDES ///////////////////////////////////
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
}
namespace pdat
{
template <int DIM, class TYPE>
class CellVariable;
}
} // namespace SAMRAI

namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
}

#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

void callTagInterfaceRefinementCellsCallbackFunction(
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too,
    void* ctx);

struct TagInterfaceRefinementCells
{
    TagInterfaceRefinementCells(SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_gradient_var,
                                double tag_min_val,
                                double tag_max_val);
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_gradient_var;
    double d_tag_min_val;
    double d_tag_max_val;
};

#endif // #ifndef included_TagInterfaceRefinementCells
