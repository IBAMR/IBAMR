// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2023 by the IBAMR developers
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

#include "ibtk/samrai_compatibility_names.h"
// SAMRAI INCLUDES
#include "SAMRAIBasePatchHierarchy.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIPointer.h"
#include <tbox/DescribedClass.h>

void callTagInterfaceRefinementCellsCallbackFunction(SAMRAIPointer<SAMRAIBasePatchHierarchy> hierarchy,
                                                     const int level_number,
                                                     const double error_data_time,
                                                     const int tag_index,
                                                     const bool initial_time,
                                                     const bool uses_richardson_extrapolation_too,
                                                     void* ctx);

struct TagInterfaceRefinementCells
{
    TagInterfaceRefinementCells(SAMRAIPointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                                SAMRAIPointer<SAMRAICellVariable<double>> scalar_var,
                                double tag_min_val,
                                double tag_max_val);
    SAMRAIPointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAIPointer<SAMRAICellVariable<double>> d_scalar_var;
    double d_tag_min_val;
    double d_tag_max_val;
};

#endif // #ifndef included_TagInterfaceRefinementCells
