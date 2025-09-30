// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_interpolation_utilities
#define included_IBTK_interpolation_utilities

#include "ibtk/config.h"

#include "ibtk/ibtk_utilities.h"

#include "Patch.h"
#include "PatchHierarchy.h"
#include "Variable.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBTK
{
/*
 * Interpolates the data stored in data_idx to the location X using the provided kernel function. The position location
 * must be synchronized across all MPI ranks prior to this function call
 *
 * Q_var must correspond to Cell, Side, Node, or Edge centered double values, or an unrecoverable error will occur.
 * Ghost cells must be filled in for the data stored in data_idx prior to this function call.
 *
 * The returned values are synchronized across all processors.
 */
std::vector<double> interpolate(const VectorNd& X,
                                int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > Q_var,
                                int Q_depth,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                std::string kernel_fcn = "IB_4");

/*
 * Interpolates the data stored in data_idx to the locations provided in X using the provided kernel function. The
 * position locations in X must be synchronized across all MPI ranks to this function call.
 *
 * Q_var must correspond to Cell, Side, Node, or Edge centered double values, or an unrecoverable error will occur.
 * Ghost cells must be filled in for the data stored in data_idx prior to this function call.
 *
 * The returned values are synchronized across all processors.
 */
std::vector<double> interpolate(const std::vector<VectorNd>& X,
                                int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > Q_var,
                                int Q_depth,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                std::string kernel_fcn = "IB_4");

} // namespace IBTK
#endif
