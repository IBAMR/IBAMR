// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_InterpolationUtilities
#define included_InterpolationUtilities

// SAMRAI INCLUDES
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIIndex.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAI_config.h"

#include <vector>

#include <ibamr/app_namespaces.h>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

namespace IBTK
{
class InterpolationUtilities
{
public:
    // Use polynomial interpolation
    static double interpolate(const std::vector<double>& X,
                              const int data_idx,
                              SAMRAIPointer<SAMRAICellVariable<double>> Q_var,
                              SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy,
                              const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
                              const double data_time,
                              const int depth = 0);
    // Use least squares interpolation
    static double interpolateL2(const std::vector<double>& X,
                                const int data_idx,
                                SAMRAIPointer<SAMRAICellVariable<double>> Q_var,
                                SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy,
                                const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
                                const double data_time,
                                const int depth = 0);
    static double weight_fcn(const std::vector<double>&, const std::vector<double>&);

private:
    static double interpolate(const double& x, const std::vector<int>& xi, const std::vector<double>& yi);

    static double interpolate_in_boxes(const SAMRAICellIndex& idx,
                                       const std::vector<double>& X,
                                       SAMRAICellData<int>& r_data,
                                       SAMRAICellData<double>& q_data,
                                       SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom,
                                       const SAMRAIBox& pbox,
                                       int dim,
                                       int cycle,
                                       std::vector<int>& completed_dims);
};

} // namespace IBTK
#endif
