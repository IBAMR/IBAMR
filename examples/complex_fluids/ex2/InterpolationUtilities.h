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

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAI_config.h"
#include "tbox/Pointer.h"

#include <Patch.h>
#include <PatchLevel.h>

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
                              IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > Q_var,
                              IBTK::SAMRAIPointer<SAMRAI::hier::PatchHierarchyNd> patch_hierarchy,
                              const std::vector<SAMRAI::solv::RobinBcCoefStrategyNd*>& bc_coefs,
                              const double data_time,
                              const int depth = 0);
    // Use least squares interpolation
    static double interpolateL2(const std::vector<double>& X,
                                const int data_idx,
                                IBTK::SAMRAIPointer<SAMRAI::pdat::CellVariableNd<double> > Q_var,
                                IBTK::SAMRAIPointer<SAMRAI::hier::PatchHierarchyNd> patch_hierarchy,
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategyNd*>& bc_coefs,
                                const double data_time,
                                const int depth = 0);
    static double weight_fcn(const std::vector<double>&, const std::vector<double>&);

private:
    static double interpolate(const double& x, const std::vector<int>& xi, const std::vector<double>& yi);

    static double interpolate_in_boxes(const CellIndexNd& idx,
                                       const std::vector<double>& X,
                                       SAMRAI::pdat::CellDataNd<int>& r_data,
                                       SAMRAI::pdat::CellDataNd<double>& q_data,
                                       SAMRAIPointer<CartesianPatchGeometryNd> pgeom,
                                       const BoxNd& pbox,
                                       int dim,
                                       int cycle,
                                       std::vector<int>& completed_dims);
};

} // namespace IBTK
#endif
