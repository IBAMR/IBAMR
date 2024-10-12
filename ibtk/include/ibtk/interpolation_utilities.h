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

#ifndef included_IBTK_interpolation_utilities
#define included_IBTK_interpolation_utilities

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

namespace IBTK
{
namespace Interpolation
{
/*
 * Interpolations the cell centered data stored in patch index data_idx to the point X. Uses bilinear interpolation.
 */
double interpolate(const VectorNd& X,
                   int data_idx,
                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                   int depth = 0);

double default_wgt_fcn(const VectorNd&, const VectorNd&);
/*
 * Interpolations the cell centered data stored in patch index data_idx to the point X. Uses moving least squares with
 * an optional weighting function.
 */
double interpolateL2(const VectorNd& X,
                     int data_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                     int stencil_width,
                     int poly_deg,
                     int depth = 0,
                     std::function<double(const VectorNd&, const VectorNd&)> wgt_fcn = default_wgt_fcn,
                     int indicator_idx = IBTK::invalid_index);

/*!
 * Evaluates the monomials up to a specified degree at the list of points provided. Returns a matrix where each row
 * consists of the monomials evaluated at a given point.
 *
 * Must specify a shift and scaling of the monomials.
 */
template <typename VectorArray>
IBTK::MatrixXd formMonomialBasis(const std::vector<VectorArray>& pts, int deg, double ds, const VectorArray& shft);

/*!
 * Returns the number of polynomials with total degree deg.
 */
int getNumberOfPolynomials(int deg);

/*!
 * Returns q^i. For the special cases when q = 0.0 or i = 0, this implementation returns 0.0 or 1.0 respectively.
 * Otherwise this returns std::pow(q,i).
 */
double pow(const double q, const int i);
} // namespace Interpolation

} // namespace IBTK

#include "ibtk/private/interpolation_utilities_inc.h"
#endif
