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

#include "ibtk/interpolation_utilities.h"

#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
IBTK_ENABLE_EXTRA_WARNINGS

#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchLevel.h>

#include <algorithm>

#include <ibtk/app_namespaces.h>

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace IBTK
{
namespace Interpolation
{
double
interpolate(const VectorNd& X,
            const int data_idx,
            Pointer<CellVariable<NDIM, double> > Q_var,
            Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            const int depth)
{
    double q_val = 0.0;
#ifndef NDEBUG
    // Check that we have enough ghost cells.
    Pointer<PatchDataFactory<NDIM> > Q_fac = patch_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_idx);
    TBOX_ASSERT(Q_fac->getGhostCellWidth().min() >= 1);
#endif
    bool done = false;
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_hierarchy->getGridGeometry();
    const double* const domain_low = grid_geom->getXLower();
    // Find this point in the patch hierarchy. Start looking on the finest level.
    for (int ln = patch_hierarchy->getFinestLevelNumber(); ln >= 0 && !done; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        // Get the cell index.
        CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, level->getGridGeometry(), level->getRatio());
        // Now search for the index.
        for (PatchLevel<NDIM>::Iterator p(level); p && !done; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* const dx = p_geom->getDx();
            const double* const x_lower = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& idx_low = patch_box.lower();
            Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
            if (patch_box.contains(idx))
            {
                // Great. This patch contains the index.
                // Interpolate to this point. Use bilinear interpolation
                // First determine "bottom left" index and the point in index space.
                hier::Index<NDIM> ll_idx;
                VectorNd X_idx;
                for (int d = 0; d < NDIM; ++d)
                {
                    X_idx[d] = (X[d] - domain_low[d]) / dx[d];
                    ll_idx(d) = std::floor((X[d] - domain_low[d]) / dx[d] - 0.5);
                }
                // Construct bilinear interpolant with Lagrange polynomial
                unsigned int row = 0;
#if (NDIM == 2)
                MatrixXd mat = MatrixXd::Ones(4, 4);
                VectorXd rhs = VectorXd::Zero(4);
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        IntVector<NDIM> ij(i, j);
                        hier::Index<NDIM> new_idx = ll_idx + ij;
                        VectorNd x_idx;
                        for (int d = 0; d < NDIM; ++d) x_idx[d] = new_idx(d) + 0.5;
                        // First column is ones
                        mat(row, 1) = x_idx(0) - X_idx[0];
                        mat(row, 2) = x_idx(1) - X_idx[1];
                        mat(row, 3) = (x_idx(0) - X_idx[0]) * (x_idx(1) - X_idx[1]);
                        rhs[row++] = (*data)(new_idx, depth);
                    }
                }
                q_val = mat.partialPivLu().solve(rhs)[0];
#endif
#if (NDIM == 3)
                MatrixXd mat = MatrixXd::Ones(8, 8);
                VectorXd rhs = VectorXd::Zero(8);
                for (int i = 0; i < 2; ++i)
                {
                    for (int j = 0; j < 2; ++j)
                    {
                        for (int k = 0; k < 2; ++k)
                        {
                            IntVector<NDIM> ijk(i, j, k);
                            hier::Index<NDIM> new_idx = ll_idx + ijk;
                            VectorNd x_idx;
                            for (int d = 0; d < NDIM; ++d) x_idx[d] = new_idx(d) + 0.5;
                            // First column is ones
                            mat(row, 1) = x_idx[0] - X_idx[0];
                            mat(row, 2) = x_idx[1] - X_idx[1];
                            mat(row, 3) = x_idx[2] - X_idx[2];
                            mat(row, 4) = (x_idx[0] - X_idx[0]) * (x_idx[1] - X_idx[1]);
                            mat(row, 5) = (x_idx[0] - X_idx[0]) * (x_idx[2] - X_idx[2]);
                            mat(row, 6) = (x_idx[1] - X_idx[1]) * (x_idx[2] - X_idx[2]);
                            mat(row, 7) = (x_idx[0] - X_idx[0]) * (x_idx[1] - X_idx[1]) * (x_idx[2] - X_idx[2]);
                            rhs[row++] = (*data)(ll_idx + ijk, depth);
                        }
                    }
                }
                q_val = mat.partialPivLu().solve(rhs)[0];
#endif
            }
        }
    }
    q_val = IBTK_MPI::sumReduction(q_val);
    return q_val;
}

double
interpolateL2(const VectorNd& X,
              const int data_idx,
              Pointer<CellVariable<NDIM, double> > Q_var,
              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
              const int stencil_width,
              const int poly_deg,
              const int depth,
              std::function<double(const VectorNd&, const VectorNd&)> wgt_fcn,
              const int indicator_idx)
{
    double q_val = 0.0;
    bool done = false;
    for (int ln = patch_hierarchy->getFinestLevelNumber(); ln >= 0 && !done; --ln)
    {
        // Start at the finest level...
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, level->getGridGeometry(), level->getRatio());
        for (PatchLevel<NDIM>::Iterator p(level); p && !done; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Pointer<CartesianPatchGeometry<NDIM> > p_geom = patch->getPatchGeometry();
            const double* const dx = p_geom->getDx();
            double dx_min = std::numeric_limits<double>::max();
            for (int d = 0; d < NDIM; ++d) dx_min = std::min(dx[d], dx_min);
            const double* const x_lower = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const hier::Index<NDIM>& idx_low = patch_box.lower();
            Pointer<CellData<NDIM, double> > data = patch->getPatchData(data_idx);
            Pointer<CellData<NDIM, double> > indicator_data =
                indicator_idx == IBTK::invalid_index ? nullptr : patch->getPatchData(indicator_idx);
            if (patch_box.contains(idx))
            {
                // Great. This patch contains the index.
                // Interpolate to this point. Use bilinear interpolation
                Box<NDIM> box(idx, idx);
                // Grow it by some number of grid cells
                box.grow(stencil_width);
                // Determine our stencil first.
                std::vector<VectorNd> interp_pts;
                std::vector<CellIndex<NDIM> > interp_idxs;
                for (CellIterator<NDIM> ci(box); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    const double indicator = indicator_data ? (*indicator_data)(idx) : -1.0;
                    if (indicator < 0.0)
                    {
                        VectorNd interp_pt;
                        for (int d = 0; d < NDIM; ++d)
                            interp_pt[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                        interp_pts.push_back(interp_pt);
                        interp_idxs.push_back(idx);
                    }
                }
                // Moving least squares with a polynomial basis. We shift the polynomials so that they are centered at
                // X. We scale the polynomials by min(dx). This shift means that the evaluation of the resulting
                // polynomial is given by the first point of the solution.
                VectorXd rhs = VectorXd::Zero(interp_pts.size());
                VectorXd wgts = VectorXd::Zero(interp_pts.size());
                for (size_t i = 0; i < interp_pts.size(); ++i)
                {
                    rhs(i) = (*data)(interp_idxs[i]);
                    wgts(i) = wgt_fcn(interp_pts[i], X);
                }
                DiagonalMatrix<double, Dynamic> W(wgts);
                MatrixXd mat = Interpolation::formMonomialBasis(interp_pts, poly_deg, dx_min, X);
                // Solve least squares system
                VectorXd soln = (W * mat).colPivHouseholderQr().solve(W * rhs);
                q_val = soln(0);
                done = true;
            }
        }
    }
    q_val = IBTK_MPI::sumReduction(q_val);
    return q_val;
}

double
default_wgt_fcn(const VectorNd&, const VectorNd&)
{
    return 1.0;
}
} // namespace Interpolation
} // namespace IBTK
