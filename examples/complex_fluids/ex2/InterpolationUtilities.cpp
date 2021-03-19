// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibtk/ibtk_macros.h"

#include "InterpolationUtilities.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <Eigen/Core>
#include <Eigen/QR>
IBTK_ENABLE_EXTRA_WARNINGS

#include <HierarchyCellDataOpsReal.h>
#include <Patch.h>
#include <PatchLevel.h>

#include <algorithm>

using Eigen::Matrix3d;
using Eigen::Vector3d;

namespace IBTK
{
double
InterpolationUtilities::interpolate(const vector<double>& X,
                                    const int data_idx,
                                    Pointer<CellVariable<NDIM, double> > Q_var,
                                    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                    const double data_time,
                                    const int depth)
{
    double q_val = 0.0;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int data_idx_temp =
        var_db->registerVariableAndContext(Q_var, var_db->getContext("Interpolation"), IntVector<NDIM>(3));
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(data_idx_temp);
    }
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy);
    hier_cc_data_ops.copyData(data_idx_temp, data_idx);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        data_idx_temp, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, bc_coefs, NULL);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
    ghost_fill_op.fillData(data_time);
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
            const double* const x_lower = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > S_data = patch->getPatchData(data_idx_temp);
            if (patch_box.contains(idx))
            {
                // Great. The patch is currently on this level
                // Let's create a box that contains this data
                Box<NDIM> box(idx, idx);
                // Grow it by some number of grid cells
                box.grow(IntVector<NDIM>(3));
                // Loop through the box, make sure the point is located
                // OUTSIDE the disk
                std::vector<double> x(NDIM);
                CellData<NDIM, int> i_data(box, 1, IntVector<NDIM>(0));
                CellData<NDIM, double> si_data(box, NDIM + 1, IntVector<NDIM>(0));
                si_data.fillAll(std::numeric_limits<double>::signaling_NaN());
                const CellIndex<NDIM> ci_l = patch_box.lower();
                int num = 0;
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
                    if (r > 1.0)
                    {
                        i_data(ci) = 1;
                        si_data(ci, NDIM) = (*S_data)(ci, depth);
                        num++;
                    }
                    else
                    {
                        i_data(ci) = -1;
                    }
                }
                // We have a box containing the point and data the says if we are inside or outside disk.
                // Find directions and interpolate. First in x, then in y.
                std::vector<int> completed_dims;
                q_val = InterpolationUtilities::interpolate_in_boxes(
                    idx, X, i_data, si_data, p_geom, patch_box, 0, 0, completed_dims);
                done = true;
            }
        }
    }
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(data_idx_temp);
    }
    q_val = IBTK_MPI::sumReduction(q_val);
    return q_val;
}

double
InterpolationUtilities::interpolateL2(const std::vector<double>& X,
                                      const int data_idx,
                                      Pointer<CellVariable<NDIM, double> > Q_var,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                      const double data_time,
                                      const int depth)
{
    double q_val = 0.0;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int data_idx_temp =
        var_db->registerVariableAndContext(Q_var, var_db->getContext("Interpolation"), IntVector<NDIM>(3));
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(data_idx_temp);
    }
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy);
    hier_cc_data_ops.copyData(data_idx_temp, data_idx);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
    ghost_cell_components[0] = InterpolationTransactionComponent(
        data_idx_temp, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, bc_coefs, NULL);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
    ghost_fill_op.fillData(data_time);
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
            const double* const x_lower = p_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > S_data = patch->getPatchData(data_idx_temp);
            if (patch_box.contains(idx))
            {
                // Great. The patch is currently on this level
                // Let's create a box that contains this data
                Box<NDIM> box(idx, idx);
                // Grow it by some number of grid cells
                box.grow(2);
                // Loop through the box, make sure the point is located
                // OUTSIDE the disk
                std::vector<double> x(NDIM);
                CellData<NDIM, int> i_data(box, 1, IntVector<NDIM>(0));
                const CellIndex<NDIM> ci_l = patch_box.lower();
                int num = 0;
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                    double r = sqrt(x[0] * x[0] + x[1] * x[1]);
                    if (r > 1.0)
                    {
                        i_data(ci) = 1;
                        num++;
                    }
                }
                // Moving least squares
                VectorXd rhs = VectorXd::Zero(6);
                VectorXd soln = VectorXd::Zero(6);
                MatrixXd mat = MatrixXd::Zero(6, 6);
                for (CellIterator<NDIM> i(box); i; i++)
                {
                    CellIndex<NDIM> ci = i();
                    if (i_data(ci) == 1)
                    {
                        for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ci(d) - ci_l(d) + 0.5);
                        // Fill in RHS and MATRIX values
                        double w = weight_fcn(X, x);
                        double f = (*S_data)(ci, depth);
                        rhs(0) += w * f;
                        rhs(1) += w * f * x[0];
                        rhs(2) += w * f * x[1];
                        rhs(3) += w * f * x[0] * x[0];
                        rhs(4) += w * f * x[1] * x[1];
                        rhs(5) += w * f * x[0] * x[1];
                        mat(0, 0) += w; /*mat(0,1) += x[0]*w; mat(0,2) += x[1]*w; mat(0,3) += x[0]*x[0]*w; mat(0,4) +=
                                           x[1]*x[1]*w; mat(0,5) += x[0]*x[1]*w;*/
                        mat(1, 0) += x[0] * w;
                        mat(1, 1) += x[0] * x[0] * w; /* mat(1,2) += x[0]*x[1]*w; mat(2,3) += x[0]*x[0]*x[0]*w; mat(2,4)
                                                         += x[0]*x[1]*x[1]*w; mat(3,5) += x[0]*x[0]*x[1]*w;*/
                        mat(2, 0) += x[1] * w;
                        mat(2, 1) += x[1] * x[0] * w;
                        mat(2, 2) += x[1] * x[1] * w;
                        mat(3, 0) += x[0] * x[0] * w;
                        mat(3, 1) += x[0] * x[0] * x[0] * w;
                        mat(3, 2) += x[0] * x[0] * x[1] * w;
                        mat(3, 3) += x[0] * x[0] * x[0] * x[0] * w;
                        mat(4, 0) += x[1] * x[1] * w;
                        mat(4, 1) += x[1] * x[1] * x[0] * w;
                        mat(4, 2) += x[1] * x[1] * x[1] * w;
                        mat(4, 3) += x[1] * x[1] * x[0] * x[0] * w;
                        mat(4, 4) += x[1] * x[1] * x[1] * x[1] * w;
                        mat(5, 0) += x[0] * x[1] * w;
                        mat(5, 1) += x[0] * x[1] * x[0] * w;
                        mat(5, 2) += x[0] * x[1] * x[1] * w;
                        mat(5, 3) += x[0] * x[1] * x[0] * x[0] * w;
                        mat(5, 4) += x[0] * x[1] * x[1] * x[1] * w;
                        mat(5, 5) += x[0] * x[1] * x[0] * x[1] * w;
                    }
                }
                //                mat(0,1) = mat(1,0); mat(0,2) = mat(2,0); mat(0,3) = mat(3,0); mat(0,4) = mat(4,0);
                //                mat(0,5) = mat(5,0); mat(1,2) = mat(2,1); mat(1,3) = mat(3,1); mat(1,4) = mat(4,1);
                //                mat(1,5) = mat(5,1); mat(2,3) = mat(3,2); mat(2,4) = mat(4,2); mat(2,5) = mat(5,2);
                //                mat(3,4) = mat(4,3); mat(3,5) = mat(5,3);
                //                mat(4,5) = mat(5,4);
                soln = mat.ldlt().solve(rhs);
                q_val = soln(0) + soln(1) * X[0] + soln(2) * X[1] + soln(3) * X[0] * X[0] + soln(4) * X[1] * X[1] +
                        soln(5) * X[0] * X[1];
                done = true;
            }
        }
    }
    q_val = IBTK_MPI::sumReduction(q_val);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(data_idx_temp);
    }
    return q_val;
}

double
InterpolationUtilities::weight_fcn(const std::vector<double>& x, const std::vector<double>& x_j)
{
    double r = sqrt((x[0] - x_j[0]) * (x[0] - x_j[0]) + (x[1] - x_j[1]) * (x[1] - x_j[1]));
    double w = exp(-10.0 * r * r);
    w = 1.0;
    return w;
}

double
InterpolationUtilities::interpolate(const double& l, const std::vector<int>& i, const std::vector<double>& y)
{
    double val = 0.0;
    TBOX_ASSERT(i.size() == 3);
    TBOX_ASSERT(y.size() == i.size());
    val = y[0] + (y[1] - y[0]) * l / static_cast<double>(i[1] - i[0]) +
          (y[2] * static_cast<double>(i[0] - i[1]) + y[0] * static_cast<double>(i[1] - i[2]) +
           y[1] * static_cast<double>(i[2] - i[0])) *
              l * (l - 1.0) / (static_cast<double>((i[0] - i[1]) * (i[0] - i[2]) * (i[1] - i[2])));
    return val;
}

double
InterpolationUtilities::interpolate_in_boxes(const CellIndex<NDIM>& idx,
                                             const std::vector<double>& X,
                                             CellData<NDIM, int>& r_data,
                                             CellData<NDIM, double>& q_data,
                                             Pointer<CartesianPatchGeometry<NDIM> > pgeom,
                                             const Box<NDIM>& pbox,
                                             int dim,
                                             int cycle,
                                             std::vector<int>& completed_dims)
{
    double q_val = 0.0;
    // form list of indices.
    std::vector<CellIndex<NDIM> > idx_list;
    std::vector<int> i_list;
    bool done = false;
    while (!done)
    {
        idx_list.clear();
        i_list.clear();
        dim = dim % NDIM;
        if (r_data(idx) == 1)
        {
            idx_list.push_back(idx);
            i_list.push_back(0);
        }
        int s = 1;
        while (idx_list.size() < 3)
        {
            IntVector<NDIM> si(0);
            si(dim) = s;
            if (r_data(idx + si) == 1)
            {
                idx_list.push_back(idx + si);
                i_list.push_back(s);
            }
            if (r_data(idx - si) == 1)
            {
                idx_list.push_back(idx - si);
                i_list.push_back(-s);
            }
            s++;
        }
        int min_val = *std::min_element(i_list.begin(), i_list.end());
        int max_val = *std::max_element(i_list.begin(), i_list.end());
        if (std::max(std::abs(min_val), std::abs(max_val)) > 3)
        {
            dim++;
            if (std::find(completed_dims.begin(), completed_dims.end(), dim) != completed_dims.end() || dim > NDIM)
            {
                TBOX_ERROR("already completed dimension, or dimension too high");
            }
        }
        else
        {
            done = true;
        }
    }
    while (idx_list.size() > 3)
    {
        idx_list.pop_back();
    }
    completed_dims.push_back(dim);
    // We start at dim = 0. Each call generates a new box to interpolate inside of, and increases the dim by 1.
    if (cycle < NDIM - 1)
    {
        std::vector<int> i_list;
        std::vector<double> y_data;
        for (std::vector<CellIndex<NDIM> >::const_iterator cit = idx_list.begin(); cit != idx_list.end(); ++cit)
        {
            const CellIndex<NDIM>& cidx = *cit;
            q_data(cidx, cycle + 1) = InterpolationUtilities::interpolate_in_boxes(
                cidx, X, r_data, q_data, pgeom, pbox, dim + 1, cycle + 1, completed_dims);
            i_list.push_back(cidx(dim));
            y_data.push_back(q_data(cidx, cycle + 1));
        }
        const double* dx = pgeom->getDx();
        const double* xlow = pgeom->getXLower();
        const CellIndex<NDIM>& idxl = pbox.lower();
        double xx = xlow[dim] + dx[dim] * (idx_list[0](dim) - idxl(dim) + 0.5);
        q_val = InterpolationUtilities::interpolate((X[dim] - xx) / dx[dim], i_list, y_data);
    }
    else if (cycle == NDIM - 1)
    {
        std::vector<int> i_list;
        std::vector<double> y_data;
        for (std::vector<CellIndex<NDIM> >::const_iterator cit = idx_list.begin(); cit != idx_list.end(); ++cit)
        {
            const CellIndex<NDIM>& cidx = *cit;
            i_list.push_back(cidx(dim));
            y_data.push_back(q_data(cidx, cycle + 1));
        }
        const double* dx = pgeom->getDx();
        const double* xlow = pgeom->getXLower();
        const CellIndex<NDIM>& idxl = pbox.lower();
        double xx = xlow[dim] + dx[dim] * (idx_list[0](dim) - idxl(dim) + 0.5);
        q_val = InterpolationUtilities::interpolate((X[dim] - xx) / dx[dim], i_list, y_data);
    }
    else
    {
        // Shouldn't get here.
        TBOX_ERROR("Invalid dimension.\n");
    }
    return q_val;
}
} // namespace IBTK
