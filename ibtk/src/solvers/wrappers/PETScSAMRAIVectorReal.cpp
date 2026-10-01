// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/CartesianCentering.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/NormOps.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>
#include <tbox/Timer.h>

#include <petscis.h>
#include <petscvec.h>

#include <ArrayData.h>
#include <Box.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>
#include <VariableDatabase.h>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_vec_duplicate;
static Timer* t_vec_dot;
static Timer* t_vec_m_dot;
static Timer* t_vec_norm;
static Timer* t_vec_t_dot;
static Timer* t_vec_m_t_dot;
static Timer* t_vec_scale;
static Timer* t_vec_copy;
static Timer* t_vec_set;
static Timer* t_vec_swap;
static Timer* t_vec_axpy;
static Timer* t_vec_axpby;
static Timer* t_vec_maxpy;
static Timer* t_vec_aypx;
static Timer* t_vec_waxpy;
static Timer* t_vec_axpbypcz;
static Timer* t_vec_pointwise_mult;
static Timer* t_vec_pointwise_divide;
static Timer* t_vec_get_size;
static Timer* t_vec_get_local_size;
static Timer* t_vec_max;
static Timer* t_vec_min;
static Timer* t_vec_set_random;
static Timer* t_vec_destroy;
static Timer* t_vec_dot_local;
static Timer* t_vec_t_dot_local;
static Timer* t_vec_norm_local;
static Timer* t_vec_m_dot_local;
static Timer* t_vec_m_t_dot_local;
static Timer* t_vec_max_pointwise_divide;
static Timer* t_vec_dot_norm2;

// Fused multi-vector operations for Krylov orthogonalization. Each follows the loop and accumulation order of the
// equivalent sequence of single-vector SAMRAI operations, so its results are bitwise identical to that sequence, but
// it reads the shared vector once per group of vectors instead of once per vector.

// Number of vectors combined in one traversal of the data; wider groups stream too many arrays at once.
constexpr int FUSED_GROUP_SIZE = 8;

// Number of vectors combined in registers while traversing one row; more independent accumulators cost more
// instructions without a benefit in cycles.
constexpr int FUSED_BLOCK = 4;

// Return the centering of the double-valued patch data allocated at the patch data index, if fused operations support
// it.
std::optional<DataCentering>
find_component_centering(const int idx)
{
    return find_data_centering<double>(
        *VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(idx));
}

// Call fn(index of the row start, offset in data_box, offset in other_box, row length) for each row of ibox, in the
// order of SAMRAI's ArrayData loops: rows in lexicographic order with the lowest remaining index fastest.
template <class RowFn>
void
for_each_row(const Box<NDIM>& ibox, const Box<NDIM>& data_box, const Box<NDIM>& other_box, RowFn&& fn)
{
    const int n0 = ibox.numberCells(0);
    Index<NDIM> idx = ibox.lower();
    while (true)
    {
        fn(idx, data_box.offset(idx), other_box.offset(idx), n0);
        int j = 1;
        for (; j < NDIM; ++j)
        {
            if (idx(j) < ibox.upper(j))
            {
                ++idx(j);
                break;
            }
            idx(j) = ibox.lower(j);
        }
        if (j == NDIM)
        {
            break;
        }
    }
}

// Add the dot products of one row of x with the rows y[0], ..., y[N - 1], weighted by the control volume cv if WITH_CV,
// to acc[0], ..., acc[N - 1], holding the sums in registers. Each acc[i] receives its terms in the same order, and
// with the same expression, as ArrayDataNormOpsReal::dot() or ArrayDataNormOpsReal::dotWithControlVolume().
template <int N, bool WITH_CV>
void
mdot_row(const double* const x, const double* const* const y, const double* const cv, const int n, double* const acc)
{
    const double* yp[N];
    double a[N];
    for (int i = 0; i < N; ++i)
    {
        yp[i] = y[i];
        a[i] = acc[i];
    }
    for (int k = 0; k < n; ++k)
    {
        const double x_val = x[k];
        if (WITH_CV)
        {
            const double cv_val = cv[k];
            for (int i = 0; i < N; ++i)
            {
                a[i] += x_val * yp[i][k] * cv_val;
            }
        }
        else
        {
            for (int i = 0; i < N; ++i)
            {
                a[i] += x_val * yp[i][k];
            }
        }
    }
    for (int i = 0; i < N; ++i)
    {
        acc[i] = a[i];
    }
}

// Call mdot_row() with the block size ny, which is at most N.
template <bool WITH_CV, int N>
void
mdot_row_tail(const double* const x,
              const double* const* const y,
              const int ny,
              const double* const cv,
              const int n,
              double* const acc)
{
    if constexpr (N > 0)
    {
        if (ny == N)
        {
            mdot_row<N, WITH_CV>(x, y, cv, n, acc);
        }
        else
        {
            mdot_row_tail<WITH_CV, N - 1>(x, y, ny, cv, n, acc);
        }
    }
}

// Add the dot products of one row of x with the rows y[0], ..., y[ny - 1] to acc, FUSED_BLOCK vectors at a time.
template <bool WITH_CV>
void
mdot_row_blocks(const double* const x,
                const double* const* const y,
                const int ny,
                const double* const cv,
                const int n,
                double* const acc)
{
    int i = 0;
    for (; i + FUSED_BLOCK <= ny; i += FUSED_BLOCK)
    {
        mdot_row<FUSED_BLOCK, WITH_CV>(x, y + i, cv, n, acc + i);
    }
    mdot_row_tail<WITH_CV, FUSED_BLOCK - 1>(x, y + i, ny - i, cv, n, acc + i);
}

// Compute dprod[i] = sum over box of x * y[i] (* cv), accumulated as in ArrayDataNormOpsReal::dot() and
// ArrayDataNormOpsReal::dotWithControlVolume(). All y arrays must share x's box and depth.
void
mdot_arrays(const ArrayData<NDIM, double>& x,
            const std::vector<const ArrayData<NDIM, double>*>& y,
            const ArrayData<NDIM, double>* cv,
            const Box<NDIM>& box,
            double* dprod)
{
    const int ny = static_cast<int>(y.size());
    for (int i = 0; i < ny; ++i)
    {
        dprod[i] = 0.0;
    }
    const Box<NDIM>& x_box = x.getBox();
    const Box<NDIM> ibox = cv ? box * x_box * x_box * cv->getBox() : box * x_box * x_box;
    if (ibox.empty())
    {
        return;
    }
    const Box<NDIM>& cv_box = cv ? cv->getBox() : x_box;
    const int cv_depth = cv ? cv->getDepth() : 1;
    std::vector<const double*> y_row(ny);
    for (int d = 0; d < x.getDepth(); ++d)
    {
        const int x_depth_offset = d * x.getOffset();
        const int cv_depth_offset = (cv && cv_depth != 1) ? d * cv->getOffset() : 0;
        for_each_row(ibox,
                     x_box,
                     cv_box,
                     [&](const Index<NDIM>& /*idx*/, const int x_offset, const int cv_offset, const int n)
                     {
                         const double* const x_row = x.getPointer() + x_depth_offset + x_offset;
                         for (int i = 0; i < ny; ++i)
                         {
                             y_row[i] = y[i]->getPointer() + x_depth_offset + x_offset;
                         }
                         if (cv)
                         {
                             mdot_row_blocks<true>(
                                 x_row, y_row.data(), ny, cv->getPointer() + cv_depth_offset + cv_offset, n, dprod);
                         }
                         else
                         {
                             mdot_row_blocks<false>(x_row, y_row.data(), ny, nullptr, n, dprod);
                         }
                     });
    }
}

// Apply y = a[i] * x[i] + y for i = 0, ..., N - 1 in order to one row, for coefficients other than 0 and +/-1.
template <int N>
void
maxpy_row_general(double* __restrict__ const y, const double* const a, const double* const* const x, const int n)
{
    const double* xp[N];
    double alpha[N];
    for (int i = 0; i < N; ++i)
    {
        xp[i] = x[i];
        alpha[i] = a[i];
    }
    for (int k = 0; k < n; ++k)
    {
        double val = y[k];
        for (int i = 0; i < N; ++i)
        {
            val = alpha[i] * xp[i][k] + val;
        }
        y[k] = val;
    }
}

// Call maxpy_row_general() with the block size nx, which is at most N.
template <int N>
void
maxpy_row_dispatch(double* const y, const double* const a, const double* const* const x, const int nx, const int n)
{
    if constexpr (N > 0)
    {
        if (nx == N)
        {
            maxpy_row_general<N>(y, a, x, n);
        }
        else
        {
            maxpy_row_dispatch<N - 1>(y, a, x, nx, n);
        }
    }
}

// Apply y = a * x + y to one row as ArrayDataBasicOps::axpy() does, including its special cases.
void
axpy_row(double* __restrict__ const y, const double a, const double* __restrict__ const x, const int n)
{
    if (a == 0.0)
    {
        return;
    }
    if (a == 1.0)
    {
        for (int k = 0; k < n; ++k)
        {
            y[k] = x[k] + y[k];
        }
    }
    else if (a == -1.0)
    {
        for (int k = 0; k < n; ++k)
        {
            y[k] = y[k] - x[k];
        }
    }
    else
    {
        for (int k = 0; k < n; ++k)
        {
            y[k] = a * x[k] + y[k];
        }
    }
}

// Compute y += sum_i alpha[i] * x[i] over box, applying the terms to each value in order as ArrayDataBasicOps::axpy()
// would. All x arrays must share y's box and depth. If norm_box and norm_sum are not null, also set *norm_sum to the
// sum of y * y (* cv) over norm_box, accumulated as in ArrayDataNormOpsReal::dot() or
// ArrayDataNormOpsReal::dotWithControlVolume() once each row holds its final values; box must contain norm_box.
void
maxpy_arrays(ArrayData<NDIM, double>& y,
             const std::vector<double>& alpha,
             const std::vector<const ArrayData<NDIM, double>*>& x,
             const Box<NDIM>& box,
             const ArrayData<NDIM, double>* const cv,
             const Box<NDIM>* const norm_box,
             double* const norm_sum)
{
    if (norm_sum)
    {
        *norm_sum = 0.0;
    }
    const int nx = static_cast<int>(x.size());
    const Box<NDIM>& y_box = y.getBox();
    const Box<NDIM> ibox = box * y_box * y_box * y_box;
    if (ibox.empty())
    {
        return;
    }
    const Box<NDIM> norm_ibox = norm_box ? (cv ? *norm_box * y_box * cv->getBox() : *norm_box * y_box) : Box<NDIM>();
    const bool accumulate_norm = norm_sum && !norm_ibox.empty();
    const int cv_depth = cv ? cv->getDepth() : 1;
    std::vector<const double*> x_row(nx);
    for (int d = 0; d < y.getDepth(); ++d)
    {
        const int depth_offset = d * y.getOffset();
        const int cv_depth_offset = (cv && cv_depth != 1) ? d * cv->getOffset() : 0;
        for_each_row(ibox,
                     y_box,
                     y_box,
                     [&](const Index<NDIM>& idx, const int offset, const int /*other_offset*/, const int n)
                     {
                         double* const y_row = y.getPointer() + depth_offset + offset;
                         for (int i = 0; i < nx; ++i)
                         {
                             x_row[i] = x[i]->getPointer() + depth_offset + offset;
                         }
                         int i = 0;
                         while (i < nx)
                         {
                             const int block = std::min(FUSED_BLOCK, nx - i);
                             bool general = true;
                             for (int j = i; j < i + block; ++j)
                             {
                                 general = general && alpha[j] != 0.0 && alpha[j] != 1.0 && alpha[j] != -1.0;
                             }
                             if (general)
                             {
                                 maxpy_row_dispatch<FUSED_BLOCK>(y_row, alpha.data() + i, x_row.data() + i, block, n);
                             }
                             else
                             {
                                 for (int j = i; j < i + block; ++j)
                                 {
                                     axpy_row(y_row, alpha[j], x_row[j], n);
                                 }
                             }
                             i += block;
                         }
                         if (!accumulate_norm)
                         {
                             return;
                         }
                         for (int j = 1; j < NDIM; ++j)
                         {
                             if (idx(j) < norm_ibox.lower(j) || idx(j) > norm_ibox.upper(j))
                             {
                                 return;
                             }
                         }
                         const int k_begin = norm_ibox.lower(0) - idx(0);
                         const int k_end = norm_ibox.upper(0) - idx(0) + 1;
                         double sum = *norm_sum;
                         if (cv)
                         {
                             Index<NDIM> cv_idx = idx;
                             cv_idx(0) = norm_ibox.lower(0);
                             const double* const cv_row =
                                 cv->getPointer() + cv_depth_offset + cv->getBox().offset(cv_idx);
                             for (int k = k_begin; k < k_end; ++k)
                             {
                                 sum += y_row[k] * y_row[k] * cv_row[k - k_begin];
                             }
                         }
                         else
                         {
                             for (int k = k_begin; k < k_end; ++k)
                             {
                                 sum += y_row[k] * y_row[k];
                             }
                         }
                         *norm_sum = sum;
                     });
    }
}

// Return whether the data of component c of each vector in vecs have the same ghost box and depth as those of ref on
// every patch, where the data are of centering C. The control volume of ref, if any, must be of centering C with depth
// one or the depth of the data, and must allocate every direction that the data allocate.
template <DataCentering C>
bool
component_conforms(const SAMRAIVectorReal<NDIM, double>& ref,
                   const std::vector<SAMRAIVectorReal<NDIM, double>*>& vecs,
                   const int c)
{
    using Traits = CartesianCentering<C>;
    using Data = typename Traits::template Data<double>;
    Pointer<PatchHierarchy<NDIM>> hierarchy = ref.getPatchHierarchy();
    const int ref_idx = ref.getComponentDescriptorIndex(c);
    const int cv_idx = ref.getControlVolumeIndex(c);
    for (int ln = ref.getCoarsestLevelNumber(); ln <= ref.getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<Data> ref_data = patch->getPatchData(ref_idx);
            if (!ref_data)
            {
                return false;
            }
            if (cv_idx >= 0)
            {
                Pointer<Data> cv_data = patch->getPatchData(cv_idx);
                if (!cv_data || (cv_data->getDepth() != 1 && cv_data->getDepth() != ref_data->getDepth()))
                {
                    return false;
                }
                for (int axis = 0; axis < Traits::num_axes(); ++axis)
                {
                    if (Traits::template has_axis<double>(*ref_data, axis) &&
                        !Traits::template has_axis<double>(*cv_data, axis))
                    {
                        return false;
                    }
                }
            }
            for (const auto* v : vecs)
            {
                Pointer<Data> data = patch->getPatchData(v->getComponentDescriptorIndex(c));
                if (!data || data->getGhostBox() != ref_data->getGhostBox() || data->getDepth() != ref_data->getDepth())
                {
                    return false;
                }
                for (int axis = 0; axis < Traits::num_axes(); ++axis)
                {
                    if (Traits::template has_axis<double>(*data, axis) !=
                        Traits::template has_axis<double>(*ref_data, axis))
                    {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// Return whether every component of each vector in vecs has double-valued data of the same supported centering as the
// corresponding component of ref, with the same ghost box and depth on every patch.
bool
can_fuse(const SAMRAIVectorReal<NDIM, double>& ref, const std::vector<SAMRAIVectorReal<NDIM, double>*>& vecs)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = ref.getPatchHierarchy();
    const int ncomp = ref.getNumberOfComponents();
    for (const auto* v : vecs)
    {
        if (v->getNumberOfComponents() != ncomp || v->getPatchHierarchy() != hierarchy ||
            v->getCoarsestLevelNumber() != ref.getCoarsestLevelNumber() ||
            v->getFinestLevelNumber() != ref.getFinestLevelNumber())
        {
            return false;
        }
    }
    for (int c = 0; c < ncomp; ++c)
    {
        const std::optional<DataCentering> centering = find_component_centering(ref.getComponentDescriptorIndex(c));
        if (!centering)
        {
            return false;
        }
        for (const auto* v : vecs)
        {
            if (find_component_centering(v->getComponentDescriptorIndex(c)) != centering)
            {
                return false;
            }
        }
        const int cv_idx = ref.getControlVolumeIndex(c);
        if (cv_idx >= 0 && find_component_centering(cv_idx) != centering)
        {
            return false;
        }
        if (!dispatch_data_centering(*centering,
                                     [&]<DataCentering C>() { return component_conforms<C>(ref, vecs, c); }))
        {
            return false;
        }
    }
    return true;
}

// Set comp_sum[i] to the local dot product of component c of x with component c of y[begin + i] for each i < ny, where
// the data are of centering C. The terms are accumulated as in SAMRAI's hierarchy, patch, and array operations: over
// levels and patches, and on each patch over the coordinate directions that the data of x allocate.
template <DataCentering C>
void
mdot_component(const SAMRAIVectorReal<NDIM, double>& x,
               const std::vector<SAMRAIVectorReal<NDIM, double>*>& y,
               const int begin,
               const int ny,
               const int c,
               double* const comp_sum)
{
    using Traits = CartesianCentering<C>;
    using Data = typename Traits::template Data<double>;
    Pointer<PatchHierarchy<NDIM>> hierarchy = x.getPatchHierarchy();
    const int x_idx = x.getComponentDescriptorIndex(c);
    const int cv_idx = x.getControlVolumeIndex(c);
    std::vector<double> patch_sum(ny), array_sum(ny);
    std::vector<const ArrayData<NDIM, double>*> y_arrays(ny);
    std::fill(comp_sum, comp_sum + ny, 0.0);
    for (int ln = x.getCoarsestLevelNumber(); ln <= x.getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<Data> x_data = patch->getPatchData(x_idx);
            Pointer<Data> cv_data = cv_idx >= 0 ? patch->getPatchData(cv_idx) : Pointer<PatchData<NDIM>>();
            const Box<NDIM> box = cv_data ? x_data->getGhostBox() : patch->getBox();
            std::fill(patch_sum.begin(), patch_sum.end(), 0.0);
            for (int axis = 0; axis < Traits::num_axes(); ++axis)
            {
                if (!Traits::template has_axis<double>(*x_data, axis))
                {
                    continue;
                }
                for (int i = 0; i < ny; ++i)
                {
                    Pointer<Data> y_data = patch->getPatchData(y[begin + i]->getComponentDescriptorIndex(c));
                    y_arrays[i] = &Traits::template array_data<double>(*y_data, axis);
                }
                mdot_arrays(Traits::template array_data<double>(*x_data, axis),
                            y_arrays,
                            cv_data ? &Traits::template array_data<double>(*cv_data, axis) : nullptr,
                            Traits::index_box(box, axis),
                            array_sum.data());
                for (int i = 0; i < ny; ++i)
                {
                    patch_sum[i] += array_sum[i];
                }
            }
            for (int i = 0; i < ny; ++i)
            {
                comp_sum[i] += patch_sum[i];
            }
        }
    }
}

// Compute val[i] = x.dot(y[i]) over the local data for all i, bitwise identical to calling SAMRAIVectorReal::dot() with
// local_only = true for each i. The data must satisfy can_fuse().
void
fused_local_mdot(const SAMRAIVectorReal<NDIM, double>& x,
                 const std::vector<SAMRAIVectorReal<NDIM, double>*>& y,
                 PetscScalar* val)
{
    const int nv = static_cast<int>(y.size());
    std::vector<double> comp_sum;
    for (int begin = 0; begin < nv; begin += FUSED_GROUP_SIZE)
    {
        const int ny = std::min(FUSED_GROUP_SIZE, nv - begin);
        comp_sum.resize(ny);
        for (int i = 0; i < ny; ++i)
        {
            val[begin + i] = 0.0;
        }
        for (int c = 0; c < x.getNumberOfComponents(); ++c)
        {
            dispatch_data_centering(*find_component_centering(x.getComponentDescriptorIndex(c)),
                                    [&]<DataCentering C>() { mdot_component<C>(x, y, begin, ny, c, comp_sum.data()); });
            for (int i = 0; i < ny; ++i)
            {
                val[begin + i] += comp_sum[i];
            }
        }
    }
}

// Apply y.axpy(alpha_group[i], x[begin + i], y) for each i in order to component c of y over its ghost boxes, where the
// data are of centering C: on each patch, for each coordinate direction that the data of y allocate. If sum_of_squares
// is not null, also add to *sum_of_squares the local dot product of component c of the result with itself, accumulated
// as in the SAMRAI operations: over levels and patches, and on each patch over the coordinate directions that the data
// allocate, in the interior of the patch.
template <DataCentering C>
void
maxpy_component(SAMRAIVectorReal<NDIM, double>& y,
                const std::vector<double>& alpha_group,
                const std::vector<SAMRAIVectorReal<NDIM, double>*>& x,
                const int begin,
                const int c,
                double* const sum_of_squares)
{
    using Traits = CartesianCentering<C>;
    using Data = typename Traits::template Data<double>;
    Pointer<PatchHierarchy<NDIM>> hierarchy = y.getPatchHierarchy();
    const int nx = static_cast<int>(alpha_group.size());
    const int y_idx = y.getComponentDescriptorIndex(c);
    const int cv_idx = y.getControlVolumeIndex(c);
    std::vector<const ArrayData<NDIM, double>*> x_arrays(nx);
    for (int ln = y.getCoarsestLevelNumber(); ln <= y.getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<Data> y_data = patch->getPatchData(y_idx);
            Pointer<Data> cv_data =
                (sum_of_squares && cv_idx >= 0) ? patch->getPatchData(cv_idx) : Pointer<PatchData<NDIM>>();
            const Box<NDIM> box = y_data->getGhostBox();
            double patch_sum = 0.0;
            for (int axis = 0; axis < Traits::num_axes(); ++axis)
            {
                if (!Traits::template has_axis<double>(*y_data, axis))
                {
                    continue;
                }
                for (int i = 0; i < nx; ++i)
                {
                    Pointer<Data> x_data = patch->getPatchData(x[begin + i]->getComponentDescriptorIndex(c));
                    x_arrays[i] = &Traits::template array_data<double>(*x_data, axis);
                }
                const Box<NDIM> interior_box = Traits::index_box(patch->getBox(), axis);
                double array_sum = 0.0;
                maxpy_arrays(Traits::template array_data<double>(*y_data, axis),
                             alpha_group,
                             x_arrays,
                             Traits::index_box(box, axis),
                             cv_data ? &Traits::template array_data<double>(*cv_data, axis) : nullptr,
                             sum_of_squares ? &interior_box : nullptr,
                             sum_of_squares ? &array_sum : nullptr);
                patch_sum += array_sum;
            }
            if (sum_of_squares)
            {
                *sum_of_squares += patch_sum;
            }
        }
    }
}

// Compute y.axpy(alpha[i], x[i], y) for each i in order over ghost boxes, bitwise identical to the sequence of calls.
// The data must satisfy can_fuse(). If sum_of_squares is not null, also set *sum_of_squares in the pass that applies
// the last group of vectors to the local sum of squares of the result, accumulated as NormOps::L2Norm() does with plain
// summation.
void
fused_maxpy(SAMRAIVectorReal<NDIM, double>& y,
            const PetscScalar* alpha,
            const std::vector<SAMRAIVectorReal<NDIM, double>*>& x,
            double* const sum_of_squares)
{
    const int nv = static_cast<int>(x.size());
    if (sum_of_squares)
    {
        *sum_of_squares = 0.0;
    }
    for (int begin = 0; begin < nv; begin += FUSED_GROUP_SIZE)
    {
        const int nx = std::min(FUSED_GROUP_SIZE, nv - begin);
        double* const group_sum_of_squares = begin + nx == nv ? sum_of_squares : nullptr;
        const std::vector<double> alpha_group(alpha + begin, alpha + begin + nx);
        for (int c = 0; c < y.getNumberOfComponents(); ++c)
        {
            dispatch_data_centering(*find_component_centering(y.getComponentDescriptorIndex(c)),
                                    [&]<DataCentering C>()
                                    { maxpy_component<C>(y, alpha_group, x, begin, c, group_sum_of_squares); });
        }
    }
}

// The way multi-vector operations are carried out, selected by the option -ibtk_vec_fusion.
enum class FusionMode
{
    EXACT,
    NONE
};

// Return the fusion mode selected by -ibtk_vec_fusion {exact,none}, reading the option on first use. The default is
// exact.
FusionMode
get_fusion_mode()
{
    static const FusionMode mode = []
    {
        char value[64] = "exact";
        PetscBool set = PETSC_FALSE;
        int ierr = PetscOptionsGetString(nullptr, nullptr, "-ibtk_vec_fusion", value, sizeof(value), &set);
        IBTK_CHKERRQ(ierr);
        const std::string mode_name(value);
        if (mode_name == "exact")
        {
            return FusionMode::EXACT;
        }
        if (mode_name == "none")
        {
            return FusionMode::NONE;
        }
        TBOX_ERROR("PETScSAMRAIVectorReal::get_fusion_mode():\n"
                   << "  unknown value " << mode_name << " of the option -ibtk_vec_fusion; expected exact or none\n");
        return FusionMode::EXACT;
    }();
    return mode;
}

// Compute val[i] = x_vec.dot(y_vecs[i]) over the local data for all i, with the fused kernel unless fusion is disabled
// or the data do not allow it, in which case SAMRAIVectorReal::dot() is used for each i.
void
compute_local_mdot(const Pointer<SAMRAIVectorReal<NDIM, double>>& x_vec,
                   const std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>& y_vecs,
                   PetscScalar* val)
{
    static const bool local_only = true;
    const PetscInt nv = static_cast<PetscInt>(y_vecs.size());
    std::vector<SAMRAIVectorReal<NDIM, double>*> y_ptrs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        y_ptrs[i] = y_vecs[i].getPointer();
    }
    if (get_fusion_mode() == FusionMode::NONE || !can_fuse(*x_vec, y_ptrs))
    {
        for (PetscInt i = 0; i < nv; ++i)
        {
            val[i] = x_vec->dot(y_vecs[i], local_only);
        }
        return;
    }
    fused_local_mdot(*x_vec, y_ptrs, val);
}

// Compute y_vec.axpy(alpha[i], x_vecs[i], y_vec) for each i in order, with the fused kernel unless fusion is disabled,
// y_vec is one of the x_vecs, or the data do not allow it, in which case SAMRAIVectorReal::axpy() is used for each i.
// If the fused kernel was used, at least one vector was added, and NormOps uses plain summation, also return the local
// sum of squares of the result that NormOps::L2Norm() computes.
std::optional<double>
compute_maxpy(const Pointer<SAMRAIVectorReal<NDIM, double>>& y_vec,
              const PetscScalar* alpha,
              const std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>& x_vecs)
{
    static const bool interior_only = false;
    const PetscInt nv = static_cast<PetscInt>(x_vecs.size());
    std::vector<SAMRAIVectorReal<NDIM, double>*> x_ptrs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        x_ptrs[i] = x_vecs[i].getPointer();
    }
    const bool y_in_x = std::find(x_ptrs.begin(), x_ptrs.end(), y_vec.getPointer()) != x_ptrs.end();
    if (get_fusion_mode() == FusionMode::NONE || y_in_x || !can_fuse(*y_vec, x_ptrs))
    {
        for (PetscInt i = 0; i < nv; ++i)
        {
            y_vec->axpy(alpha[i], x_vecs[i], y_vec, interior_only);
        }
        return std::nullopt;
    }
    // Without vectors to add, no pass of the data computes the sum of squares.
    const bool compute_norm = nv > 0 && !NormOps::getSortedSummation();
    double sum_of_squares = 0.0;
    fused_maxpy(*y_vec, alpha, x_ptrs, compute_norm ? &sum_of_squares : nullptr);
    return compute_norm ? std::optional<double>(sum_of_squares) : std::nullopt;
}

#define PSVR_CAST1(v) (static_cast<PETScSAMRAIVectorReal*>(v->data))
#define PSVR_CAST2(v) (static_cast<PETScSAMRAIVectorReal*>(v->data)->d_samrai_vector)

#if !defined(NDEBUG)
#define PSVR_CHECK1(v)                                                                                                 \
    TBOX_ASSERT((v));                                                                                                  \
    TBOX_ASSERT(!PSVR_CAST1((v))->d_vector_checked_out_read);
#define PSVR_CHECK2(v1, v2)                                                                                            \
    PSVR_CHECK1((v1));                                                                                                 \
    PSVR_CHECK1((v2));
#define PSVR_CHECK3(v1, v2, v3)                                                                                        \
    PSVR_CHECK1((v1));                                                                                                 \
    PSVR_CHECK1((v2));                                                                                                 \
    PSVR_CHECK1((v3));
#define PSVR_CHECKN(v, N)                                                                                              \
    for (int i = 0; i < static_cast<int>(N); ++i)                                                                      \
    {                                                                                                                  \
        PSVR_CHECK1((v)[i]);                                                                                           \
    }
#else
#define PSVR_CHECK1(v)
#define PSVR_CHECK2(v1, v2)
#define PSVR_CHECK3(v1, v2, v3)
#define PSVR_CHECKN(v, N)
#endif
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

PETScSAMRAIVectorReal::PETScSAMRAIVectorReal(Pointer<SAMRAIVectorReal<NDIM, PetscScalar>> samrai_vector,
                                             bool vector_created_via_duplicate,
                                             MPI_Comm comm)
    : d_samrai_vector(samrai_vector), d_vector_created_via_duplicate(vector_created_via_duplicate)
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_vec_duplicate = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDuplicate()");
        t_vec_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDot()");
        t_vec_m_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMDot()");
        t_vec_norm = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecNorm()");
        t_vec_t_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecTDot()");
        t_vec_m_t_dot = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMTDot()");
        t_vec_scale = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecScale()");
        t_vec_copy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecCopy()");
        t_vec_set = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSet()");
        t_vec_swap = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSwap()");
        t_vec_axpy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPY()");
        t_vec_axpby = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPBY()");
        t_vec_maxpy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMAXPY()");
        t_vec_aypx = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAYPX()");
        t_vec_waxpy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecWAXPY()");
        t_vec_axpbypcz = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecAXPBYPCZ()");
        t_vec_pointwise_mult = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecPointwiseMult()");
        t_vec_pointwise_divide =
            TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecPointwiseDivide()");
        t_vec_get_size = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecGetSize()");
        t_vec_get_local_size = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecGetLocalSize()");
        t_vec_max = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMax()");
        t_vec_min = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMin()");
        t_vec_set_random = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecSetRandom()");
        t_vec_destroy = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDestroy()");
        t_vec_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDot_local()");
        t_vec_t_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecTDot_local()");
        t_vec_norm_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecNorm_local()");
        t_vec_m_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMDot_local()");
        t_vec_m_t_dot_local = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMTDot_local()");
        t_vec_max_pointwise_divide =
            TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecMaxPointwiseDivide()");
        t_vec_dot_norm2 = TimerManager::getManager()->getTimer("IBTK::PETScSAMRAIVectorReal::VecDotNorm2()"););

    int ierr;
    ierr = VecCreate(comm, &d_petsc_vector);
    IBTK_CHKERRQ(ierr);

    // Assign vector operations to PETSc vector object.
    static struct _VecOps DvOps;
    IBTK_DO_ONCE(DvOps.duplicate = PETScSAMRAIVectorReal::VecDuplicate_SAMRAI;
                 DvOps.duplicatevecs = VecDuplicateVecs_SAMRAI;
                 DvOps.destroyvecs = VecDestroyVecs_SAMRAI;
                 DvOps.dot = VecDot_SAMRAI;
                 DvOps.mdot = VecMDot_SAMRAI;
                 DvOps.norm = VecNorm_SAMRAI;
                 DvOps.tdot = VecTDot_SAMRAI;
                 DvOps.mtdot = VecMTDot_SAMRAI;
                 DvOps.scale = VecScale_SAMRAI;
                 DvOps.copy = VecCopy_SAMRAI;
                 DvOps.set = VecSet_SAMRAI;
                 DvOps.swap = VecSwap_SAMRAI;
                 DvOps.axpy = VecAXPY_SAMRAI;
                 DvOps.axpby = VecAXPBY_SAMRAI;
                 DvOps.maxpy = VecMAXPY_SAMRAI;
                 DvOps.aypx = VecAYPX_SAMRAI;
                 DvOps.waxpy = VecWAXPY_SAMRAI;
                 DvOps.axpbypcz = VecAXPBYPCZ_SAMRAI;
                 DvOps.pointwisemult = VecPointwiseMult_SAMRAI;
                 DvOps.pointwisedivide = VecPointwiseDivide_SAMRAI;
                 DvOps.getsize = VecGetSize_SAMRAI;
                 DvOps.getlocalsize = VecGetLocalSize_SAMRAI;
                 DvOps.max = VecMax_SAMRAI;
                 DvOps.min = VecMin_SAMRAI;
                 DvOps.setrandom = VecSetRandom_SAMRAI;
                 DvOps.destroy = PETScSAMRAIVectorReal::VecDestroy_SAMRAI;
                 DvOps.dot_local = VecDot_local_SAMRAI;
                 DvOps.tdot_local = VecTDot_local_SAMRAI;
                 DvOps.norm_local = VecNorm_local_SAMRAI;
                 DvOps.mdot_local = VecMDot_local_SAMRAI;
                 DvOps.mtdot_local = VecMTDot_local_SAMRAI;
                 DvOps.maxpointwisedivide = VecMaxPointwiseDivide_SAMRAI;
                 DvOps.dotnorm2 = VecDotNorm2_SAMRAI;);
    ierr = PetscMemcpy(d_petsc_vector->ops, &DvOps, sizeof(DvOps));
    IBTK_CHKERRQ(ierr);

    // Set PETSc vector data.
    d_petsc_vector->data = this;
    d_petsc_vector->petscnative = PETSC_FALSE;
    int size;
    MPI_Comm_size(comm, &size);
    d_petsc_vector->map->n = 1;    // NOTE: Here we are giving a bogus local  size.
    d_petsc_vector->map->N = size; // NOTE: Here we are giving a bogus global size.
    d_petsc_vector->map->bs = 1;   // NOTE: Here we are giving a bogus block  size.

    // Set the PETSc vector type name.
    ierr = PetscObjectChangeTypeName(reinterpret_cast<PetscObject>(d_petsc_vector), "Vec_SAMRAI");
    IBTK_CHKERRQ(ierr);

    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(d_petsc_vector));
    IBTK_CHKERRQ(ierr);
}

PETScSAMRAIVectorReal::~PETScSAMRAIVectorReal()
{
    if (!d_vector_created_via_duplicate)
    {
        d_petsc_vector->ops->destroy = nullptr;
        int ierr = VecDestroy(&d_petsc_vector);
        IBTK_CHKERRQ(ierr);
    }
}

/////////////////////////////// PRIVATE //////////////////////////////////////

PetscErrorCode
PETScSAMRAIVectorReal::VecDuplicate_SAMRAI(Vec v, Vec* newv)
{
    IBTK_TIMER_START(t_vec_duplicate);
    PetscFunctionBeginUser;
    PSVR_CHECK1(v);
    PetscErrorCode ierr;
    Pointer<SAMRAIVectorReal<NDIM, PetscScalar>> samrai_vec = PSVR_CAST2(v)->cloneVector(PSVR_CAST2(v)->getName());
    samrai_vec->allocateVectorData();
    static const bool vector_created_via_duplicate = true;
    MPI_Comm comm;
    ierr = PetscObjectGetComm(reinterpret_cast<PetscObject>(v), &comm);
    CHKERRQ(ierr);
    PETScSAMRAIVectorReal* new_psv = new PETScSAMRAIVectorReal(samrai_vec, vector_created_via_duplicate, comm);
    *newv = new_psv->d_petsc_vector;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(*newv));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_duplicate);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDestroy_SAMRAI(Vec v)
{
    IBTK_TIMER_START(t_vec_destroy);
    PetscFunctionBeginUser;
    PSVR_CHECK1(v);
    if (PSVR_CAST1(v)->d_vector_created_via_duplicate)
    {
        PSVR_CAST2(v)->resetLevels(0,
                                   std::min(PSVR_CAST2(v)->getFinestLevelNumber(),
                                            PSVR_CAST2(v)->getPatchHierarchy()->getFinestLevelNumber()));
        deallocate_vector_data(*PSVR_CAST2(v));
        free_vector_components(*PSVR_CAST2(v));
        PSVR_CAST2(v).setNull();
        destroyPETScVector(PSVR_CAST1(v)->d_petsc_vector);
    }
    IBTK_TIMER_STOP(t_vec_destroy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDuplicateVecs_SAMRAI(Vec v, PetscInt m, Vec* V[])
{
    PetscFunctionBeginUser;
    PSVR_CHECK1(v);
    PetscErrorCode ierr;
    ierr = PetscMalloc1(m, V);
    CHKERRQ(ierr);
    for (PetscInt i = 0; i < m; ++i)
    {
        ierr = VecDuplicate(v, *V + i);
        CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDestroyVecs_SAMRAI(PetscInt m, Vec vv[])
{
    PetscFunctionBeginUser;
    PSVR_CHECKN(vv, m);
    PetscErrorCode ierr;
    for (PetscInt i = 0; i < m; ++i)
    {
        ierr = VecDestroy(&vv[i]);
        CHKERRQ(ierr);
    }
    ierr = PetscFree(vv);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDot_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_dot);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y));
    IBTK_TIMER_STOP(t_vec_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMDot_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_dot);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> y_vecs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        y_vecs[i] = PSVR_CAST2(y[i]);
    }
    compute_local_mdot(PSVR_CAST2(x), y_vecs, val);
    IBTK_MPI::sumReduction(val, nv);
    IBTK_TIMER_STOP(t_vec_m_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecNorm_SAMRAI(Vec x, NormType type, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_norm);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    if (type == NORM_1)
    {
        *val = NormOps::L1Norm(PSVR_CAST2(x));
    }
    else if (type == NORM_2)
    {
        const bool local_only = false;
        PetscErrorCode ierr = L2Norm_SAMRAI(x, local_only, val);
        CHKERRQ(ierr);
    }
    else if (type == NORM_INFINITY)
    {
        *val = NormOps::maxNorm(PSVR_CAST2(x));
    }
    else if (type == NORM_1_AND_2)
    {
        static const bool local_only = true;
        val[0] = NormOps::L1Norm(PSVR_CAST2(x), local_only);
        val[1] = NormOps::L2Norm(PSVR_CAST2(x), local_only);
        val[1] = val[1] * val[1];
        IBTK_MPI::sumReduction(val, 2);
        val[1] = std::sqrt(val[1]);
    }
    else
    {
        TBOX_ERROR("PETScSAMRAIVectorReal::norm()\n"
                   << "  vector norm type " << static_cast<int>(type) << " unsupported" << std::endl);
    }
    IBTK_TIMER_STOP(t_vec_norm);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::L2Norm_SAMRAI(Vec x, const bool local_only, PetscScalar* val)
{
    PetscFunctionBeginUser;
    PetscObjectState state;
    int ierr = PetscObjectStateGet(reinterpret_cast<PetscObject>(x), &state);
    CHKERRQ(ierr);
    const PETScSAMRAIVectorReal* const x_wrapper = PSVR_CAST1(x);
    if (x_wrapper->d_norm_cache_state == state && !NormOps::getSortedSummation())
    {
        *val = NormOps::L2NormFromSumOfSquares(x_wrapper->d_norm_cache_sum_of_squares, local_only);
#if !defined(NDEBUG)
        // The comparison is local so that every process performs the same reductions.
        static const bool check_local_only = true;
        TBOX_ASSERT(NormOps::L2NormFromSumOfSquares(x_wrapper->d_norm_cache_sum_of_squares, check_local_only) ==
                    NormOps::L2Norm(PSVR_CAST2(x), check_local_only));
#endif
    }
    else
    {
        *val = NormOps::L2Norm(PSVR_CAST2(x), local_only);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecTDot_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_t_dot);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y));
    IBTK_TIMER_STOP(t_vec_t_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMTDot_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_t_dot);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> y_vecs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        y_vecs[i] = PSVR_CAST2(y[i]);
    }
    compute_local_mdot(PSVR_CAST2(x), y_vecs, val);
    IBTK_MPI::sumReduction(val, nv);
    IBTK_TIMER_STOP(t_vec_m_t_dot);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecScale_SAMRAI(Vec x, PetscScalar alpha)
{
    IBTK_TIMER_START(t_vec_scale);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    static const bool interior_only = false;
    PSVR_CAST2(x)->scale(alpha, PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_scale);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecCopy_SAMRAI(Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_copy);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    PSVR_CAST2(y)->copyVector(PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_copy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecSet_SAMRAI(Vec x, PetscScalar alpha)
{
    IBTK_TIMER_START(t_vec_set);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    static const bool interior_only = false;
    PSVR_CAST2(x)->setToScalar(alpha, interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_set);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecSwap_SAMRAI(Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_swap);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    PSVR_CAST2(x)->swapVectors(PSVR_CAST2(y));
    int ierr;
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_swap);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAXPY_SAMRAI(Vec y, PetscScalar alpha, Vec x)
{
    IBTK_TIMER_START(t_vec_axpy);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_axpy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAXPBY_SAMRAI(Vec y, PetscScalar alpha, PetscScalar beta, Vec x)
{
    IBTK_TIMER_START(t_vec_axpby);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    PSVR_CAST2(y)->linearSum(alpha, PSVR_CAST2(x), beta, PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_axpby);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMAXPY_SAMRAI(Vec y, PetscInt nv, const PetscScalar* alpha, Vec* x)
{
    IBTK_TIMER_START(t_vec_maxpy);
    PetscFunctionBeginUser;
    PSVR_CHECK1(y);
    PSVR_CHECKN(x, nv);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> x_vecs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        x_vecs[i] = PSVR_CAST2(x[i]);
    }
    const std::optional<double> sum_of_squares = compute_maxpy(PSVR_CAST2(y), alpha, x_vecs);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    PETScSAMRAIVectorReal* const y_wrapper = PSVR_CAST1(y);
    y_wrapper->d_norm_cache_state = -1;
    if (sum_of_squares)
    {
        // PETSc's VecMAXPY() increments the object state once more after this operation returns, so the cached sum
        // of squares corresponds to the next state.
        ierr = PetscObjectStateGet(reinterpret_cast<PetscObject>(y), &y_wrapper->d_norm_cache_state);
        CHKERRQ(ierr);
        ++y_wrapper->d_norm_cache_state;
        y_wrapper->d_norm_cache_sum_of_squares = *sum_of_squares;
    }
    IBTK_TIMER_STOP(t_vec_maxpy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAYPX_SAMRAI(Vec y, const PetscScalar alpha, Vec x)
{
    IBTK_TIMER_START(t_vec_aypx);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    static const bool interior_only = false;
    PSVR_CAST2(y)->axpy(alpha, PSVR_CAST2(y), PSVR_CAST2(x), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_aypx);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecWAXPY_SAMRAI(Vec w, PetscScalar alpha, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_waxpy);
    PetscFunctionBeginUser;
    PSVR_CHECK3(w, x, y);
    static const bool interior_only = false;
    PSVR_CAST2(w)->axpy(alpha, PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_waxpy);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecAXPBYPCZ_SAMRAI(Vec z, PetscScalar alpha, PetscScalar beta, PetscScalar gamma, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_axpbypcz);
    PetscFunctionBeginUser;
    PSVR_CHECK3(x, y, z);
    static const bool interior_only = false;
    PSVR_CAST2(z)->linearSum(alpha, PSVR_CAST2(x), gamma, PSVR_CAST2(z), interior_only);
    PSVR_CAST2(z)->axpy(beta, PSVR_CAST2(y), PSVR_CAST2(z), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(z));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_axpbypcz);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecPointwiseMult_SAMRAI(Vec w, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_pointwise_mult);
    PetscFunctionBeginUser;
    PSVR_CHECK3(w, x, y);
    static const bool interior_only = false;
    PSVR_CAST2(w)->multiply(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_pointwise_mult);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecPointwiseDivide_SAMRAI(Vec w, Vec x, Vec y)
{
    IBTK_TIMER_START(t_vec_pointwise_divide);
    PetscFunctionBeginUser;
    PSVR_CHECK3(w, x, y);
    static const bool interior_only = false;
    PSVR_CAST2(w)->divide(PSVR_CAST2(x), PSVR_CAST2(y), interior_only);
    int ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(w));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_pointwise_divide);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecGetSize_SAMRAI(Vec v, PetscInt* size)
{
    IBTK_TIMER_START(t_vec_get_size);
    PetscFunctionBeginUser;
    PSVR_CHECK1(v);
    *size = v->map->N;
    IBTK_TIMER_STOP(t_vec_get_size);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecGetLocalSize_SAMRAI(Vec v, PetscInt* size)
{
    IBTK_TIMER_START(t_vec_get_local_size);
    PetscFunctionBeginUser;
    PSVR_CHECK1(v);
    *size = v->map->n;
    IBTK_TIMER_STOP(t_vec_get_local_size);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMax_SAMRAI(Vec x, PetscInt* p, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_max);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    *p = -1;
    *val = PSVR_CAST2(x)->max();
    IBTK_TIMER_STOP(t_vec_max);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMin_SAMRAI(Vec x, PetscInt* p, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_min);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    *p = -1;
    *val = PSVR_CAST2(x)->min();
    IBTK_TIMER_STOP(t_vec_min);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecSetRandom_SAMRAI(Vec x, PetscRandom rctx)
{
    IBTK_TIMER_START(t_vec_set_random);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    PetscScalar lo, hi;
    int ierr;
    ierr = PetscRandomGetInterval(rctx, &lo, &hi);
    CHKERRQ(ierr);
    PSVR_CAST2(x)->setRandomValues(hi - lo, lo);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(x));
    CHKERRQ(ierr);
    IBTK_TIMER_STOP(t_vec_set_random);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDot_local_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_dot_local);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    static const bool local_only = true;
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y), local_only);
    IBTK_TIMER_STOP(t_vec_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecTDot_local_SAMRAI(Vec x, Vec y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_t_dot_local);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    static const bool local_only = true;
    *val = PSVR_CAST2(x)->dot(PSVR_CAST2(y), local_only);
    IBTK_TIMER_STOP(t_vec_t_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecNorm_local_SAMRAI(Vec x, NormType type, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_norm_local);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    static const bool local_only = true;
    if (type == NORM_1)
    {
        *val = NormOps::L1Norm(PSVR_CAST2(x), local_only);
    }
    else if (type == NORM_2)
    {
        PetscErrorCode ierr = L2Norm_SAMRAI(x, local_only, val);
        CHKERRQ(ierr);
    }
    else if (type == NORM_INFINITY)
    {
        *val = NormOps::maxNorm(PSVR_CAST2(x), local_only);
    }
    else if (type == NORM_1_AND_2)
    {
        val[0] = NormOps::L1Norm(PSVR_CAST2(x), local_only);
        val[1] = NormOps::L2Norm(PSVR_CAST2(x), local_only);
    }
    else
    {
        TBOX_ERROR("PETScSAMRAIVectorReal::norm()\n"
                   << "  vector norm type " << static_cast<int>(type) << " unsupported" << std::endl);
    }
    IBTK_TIMER_STOP(t_vec_norm_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMDot_local_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_dot_local);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> y_vecs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        y_vecs[i] = PSVR_CAST2(y[i]);
    }
    compute_local_mdot(PSVR_CAST2(x), y_vecs, val);
    IBTK_TIMER_STOP(t_vec_m_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMTDot_local_SAMRAI(Vec x, PetscInt nv, const Vec* y, PetscScalar* val)
{
    IBTK_TIMER_START(t_vec_m_t_dot_local);
    PetscFunctionBeginUser;
    PSVR_CHECK1(x);
    PSVR_CHECKN(y, nv);
    std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>> y_vecs(nv);
    for (PetscInt i = 0; i < nv; ++i)
    {
        y_vecs[i] = PSVR_CAST2(y[i]);
    }
    compute_local_mdot(PSVR_CAST2(x), y_vecs, val);
    IBTK_TIMER_STOP(t_vec_m_t_dot_local);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecMaxPointwiseDivide_SAMRAI(Vec x, Vec y, PetscScalar* max)
{
    IBTK_TIMER_START(t_vec_max_pointwise_divide);
    PetscFunctionBeginUser;
    PSVR_CHECK2(x, y);
    *max = PSVR_CAST2(x)->maxPointwiseDivide(PSVR_CAST2(y));
    IBTK_TIMER_STOP(t_vec_max_pointwise_divide);
    PetscFunctionReturn(0);
}

PetscErrorCode
PETScSAMRAIVectorReal::VecDotNorm2_SAMRAI(Vec s, Vec t, PetscScalar* dp, PetscScalar* nm)
{
    IBTK_TIMER_START(t_vec_dot_norm2);
    PetscFunctionBeginUser;
    PSVR_CHECK2(s, t);
    *dp = PSVR_CAST2(s)->dot(PSVR_CAST2(t));
    *nm = PSVR_CAST2(t)->dot(PSVR_CAST2(t));
    IBTK_TIMER_STOP(t_vec_dot_norm2);
    PetscFunctionReturn(0);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
