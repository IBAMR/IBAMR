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

#include "ibtk/IBTK_MPI.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/interpolation_utilities.h"

#include <ibtk/app_namespaces.h>

namespace IBTK
{
int
determine_depth(Pointer<hier::Variable<NDIM> > var, int depth)
{
    Pointer<SideVariable<NDIM, double> > sc_var = var;
    Pointer<FaceVariable<NDIM, double> > fc_var = var;
    if (sc_var || fc_var)
        return depth * NDIM;
    else
        return depth;
}

std::vector<double>
flatten_eig_vec(const std::vector<IBTK::VectorNd>& eig_vec)
{
    std::vector<double> X_vec(eig_vec.size() * NDIM);
    int i = 0;
    for (const auto& X : eig_vec)
    {
        for (int d = 0; d < NDIM; ++d) X_vec[i++] = X[d];
    }
    return X_vec;
}

// Checks if a vector is the same on all ranks. If the min and max reduction are equal, the elements are the same.
template <typename T>
bool
check_consistent_across_ranks(std::vector<T> data)
{
    std::vector<T> copied = data;
    // Do a max reduction
    IBTK_MPI::maxReduction(data.data(), data.size());
    IBTK_MPI::minReduction(copied.data(), copied.size());
    return std::equal(
        data.begin(), data.begin() + data.size(), copied.data(), [](T v1, T v2) -> bool { return v1 == v2; });
}

std::vector<double>
interpolate(const VectorNd& X,
            const int data_idx,
            Pointer<hier::Variable<NDIM> > Q_var,
            int Q_depth,
            Pointer<PatchHierarchy<NDIM> > hierarchy,
            std::string interp_fcn)
{
    std::vector<VectorNd> X_vec = { X };
    return interpolate(X_vec, data_idx, Q_var, Q_depth, hierarchy, std::move(interp_fcn));
}

std::vector<double>
interpolate(const std::vector<VectorNd>& X,
            const int data_idx,
            Pointer<hier::Variable<NDIM> > Q_var,
            int Q_depth,
            Pointer<PatchHierarchy<NDIM> > hierarchy,
            std::string interp_fcn)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const int coarsest_ln = 0;
    // Determine the actual depth of the Q_var, taking into account data layout.
    int actual_depth = determine_depth(Q_var, Q_depth);
    // We store interpolated data separately on each level, so we can correctly reduce it later.
    std::vector<std::vector<double> > Q_data_ln_vec(finest_ln + 1, std::vector<double>(X.size() * actual_depth, 0.0));
    // Flatten the vector of VectorNd to a single list of doubles.
    std::vector<double> X_data = flatten_eig_vec(X);
#ifndef NDEBUG
    if (!check_consistent_across_ranks(X_data))
    {
        std::ostringstream msg;
        msg << "The input data is not consistent across ranks. The function IBTK::interpolate currently requires that "
               "the points are synchronized across all MPI ranks.";
        TBOX_ERROR(msg.str());
    }
#endif
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        std::vector<double>& Q_data = Q_data_ln_vec[ln];
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();
            // Note that LEInteractor currently only interpolates cell, side, node, and edge data.
            Pointer<CellData<NDIM, double> > cc_data = patch->getPatchData(data_idx);
            Pointer<SideData<NDIM, double> > sc_data = patch->getPatchData(data_idx);
            Pointer<NodeData<NDIM, double> > nc_data = patch->getPatchData(data_idx);
            Pointer<EdgeData<NDIM, double> > ec_data = patch->getPatchData(data_idx);
            if (cc_data)
                LEInteractor::interpolate(Q_data.data(),
                                          Q_data.size(),
                                          Q_depth,
                                          X_data.data(),
                                          X_data.size(),
                                          NDIM,
                                          cc_data,
                                          patch,
                                          box,
                                          std::move(interp_fcn));
            else if (sc_data)
                LEInteractor::interpolate(Q_data.data(),
                                          Q_data.size(),
                                          NDIM,
                                          X_data.data(),
                                          X_data.size(),
                                          NDIM,
                                          sc_data,
                                          patch,
                                          box,
                                          std::move(interp_fcn));
            else if (nc_data)
                LEInteractor::interpolate(Q_data.data(),
                                          Q_data.size(),
                                          Q_depth,
                                          X_data.data(),
                                          X_data.size(),
                                          NDIM,
                                          nc_data,
                                          patch,
                                          box,
                                          std::move(interp_fcn));
            else if (ec_data)
                LEInteractor::interpolate(Q_data.data(),
                                          Q_data.size(),
                                          Q_depth,
                                          X_data.data(),
                                          X_data.size(),
                                          NDIM,
                                          ec_data,
                                          patch,
                                          box,
                                          std::move(interp_fcn));
        }
    }
    // Different processors may have interpolated to values on different levels. So we need to do a reduction to make
    // sure the correct value is located on all processors, while prefering the value on the finest level.
    std::vector<double> final_Q_vals(X.size() * actual_depth, 0.0);
    std::vector<int> Q_val_found(X.size() * actual_depth, 0);
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        std::vector<double>& Q_vals = Q_data_ln_vec[ln];
        for (size_t i = 0; i < final_Q_vals.size(); ++i)
        {
            if (Q_val_found[i] == 0 && Q_vals[i] != 0.0)
            {
                final_Q_vals[i] = Q_vals[i];
                Q_val_found[i] = 1;
            }
        }
        IBTK_MPI::maxReduction(Q_val_found.data(), Q_val_found.size());
    }
    IBTK_MPI::sumReduction(final_Q_vals.data(), final_Q_vals.size());
    return final_Q_vals;
}

} // namespace IBTK
