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

#include <tbox/Pointer.h>

#include <petscsys.h>

#include <Box.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>
#include <Variable.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// WARNING: This function will sort the input vector in ascending order.
inline double
accurate_sum(std::vector<double>& vec)
{
    if (vec.size() == 1) return vec[0];
    std::sort(vec.begin(), vec.end(), std::less<double>());
    return std::accumulate(vec.begin(), vec.end(), 0.0);
} // accurate_sum

// WARNING: This function will sort the input vector in ascending order.
inline double
accurate_sum_of_squares(std::vector<double>& vec)
{
    if (vec.size() == 1) return vec[0] * vec[0];
    std::sort(vec.begin(), vec.end(), std::less<double>());
    return std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
} // accurate_sum_of_squares

// Call function(patch_ops, data, patch_box, cvol) for every patch of one component of the vector, which has the data
// centering C, where patch_ops is the norm operations of C, data and cvol are the patch data of the component and of
// its control volume (null if the component has none), and patch_box is the box of the patch. The calls are ordered by
// level, from the coarsest to the finest level of the vector, and for each level by the iteration order of its patches.
template <DataCentering C, typename Function>
void
for_each_patch(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const int comp, Function&& function)
{
    using Data = typename CartesianCentering<C>::template Data<double>;
    typename CartesianCentering<C>::template PatchNormOps<double> patch_ops;
    Pointer<PatchHierarchy<NDIM>> hierarchy = samrai_vector->getPatchHierarchy();
    const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
    const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
    const bool has_cvol = cvol_idx >= 0;
    for (int ln = samrai_vector->getCoarsestLevelNumber(); ln <= samrai_vector->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<Data> data = patch->getPatchData(comp_idx);
            Pointer<Data> cvol = has_cvol ? patch->getPatchData(cvol_idx) : Pointer<PatchData<NDIM>>(nullptr);
            function(patch_ops, data, patch->getBox(), cvol);
        }
    }
}

// Call function(patch_ops, data, patch_box, cvol), as described for for_each_patch(), for every patch of every
// component of the vector, in the order of the components. It is a fatal error if a component is not cell-, node-,
// side-, face-, or edge-centered double-precision data.
template <typename Function>
void
for_each_component_patch(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, Function&& function)
{
    for (int comp = 0; comp < samrai_vector->getNumberOfComponents(); ++comp)
    {
        const DataCentering centering =
            get_data_centering<double>(*samrai_vector->getComponentVariable(comp)->getPatchDataFactory());
        dispatch_data_centering(centering,
                                [&]<DataCentering C>() { for_each_patch<C>(samrai_vector, comp, function); });
    }
}
} // namespace

std::optional<bool> NormOps::s_sorted_summation;

/////////////////////////////// PUBLIC ///////////////////////////////////////

double
NormOps::L1Norm(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const bool local_only)
{
    if (!getSortedSummation())
    {
        const double L1_norm_local = L1Norm_local_unsorted(samrai_vector);
        return local_only ? L1_norm_local : IBTK_MPI::sumReduction(L1_norm_local);
    }
    const double L1_norm_local = L1Norm_local(samrai_vector);
    if (local_only) return L1_norm_local;

    const int nprocs = IBTK_MPI::getNodes();
    std::vector<double> L1_norm_proc(nprocs, 0.0);
    IBTK_MPI::allGather(L1_norm_local, &L1_norm_proc[0]);
    const double ret_val = accurate_sum(L1_norm_proc);
    return ret_val;
} // L1Norm

double
NormOps::L2Norm(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const bool local_only)
{
    if (!getSortedSummation())
    {
        return L2NormFromSumOfSquares(L2NormSquared_local_unsorted(samrai_vector), local_only);
    }
    const double L2_norm_local = L2Norm_local(samrai_vector);
    if (local_only) return L2_norm_local;

    const int nprocs = IBTK_MPI::getNodes();
    std::vector<double> L2_norm_proc(nprocs, 0.0);
    IBTK_MPI::allGather(L2_norm_local, &L2_norm_proc[0]);
    const double ret_val = std::sqrt(accurate_sum_of_squares(L2_norm_proc));
    return ret_val;
} // L2Norm

double
NormOps::maxNorm(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const bool local_only)
{
    return samrai_vector->maxNorm(local_only);
} // maxNorm

void
NormOps::setSortedSummation(const bool sorted_summation)
{
    s_sorted_summation = sorted_summation;
} // setSortedSummation

bool
NormOps::getSortedSummation()
{
    if (!s_sorted_summation.has_value())
    {
        // The option cannot be read before PETSc is initialized; use the default and read the option on a later call.
        PetscBool petsc_initialized = PETSC_FALSE;
        int ierr = PetscInitialized(&petsc_initialized);
        IBTK_CHKERRQ(ierr);
        if (!petsc_initialized)
        {
            return false;
        }
        PetscBool sorted_summation = PETSC_FALSE;
        ierr = PetscOptionsGetBool(nullptr, nullptr, "-ibtk_sorted_norm_summation", &sorted_summation, nullptr);
        IBTK_CHKERRQ(ierr);
        s_sorted_summation = sorted_summation == PETSC_TRUE;
    }
    return *s_sorted_summation;
} // getSortedSummation

double
NormOps::L2NormFromSumOfSquares(const double local_sum_of_squares, const bool local_only)
{
    return std::sqrt(local_only ? local_sum_of_squares : IBTK_MPI::sumReduction(local_sum_of_squares));
} // L2NormFromSumOfSquares

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

double
NormOps::L1Norm_local(const SAMRAIVectorReal<NDIM, double>* const samrai_vector)
{
    std::vector<double> L1_norm_local_patch;
    for_each_component_patch(samrai_vector,
                             [&](const auto& patch_ops, const auto& data, const Box<NDIM>& patch_box, const auto& cvol)
                             { L1_norm_local_patch.push_back(patch_ops.L1Norm(data, patch_box, cvol)); });
    return accurate_sum(L1_norm_local_patch);
} // L1Norm_local

double
NormOps::L2Norm_local(const SAMRAIVectorReal<NDIM, double>* const samrai_vector)
{
    std::vector<double> L2_norm_local_patch;
    for_each_component_patch(samrai_vector,
                             [&](const auto& patch_ops, const auto& data, const Box<NDIM>& patch_box, const auto& cvol)
                             { L2_norm_local_patch.push_back(patch_ops.L2Norm(data, patch_box, cvol)); });
    return std::sqrt(accurate_sum_of_squares(L2_norm_local_patch));
} // L2Norm_local

double
NormOps::L1Norm_local_unsorted(const SAMRAIVectorReal<NDIM, double>* const samrai_vector)
{
    double L1_norm = 0.0;
    for_each_component_patch(samrai_vector,
                             [&](const auto& patch_ops, const auto& data, const Box<NDIM>& patch_box, const auto& cvol)
                             { L1_norm += patch_ops.L1Norm(data, patch_box, cvol); });
    return L1_norm;
} // L1Norm_local_unsorted

double
NormOps::L2NormSquared_local_unsorted(const SAMRAIVectorReal<NDIM, double>* const samrai_vector)
{
    double sum_of_squares = 0.0;
    for_each_component_patch(samrai_vector,
                             [&](const auto& patch_ops, const auto& data, const Box<NDIM>& patch_box, const auto& cvol)
                             { sum_of_squares += patch_ops.dot(data, data, patch_box, cvol); });
    return sum_of_squares;
} // L2NormSquared_local_unsorted

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
