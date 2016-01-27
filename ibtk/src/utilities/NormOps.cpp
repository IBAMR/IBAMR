// Filename: NormOps.cpp
// Created on 08 Dec 2008 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>

#include "CellData.h"
#include "CellVariable.h"
#include "IntVector.h"
#include "ibtk/NormOps.h"
#include "Patch.h"
#include "PatchCellDataNormOpsReal.h"
#include "PatchData.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataNormOpsReal.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideVariable.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

double
NormOps::L1Norm(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const bool local_only)
{
    const double L1_norm_local = L1Norm_local(samrai_vector);
    if (local_only) return L1_norm_local;

    const int nprocs = SAMRAI_MPI::getNodes();
    std::vector<double> L1_norm_proc(nprocs, 0.0);
    SAMRAI_MPI::allGather(L1_norm_local, &L1_norm_proc[0]);
    const double ret_val = accurate_sum(L1_norm_proc);
    return ret_val;
} // L1Norm

double
NormOps::L2Norm(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const bool local_only)
{
    const double L2_norm_local = L2Norm_local(samrai_vector);
    if (local_only) return L2_norm_local;

    const int nprocs = SAMRAI_MPI::getNodes();
    std::vector<double> L2_norm_proc(nprocs, 0.0);
    SAMRAI_MPI::allGather(L2_norm_local, &L2_norm_proc[0]);
    const double ret_val = std::sqrt(accurate_sum_of_squares(L2_norm_proc));
    return ret_val;
} // L2Norm

double
NormOps::maxNorm(const SAMRAIVectorReal<NDIM, double>* const samrai_vector, const bool local_only)
{
    return samrai_vector->maxNorm(local_only);
} // maxNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

double
NormOps::L1Norm_local(const SAMRAIVectorReal<NDIM, double>* const samrai_vector)
{
    std::vector<double> L1_norm_local_patch;
    Pointer<PatchHierarchy<NDIM> > hierarchy = samrai_vector->getPatchHierarchy();
    const int coarsest_ln = samrai_vector->getCoarsestLevelNumber();
    const int finest_ln = samrai_vector->getFinestLevelNumber();
    const int ncomp = samrai_vector->getNumberOfComponents();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const Pointer<Variable<NDIM> >& comp_var = samrai_vector->getComponentVariable(comp);
        const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
        const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
        const bool has_cvol = cvol_idx >= 0;

        Pointer<CellVariable<NDIM, double> > comp_cc_var = comp_var;
        if (comp_cc_var)
        {
            PatchCellDataNormOpsReal<NDIM, double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    Pointer<CellData<NDIM, double> > comp_data = patch->getPatchData(comp_idx);
                    Pointer<CellData<NDIM, double> > cvol_data =
                        (has_cvol ? patch->getPatchData(cvol_idx) : Pointer<PatchData<NDIM> >(NULL));
                    L1_norm_local_patch.push_back(patch_ops.L1Norm(comp_data, patch_box, cvol_data));
                }
            }
        }

        Pointer<SideVariable<NDIM, double> > comp_sc_var = comp_var;
        if (comp_sc_var)
        {
            PatchSideDataNormOpsReal<NDIM, double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    Pointer<SideData<NDIM, double> > comp_data = patch->getPatchData(comp_idx);
                    Pointer<SideData<NDIM, double> > cvol_data =
                        (has_cvol ? patch->getPatchData(cvol_idx) : Pointer<PatchData<NDIM> >(NULL));
                    L1_norm_local_patch.push_back(patch_ops.L1Norm(comp_data, patch_box, cvol_data));
                }
            }
        }
    }
    return accurate_sum(L1_norm_local_patch);
} // L1Norm_local

double
NormOps::L2Norm_local(const SAMRAIVectorReal<NDIM, double>* const samrai_vector)
{
    std::vector<double> L2_norm_local_patch;
    Pointer<PatchHierarchy<NDIM> > hierarchy = samrai_vector->getPatchHierarchy();
    const int coarsest_ln = samrai_vector->getCoarsestLevelNumber();
    const int finest_ln = samrai_vector->getFinestLevelNumber();
    const int ncomp = samrai_vector->getNumberOfComponents();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        const Pointer<Variable<NDIM> >& comp_var = samrai_vector->getComponentVariable(comp);
        const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
        const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
        const bool has_cvol = cvol_idx >= 0;

        Pointer<CellVariable<NDIM, double> > comp_cc_var = comp_var;
        if (comp_cc_var)
        {
            PatchCellDataNormOpsReal<NDIM, double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    Pointer<CellData<NDIM, double> > comp_data = patch->getPatchData(comp_idx);
                    Pointer<CellData<NDIM, double> > cvol_data =
                        (has_cvol ? patch->getPatchData(cvol_idx) : Pointer<PatchData<NDIM> >(NULL));
                    L2_norm_local_patch.push_back(patch_ops.L2Norm(comp_data, patch_box, cvol_data));
                }
            }
        }

        Pointer<SideVariable<NDIM, double> > comp_sc_var = comp_var;
        if (comp_sc_var)
        {
            PatchSideDataNormOpsReal<NDIM, double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    Pointer<SideData<NDIM, double> > comp_data = patch->getPatchData(comp_idx);
                    Pointer<SideData<NDIM, double> > cvol_data =
                        (has_cvol ? patch->getPatchData(cvol_idx) : Pointer<PatchData<NDIM> >(NULL));
                    L2_norm_local_patch.push_back(patch_ops.L2Norm(comp_data, patch_box, cvol_data));
                }
            }
        }
    }
    return std::sqrt(accurate_sum_of_squares(L2_norm_local_patch));
} // L2Norm_local

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
