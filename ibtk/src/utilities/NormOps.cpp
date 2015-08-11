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

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "ibtk/NormOps.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/math/PatchSideDataNormOpsReal.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "SAMRAI/tbox/SAMRAI_MPI.h"

namespace SAMRAI
{
namespace hier
{

class Variable;

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
inline double accurate_sum(std::vector<double>& vec)
{
    if (vec.size() == 1) return vec[0];
    std::sort(vec.begin(), vec.end(), std::less<double>());
    return std::accumulate(vec.begin(), vec.end(), 0.0);
}

// WARNING: This function will sort the input vector in ascending order.
inline double accurate_sum_of_squares(std::vector<double>& vec)
{
    if (vec.size() == 1) return vec[0] * vec[0];
    std::sort(vec.begin(), vec.end(), std::less<double>());
    return std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0);
}
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

double NormOps::L1Norm(const boost::shared_ptr<SAMRAIVectorReal<double>>& samrai_vector, const bool local_only)
{
    double L1_norm_local = L1Norm_local(samrai_vector.get());
    if (local_only) return L1_norm_local;

    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    const int nprocs = comm.getSize();
    std::vector<double> L1_norm_proc(nprocs, 0.0);
    comm.Allgather(&L1_norm_local, 1, MPI_DOUBLE, &L1_norm_proc[0], nprocs, MPI_DOUBLE);
    const double ret_val = accurate_sum(L1_norm_proc);
    return ret_val;
}

double NormOps::L2Norm(const boost::shared_ptr<SAMRAIVectorReal<double>>& samrai_vector, const bool local_only)
{
    double L2_norm_local = L2Norm_local(samrai_vector.get());
    if (local_only) return L2_norm_local;

    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    const int nprocs = comm.getSize();
    std::vector<double> L2_norm_proc(nprocs, 0.0);
    comm.Allgather(&L2_norm_local, 1, MPI_DOUBLE, &L2_norm_proc[0], nprocs, MPI_DOUBLE);
    const double ret_val = std::sqrt(accurate_sum_of_squares(L2_norm_proc));
    return ret_val;
}

double NormOps::maxNorm(const boost::shared_ptr<SAMRAIVectorReal<double>>& samrai_vector, const bool local_only)
{
    return samrai_vector->maxNorm(local_only);
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

double NormOps::L1Norm_local(const SAMRAIVectorReal<double>* const samrai_vector)
{
    std::vector<double> L1_norm_local_patch;
    auto hierarchy = samrai_vector->getPatchHierarchy();
    const int coarsest_ln = samrai_vector->getCoarsestLevelNumber();
    const int finest_ln = samrai_vector->getFinestLevelNumber();
    const int ncomp = samrai_vector->getNumberOfComponents();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        auto comp_var = samrai_vector->getComponentVariable(comp);
        const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
        const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
        const bool has_cvol = cvol_idx >= 0;

        auto comp_cc_var = boost::dynamic_pointer_cast<CellVariable<double>>(comp_var);
        if (comp_cc_var)
        {
            PatchCellDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                auto level = hierarchy->getPatchLevel(ln);
                for (auto p = level->begin(); p != level->end(); ++p)
                {
                    auto patch = *p;
                    const Box& patch_box = patch->getBox();
                    auto comp_data = BOOST_CAST<CellData<double>>(patch->getPatchData(comp_idx));
                    auto cvol_data = BOOST_CAST<CellData<double>>((has_cvol ? patch->getPatchData(cvol_idx) : NULL));
                    L1_norm_local_patch.push_back(patch_ops.L1Norm(comp_data, patch_box, cvol_data));
                }
            }
        }

        auto comp_sc_var = boost::dynamic_pointer_cast<SideVariable<double>>(comp_var);
        if (comp_sc_var)
        {
            PatchSideDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                auto level = hierarchy->getPatchLevel(ln);
                for (auto p = level->begin(); p != level->end(); ++p)
                {
                    auto patch = *p;
                    const Box& patch_box = patch->getBox();
                    auto comp_data = BOOST_CAST<SideData<double>>(patch->getPatchData(comp_idx));
                    auto cvol_data = BOOST_CAST<SideData<double>>((has_cvol ? patch->getPatchData(cvol_idx) : NULL));
                    L1_norm_local_patch.push_back(patch_ops.L1Norm(comp_data, patch_box, cvol_data));
                }
            }
        }
    }
    return accurate_sum(L1_norm_local_patch);
}

double NormOps::L2Norm_local(const SAMRAIVectorReal<double>* const samrai_vector)
{
    std::vector<double> L2_norm_local_patch;
    auto hierarchy = samrai_vector->getPatchHierarchy();
    const int coarsest_ln = samrai_vector->getCoarsestLevelNumber();
    const int finest_ln = samrai_vector->getFinestLevelNumber();
    const int ncomp = samrai_vector->getNumberOfComponents();
    for (int comp = 0; comp < ncomp; ++comp)
    {
        auto comp_var = samrai_vector->getComponentVariable(comp);
        const int comp_idx = samrai_vector->getComponentDescriptorIndex(comp);
        const int cvol_idx = samrai_vector->getControlVolumeIndex(comp);
        const bool has_cvol = cvol_idx >= 0;

        auto comp_cc_var = boost::dynamic_pointer_cast<CellVariable<double>>(comp_var);
        if (comp_cc_var)
        {
            PatchCellDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                auto level = hierarchy->getPatchLevel(ln);
                for (auto p = level->begin(); p != level->end(); ++p)
                {
                    auto patch = *p;
                    const Box& patch_box = patch->getBox();
                    auto comp_data = BOOST_CAST<CellData<double>>(patch->getPatchData(comp_idx));
                    auto cvol_data = BOOST_CAST<CellData<double>>((has_cvol ? patch->getPatchData(cvol_idx) : NULL));
                    L2_norm_local_patch.push_back(patch_ops.L2Norm(comp_data, patch_box, cvol_data));
                }
            }
        }

        auto comp_sc_var = boost::dynamic_pointer_cast<SideVariable<double>>(comp_var);
        if (comp_sc_var)
        {
            PatchSideDataNormOpsReal<double> patch_ops;
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                auto level = hierarchy->getPatchLevel(ln);
                for (auto p = level->begin(); p != level->end(); ++p)
                {
                    auto patch = *p;
                    const Box& patch_box = patch->getBox();
                    auto comp_data = BOOST_CAST<SideData<double>>(patch->getPatchData(comp_idx));
                    auto cvol_data = BOOST_CAST<SideData<double>>((has_cvol ? patch->getPatchData(cvol_idx) : NULL));
                    L2_norm_local_patch.push_back(patch_ops.L2Norm(comp_data, patch_box, cvol_data));
                }
            }
        }
    }
    return std::sqrt(accurate_sum_of_squares(L2_norm_local_patch));
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
