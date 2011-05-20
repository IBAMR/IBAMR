// Filename: IBStandardSourceGen.C
// Created on 28 Apr 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "IBStandardSourceGen.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBSourceSpec.h>
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardSourceGen::IBStandardSourceGen()
    : d_n_src(),
      d_source_names(),
      d_r_src(),
      d_num_perimeter_nodes(),
      d_Q_src(),
      d_P_src()
{
    // intentionally blank
    return;
}// IBStandardSourceGen

IBStandardSourceGen::~IBStandardSourceGen()
{
    // intentionally blank
    return;
}// ~IBStandardSourceGen

std::vector<double>&
IBStandardSourceGen::getSourceStrengths(
    const int ln)
{
    return d_Q_src[ln];
}// getSourceStrengths

const std::vector<double>&
IBStandardSourceGen::getSourceStrengths(
    const int ln) const
{
    return d_Q_src[ln];
}// getSourceStrengths

const std::vector<double>&
IBStandardSourceGen::getSourcePressures(
    const int ln) const
{
    return d_P_src[ln];
}// getSourcePressures

unsigned int
IBStandardSourceGen::getNumSources(
    const int ln) const
{
    return d_n_src[ln];
}// getNumSources

std::vector<std::string>&
IBStandardSourceGen::getSourceNames(
    const int ln)
{
    return d_source_names[ln];
}// getSourceNames

const std::vector<std::string>&
IBStandardSourceGen::getSourceNames(
    const int ln) const
{
    return d_source_names[ln];
}// getSourceNames

std::vector<double>&
IBStandardSourceGen::getSourceRadii(
    const int ln)
{
    return d_r_src[ln];
}// getSourceRadii

const std::vector<double>&
IBStandardSourceGen::getSourceRadii(
    const int ln) const
{
    return d_r_src[ln];
}// getSourceRadii

void
IBStandardSourceGen::initializeLevelData(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const l_data_manager)
{
    d_n_src              .resize(std::max(level_number+1,static_cast<int>(d_n_src.size())),0);
    d_source_names       .resize(std::max(level_number+1,static_cast<int>(d_source_names.size())));
    d_r_src              .resize(std::max(level_number+1,static_cast<int>(d_r_src.size())));
    d_num_perimeter_nodes.resize(std::max(level_number+1,static_cast<int>(d_num_perimeter_nodes.size())));
    d_Q_src              .resize(std::max(level_number+1,static_cast<int>(d_Q_src.size())));
    d_P_src              .resize(std::max(level_number+1,static_cast<int>(d_P_src.size())));

    d_n_src[level_number] = IBSourceSpec::getNumSources(level_number);
    if (d_n_src[level_number] == 0) return;

    d_source_names[level_number] = IBSourceSpec::getSourceNames(level_number);
    d_r_src[level_number] = IBSourceSpec::getSourceRadii(level_number);

    d_num_perimeter_nodes[level_number].resize(d_n_src[level_number],0);
    d_Q_src[level_number].resize(d_n_src[level_number],0.0);
    d_P_src[level_number].resize(d_n_src[level_number],0.0);

    std::fill(d_num_perimeter_nodes[level_number].begin(),d_num_perimeter_nodes[level_number].end(),0);
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBSourceSpec* const spec = node_idx->getNodeDataItem<IBSourceSpec>();
        if (spec != NULL)
        {
            const int source_idx = spec->getSourceIndex();
            ++d_num_perimeter_nodes[level_number][source_idx];
        }
    }
    SAMRAI_MPI::sumReduction(&d_num_perimeter_nodes[level_number][0],d_num_perimeter_nodes[level_number].size());
    return;
}// initializeLevelData

unsigned int
IBStandardSourceGen::getNumSources(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const l_data_manager)
{
    return d_n_src[level_number];
}// getNumSources

void
IBStandardSourceGen::getSourceLocations(
    std::vector<blitz::TinyVector<double,NDIM> >& X_src,
    std::vector<double>& r_src,
    Pointer<LData> X_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const l_data_manager)
{
    if (d_n_src[level_number] == 0) return;

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(X_src.size() == d_n_src[level_number]);
    TBOX_ASSERT(r_src.size() == d_n_src[level_number]);
#endif

    // Set the radii of the sources.
    r_src = d_r_src[level_number];

    // Determine the positions of the sources.
    for (unsigned int k = 0; k < X_src.size(); ++k) X_src[k] = 0.0;
    const double* const restrict X_node = X_data->getGhostedLocalFormVecArray()->data();
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBSourceSpec* const spec = node_idx->getNodeDataItem<IBSourceSpec>();
        if (spec != NULL)
        {
            const int& petsc_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_node[NDIM*petsc_idx];
            const int source_idx = spec->getSourceIndex();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_src[source_idx][d] += X[d]/static_cast<double>(d_num_perimeter_nodes[level_number][source_idx]);
            }
        }
    }
    X_data->restoreArrays();

    std::vector<double> X_src_flattened(NDIM*X_src.size());
    for (unsigned int k = 0; k < X_src.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_src_flattened[NDIM*k+d] = X_src[k][d];
        }
    }
    SAMRAI_MPI::sumReduction(&X_src_flattened[0],X_src_flattened.size());
    for (unsigned int k = 0; k < X_src.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_src[k][d] = X_src_flattened[NDIM*k+d];
        }
    }
    return;
}// getSourceLocations

void
IBStandardSourceGen::setSourcePressures(
    const std::vector<double>& P_src,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const l_data_manager)
{
    d_P_src[level_number] = P_src;
    return;
}// computeSourceStrengths

void
IBStandardSourceGen::computeSourceStrengths(
    std::vector<double>& Q_src,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const l_data_manager)
{
    Q_src = d_Q_src[level_number];
    return;
}// computeSourceStrengths

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBStandardSourceGen>;

//////////////////////////////////////////////////////////////////////////////
