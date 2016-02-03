// Filename: IBStandardSourceGen.cpp
// Created on 28 Apr 2011 by Boyce Griffith
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

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "boost/multi_array.hpp"
#include "ibamr/IBSourceSpec.h"
#include "ibamr/IBStandardSourceGen.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

std::vector<int> IBStandardSourceGen::s_num_sources;
std::vector<std::vector<std::string> > IBStandardSourceGen::s_source_names;
std::vector<std::vector<double> > IBStandardSourceGen::s_source_radii;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardSourceGen::IBStandardSourceGen()
    : d_n_src(), d_source_names(), d_r_src(), d_num_perimeter_nodes(), d_Q_src(), d_P_src()
{
    RestartManager::getManager()->registerRestartItem("IBStandardSourceGen", this);
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBStandardSourceGen

IBStandardSourceGen::~IBStandardSourceGen()
{
    // intentionally blank
    return;
} // ~IBStandardSourceGen

void
IBStandardSourceGen::setNumSources(const int ln, const unsigned int num_sources)
{
    s_num_sources.resize(std::max(static_cast<int>(s_num_sources.size()), ln + 1), 0);
    s_num_sources[ln] = num_sources;
    return;
} // getNumSources

unsigned int
IBStandardSourceGen::getNumSources(const int ln)
{
    return s_num_sources[ln];
} // getNumSources

const std::vector<std::string>&
IBStandardSourceGen::getSourceNames(const int ln)
{
    return s_source_names[ln];
} // getSourceNames

void
IBStandardSourceGen::setSourceNames(const int ln, const std::vector<std::string>& names)
{
    s_source_names.resize(std::max(static_cast<int>(s_source_names.size()), ln + 1));
    s_source_names[ln] = names;
    return;
} // getSourceNames

void
IBStandardSourceGen::setSourceRadii(const int ln, const std::vector<double>& radii)
{
    s_source_radii.resize(std::max(static_cast<int>(s_source_radii.size()), ln + 1));
    s_source_radii[ln] = radii;
    return;
} // getSourceRadii

const std::vector<double>&
IBStandardSourceGen::getSourceRadii(const int ln)
{
    return s_source_radii[ln];
} // getSourceRadii

std::vector<double>&
IBStandardSourceGen::getSourceStrengths(const int ln)
{
    return d_Q_src[ln];
} // getSourceStrengths

const std::vector<double>&
IBStandardSourceGen::getSourceStrengths(const int ln) const
{
    return d_Q_src[ln];
} // getSourceStrengths

const std::vector<double>&
IBStandardSourceGen::getSourcePressures(const int ln) const
{
    return d_P_src[ln];
} // getSourcePressures

void
IBStandardSourceGen::initializeLevelData(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                         const int level_number,
                                         const double /*init_data_time*/,
                                         const bool /*initial_time*/,
                                         IBTK::LDataManager* const l_data_manager)
{
    d_n_src.resize(std::max(level_number + 1, static_cast<int>(d_n_src.size())), 0);
    d_source_names.resize(std::max(level_number + 1, static_cast<int>(d_source_names.size())));
    d_r_src.resize(std::max(level_number + 1, static_cast<int>(d_r_src.size())));
    d_num_perimeter_nodes.resize(std::max(level_number + 1, static_cast<int>(d_num_perimeter_nodes.size())));
    d_Q_src.resize(std::max(level_number + 1, static_cast<int>(d_Q_src.size())));
    d_P_src.resize(std::max(level_number + 1, static_cast<int>(d_P_src.size())));

    d_n_src[level_number] = getNumSources(level_number);
    if (d_n_src[level_number] == 0) return;
#if !defined(NDEBUG)
    TBOX_ASSERT(l_data_manager->levelContainsLagrangianData(level_number));
#endif
    d_source_names[level_number] = getSourceNames(level_number);
    d_r_src[level_number] = getSourceRadii(level_number);

    d_num_perimeter_nodes[level_number].resize(d_n_src[level_number], 0);
    d_Q_src[level_number].resize(d_n_src[level_number], 0.0);
    d_P_src[level_number].resize(d_n_src[level_number], 0.0);

    std::fill(d_num_perimeter_nodes[level_number].begin(), d_num_perimeter_nodes[level_number].end(), 0);
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBSourceSpec* const spec = node_idx->getNodeDataItem<IBSourceSpec>();
        if (!spec) continue;
        const int source_idx = spec->getSourceIndex();
        ++d_num_perimeter_nodes[level_number][source_idx];
    }
    SAMRAI_MPI::sumReduction(&d_num_perimeter_nodes[level_number][0],
                             static_cast<int>(d_num_perimeter_nodes[level_number].size()));
    return;
} // initializeLevelData

unsigned int
IBStandardSourceGen::getNumSources(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                   const int level_number,
                                   const double /*data_time*/,
                                   LDataManager* const /*l_data_manager*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_n_src[level_number] >= 0);
#endif
    return d_n_src[level_number];
} // getNumSources

void
IBStandardSourceGen::getSourceLocations(std::vector<Point>& X_src,
                                        std::vector<double>& r_src,
                                        Pointer<LData> X_data,
                                        const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                        const int level_number,
                                        const double /*data_time*/,
                                        LDataManager* const l_data_manager)
{
    if (d_n_src[level_number] == 0) return;

#if !defined(NDEBUG)
    TBOX_ASSERT(X_src.size() == static_cast<unsigned int>(d_n_src[level_number]));
    TBOX_ASSERT(r_src.size() == static_cast<unsigned int>(d_n_src[level_number]));
#endif

    // Set the radii of the sources.
    r_src = d_r_src[level_number];

    // Determine the positions of the sources.
    std::fill(X_src.begin(), X_src.end(), Point::Zero());
    const double* const X_node = X_data->getLocalFormVecArray()->data();
    const Pointer<LMesh> mesh = l_data_manager->getLMesh(level_number);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
        const LNode* const node_idx = *cit;
        const IBSourceSpec* const spec = node_idx->getNodeDataItem<IBSourceSpec>();
        if (!spec) continue;
        const int& petsc_idx = node_idx->getLocalPETScIndex();
        const double* const X = &X_node[NDIM * petsc_idx];
        const int source_idx = spec->getSourceIndex();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_src[source_idx][d] += X[d] / static_cast<double>(d_num_perimeter_nodes[level_number][source_idx]);
        }
    }
    X_data->restoreArrays();

    std::vector<double> X_src_flattened(NDIM * X_src.size());
    for (unsigned int k = 0; k < X_src.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_src_flattened[NDIM * k + d] = X_src[k][d];
        }
    }
    SAMRAI_MPI::sumReduction(&X_src_flattened[0], static_cast<int>(X_src_flattened.size()));
    for (unsigned int k = 0; k < X_src.size(); ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_src[k][d] = X_src_flattened[NDIM * k + d];
        }
    }
    return;
} // getSourceLocations

void
IBStandardSourceGen::setSourcePressures(const std::vector<double>& P_src,
                                        const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                        const int level_number,
                                        const double /*data_time*/,
                                        LDataManager* const /*l_data_manager*/)
{
    d_P_src[level_number] = P_src;
    return;
} // setSourcePressures

void
IBStandardSourceGen::computeSourceStrengths(std::vector<double>& Q_src,
                                            const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                            const int level_number,
                                            const double /*data_time*/,
                                            LDataManager* const /*l_data_manager*/)
{
    Q_src = d_Q_src[level_number];
    return;
} // computeSourceStrengths

void
IBStandardSourceGen::putToDatabase(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    const int s_num_sources_sz = static_cast<int>(s_num_sources.size());
    db->putInteger("s_num_sources.size()", s_num_sources_sz);
    db->putIntegerArray("s_num_sources", &s_num_sources[0], s_num_sources_sz);
    for (unsigned int ln = 0; ln < s_num_sources.size(); ++ln)
    {
        for (int n = 0; n < s_num_sources[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->putString("s_source_names_" + id_string, s_source_names[ln][n]);
            db->putDouble("s_source_radii_" + id_string, s_source_radii[ln][n]);
        }
    }

    const int d_n_src_sz = static_cast<int>(d_n_src.size());
    db->putInteger("finest_hier_level", d_n_src_sz - 1);
    db->putIntegerArray("d_n_src", &d_n_src[0], d_n_src_sz);
    for (unsigned int ln = 0; ln < d_n_src.size(); ++ln)
    {
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            db->putString("d_source_names_" + id_string, d_source_names[ln][n]);
            db->putDouble("d_r_src_" + id_string, d_r_src[ln][n]);
            db->putInteger("d_num_perimeter_nodes_" + id_string, d_num_perimeter_nodes[ln][n]);
            db->putDouble("d_Q_src_" + id_string, d_Q_src[ln][n]);
            db->putDouble("d_P_src_" + id_string, d_P_src[ln][n]);
        }
    }
    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBStandardSourceGen::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase("IBStandardSourceGen")) // TODO: Make this ID string a variable.
    {
        db = restart_db->getDatabase("IBStandardSourceGen");
    }
    else
    {
        TBOX_ERROR("Restart database corresponding to "
                   << "IBStandardSourceGen"
                   << " not found in restart file.");
    }

    const int s_num_sources_size = db->getInteger("s_num_sources.size()");
    s_num_sources.resize(s_num_sources_size);
    s_source_names.resize(s_num_sources_size);
    s_source_radii.resize(s_num_sources_size);
    db->getIntegerArray("s_num_sources", &s_num_sources[0], s_num_sources_size);
    for (unsigned int ln = 0; ln < s_num_sources.size(); ++ln)
    {
        s_source_names[ln].resize(s_num_sources[ln]);
        s_source_radii[ln].resize(s_num_sources[ln]);
        for (int n = 0; n < s_num_sources[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            s_source_names[ln][n] = db->getString("s_source_names_" + id_string);
            s_source_radii[ln][n] = db->getDouble("s_source_radii_" + id_string);
        }
    }

    const int finest_hier_level = db->getInteger("finest_hier_level");
    d_n_src.resize(finest_hier_level + 1, 0);
    d_source_names.resize(finest_hier_level + 1);
    d_r_src.resize(finest_hier_level + 1);
    d_num_perimeter_nodes.resize(finest_hier_level + 1);
    d_Q_src.resize(finest_hier_level + 1);
    d_P_src.resize(finest_hier_level + 1);
    db->getIntegerArray("d_n_src", &d_n_src[0], finest_hier_level + 1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        d_source_names[ln].resize(d_n_src[ln]);
        d_r_src[ln].resize(d_n_src[ln], std::numeric_limits<double>::quiet_NaN());
        d_num_perimeter_nodes[ln].resize(d_n_src[ln], -1);
        d_Q_src[ln].resize(d_n_src[ln], std::numeric_limits<double>::quiet_NaN());
        d_P_src[ln].resize(d_n_src[ln], std::numeric_limits<double>::quiet_NaN());
        for (int n = 0; n < d_n_src[ln]; ++n)
        {
            std::ostringstream id_stream;
            id_stream << ln << "_" << n;
            const std::string id_string = id_stream.str();
            d_source_names[ln][n] = db->getString("d_source_names_" + id_string);
            d_r_src[ln][n] = db->getDouble("d_r_src_" + id_string);
            d_num_perimeter_nodes[ln][n] = db->getInteger("d_num_perimeter_nodes_" + id_string);
            d_Q_src[ln][n] = db->getDouble("d_Q_src_" + id_string);
            d_P_src[ln][n] = db->getDouble("d_P_src_" + id_string);
        }
    }
    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
