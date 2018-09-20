// Filename: IBRedundantInitializer.cpp
// Created on 24 May 2018 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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
#include <cmath>
#include <ios>
#include <iosfwd>
#include <istream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "boost/array.hpp"
#include "boost/math/special_functions/round.hpp"
#include "boost/multi_array.hpp"
#include "ibamr/IBAnchorPointSpec.h"
#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/IBInstrumentationSpec.h"
#include "ibamr/IBRedundantInitializer.h"
#include "ibamr/IBRodForceSpec.h"
#include "ibamr/IBSourceSpec.h"
#include "ibamr/IBSpringForceSpec.h"
#include "ibamr/IBStandardSourceGen.h"
#include "ibamr/IBTargetPointForceSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/Streamable.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include <../contrib/eigen/Eigen/src/Core/products/GeneralBlockPanelKernel.h>

namespace IBTK
{
class LDataManager;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBRedundantInitializer::IBRedundantInitializer(const std::string& object_name, Pointer<Database> input_db)
    : d_object_name(object_name),
      d_max_levels(-1),
      d_level_is_initialized(),
      d_silo_writer(NULL),
      d_base_filename(),
      d_length_scale_factor(1.0),
      d_posn_shift(Vector::Zero()),
      d_num_vertex(),
      d_vertex_offset(),
      d_vertex_posn(),
      d_spring_edge_map(),
      d_spring_spec_data(),
      d_xspring_edge_map(),
      d_xspring_spec_data(),
      d_beam_spec_data(),
      d_rod_edge_map(),
      d_rod_spec_data(),
      d_target_spec_data(),
      d_anchor_spec_data(),
      d_bdry_mass_spec_data(),
      d_directors(),
      d_instrument_idx(),
      d_source_idx(),
      d_global_index_offset(),
      d_data_processed(false),
      d_init_structure_on_level_fcn(NULL),
      d_init_spring_on_level_fcn(NULL),
      d_init_beam_on_level_fcn(NULL),
      d_init_director_and_rod_on_level_fcn(NULL),
      d_init_boundary_mass_on_level_fcn(NULL),
      d_init_target_pt_on_level_fcn(NULL),
      d_init_anchor_pt_on_level_fcn(NULL),
      d_init_instrumentation_on_level_fcn(NULL),
      d_init_source_on_level_fcn(NULL)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Register the specification objects with the StreamableManager class.
    IBAnchorPointSpec::registerWithStreamableManager();
    IBBeamForceSpec::registerWithStreamableManager();
    IBInstrumentationSpec::registerWithStreamableManager();
    IBRodForceSpec::registerWithStreamableManager();
    IBSourceSpec::registerWithStreamableManager();
    IBSpringForceSpec::registerWithStreamableManager();
    IBTargetPointForceSpec::registerWithStreamableManager();

    getFromInput(input_db);
    return;
} // IBRedundantInitializer

IBRedundantInitializer::~IBRedundantInitializer()
{
    pout << d_object_name << ":  Deallocating initialization data.\n";
    return;
} // ~IBRedundantInitializer

void
IBRedundantInitializer::registerLSiloDataWriter(Pointer<LSiloDataWriter> silo_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(silo_writer);
#endif

    // Cache a pointer to the data writer.
    d_silo_writer = silo_writer;

    // Check to see if we are starting from a restart file.
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();

    // Initialize the Silo data writer only if we are not starting from a
    // restart file.
    if (!is_from_restart)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_data_processed);
#endif

        for (int ln = 0; ln < d_max_levels; ++ln)
        {
            if (d_level_is_initialized[ln])
            {
                initializeLSiloDataWriter(ln);
            }
        }
    }
    return;
} // registerLSiloDataWriter

bool
IBRedundantInitializer::getLevelHasLagrangianData(const int level_number, const bool /*can_be_refined*/) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    return !d_num_vertex[level_number].empty();
} // getLevelHasLagrangianData

unsigned int
IBRedundantInitializer::computeGlobalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                           const int level_number,
                                                           const double /*init_data_time*/,
                                                           const bool /*can_be_refined*/,
                                                           const bool /*initial_time*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    return std::accumulate(d_num_vertex[level_number].begin(), d_num_vertex[level_number].end(), 0);
}

unsigned int
IBRedundantInitializer::computeLocalNodeCountOnPatchLevel(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                          const int level_number,
                                                          const double /*init_data_time*/,
                                                          const bool /*can_be_refined*/,
                                                          const bool /*initial_time*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif
    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Loop over all patches in the specified level of the patch level and count
    // the number of local vertices.
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of vertices whose initial locations will be within
        // the given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
    }
    return local_node_count;
} // computeLocalNodeCountOnPatchLevel

void
IBRedundantInitializer::initializeStructureIndexingOnPatchLevel(
    std::map<int, std::string>& strct_id_to_strct_name_map,
    std::map<int, std::pair<int, int> >& strct_id_to_lag_idx_range_map,
    const int level_number,
    const double /*init_data_time*/,
    const bool /*can_be_refined*/,
    const bool /*initial_time*/,
    LDataManager* const /*l_data_manager*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    int offset = 0;
    for (int j = 0; j < static_cast<int>(d_base_filename[level_number].size()); ++j)
    {
        strct_id_to_strct_name_map[j] = d_base_filename[level_number][j];
        strct_id_to_lag_idx_range_map[j] = std::make_pair(offset, offset + d_num_vertex[level_number][j]);
        offset += d_num_vertex[level_number][j];
    }
    return;
} // initializeStructureIndexingOnPatchLevel

void
IBRedundantInitializer::registerInitStructureFunction(InitStructureOnLevel fcn)
{
    d_init_structure_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitSpringDataFunction(InitSpringDataOnLevel fcn)
{
    d_init_spring_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitXSpringDataFunction(InitXSpringDataOnLevel fcn)
{
    d_init_xspring_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitBeamDataFunction(InitBeamDataOnLevel fcn)
{
    d_init_beam_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitDirectorAndRodFunction(InitDirectorAndRodOnLevel fcn)
{
    d_init_director_and_rod_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitBoundaryMassFunction(InitBoundaryMassOnLevel fcn)
{
    d_init_boundary_mass_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitTargetPtFunction(InitTargetPtOnLevel fcn)
{
    d_init_target_pt_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitAnchorPtFunction(InitAnchorPtOnLevel fcn)
{
    d_init_anchor_pt_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitInstrumentationFunction(InitInstrumentationOnLevel fcn)
{
    d_init_instrumentation_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::registerInitSourceFunction(InitSourceOnLevel fcn)
{
    d_init_source_on_level_fcn = fcn;
    return;
}

void
IBRedundantInitializer::initializeStructurePosition()
{
    if (!d_init_structure_on_level_fcn)
    {
        TBOX_ERROR("IBRedundantInitializer::initializeStructurePosition()\n"
                   << " no function registered to initialize structure.\n");
    }
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_num_vertex[ln].resize(num_base_filename, 0);
        d_vertex_offset[ln].resize(num_base_filename, std::numeric_limits<int>::max());
        d_vertex_posn[ln].resize(num_base_filename);
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            if (j == 0)
            {
                d_vertex_offset[ln][j] = 0;
            }
            else
            {
                d_vertex_offset[ln][j] = d_vertex_offset[ln][j - 1] + d_num_vertex[ln][j - 1];
            }

            d_init_structure_on_level_fcn(j, ln, d_num_vertex[ln][j], d_vertex_posn[ln][j]);
#if !defined(NDEBUG)
            if (d_vertex_posn[ln][j].size() != std::size_t(d_num_vertex[ln][j]))
            {
                TBOX_ERROR(d_object_name << ":\n Invalid number of vertices " << d_vertex_posn[ln][j].size() << " of structure "
                                         << j << " on level " << ln << ".\n"
                                         << "Expected " << d_num_vertex[ln][j] << " vertices.");
            }
#endif

            // Shift and scale the position of structures
            for (unsigned int k = 0; k < unsigned(d_num_vertex[ln][j]); ++k)
            {
                Point& X = d_vertex_posn[ln][j][k];
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = d_length_scale_factor * (X[d] + d_posn_shift[d]);
                }
            }
        }
    }

    return;
} // initializeStructurePosition

void
IBRedundantInitializer::initializeSprings()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_spring_edge_map[ln].resize(num_base_filename);
        d_spring_spec_data[ln].resize(num_base_filename);
        if (d_init_spring_on_level_fcn)
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_spring_on_level_fcn(j, ln, d_spring_edge_map[ln][j], d_spring_spec_data[ln][j]);

                int min_idx = 0;
                int max_idx = d_num_vertex[ln][j];
                for (std::multimap<int, Edge>::iterator it = d_spring_edge_map[ln][j].begin();
                     it != d_spring_edge_map[ln][j].end();
                     ++it)
                {
                    Edge& e = it->second;
                    const SpringSpec& spec = d_spring_spec_data[ln][j][e];
                    if ((e.first < min_idx) || (e.first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid spring edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.first << " is not a valid index.");
                    }
                    if ((e.second < min_idx) || (e.second > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid spring edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.second << " is not a valid index.");
                    }
                    if (it->first > e.second)
                    {
                        TBOX_ERROR(d_object_name << ":\n Error on level " << ln << " and structure number " << j
                                                 << ".\n Master index must be lower than the slave index for springs.");
                    }
                    if (spec.parameters[0] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid spring constant encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << spec.parameters[0] << " for index " << e.first << " is negative.");
                    }
                    if (spec.parameters[1] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid resting length encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << spec.parameters[1] << " for index " << e.first << " is negative.");
                    }
                }
            }
        }
    }
    return;
} // initializeSprings

void
IBRedundantInitializer::initializeXSprings()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_xspring_edge_map[ln].resize(num_base_filename);
        d_xspring_spec_data[ln].resize(num_base_filename);
        if (d_init_xspring_on_level_fcn)
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_xspring_on_level_fcn(j, ln, d_xspring_edge_map[ln][j], d_xspring_spec_data[ln][j]);
                const int min_idx = 0;
                const int max_idx = std::accumulate(d_num_vertex[ln].begin(), d_num_vertex[ln].end(), 0);
                for (std::multimap<int, Edge>::iterator it = d_xspring_edge_map[ln][j].begin();
                     it != d_xspring_edge_map[ln][j].end();
                     ++it)
                {
                    Edge& e = it->second;
                    const XSpringSpec& spec = d_xspring_spec_data[ln][j][e];
                    if ((e.first < min_idx) || (e.first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid xspring edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.first << " is not a valid index.");
                    }
                    if ((e.second < min_idx) || (e.second > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid xspring edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.second << " is not a valid index.");
                    }
                    if (it->first > e.second)
                    {
                        TBOX_ERROR(d_object_name
                                   << ":\n Error on level " << ln << " and structure number " << j
                                   << ".\n Master index must be lower than the slave index for xsprings.");
                    }
                    if (spec.parameters[0] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid spring constant encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << spec.parameters[0] << " for index " << e.first << " is negative.");
                    }
                    if (spec.parameters[1] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid resting length encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << spec.parameters[1] << " for index " << e.first << " is negative.");
                    }
                }
            }
        }
    }
    return;
} // initializeXSprings

void
IBRedundantInitializer::initializeBeams()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_beam_spec_data[ln].resize(num_base_filename);
        if (d_init_beam_on_level_fcn)
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_beam_on_level_fcn(j, ln, d_beam_spec_data[ln][j]);

                const int min_idx = 0;
                const int max_idx = d_num_vertex[ln][j];
                for (std::multimap<int, BeamSpec>::const_iterator it = d_beam_spec_data[ln][j].begin();
                     it != d_beam_spec_data[ln][j].end();
                     ++it)
                {
                    const BeamSpec& e = it->second;
                    const std::pair<int, int>& idxs = e.neighbor_idxs;
                    if ((idxs.first < min_idx) || (idxs.first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid beam edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << idxs.first << " is not a valid index");
                    }
                    if ((idxs.second < min_idx) || (idxs.second > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid beam edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << idxs.second << " is not a valid index");
                    }
                    if (e.bend_rigidity < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid bending rigidity encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.bend_rigidity << " for index " << it->first << " is negative");
                    }
                }
            }
        }
    }
    return;
} // initializeBeams

void
IBRedundantInitializer::initializeTargetPts()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        TargetSpec default_spec;
        default_spec.stiffness = 0.0;
        default_spec.damping = 0.0;
        d_target_spec_data[ln].resize(num_base_filename);
        std::multimap<int, TargetSpec> tg_pt_spec;
        if (d_init_target_pt_on_level_fcn)
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_target_pt_on_level_fcn(j, ln, tg_pt_spec);

                int min_idx = 0;
                int max_idx = d_num_vertex[ln][j];
                d_target_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);
                for (std::multimap<int, TargetSpec>::const_iterator it = tg_pt_spec.begin(); it != tg_pt_spec.end();
                     ++it)
                {
                    if ((it->first < min_idx) || (it->first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid target point index on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << it->first);
                    }
                    const TargetSpec& tg_spec = it->second;
                    if (tg_spec.stiffness < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid target point stiffness encountered on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << tg_spec.stiffness << " is negative");
                    }
                    if (tg_spec.damping < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid target point damping encountered on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << tg_spec.damping << " is negative");
                    }
                    d_target_spec_data[ln][j][it->first] = tg_spec;
                }
            }
        }
        else
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_target_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);
            }
        }
    }
    return;
} // initializeTargetPts

void
IBRedundantInitializer::initializeDirectorAndRods()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_directors[ln].resize(num_base_filename);
        d_rod_edge_map[ln].resize(num_base_filename);
        d_rod_spec_data[ln].resize(num_base_filename);
        if (d_init_director_and_rod_on_level_fcn)
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_director_and_rod_on_level_fcn(
                    j, ln, d_directors[ln][j], d_rod_edge_map[ln][j], d_rod_spec_data[ln][j]);

                const int min_idx = 0;
                const int max_idx = d_num_vertex[ln][j];
                if (d_directors[ln][j].size() != unsigned(max_idx))
                {
                    TBOX_ERROR(d_object_name << "\n Not enough director vectors supplied for structure " << j
                                             << "on level " << ln << ".");
                }
                for (unsigned int k = 0; k < d_directors[ln][j].size(); ++k)
                {
                    if (d_directors[ln][j][k].size() != 9)
                    {
                        TBOX_ERROR(d_object_name << "\n Not enough director vectors supplied on index " << k
                                                 << " for structure " << j << "on level " << ln << ".\n");
                    }
                    for (int n = 0; n < 3; ++n)
                    {
                        double D_norm_squared = 0.0;
                        for (int d = 0; d < 3; ++d)
                        {
                            D_norm_squared += d_directors[ln][j][k][3 * n + d] * d_directors[ln][j][k][3 * n + d];
                        }
                        const double D_norm = sqrt(D_norm_squared);
                        if (!MathUtilities<double>::equalEps(D_norm, 1.0))
                        {
                            TBOX_WARNING(d_object_name << ":\n  Director vector for index " << k << " of structure "
                                                       << j << " on level " << ln
                                                       << " is not normalized; norm = " << D_norm << std::endl);
                            for (int d = 0; d < 3; ++d)
                            {
                                d_directors[ln][j][k][3 * n + d] /= D_norm;
                            }
                        }
                    }
                }
                for (std::multimap<int, Edge>::const_iterator it = d_rod_edge_map[ln][j].begin();
                     it != d_rod_edge_map[ln][j].end();
                     ++it)
                {
                    const Edge& e = it->second;
                    const RodSpec& rod_spec = d_rod_spec_data[ln][j][e];
                    const boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> parameters = rod_spec.properties;
                    if ((e.first < min_idx) || (e.first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid rod edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.first << " is not a valid index");
                    }
                    if ((e.second < min_idx) || (e.second > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid rod edge encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << e.second << " is not a valid index");
                    }
                    if (parameters[0] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid rod length encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[0] << "for index " << e.first << " is negative");
                    }
                    if (parameters[1] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid parameter a1 encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[1] << "for index " << e.first << " is negative");
                    }
                    if (parameters[2] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid parameter a2 encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[2] << "for index " << e.first << " is negative");
                    }
                    if (parameters[3] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid a3 length encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[3] << "for index " << e.first << " is negative");
                    }
                    if (parameters[4] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid b1 encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[4] << "for index " << e.first << " is negative");
                    }
                    if (parameters[5] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid b2 encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[5] << "for index " << e.first << " is negative");
                    }
                    if (parameters[6] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid b3 encountered on level " << ln
                                                 << " and structure number " << j << ":\n"
                                                 << parameters[6] << "for index " << e.first << " is negative");
                    }
                }
            }
        }
    }
    return;
} // initializeDirectorAndRods

void
IBRedundantInitializer::initializeBoundaryMass()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        BdryMassSpec default_spec;
        default_spec.bdry_mass = 0.0;
        default_spec.stiffness = 0.0;
        d_bdry_mass_spec_data[ln].resize(num_base_filename);
        if (d_init_boundary_mass_on_level_fcn)
        {
            std::multimap<int, BdryMassSpec> bdry_map;
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_boundary_mass_on_level_fcn(j, ln, bdry_map);
                d_bdry_mass_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);
                int min_idx = 0;
                int max_idx = d_num_vertex[ln][j];
                for (std::multimap<int, BdryMassSpec>::const_iterator it = bdry_map.begin(); it != bdry_map.end(); ++it)
                {
                    if ((it->first < min_idx) || (it->first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid massive point index on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << it->first);
                    }
                    const BdryMassSpec& bdry_mass_spec = it->second;
                    if (bdry_mass_spec.bdry_mass < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid boundary mass encountered on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << bdry_mass_spec.bdry_mass << " is negative");
                    }
                    if (bdry_mass_spec.stiffness < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid boundary mass stiffness encountered on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << bdry_mass_spec.stiffness << " is negative");
                    }
                    d_bdry_mass_spec_data[ln][j][it->first] = bdry_mass_spec;
                }
            }
        }
        else
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_bdry_mass_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);
            }
        }
    }
    return;
}

void
IBRedundantInitializer::initializeAnchorPts()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_anchor_spec_data[ln].resize(num_base_filename);
        AnchorSpec default_spec;
        default_spec.is_anchor_point = false;
        if (d_init_anchor_pt_on_level_fcn)
        {
            std::multimap<int, AnchorSpec> anchor_map;
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_anchor_pt_on_level_fcn(j, ln, anchor_map);
                d_anchor_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);
                int min_idx = 0;
                int max_idx = d_num_vertex[ln][j];
                for (std::multimap<int, AnchorSpec>::const_iterator it = anchor_map.begin(); it != anchor_map.end();
                     ++it)
                {
                    if ((it->first < min_idx) || (it->first > max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid anchor point index on level " << ln
                                                 << " and structure number " << j << ": \n"
                                                 << it->first);
                    }
                    d_anchor_spec_data[ln][j][it->first] = it->second;
                }
            }
        }
        else
        {
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_anchor_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);
            }
        }
    }
    return;
}

void
IBRedundantInitializer::initializeInstrumentationData()
{
    std::vector<std::string> instrument_names;
    int instrument_offset = 0;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const size_t num_base_filename = d_base_filename[ln].size();
        d_instrument_idx[ln].resize(num_base_filename);
        if (d_init_instrumentation_on_level_fcn)
        {
            std::vector<std::string> new_names;
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                d_init_instrumentation_on_level_fcn(j, ln, new_names, d_instrument_idx[ln][j]);
                std::vector<bool> encountered_instrument_idx;
                std::map<int, std::vector<bool> > encountered_node_idxs;
                for (std::vector<std::string>::iterator i = new_names.begin(); i != new_names.end(); ++i)
                    instrument_names.push_back(*i);
                const int min_idx = 0;
                const int max_idx = d_num_vertex[ln][j];
                for (std::map<int, std::pair<int, int> >::iterator it = d_instrument_idx[ln][j].begin();
                     it != d_instrument_idx[ln][j].end();
                     ++it)
                {
                    if ((it->first < min_idx) || (it->first >= max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid instrument master index on level " << ln
                                                 << " and structure number " << j << ".\n Vertex " << it->first
                                                 << " is out of range.\n");
                    }
                    std::pair<int, int>& meter_map = it->second;
                    if (meter_map.first < 0 || unsigned(meter_map.first) >= instrument_names.size())
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid meter number on level " << ln
                                                 << " and structure number " << j << ".\n Meter index "
                                                 << meter_map.first << " is out of range.\n");
                    }

                    if (meter_map.second < 0)
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid node index on level " << ln << " and structure number "
                                                 << j << ".\n Node Index " << meter_map.second << " is invalid.\n");
                    }

                    if (meter_map.first >= static_cast<int>(encountered_instrument_idx.size()))
                    {
                        encountered_instrument_idx.resize(meter_map.first + 1, false);
                    }
                    encountered_instrument_idx[meter_map.first] = true;

                    if (meter_map.second >= static_cast<int>(encountered_node_idxs[meter_map.first].size()))
                    {
                        encountered_node_idxs[meter_map.first].resize(meter_map.second + 1, false);
                    }
                    encountered_node_idxs[meter_map.first][meter_map.second] = true;

                    meter_map.first += instrument_offset;
                }
                for (std::vector<bool>::iterator meter_it = encountered_instrument_idx.begin();
                     meter_it != encountered_instrument_idx.end();
                     ++meter_it)
                {
                    const size_t meter_idx = std::distance(encountered_instrument_idx.begin(), meter_it);
                    if ((*meter_it) == false)
                    {
                        TBOX_ERROR(d_object_name << ":\n Instrument index " << meter_idx << " not found on level " << ln
                                                 << " and structure number " << j);
                    }

                    std::vector<bool>& meter_node_idxs = encountered_node_idxs[meter_idx];
                    for (std::vector<bool>::iterator node_it = meter_node_idxs.begin();
                         node_it != meter_node_idxs.end();
                         ++node_it)
                    {
                        const size_t node_idx = std::distance(meter_node_idxs.begin(), node_it);
                        if ((*node_it) == false)
                        {
                            TBOX_ERROR(d_object_name << ":\n Node index " << node_idx << " associated with meter index "
                                                     << meter_idx << "not found on level " << ln
                                                     << " and structure number " << j);
                        }
                    }
                }
                if (encountered_instrument_idx.size() != new_names.size())
                {
                    TBOX_ERROR(d_object_name << ":\n Not all anticipated instrument indices were found on level " << ln
                                             << " and structure number " << j << ". Expected to find "
                                             << new_names.size() << " distinct meter indices.");
                }
            }
        }
    }
    IBInstrumentationSpec::setInstrumentNames(instrument_names);
    return;
}

void
IBRedundantInitializer::initializeSourceData()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        std::vector<std::string> source_names;
        std::vector<double> source_radii;
        int source_offset = 0;
        const size_t num_base_filename = d_base_filename[ln].size();
        d_source_idx[ln].resize(num_base_filename);
        if (d_init_source_on_level_fcn)
        {
            std::vector<std::string> new_names;
            std::vector<double> new_radii;
            for (unsigned int j = 0; j < num_base_filename; ++j)
            {
                const int min_idx = 0;
                const int max_idx = d_num_vertex[ln][j];
                int num_source;
                d_init_source_on_level_fcn(j, ln, d_source_idx[ln][j], new_names, new_radii);
                if (source_names.size() != source_radii.size())
                {
                    TBOX_ERROR(d_object_name << ":\n Invalid number of sources/sinks. Number of sources "
                                             << source_names.size() << " is not equal to number of radii "
                                             << source_radii.size() << ".\n");
                }
                for (std::vector<std::string>::iterator i = new_names.begin(); i != new_names.end(); ++i)
                    source_names.push_back(*i);
                for (std::vector<double>::iterator i = new_radii.begin(); i != new_radii.end(); ++i)
                    source_radii.push_back(*i);
                num_source = new_names.size();
                for (std::map<int, int>::iterator it = d_source_idx[ln][j].begin(); it != d_source_idx[ln][j].end();
                     ++it)
                {
                    int& src_num = it->second;
                    if ((it->first < min_idx) || (it->first >= max_idx))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid source vertex on level " << ln
                                                 << " and structure number " << j << ".\n Source vertex " << it->first
                                                 << ".\n");
                    }
                    if ((src_num < 0) || (src_num < num_source))
                    {
                        TBOX_ERROR(d_object_name << ":\n Invalid source number on level " << ln
                                                 << " and structure number " << j << ".\n Source number " << src_num
                                                 << " but there are " << num_source << " Sources.\n");
                    }
                    src_num += source_offset;
                }
                source_offset += num_source;
            }
            IBStandardSourceGen::setNumSources(ln, source_offset);
            IBStandardSourceGen::setSourceNames(ln, source_names);
            IBStandardSourceGen::setSourceRadii(ln, source_radii);
        }
    }
    return;
}

unsigned int
IBRedundantInitializer::initializeDataOnPatchLevel(const int lag_node_index_idx,
                                                   const unsigned int global_index_offset,
                                                   const unsigned int local_index_offset,
                                                   Pointer<LData> X_data,
                                                   Pointer<LData> U_data,
                                                   const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                   const int level_number,
                                                   const double /*init_data_time*/,
                                                   const bool /*can_be_refined*/,
                                                   const bool /*initial_time*/,
                                                   LDataManager* const /*l_data_manager*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    Vector domain_length;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        domain_length[d] = domain_x_upper[d] - domain_x_lower[d];
    }

    // Set the global index offset.  This is equal to the number of Lagrangian
    // indices that have already been initialized on the specified level.
    d_global_index_offset[level_number] = global_index_offset;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
    boost::multi_array_ref<double, 2>& U_array = *U_data->getLocalFormVecArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(ratio);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();

        Pointer<LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int lagrangian_idx = getCanonicalLagrangianIndex(point_idx, level_number) + global_index_offset;
            const int local_petsc_idx = ++local_idx + local_index_offset;
            const int global_petsc_idx = local_petsc_idx + global_index_offset;

            // Get the coordinates and periodic shifters of the present vertex.
            Point X_real = getVertexPosn(point_idx, level_number);
            Point X = getShiftedVertexPosn(point_idx, level_number, domain_x_lower, domain_x_upper, periodic_shift);
            Vector periodic_displacement = X_real - X;
            IntVector<NDIM> periodic_offset;
            for (int d = 0; d < NDIM; ++d)
            {
                periodic_offset[d] = boost::math::round(periodic_displacement[d] / patch_dx[d]);
            }

            // Ensure that all points are initially within the computational
            // domain.
            for (int d = 0; d < NDIM; ++d)
            {
                if (!periodic_shift[d] && X[d] < domain_x_lower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                                             << "  encountered node below lower physical boundary\n"
                                             << "  please ensure that all nodes are within the "
                                                "computational domain."
                                             << std::endl);
                }

                if (!periodic_shift[d] && X[d] >= domain_x_upper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                                             << "  encountered node above upper physical boundary\n"
                                             << "  please ensure that all nodes are within the "
                                                "computational domain."
                                             << std::endl);
                }
            }

            // Set X_array.
            for (int d = 0; d < NDIM; ++d)
            {
                X_array[local_petsc_idx][d] = X[d];
            }

            // Get the index of the cell in which the present vertex is
            // initially located.
            const CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, grid_geom, ratio);

            // Initialize the specification objects associated with the present
            // vertex.
            std::vector<Pointer<Streamable> > node_data =
                initializeNodeData(point_idx, global_index_offset, level_number);
            for (std::vector<Pointer<Streamable> >::iterator it = node_data.begin(); it != node_data.end(); ++it)
            {
                (*it)->registerPeriodicShift(periodic_offset, periodic_displacement);
            }

            // Create or retrieve a pointer to the LNodeSet associated with the
            // current Cartesian grid cell.
            if (!index_data->isElement(idx))
            {
                index_data->appendItemPointer(idx, new LNodeSet());
            }
            LNodeSet* const node_set = index_data->getItem(idx);
            node_set->push_back(new LNode(lagrangian_idx,
                                          global_petsc_idx,
                                          local_petsc_idx,
                                          /*initial*/ periodic_offset,
                                          /*current*/ periodic_offset,
                                          /*initial*/ periodic_displacement,
                                          /*current*/ periodic_displacement,
                                          node_data));

            // Initialize the velocity of the present vertex.
            std::fill(&U_array[local_petsc_idx][0], &U_array[local_petsc_idx][0] + NDIM, 0.0);
        }
    }
    X_data->restoreArrays();
    U_data->restoreArrays();

    d_level_is_initialized[level_number] = true;

    // If a Lagrangian Silo data writer is registered with the initializer,
    // setup the visualization data corresponding to the present level of the
    // locally refined grid.
    if (d_silo_writer)
    {
        initializeLSiloDataWriter(level_number);
    }
    return local_node_count;
} // initializeDataOnPatchLevel

unsigned int
IBRedundantInitializer::initializeMassDataOnPatchLevel(const unsigned int /*global_index_offset*/,
                                                       const unsigned int local_index_offset,
                                                       Pointer<LData> M_data,
                                                       Pointer<LData> K_data,
                                                       const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                       const int level_number,
                                                       const double /*init_data_time*/,
                                                       const bool /*can_be_refined*/,
                                                       const bool /*initial_time*/,
                                                       LDataManager* const /*l_data_manager*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 1>& M_array = *M_data->getLocalFormArray();
    boost::multi_array_ref<double, 1>& K_array = *K_data->getLocalFormArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int local_petsc_idx = ++local_idx + local_index_offset;

            // Initialize the mass and penalty stiffness coefficient
            // corresponding to the present vertex.
            const BdryMassSpec& spec = getVertexBdryMassSpec(point_idx, level_number);
            const double M = spec.bdry_mass;
            const double K = spec.stiffness;

            // Avoid division by zero at massless nodes.
            if (MathUtilities<double>::equalEps(M, 0.0))
            {
                M_array[local_petsc_idx] = std::numeric_limits<double>::epsilon();
                K_array[local_petsc_idx] = 0.0;
            }
            else
            {
                M_array[local_petsc_idx] = M;
                K_array[local_petsc_idx] = K;
            }
        }
    }
    M_data->restoreArrays();
    K_data->restoreArrays();
    return local_node_count;
} // initializeMassOnPatchLevel

unsigned int
IBRedundantInitializer::initializeDirectorDataOnPatchLevel(const unsigned int /*global_index_offset*/,
                                                           const unsigned int local_index_offset,
                                                           Pointer<LData> D_data,
                                                           const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                           const int level_number,
                                                           const double /*init_data_time*/,
                                                           const bool /*can_be_refined*/,
                                                           const bool /*initial_time*/,
                                                           LDataManager* const /*l_data_manager*/)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    boost::multi_array_ref<double, 2>& D_array = *D_data->getLocalFormVecArray();
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int, int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, hierarchy);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin(); it != patch_vertices.end();
             ++it)
        {
            const std::pair<int, int>& point_idx = (*it);
            const int local_petsc_idx = ++local_idx + local_index_offset;

            // Initialize the director corresponding to the present vertex.
            const std::vector<double>& D = getVertexDirectors(point_idx, level_number);
            for (int d = 0; d < 3 * 3; ++d)
            {
                D_array[local_petsc_idx][d] = D[d];
            }
        }
    }
    D_data->restoreArrays();
    return local_node_count;
} // initializeDirectorOnPatchLevel

void
IBRedundantInitializer::tagCellsForInitialRefinement(const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                     const int level_number,
                                                     const double /*error_data_time*/,
                                                     const int tag_index)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_data_processed);
#endif

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();

    // Loop over all patches in the specified level of the patch level and tag
    // cells for refinement wherever there are vertices assigned to a finer
    // level of the Cartesian grid.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(ratio);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();

        Pointer<CellData<NDIM, int> > tag_data = patch->getPatchData(tag_index);

        // Tag cells for refinement whenever there are vertices whose initial
        // locations will be within the index space of the given patch, but on
        // the finer levels of the AMR patch hierarchy.
        for (int ln = level_number + 1; ln < d_max_levels; ++ln)
        {
            std::vector<std::pair<int, int> > patch_vertices;
            getPatchVerticesAtLevel(patch_vertices, patch, hierarchy, ln);
            for (std::vector<std::pair<int, int> >::const_iterator it = patch_vertices.begin();
                 it != patch_vertices.end();
                 ++it)
            {
                const std::pair<int, int>& point_idx = (*it);

                // Get the coordinates of the present vertex.
                const Point& X = getShiftedVertexPosn(point_idx, ln, domain_x_lower, domain_x_upper, periodic_shift);

                // Get the index of the cell in which the present vertex is
                // initially located.
                const CellIndex<NDIM> i = IndexUtilities::getCellIndex(X, grid_geom, ratio);

                // Tag the cell for refinement.
                if (patch_box.contains(i)) (*tag_data)(i) = 1;
            }
        }
    }
    return;
} // tagCellsForInitialRefinement

void
IBRedundantInitializer::setStructureNamesOnLevel(const int& level_num, const std::vector<std::string>& strct_names)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level_num >= 0);
    TBOX_ASSERT(level_num < d_max_levels);
#endif
    d_base_filename[level_num] = strct_names;
    return;
}

void
IBRedundantInitializer::init()
{
    if (d_data_processed)
    {
        return;
    }
    else
    {
        // Process structure information.
        initializeStructurePosition();
        initializeSprings();
        initializeXSprings();
        initializeBeams();
        initializeDirectorAndRods();
        initializeBoundaryMass();
        initializeTargetPts();
        initializeAnchorPts();
        initializeInstrumentationData();
        initializeSourceData();
    }

    // Indicate that we have processed data.
    d_data_processed = true;

    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBRedundantInitializer::initializeLSiloDataWriter(const int level_number)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(level_number < d_max_levels);
    TBOX_ASSERT(d_level_is_initialized[level_number]);
#endif

    // WARNING: This code does not work if the global node offset is nonzero on
    // any of the levels of the locally refined Cartesian grid.
    if (d_global_index_offset[level_number] != 0)
    {
        TBOX_ERROR("This is broken --- please submit a bug report if you encounter this error.\n");
    }

    // WARNING: For now, we just register the visualization data on MPI process
    // 0.  This will fail if the structure is too large to be stored in the
    // memory available to a single MPI process.
    if (SAMRAI_MPI::getRank() == 0)
    {
        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_num_vertex[level_number][j] > 0)
            {
                const std::string postfix = "_vertices";
                d_silo_writer->registerMarkerCloud(d_base_filename[level_number][j] + postfix,
                                                   d_num_vertex[level_number][j],
                                                   d_vertex_offset[level_number][j],
                                                   level_number);
            }
        }

        bool registered_spring_edge_map = false;
        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_spring_edge_map[level_number][j].size() > 0)
            {
                registered_spring_edge_map = true;
                const std::string postfix = "_mesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix, d_spring_edge_map[level_number][j], level_number);
            }
        }

        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_xspring_edge_map[level_number][j].size() > 0)
            {
                const std::string postfix = "_xmesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix, d_xspring_edge_map[level_number][j], level_number);
            }
        }

        for (unsigned int j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_rod_edge_map[level_number][j].size() > 0)
            {
                const std::string postfix = registered_spring_edge_map ? "_rod_mesh" : "_mesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix, d_rod_edge_map[level_number][j], level_number);
            }
        }
    }
    return;
} // initializeLSiloDataWriter

void
IBRedundantInitializer::getPatchVertices(std::vector<std::pair<int, int> >& patch_vertices,
                                         const Pointer<Patch<NDIM> > patch,
                                         const Pointer<PatchHierarchy<NDIM> > hierarchy) const
{
#if !defined(NDEBUG)
    TBOX_ASSERT(patch->inHierarchy());
#endif
    const int level_number = patch->getPatchLevelNumber();
    getPatchVerticesAtLevel(patch_vertices, patch, hierarchy, level_number);
    return;
} // getPatchVertices

void
IBRedundantInitializer::getPatchVerticesAtLevel(std::vector<std::pair<int, int> >& patch_vertices,
                                                const Pointer<Patch<NDIM> > patch,
                                                const Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                const int vertex_level_number) const
{
    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
#if !defined(NDEBUG)
    TBOX_ASSERT(patch->inHierarchy());
#endif
    const int level_number = patch->getPatchLevelNumber();
    const Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    const IntVector<NDIM>& ratio = level->getRatio();
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift(ratio);

    // Loop over all of the vertices to determine the indices of those vertices
    // within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will work for
    // now.
    const Box<NDIM>& patch_box = patch->getBox();
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
    for (unsigned int j = 0; j < d_num_vertex[vertex_level_number].size(); ++j)
    {
        for (int k = 0; k < d_num_vertex[vertex_level_number][j]; ++k)
        {
            std::pair<int, int> point_index(j, k);
            const Point& X =
                getShiftedVertexPosn(point_index, vertex_level_number, domain_x_lower, domain_x_upper, periodic_shift);
            const CellIndex<NDIM> idx = IndexUtilities::getCellIndex(X, grid_geom, ratio);
            if (patch_box.contains(idx)) patch_vertices.push_back(point_index);
        }
    }
    return;
} // getPatchVerticesAtLevel

int
IBRedundantInitializer::getCanonicalLagrangianIndex(const std::pair<int, int>& point_index,
                                                    const int level_number) const
{
    return d_vertex_offset[level_number][point_index.first] + point_index.second;
} // getCanonicalLagrangianIndex

Point
IBRedundantInitializer::getVertexPosn(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_vertex_posn[level_number][point_index.first][point_index.second];
} // getVertexPosn

Point
IBRedundantInitializer::getShiftedVertexPosn(const std::pair<int, int>& point_index,
                                             const int level_number,
                                             const double* const domain_x_lower,
                                             const double* const domain_x_upper,
                                             const IntVector<NDIM>& periodic_shift) const
{
    Point X = getVertexPosn(point_index, level_number);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (periodic_shift[d])
        {
            double domain_length = domain_x_upper[d] - domain_x_lower[d];
            while (X[d] < domain_x_lower[d]) X[d] += domain_length;
            while (X[d] >= domain_x_upper[d]) X[d] -= domain_length;
            TBOX_ASSERT(X[d] >= domain_x_lower[d] && X[d] < domain_x_upper[d]);
            X[d] = std::max(X[d], domain_x_lower[d]);
            X[d] = std::min(X[d], domain_x_upper[d] - std::numeric_limits<double>::epsilon());
        }
    }
    return X;
} // getShiftedVertexPosn

const IBRedundantInitializer::TargetSpec&
IBRedundantInitializer::getVertexTargetSpec(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_target_spec_data[level_number][point_index.first][point_index.second];
} // getVertexTargetSpec

const IBRedundantInitializer::AnchorSpec&
IBRedundantInitializer::getVertexAnchorSpec(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_anchor_spec_data[level_number][point_index.first][point_index.second];
} // getVertexAnchorSpec

const IBRedundantInitializer::BdryMassSpec&
IBRedundantInitializer::getVertexBdryMassSpec(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_bdry_mass_spec_data[level_number][point_index.first][point_index.second];
} // getVertexBdryMassSpec

const std::vector<double>&
IBRedundantInitializer::getVertexDirectors(const std::pair<int, int>& point_index, const int level_number) const
{
    return d_directors[level_number][point_index.first][point_index.second];
} // getVertexDirectors

std::pair<int, int>
IBRedundantInitializer::getVertexInstrumentationIndices(const std::pair<int, int>& point_index,
                                                        const int level_number) const
{
    std::map<int, std::pair<int, int> >::const_iterator it =
        d_instrument_idx[level_number][point_index.first].find(point_index.second);
    if (it != d_instrument_idx[level_number][point_index.first].end())
    {
        return it->second;
    }
    else
    {
        return std::make_pair(-1, -1);
    }
} // getVertexInstrumentationIndices

int
IBRedundantInitializer::getVertexSourceIndices(const std::pair<int, int>& point_index, const int level_number) const
{
    std::map<int, int>::const_iterator it = d_source_idx[level_number][point_index.first].find(point_index.second);
    if (it != d_source_idx[level_number][point_index.first].end())
    {
        return it->second;
    }
    else
    {
        return -1;
    }
} // getVertexSourceIndices

std::vector<Pointer<Streamable> >
IBRedundantInitializer::initializeNodeData(const std::pair<int, int>& point_index,
                                           const unsigned int global_index_offset,
                                           const int level_number) const
{
    std::vector<Pointer<Streamable> > node_data;

    const int j = point_index.first;
    const int mastr_idx = getCanonicalLagrangianIndex(point_index, level_number);

    // Initialize any spring specifications associated with the present vertex.
    {
        std::vector<int> slave_idxs, force_fcn_idxs;
        std::vector<std::vector<double> > parameters;
        for (std::multimap<int, Edge>::const_iterator it = d_spring_edge_map[level_number][j].lower_bound(mastr_idx);
             it != d_spring_edge_map[level_number][j].upper_bound(mastr_idx);
             ++it)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mastr_idx == it->first);
#endif
            // The connectivity information.d(
            const Edge& e = it->second;
            if (e.first == mastr_idx)
            {
                slave_idxs.push_back(e.second + global_index_offset);
            }
            else
            {
                slave_idxs.push_back(e.first + global_index_offset);
            }

            // The material properties.
            const SpringSpec& spec_data = d_spring_spec_data[level_number][j].find(e)->second;
            parameters.push_back(spec_data.parameters);
            force_fcn_idxs.push_back(spec_data.force_fcn_idx);
        }
        const size_t num_base_filename = d_base_filename[level_number].size();
        for (unsigned int j = 0; j < num_base_filename; ++j)
        {
            for (std::multimap<int, Edge>::const_iterator it =
                     d_xspring_edge_map[level_number][j].lower_bound(mastr_idx);
                 it != d_xspring_edge_map[level_number][j].upper_bound(mastr_idx);
                 ++it)
            {
#if !defined(NDEBUG)
                TBOX_ASSERT(mastr_idx == it->first);
#endif
                // The connectivity information.
                const Edge& e = it->second;
                if (e.first == mastr_idx)
                {
                    slave_idxs.push_back(e.second + global_index_offset);
                }
                else
                {
                    slave_idxs.push_back(e.first + global_index_offset);
                }

                // The material properties.
                const XSpringSpec& spec_data = d_xspring_spec_data[level_number][j].find(e)->second;
                parameters.push_back(spec_data.parameters);
                force_fcn_idxs.push_back(spec_data.force_fcn_idx);
            }
        }
        if (slave_idxs.size() > 0)
        {
            node_data.push_back(new IBSpringForceSpec(mastr_idx, slave_idxs, force_fcn_idxs, parameters));
        }
    }

    // Initialize any beam specifications associated with the present vertex.
    {
        std::vector<std::pair<int, int> > beam_neighbor_idxs;
        std::vector<double> beam_bend_rigidity;
        std::vector<Vector> beam_mesh_dependent_curvature;
        for (std::multimap<int, BeamSpec>::const_iterator it = d_beam_spec_data[level_number][j].lower_bound(mastr_idx);
             it != d_beam_spec_data[level_number][j].upper_bound(mastr_idx);
             ++it)
        {
            const BeamSpec& spec_data = it->second;
            beam_neighbor_idxs.push_back(spec_data.neighbor_idxs);
            beam_bend_rigidity.push_back(spec_data.bend_rigidity);
            beam_mesh_dependent_curvature.push_back(spec_data.curvature);
        }
        if (!beam_neighbor_idxs.empty())
        {
            node_data.push_back(
                new IBBeamForceSpec(mastr_idx, beam_neighbor_idxs, beam_bend_rigidity, beam_mesh_dependent_curvature));
        }
    }

    // Initialize any rod specifications associated with the present vertex.
    {
        std::vector<int> rod_next_idxs;
        std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> > rod_material_params;
        for (std::multimap<int, Edge>::const_iterator it = d_rod_edge_map[level_number][j].lower_bound(mastr_idx);
             it != d_rod_edge_map[level_number][j].upper_bound(mastr_idx);
             ++it)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(mastr_idx == it->first);
#endif
            // The connectivity information.
            const Edge& e = it->second;
            if (e.first == mastr_idx)
            {
                rod_next_idxs.push_back(e.second + global_index_offset);
            }
            else
            {
                rod_next_idxs.push_back(e.first + global_index_offset);
            }

            // The material properties.
            const RodSpec& spec_data = d_rod_spec_data[level_number][j].find(e)->second;
            rod_material_params.push_back(spec_data.properties);
        }
        if (!rod_next_idxs.empty())
        {
            node_data.push_back(new IBRodForceSpec(mastr_idx, rod_next_idxs, rod_material_params));
        }
    }

    // Initialize any target point specifications associated with the present
    // vertex.
    {
        const TargetSpec& spec_data = getVertexTargetSpec(point_index, level_number);
        const double kappa_target = spec_data.stiffness;
        const double eta_target = spec_data.damping;
        const Point& X_target = getVertexPosn(point_index, level_number);
        node_data.push_back(new IBTargetPointForceSpec(mastr_idx, kappa_target, eta_target, X_target));
    }

    // Initialize any anchor point specifications associated with the present
    // vertex.
    {
        const AnchorSpec& spec_data = getVertexAnchorSpec(point_index, level_number);
        const bool is_anchor_point = spec_data.is_anchor_point;
        if (is_anchor_point)
        {
            node_data.push_back(new IBAnchorPointSpec(mastr_idx));
        }
    }

    // Initialize any instrumentation specifications associated with the present
    // vertex.
    {
        const std::pair<int, int> inst_idx = getVertexInstrumentationIndices(point_index, level_number);
        if (inst_idx.first != -1 && inst_idx.second != -1)
        {
            node_data.push_back(new IBInstrumentationSpec(mastr_idx, inst_idx.first, inst_idx.second));
        }
    }

    // Initialize any source specifications associated with the present
    // vertex.
    {
        const int source_idx = getVertexSourceIndices(point_index, level_number);
        if (source_idx != -1)
        {
            node_data.push_back(new IBSourceSpec(mastr_idx, source_idx));
        }
    }

    return node_data;
} // initializeNodeData

void
IBRedundantInitializer::getFromInput(Pointer<Database> db)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(db);
#endif
    if (db->keyExists("max_levels"))
    {
        d_max_levels = db->getInteger("max_levels");
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Key data `max_levels' not found in input.");
    }

    if (d_max_levels < 1)
    {
        TBOX_ERROR(d_object_name << ":  "
                                 << "Key data `max_levels' found in input is < 1.");
    }

    d_level_is_initialized.resize(d_max_levels, false);
    d_base_filename.resize(d_max_levels);
    d_num_vertex.resize(d_max_levels);
    d_vertex_offset.resize(d_max_levels);
    d_vertex_posn.resize(d_max_levels);
    d_spring_edge_map.resize(d_max_levels);
    d_spring_spec_data.resize(d_max_levels);
    d_xspring_edge_map.resize(d_max_levels);
    d_xspring_spec_data.resize(d_max_levels);
    d_beam_spec_data.resize(d_max_levels);
    d_rod_edge_map.resize(d_max_levels);
    d_rod_spec_data.resize(d_max_levels);
    d_target_spec_data.resize(d_max_levels);
    d_anchor_spec_data.resize(d_max_levels);
    d_bdry_mass_spec_data.resize(d_max_levels);
    d_directors.resize(d_max_levels);
    d_instrument_idx.resize(d_max_levels);
    d_source_idx.resize(d_max_levels);

    d_global_index_offset.resize(d_max_levels);

    // Determine the various input file names.
    //
    // Prefer to use the new ``structure_names'' key, but revert to the
    // level-by-level ``base_filenames'' keys if necessary.
    if (db->keyExists("structure_names"))
    {
        const int num_strcts = db->getArraySize("structure_names");
        std::vector<std::string> structure_names(num_strcts);
        db->getStringArray("structure_names", &structure_names[0], num_strcts);
        for (int n = 0; n < num_strcts; ++n)
        {
            const std::string& strct_name = structure_names[n];
            if (db->keyExists(strct_name))
            {
                Pointer<Database> sub_db = db->getDatabase((strct_name));
                if (sub_db->keyExists("level_number"))
                {
                    const int ln = sub_db->getInteger("level_number");
                    if (ln < 0)
                    {
                        TBOX_ERROR(d_object_name << ":  "
                                                 << "Key data `level_number' associated with structure `" << strct_name
                                                 << "' is negative.");
                    }
                    else if (ln > d_max_levels)
                    {
                        TBOX_ERROR(d_object_name << ":  "
                                                 << "Key data `level_number' associated with structure `" << strct_name
                                                 << "' is greater than the expected maximum level number "
                                                 << d_max_levels << ".");
                    }
                    d_base_filename[ln].push_back(strct_name);
                }
                else
                {
                    TBOX_ERROR(d_object_name << ":  "
                                             << "Key data `level_number' not found in structure `" << strct_name
                                             << "' input.");
                }
            }
            else
            {
                TBOX_ERROR(d_object_name << ":  "
                                         << "Key data `" << strct_name << "' not found in input.");
            }
        }
    }
    else if (db->keyExists("structure_levels"))
    {
        const int num_levels = db->getArraySize("structure_levels");
        std::vector<int> strct_levels(num_levels);
        db->getIntegerArray("structure_levels", &strct_levels[0], num_levels);
        std::vector<int> num_strcts(num_levels);
        db->getIntegerArray("num_structures_levels", &num_strcts[0], num_levels);

        for (int k = 0; k < num_levels; ++k)
        {
            const int ln = strct_levels[k];
            const int num_strcts_on_ln = num_strcts[k];
            d_base_filename[ln].resize(num_strcts_on_ln);

            for (int s = 0; s < num_strcts_on_ln; ++s)
            {
                std::ostringstream strct_name_stream;
                strct_name_stream << "body_" << s << "_ln_" << ln;
                d_base_filename[ln][s] = strct_name_stream.str();
            }
        }
    }
    else
    {
        for (int ln = 0; ln < d_max_levels; ++ln)
        {
            std::ostringstream db_key_name_stream;
            db_key_name_stream << "base_filenames_" << ln;
            const std::string db_key_name = db_key_name_stream.str();
            if (db->keyExists(db_key_name))
            {
                const int num_files = db->getArraySize(db_key_name);
                d_base_filename[ln].resize(num_files);
                db->getStringArray(db_key_name, &d_base_filename[ln][0], num_files);
            }
            else
            {
                TBOX_WARNING(d_object_name << ":  "
                                           << "Key data `" << db_key_name << "' not found in input.");
            }
        }
    }

    return;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
