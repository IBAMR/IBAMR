// Filename: IBStandardInitializer.C
// Created on 22 Nov 2006 by Boyce Griffith
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

#include "IBStandardInitializer.h"

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
#include <ibamr/IBAnchorPointSpec.h>
#include <ibamr/IBBeamForceSpec.h>
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/IBRodForceSpec.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>
#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <fstream>
#include <iostream>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int NUM_ROD_PARAMS = 10;

inline std::string
discard_comments(
    const std::string& input_string)
{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    std::string output_string = input_string;
    std::istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '%');
    string_stream.clear();
    return output_string;
}// discard_comments
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardInitializer::IBStandardInitializer(
    const std::string& object_name,
    Pointer<Database> input_db)
    : d_object_name(object_name),
      d_use_file_batons(true),
      d_max_levels(-1),
      d_level_is_initialized(),
      d_silo_writer(NULL),
      d_base_filename(),
      d_length_scale_factor(1.0),
      d_posn_shift(NDIM,0.0),
      d_num_vertex(),
      d_vertex_offset(),
      d_vertex_posn(),
      d_enable_springs(),
      d_spring_edge_map(),
      d_spring_spec_data(),
      d_using_uniform_spring_stiffness(),
      d_uniform_spring_stiffness(),
      d_using_uniform_spring_rest_length(),
      d_uniform_spring_rest_length(),
      d_using_uniform_spring_force_fcn_idx(),
      d_uniform_spring_force_fcn_idx(),
      d_using_uniform_spring_subdomain_idx(),
      d_uniform_spring_subdomain_idx(),
      d_enable_beams(),
      d_beam_spec_data(),
      d_using_uniform_beam_bend_rigidity(),
      d_uniform_beam_bend_rigidity(),
      d_using_uniform_beam_curvature(),
      d_uniform_beam_curvature(),
      d_using_uniform_beam_subdomain_idx(),
      d_uniform_beam_subdomain_idx(),
      d_enable_rods(),
      d_rod_edge_map(),
      d_rod_spec_data(),
      d_using_uniform_rod_properties(),
      d_uniform_rod_properties(),
      d_using_uniform_rod_subdomain_idx(),
      d_uniform_rod_subdomain_idx(),
      d_enable_target_points(),
      d_target_spec_data(),
      d_using_uniform_target_stiffness(),
      d_uniform_target_stiffness(),
      d_using_uniform_target_damping(),
      d_uniform_target_damping(),
      d_using_uniform_target_subdomain_idx(),
      d_uniform_target_subdomain_idx(),
      d_enable_anchor_points(),
      d_anchor_spec_data(),
      d_using_uniform_anchor_subdomain_idx(),
      d_uniform_anchor_subdomain_idx(),
      d_enable_bdry_mass(),
      d_bdry_mass_spec_data(),
      d_using_uniform_bdry_mass(),
      d_uniform_bdry_mass(),
      d_using_uniform_bdry_mass_stiffness(),
      d_uniform_bdry_mass_stiffness(),
      d_directors(),
      d_enable_instrumentation(),
      d_instrument_idx(),
      d_global_index_offset()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif

    // Register the specification objects with the StreamableManager class.
    IBAnchorPointSpec::registerWithStreamableManager();
    IBBeamForceSpec::registerWithStreamableManager();
    IBInstrumentationSpec::registerWithStreamableManager();
    IBRodForceSpec::registerWithStreamableManager();
    IBSpringForceSpec::registerWithStreamableManager();
    IBTargetPointForceSpec::registerWithStreamableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Check to see if we are starting from a restart file.
    RestartManager* restart_manager = RestartManager::getManager();
    const bool is_from_restart = restart_manager->isFromRestart();

    // Process the input files only if we are not starting from a restart file.
    if (!is_from_restart)
    {
        // Process the vertex information.
        readVertexFiles();

        // Process the (optional) spring information.
        readSpringFiles();

        // Process the (optional) beam information.
        readBeamFiles();

        // Process the (optional) rod information.
        readRodFiles();

        // Process the (optional) target point information.
        readTargetPointFiles();

        // Process the (optional) anchor point information.
        readAnchorPointFiles();

        // Process the (optional) mass information.
        readBoundaryMassFiles();

        // Process the (optional) directors information.
        readDirectorFiles();

        // Process the (optional) instrumentation information.
        readInstrumentationFiles();

        // Wait for all processes to finish.
        SAMRAI_MPI::barrier();
    }
    return;
}// IBStandardInitializer

IBStandardInitializer::~IBStandardInitializer()
{
    // intentionally blank
    return;
}// ~IBStandardInitializer

void
IBStandardInitializer::registerLagSiloDataWriter(
    Pointer<LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
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
        for (int ln = 0; ln < d_max_levels; ++ln)
        {
            if (d_level_is_initialized[ln])
            {
                initializeLagSiloDataWriter(ln);
            }
        }
    }
    return;
}// registerLagSiloDataWriter

bool
IBStandardInitializer::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    return !d_num_vertex[level_number].empty();
}// getLevelHasLagrangianData

int
IBStandardInitializer::getLocalNodeCountOnPatchLevel(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Loop over all patches in the specified level of the patch level and count
    // the number of local vertices.
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of vertices whose initial locations will be within
        // the given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
    }
    return local_node_count;
}// getLocalNodeCountOnPatchLevel

void
IBStandardInitializer::initializeStructureIndexingOnPatchLevel(
    std::map<int,std::string>& strct_id_to_strct_name_map,
    std::map<int,std::pair<int,int> >& strct_id_to_lag_idx_range_map,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    (void) lag_manager;
    int offset = 0;
    for (int j = 0; j < int(d_base_filename[level_number].size()); ++j)
    {
        strct_id_to_strct_name_map   [j] = d_base_filename[level_number][j];
        strct_id_to_lag_idx_range_map[j] = std::make_pair(offset,offset+d_num_vertex[level_number][j]);
        offset += d_num_vertex[level_number][j];
    }
    return;
}// initializeStructureIndexingOnPatchLevel

int
IBStandardInitializer::initializeDataOnPatchLevel(
    const int lag_node_index_idx,
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LNodeLevelData>& X_data,
    Pointer<LNodeLevelData>& U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    (void) lag_manager;

    // Determine the extents of the physical domain.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const XLower = grid_geom->getXLower();
    const double* const XUpper = grid_geom->getXUpper();

    // Set the global index offset.  This is equal to the number of Lagrangian
    // indices that have already been initialized on the specified level.
    d_global_index_offset[level_number] = global_index_offset;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const CellIndex<NDIM>& patch_upper = patch_box.upper();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        Pointer<LNodeIndexData> index_data = patch->getPatchData(lag_node_index_idx);

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int>& point_idx = (*it);
            const int current_global_idx = getCanonicalLagrangianIndex(
                point_idx, level_number) + global_index_offset;
            const int current_local_idx = ++local_idx + local_index_offset;

            // Get the coordinates of the present vertex.
            const std::vector<double> X = getVertexPosn(point_idx, level_number);

            // Initialize the location of the present vertex.
            double* const node_X = &(*X_data)(current_local_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                node_X[d] = X[d];

                if (X[d] <= XLower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node below lower physical boundary\n"
                               << "  please ensure that all nodes are within the computational domain."<< std::endl);
                }

                if (X[d] >= XUpper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node above upper physical boundary\n"
                               << "  please ensure that all nodes are within the computational domain."<< std::endl);
                }
            }

            // Get the index of the cell in which the present vertex is
            // initially located.
            const CellIndex<NDIM> idx = IndexUtilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);

            // Initialize the specification objects associated with the present
            // vertex.
            std::vector<Pointer<Streamable> > specs = initializeSpecs(
                point_idx, global_index_offset, level_number);

            // Create or retrieve a pointer to the LNodeIndexSet associated with
            // the current Cartesian grid cell.
            if (!index_data->isElement(idx))
            {
                index_data->appendItemPointer(idx, new LNodeIndexSet());
            }
            LNodeIndexSet* const node_set = index_data->getItem(idx);
            static const IntVector<NDIM> periodic_offset(0);
            static const std::vector<double> periodic_displacement(NDIM,0.0);
            node_set->push_back(new LNodeIndex(current_global_idx, current_local_idx,
                                               &(*X_data)(current_local_idx),
                                               periodic_offset, periodic_displacement,
                                               specs));

            // Initialize the velocity of the present vertex.
            double* const node_U = &(*U_data)(current_local_idx);
            std::fill(node_U,node_U+NDIM,0.0);
        }
    }

    d_level_is_initialized[level_number] = true;

    // If a Lagrangian Silo data writer is registered with the initializer,
    // setup the visualization data corresponding to the present level of the
    // locally refined grid.
    if (!d_silo_writer.isNull())
    {
        initializeLagSiloDataWriter(level_number);
    }
    return local_node_count;
}// initializeDataOnPatchLevel

int
IBStandardInitializer::initializeMassDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LNodeLevelData>& M_data,
    Pointer<LNodeLevelData>& K_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    (void) lag_manager;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int>& point_idx = (*it);
            const int current_local_idx = ++local_idx + local_index_offset;

            // Initialize the mass and penalty stiffness coefficient
            // corresponding to the present vertex.
            const BdryMassSpec& spec = getVertexBdryMassSpec(point_idx, level_number);
            const double M = spec.bdry_mass;
            const double K = spec.stiffness;

            // Avoid division by zero at massless nodes.
            if (MathUtilities<double>::equalEps(M,0.0))
            {
                (*M_data)(current_local_idx) = std::numeric_limits<double>::epsilon();
                (*K_data)(current_local_idx) = 0.0;
            }
            else
            {
                (*M_data)(current_local_idx) = M;
                (*K_data)(current_local_idx) = K;
            }
        }
    }
    return local_node_count;
}// initializeMassOnPatchLevel

int
IBStandardInitializer::initializeDirectorDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    Pointer<LNodeLevelData>& D_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    (void) lag_manager;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    int local_idx = -1;
    int local_node_count = 0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int>& point_idx = (*it);
            const int current_local_idx = ++local_idx + local_index_offset;

            // Initialize the director corresponding to the present vertex.
            const std::vector<double>& D = getVertexDirectors(point_idx, level_number);
            for (int d = 0; d < 3*3; ++d)
            {
                (*D_data)(current_local_idx,d)= D[d];
            }
        }
    }
    return local_node_count;
}// initializeDirectorOnPatchLevel

void
IBStandardInitializer::tagCellsForInitialRefinement(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    // Loop over all patches in the specified level of the patch level and tag
    // cells for refinement wherever there are vertices assigned to a finer
    // level of the Cartesian grid.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const Box<NDIM>& patch_box = patch->getBox();
        const CellIndex<NDIM>& patch_lower = patch_box.lower();
        const CellIndex<NDIM>& patch_upper = patch_box.upper();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        Pointer<CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

        // Tag cells for refinement whenever there are vertices whose initial
        // locations will be within the index space of the given patch, but on
        // the finer levels of the AMR patch hierarchy.
        const bool can_be_refined = level_number+2 < d_max_levels;
        for (int ln = level_number+1; ln < d_max_levels; ++ln)
        {
            std::vector<std::pair<int,int> > patch_vertices;
            getPatchVertices(patch_vertices, patch, ln, can_be_refined);
            for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
                 it != patch_vertices.end(); ++it)
            {
                const std::pair<int,int>& point_idx = (*it);

                // Get the coordinates of the present vertex.
                const std::vector<double> X = getVertexPosn(point_idx, ln);

                // Get the index of the cell in which the present vertex is
                // initially located.
                const CellIndex<NDIM> i = IndexUtilities::getCellIndex(
                    X, xLower, xUpper, dx, patch_lower, patch_upper);

                // Tag the cell for refinement.
                if (patch_box.contains(i)) (*tag_data)(i) = 1;
            }
        }
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBStandardInitializer::initializeLagSiloDataWriter(
    const int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
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
        for (unsigned j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            const std::string postfix = "_vertices";
            d_silo_writer->registerMarkerCloud(
                d_base_filename[level_number][j] + postfix,
                d_num_vertex[level_number][j], d_vertex_offset[level_number][j], level_number);
        }

        bool registered_spring_edge_map = false;
        for (unsigned j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_spring_edge_map[level_number][j].size() > 0)
            {
                registered_spring_edge_map = true;
                const std::string postfix = "_mesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix,
                    d_spring_edge_map[level_number][j], level_number);
            }
        }

        for (unsigned j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            if (d_rod_edge_map[level_number][j].size() > 0)
            {
                const std::string postfix = registered_spring_edge_map ? "_rod_mesh" : "_mesh";
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + postfix,
                    d_rod_edge_map[level_number][j], level_number);
            }
        }
    }
    return;
}// initializeLagSiloDataWriter

void
IBStandardInitializer::readVertexFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_num_vertex[ln].resize(num_base_filename,std::numeric_limits<int>::max());
        d_vertex_offset[ln].resize(num_base_filename,std::numeric_limits<int>::max());
        d_vertex_posn[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            if (j == 0)
            {
                d_vertex_offset[ln][j] = 0;
            }
            else
            {
                d_vertex_offset[ln][j] = d_vertex_offset[ln][j-1]+d_num_vertex[ln][j-1];
            }

            // Ensure that the file exists.
            const std::string vertex_filename = d_base_filename[ln][j] + ".vertex";
            std::ifstream file_stream;
            file_stream.open(vertex_filename.c_str(), std::ios::in);
            if (!file_stream.is_open()) TBOX_ERROR(d_object_name << ":\n  Unable to open input file " << vertex_filename << std::endl);

            plog << d_object_name << ":  "
                 << "processing vertex data from ASCII input file named " << vertex_filename << std::endl
                 << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

            // The first entry in the file is the number of vertices.
            if (!std::getline(file_stream, line_string))
            {
                TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << vertex_filename << std::endl);
            }
            else
            {
                line_string = discard_comments(line_string);
                std::istringstream line_stream(line_string);
                if (!(line_stream >> d_num_vertex[ln][j]))
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << vertex_filename << std::endl);
                }
            }

            if (d_num_vertex[ln][j] <= 0)
            {
                TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << vertex_filename << std::endl);
            }

            // Each successive line provides the initial position of each vertex
            // in the input file.
            d_vertex_posn[ln][j].resize(d_num_vertex[ln][j]*NDIM);
            for (int k = 0; k < d_num_vertex[ln][j]; ++k)
            {
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << vertex_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (!(line_stream >> d_vertex_posn[ln][j][k*NDIM+d]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << vertex_filename << std::endl);
                        }
                        d_vertex_posn[ln][j][k*NDIM+d] = d_length_scale_factor*(d_vertex_posn[ln][j][k*NDIM+d] + d_posn_shift[d]);
                    }
                }
            }

            // Close the input file.
            file_stream.close();

            plog << d_object_name << ":  "
                 << "read " << d_num_vertex[ln][j] << " vertices from ASCII input file named " << vertex_filename << std::endl
                 << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
}// readVertexFiles

void
IBStandardInitializer::readSpringFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_spring_edge_map[ln].resize(num_base_filename);
        d_spring_spec_data[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            bool warned = false;

            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            // Ensure that the file exists.
            const std::string spring_filename = d_base_filename[ln][j] + ".spring";
            std::ifstream file_stream;
            file_stream.open(spring_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing spring data from ASCII input file named " << spring_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of edges in the input
                // file.
                int num_edges;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << spring_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_edges))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << spring_filename << std::endl);
                    }
                }

                if (num_edges <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << spring_filename << std::endl);
                }

                // Each successive line provides the connectivity and material parameter
                // information for each spring in the structure.
                for (int k = 0; k < num_edges; ++k)
                {
                    Edge e;
                    double kappa, length;
                    int force_fcn_idx, subdomain_idx;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << spring_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> e.first))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl);
                        }
                        else if ((e.first < 0) || (e.first >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl
                                       << "  vertex index " << e.first << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> e.second))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl);
                        }
                        else if ((e.second < 0) || (e.second >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl
                                       << "  vertex index " << e.second << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> kappa))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl);
                        }
                        else if (kappa < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl
                                       << "  spring constant is negative" << std::endl);
                        }

                        if (!(line_stream >> length))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl);
                        }
                        else if (length < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl
                                       << "  spring resting length is negative" << std::endl);
                        }
                        length *= d_length_scale_factor;

                        if (!(line_stream >> force_fcn_idx))
                        {
                            force_fcn_idx = 0;  // default force function specification.
                        }

                        if (!(line_stream >> subdomain_idx))
                        {
                            subdomain_idx = -1;  // default subdomain index.
                        }
                    }

                    // Modify kappa, length, and subdomain_idx according to
                    // whether uniform values are to be employed for this
                    // particular structure.
                    if (d_using_uniform_spring_stiffness[ln][j])
                    {
                        kappa = d_uniform_spring_stiffness[ln][j];
                    }
                    if (d_using_uniform_spring_rest_length[ln][j])
                    {
                        length = d_uniform_spring_rest_length[ln][j];
                    }
                    if (d_using_uniform_spring_force_fcn_idx[ln][j])
                    {
                        force_fcn_idx = d_uniform_spring_force_fcn_idx[ln][j];
                    }
                    if (d_using_uniform_spring_subdomain_idx[ln][j])
                    {
                        subdomain_idx = d_uniform_spring_subdomain_idx[ln][j];
                    }

                    // Correct the edge numbers to be in the global Lagrangian indexing
                    // scheme.
                    e.first  += d_vertex_offset[ln][j];
                    e.second += d_vertex_offset[ln][j];

                    // Always place the lower index first.
                    if (e.first > e.second)
                    {
                        std::swap<int>(e.first, e.second);
                    }

                    // Check to see if the edge has already been inserted in the edge map.
                    bool duplicate_edge = false;
                    for (std::multimap<int,Edge>::const_iterator it =
                             d_spring_edge_map[ln][j].lower_bound(e.first);
                         it != d_spring_edge_map[ln][j].upper_bound(e.first); ++it)
                    {
                        const Edge& other_e = it->second;
                        if (e.first  == other_e.first &&
                            e.second == other_e.second)
                        {
                            // This is a duplicate edge and should not be inserted into the
                            // edge map.
                            duplicate_edge = true;

                            // Ensure that the link information is consistent.
                            if (!MathUtilities<double>::equalEps(d_spring_spec_data[ln][j].find(e)->second.stiffness  , kappa ) ||
                                !MathUtilities<double>::equalEps(d_spring_spec_data[ln][j].find(e)->second.rest_length, length) ||
                                (d_spring_spec_data[ln][j].find(e)->second.force_fcn_idx != force_fcn_idx) ||
                                (d_spring_spec_data[ln][j].find(e)->second.subdomain_idx != subdomain_idx))
                            {
                                TBOX_ERROR(d_object_name << ":\n  Inconsistent duplicate edges in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl
                                           << "  first vertex = " << e.first-d_vertex_offset[ln][j] << " second vertex = " << e.second-d_vertex_offset[ln][j] << std::endl
                                           << "  original spring constant      = " << d_spring_spec_data[ln][j].find(e)->second.stiffness     << std::endl
                                           << "  original resting length       = " << d_spring_spec_data[ln][j].find(e)->second.rest_length   << std::endl
                                           << "  original force function index = " << d_spring_spec_data[ln][j].find(e)->second.force_fcn_idx << std::endl
                                           << "  original subdomain index      = " << d_spring_spec_data[ln][j].find(e)->second.subdomain_idx << std::endl);
                            }
                        }
                    }

                    // Initialize the map data corresponding to the present edge.
                    //
                    // Note that in the edge map, each edge is associated with only the
                    // first vertex.
                    if (!duplicate_edge)
                    {
                        d_spring_edge_map[ln][j].insert(std::make_pair(e.first,e));
                        SpringSpec& spec_data = d_spring_spec_data[ln][j][e];
                        spec_data.stiffness     = kappa;
                        spec_data.rest_length   = length;
                        spec_data.force_fcn_idx = force_fcn_idx;
                        spec_data.subdomain_idx = subdomain_idx;
                    }

                    // Check to see if the spring constant is zero and, if so,
                    // emit a warning.
                    if (!warned && d_enable_springs[ln][j] &&
                        (kappa == 0.0 || MathUtilities<double>::equalEps(kappa,0.0)))
                    {
                        TBOX_WARNING(d_object_name << ":\n  Spring with zero spring constant encountered in ASCII input file named " << spring_filename << "." << std::endl);
                        warned = true;
                    }
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_edges << " edges from ASCII input file named " << spring_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
}// readSpringFiles

void
IBStandardInitializer::readBeamFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_beam_spec_data[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            bool warned = false;

            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            const std::string beam_filename = d_base_filename[ln][j] + ".beam";
            std::ifstream file_stream;
            file_stream.open(beam_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing beam data from ASCII input file named " << beam_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of beams in
                // the input file.
                int num_beams;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << beam_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_beams))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << beam_filename << std::endl);
                    }
                }

                if (num_beams <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << beam_filename << std::endl);
                }

                // Each successive line provides the connectivity and material
                // parameter information for each beam in the structure.
                for (int k = 0; k < num_beams; ++k)
                {
                    int prev_idx, curr_idx, next_idx;
                    double bend;
                    std::vector<double> curv(NDIM,0.0);
                    int subdomain_idx;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << beam_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> prev_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl);
                        }
                        else if ((prev_idx < 0) || (prev_idx >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl
                                       << "  vertex index " << prev_idx << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> curr_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl);
                        }
                        else if ((curr_idx < 0) || (curr_idx >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl
                                       << "  vertex index " << curr_idx << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> next_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl);
                        }
                        else if ((next_idx < 0) || (next_idx >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl
                                       << "  vertex index " << next_idx << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> bend))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl);
                        }
                        else if (bend < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl
                                       << "  beam constant is negative" << std::endl);
                        }

                        bool curv_found_in_input = false;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            double c;
                            if (!(line_stream >> c))
                            {
                                if (curv_found_in_input)
                                {
                                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl
                                               << "  incomplete beam curvature specification" << std::endl);
                                }
                            }
                            else
                            {
                                curv_found_in_input = true;
                                curv[d] = c;
                            }
                        }

                        if (!(line_stream >> subdomain_idx))
                        {
                            subdomain_idx = -1;  // default subdomain index.
                        }
                    }

                    // Modify bend, curvature, and subdomain_idx according to
                    // whether uniform values are to be employed for this
                    // particular structure.
                    if (d_using_uniform_beam_bend_rigidity[ln][j])
                    {
                        bend = d_uniform_beam_bend_rigidity[ln][j];
                    }
                    if (d_using_uniform_beam_curvature[ln][j])
                    {
                        curv = d_uniform_beam_curvature[ln][j];
                    }
                    if (d_using_uniform_beam_subdomain_idx[ln][j])
                    {
                        subdomain_idx = d_uniform_beam_subdomain_idx[ln][j];
                    }

                    // Correct the node numbers to be in the global Lagrangian
                    // indexing scheme.
                    prev_idx += d_vertex_offset[ln][j];
                    curr_idx += d_vertex_offset[ln][j];
                    next_idx += d_vertex_offset[ln][j];

                    // Initialize the map data corresponding to the present
                    // beam.
                    //
                    // Note that in the beam property map, each edge is
                    // associated with only the "current" vertex.
                    BeamSpec spec_data;
                    spec_data.neighbor_idxs = std::make_pair(next_idx,prev_idx);
                    spec_data.bend_rigidity = bend;
                    spec_data.curvature     = curv;
                    spec_data.subdomain_idx = subdomain_idx;
                    d_beam_spec_data[ln][j].insert(std::make_pair(curr_idx,spec_data));

                    // Check to see if the bending rigidity is zero and, if so,
                    // emit a warning.
                    if (!warned && d_enable_beams[ln][j] &&
                        (bend == 0.0 || MathUtilities<double>::equalEps(bend,0.0)))
                    {
                        TBOX_WARNING(d_object_name << ":\n  Beam with zero bending rigidity encountered in ASCII input file named " << beam_filename << "." << std::endl);
                        warned = true;
                    }
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_beams << " beams from ASCII input file named " << beam_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
}// readBeamFiles

void
IBStandardInitializer::readRodFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_rod_edge_map[ln].resize(num_base_filename);
        d_rod_spec_data[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            const std::string rod_filename = d_base_filename[ln][j] + ".rod";
            std::ifstream file_stream;
            file_stream.open(rod_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing rod data from ASCII input file named " << rod_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of rods in
                // the input file.
                int num_rods;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << rod_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_rods))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << rod_filename << std::endl);
                    }
                }

                if (num_rods <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << rod_filename << std::endl);
                }

                // Each successive line provides the connectivity and material
                // parameter information for each rod in the structure.
                for (int k = 0; k < num_rods; ++k)
                {
                    int curr_idx, next_idx;
                    std::vector<double> properties(NUM_ROD_PARAMS);
                    double& ds = properties[0];
                    double& a1 = properties[1];
                    double& a2 = properties[2];
                    double& a3 = properties[3];
                    double& b1 = properties[4];
                    double& b2 = properties[5];
                    double& b3 = properties[6];
                    double& kappa1 = properties[7];
                    double& kappa2 = properties[8];
                    double& tau = properties[9];
                    int subdomain_idx;

                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << rod_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);

                        if (!(line_stream >> curr_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if ((curr_idx < 0) || (curr_idx >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  vertex index " << curr_idx << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> next_idx))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if ((next_idx < 0) || (next_idx >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  vertex index " << next_idx << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> ds))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (ds < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant ds is negative" << std::endl);
                        }

                        if (!(line_stream >> a1))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (a1 < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant a1 is negative" << std::endl);
                        }

                        if (!(line_stream >> a2))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (a2 < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant a2 is negative" << std::endl);
                        }

                        if (!(line_stream >> a3))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (a3 < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant a3 is negative" << std::endl);
                        }

                        if (!(line_stream >> b1))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (b1 < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant b1 is negative" << std::endl);
                        }

                        if (!(line_stream >> b2))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (b2 < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant b2 is negative" << std::endl);
                        }

                        if (!(line_stream >> b3))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl);
                        }
                        else if (b3 < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                       << "  rod material constant b3 is negative" << std::endl);
                        }

                        bool curvature_data_found_in_input = false;

                        if (!(line_stream >> kappa1))
                        {
                            kappa1 = 0.0;
                        }
                        else
                        {
                            curvature_data_found_in_input = true;
                        }

                        if (!(line_stream >> kappa2))
                        {
                            kappa2 = 0.0;
                            if (curvature_data_found_in_input)
                            {
                                TBOX_WARNING(d_object_name << ":\n  Potentially invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                             << "  intrinsic curvature kappa1 was specified but kappa2 was not" << std::endl);
                            }
                        }
                        else
                        {
                            curvature_data_found_in_input = true;
                        }

                        if (!(line_stream >> tau))
                        {
                            tau = 0.0;
                            if (curvature_data_found_in_input)
                            {
                                TBOX_WARNING(d_object_name << ":\n  Potentially invalid entry in input file encountered on line " << k+2 << " of file " << rod_filename << std::endl
                                             << "  intrinsic curvatures kappa1 and kappa2 were specified but intrinsic twist tau was not" << std::endl);
                            }
                        }
                        else
                        {
                            curvature_data_found_in_input = true;
                        }

                        if (!(line_stream >> subdomain_idx))
                        {
                            subdomain_idx = -1;  // default subdomain index.
                        }
                    }

                    // Modify properties and subdomain_idx according to whether
                    // uniform values are to be employed for this particular
                    // structure.
                    if (d_using_uniform_rod_properties[ln][j])
                    {
                        properties = d_uniform_rod_properties[ln][j];
                    }
                    if (d_using_uniform_rod_subdomain_idx[ln][j])
                    {
                        subdomain_idx = d_uniform_rod_subdomain_idx[ln][j];
                    }

                    // Correct the node numbers to be in the global Lagrangian
                    // indexing scheme.
                    curr_idx += d_vertex_offset[ln][j];
                    next_idx += d_vertex_offset[ln][j];

                    // Initialize the map data corresponding to the present
                    // edge.
                    Edge e;
                    e.first  = curr_idx;
                    e.second = next_idx;
                    d_rod_edge_map[ln][j].insert(std::make_pair(e.first,e));

                    // Initialize the map data corresponding to the present rod.
                    //
                    // Note that in the rod property map, each edge is
                    // associated with only the "cur1rent" vertex.
                    RodSpec& rod_spec = d_rod_spec_data[ln][j][e];
                    rod_spec.properties = properties;
                    rod_spec.subdomain_idx = subdomain_idx;
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_rods << " rods from ASCII input file named " << rod_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
}// readRodFiles

void
IBStandardInitializer::readTargetPointFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_target_spec_data[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            bool warned = false;

            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            TargetSpec default_spec;
            default_spec.stiffness = 0.0;
            default_spec.damping = 0.0;
            default_spec.subdomain_idx = -1;
            d_target_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);

            const std::string target_point_stiffness_filename = d_base_filename[ln][j] + ".target";
            std::ifstream file_stream;
            file_stream.open(target_point_stiffness_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing target point data from ASCII input file named " << target_point_stiffness_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of target
                // point specifications in the input file.
                int num_target_points;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << target_point_stiffness_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_target_points))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << target_point_stiffness_filename << std::endl);
                    }
                }

                if (num_target_points <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << target_point_stiffness_filename << std::endl);
                }

                // Each successive line indicates the vertex number and spring
                // constant associated with any target points.
                for (int k = 0; k < num_target_points; ++k)
                {
                    int n;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << target_point_stiffness_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> n))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl);
                        }
                        else if ((n < 0) || (n >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl
                                       << "  vertex index " << n << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> d_target_spec_data[ln][j][n].stiffness))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl);
                        }
                        else if (d_target_spec_data[ln][j][n].stiffness < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl
                                       << "  target point spring constant is negative" << std::endl);
                        }

                        if (!(line_stream >> d_target_spec_data[ln][j][n].damping))
                        {
                            d_target_spec_data[ln][j][n].damping = 0.0;
                        }
                        else if (d_target_spec_data[ln][j][n].damping < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl
                                       << "  target point damping coefficient is negative" << std::endl);
                        }

                        if (!(line_stream >> d_target_spec_data[ln][j][n].subdomain_idx))
                        {
                            d_target_spec_data[ln][j][n].subdomain_idx = -1;  // default subdomain index.
                        }
                    }

                    // Check to see if the penalty spring constant is zero and,
                    // if so, emit a warning.
                    const double kappa = d_target_spec_data[ln][j][n].stiffness;
                    if (!warned && d_enable_target_points[ln][j] &&
                        (kappa == 0.0 || MathUtilities<double>::equalEps(kappa,0.0)))
                    {
                        TBOX_WARNING(d_object_name << ":\n  Target point with zero penalty spring constant encountered in ASCII input file named " << target_point_stiffness_filename << "." << std::endl);
                        warned = true;
                    }
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_target_points << " target points from ASCII input file named " << target_point_stiffness_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Modify the target point stiffness constants according to whether
            // target point penalty forces are enabled, or whether uniform
            // values are to be employed, for this particular structure.
            if (!d_enable_target_points[ln][j])
            {
                for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                {
                    d_target_spec_data[ln][j][k].stiffness = 0.0;
                    d_target_spec_data[ln][j][k].damping   = 0.0;
                }
            }
            else
            {
                if (d_using_uniform_target_stiffness[ln][j])
                {
                    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                    {
                        d_target_spec_data[ln][j][k].stiffness = d_uniform_target_stiffness[ln][j];
                    }
                }
                if (d_using_uniform_target_damping[ln][j])
                {
                    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                    {
                        d_target_spec_data[ln][j][k].damping = d_uniform_target_damping[ln][j];
                    }
                }
                if (d_using_uniform_target_subdomain_idx[ln][j])
                {
                    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                    {
                        d_target_spec_data[ln][j][k].subdomain_idx = d_uniform_target_subdomain_idx[ln][j];
                    }
                }
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI_MPI::barrier();
    return;
}// readTargetPointFiles

void
IBStandardInitializer::readAnchorPointFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_anchor_spec_data[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            AnchorSpec default_spec;
            default_spec.is_anchor_point = false;
            default_spec.subdomain_idx = -1;
            d_anchor_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);

            const std::string anchor_point_filename = d_base_filename[ln][j] + ".anchor";
            std::ifstream file_stream;
            file_stream.open(anchor_point_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing anchor point data from ASCII input file named " << anchor_point_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of anchor
                // points in the input file.
                int num_anchor_pts;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << anchor_point_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_anchor_pts))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << anchor_point_filename << std::endl);
                    }
                }

                if (num_anchor_pts <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << anchor_point_filename << std::endl);
                }

                // Each successive line indicates the vertex number of the
                // anchor points in the input file.
                for (int k = 0; k < num_anchor_pts; ++k)
                {
                    int n;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << anchor_point_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> n))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << anchor_point_filename << std::endl);
                        }
                        else if ((n < 0) || (n >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << anchor_point_filename << std::endl
                                       << "  vertex index " << n << " is out of range" << std::endl);
                        }

                        d_anchor_spec_data[ln][j][n].is_anchor_point = true;

                        if (!(line_stream >> d_anchor_spec_data[ln][j][n].subdomain_idx))
                        {
                            d_anchor_spec_data[ln][j][n].subdomain_idx = -1;  // default subdomain index.
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_anchor_pts << " anchor points from ASCII input file named " << anchor_point_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Modify the anchor point properties according to whether uniform
            // values are to be employed for this particular structure.
            if (d_using_uniform_anchor_subdomain_idx[ln][j])
            {
                for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                {
                    d_anchor_spec_data[ln][j][k].subdomain_idx = d_uniform_anchor_subdomain_idx[ln][j];
                }
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }
    return;
}// readAnchorPointFiles

void
IBStandardInitializer::readBoundaryMassFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_bdry_mass_spec_data[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            BdryMassSpec default_spec;
            default_spec.bdry_mass = 0.0;
            default_spec.stiffness = 0.0;
            d_bdry_mass_spec_data[ln][j].resize(d_num_vertex[ln][j], default_spec);

            const std::string bdry_mass_filename = d_base_filename[ln][j] + ".mass";
            std::ifstream file_stream;
            file_stream.open(bdry_mass_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing boundary mass data from ASCII input file named " << bdry_mass_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of massive IB
                // points in the input file.
                int num_bdry_mass_pts;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << bdry_mass_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_bdry_mass_pts))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << bdry_mass_filename << std::endl);
                    }
                }

                if (num_bdry_mass_pts <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << bdry_mass_filename << std::endl);
                }

                // Each successive line indicates the vertex number, mass, and
                // penalty spring constant associated with any massive IB
                // points.
                for (int k = 0; k < num_bdry_mass_pts; ++k)
                {
                    int n;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << bdry_mass_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> n))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl);
                        }
                        else if ((n < 0) || (n >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl
                                       << "  vertex index " << n << " is out of range" << std::endl);
                        }

                        if (!(line_stream >> d_bdry_mass_spec_data[ln][j][n].bdry_mass))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl);
                        }
                        else if (d_bdry_mass_spec_data[ln][j][n].bdry_mass < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl
                                       << "  boundary mass is negative" << std::endl);
                        }

                        if (!(line_stream >> d_bdry_mass_spec_data[ln][j][n].stiffness))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl);
                        }
                        else if (d_bdry_mass_spec_data[ln][j][n].stiffness < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl
                                       << "  boundary mass spring constant is negative" << std::endl);
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_bdry_mass_pts << " boundary mass points from ASCII input file named " << bdry_mass_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Modify the boundary masses and boundary mass stiffness constants
            // according to whether boundary mass is enabled, or whether uniform
            // values are to be employed, for this particular structure.
            if (!d_enable_bdry_mass[ln][j])
            {
                for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                {
                    d_bdry_mass_spec_data[ln][j][k].bdry_mass = 0.0;
                    d_bdry_mass_spec_data[ln][j][k].stiffness = 0.0;
                }
            }
            else
            {
                if (d_using_uniform_bdry_mass[ln][j])
                {
                    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                    {
                        d_bdry_mass_spec_data[ln][j][k].bdry_mass = d_uniform_bdry_mass[ln][j];
                    }
                }
                if (d_using_uniform_bdry_mass_stiffness[ln][j])
                {
                    for (int k = 0; k < d_num_vertex[ln][j]; ++k)
                    {
                        d_bdry_mass_spec_data[ln][j][k].stiffness = d_uniform_bdry_mass_stiffness[ln][j];
                    }
                }
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }
    return;
}// readBoundaryMassFiles

void
IBStandardInitializer::readDirectorFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_directors[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            d_directors[ln][j].resize(d_num_vertex[ln][j], std::vector<double>(3*3,0.0));

            const std::string directors_filename = d_base_filename[ln][j] + ".director";
            std::ifstream file_stream;
            file_stream.open(directors_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                plog << d_object_name << ":  "
                     << "processing director data from ASCII input file named " << directors_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of sets of
                // directors in the input file.
                int num_directors_pts;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << directors_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_directors_pts))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << directors_filename << std::endl);
                    }
                }

                if (num_directors_pts != d_num_vertex[ln][j])
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << directors_filename << std::endl);
                }

                // Each successive set of three lines indicates the initial
                // configuration of a triad.
                for (int k = 0; k < num_directors_pts; ++k)
                {
                    for (int n = 0; n < 3; ++n)
                    {
                        if (!std::getline(file_stream, line_string))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << 3*k+n+2 << " of file " << directors_filename << std::endl);
                        }
                        else
                        {
                            line_string = discard_comments(line_string);
                            std::istringstream line_stream(line_string);
                            double D_norm_squared = 0.0;
                            for (int d = 0; d < 3; ++d)
                            {
                                if (!(line_stream >> d_directors[ln][j][k][3*n+d]))
                                {
                                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << 3*k+n+2 << " of file " << directors_filename << std::endl);
                                }
                                D_norm_squared += d_directors[ln][j][k][3*n+d]*d_directors[ln][j][k][3*n+d];
                            }
                            const double D_norm = sqrt(D_norm_squared);
                            if (!MathUtilities<double>::equalEps(D_norm,1.0))
                            {
                                TBOX_WARNING(d_object_name << ":\n  Director vector on line " << 3*k+n+2 << " of file " << directors_filename << " is not normalized; norm = " << D_norm << std::endl);
                                for (int d = 0; d < 3; ++d)
                                {
                                    d_directors[ln][j][k][3*n+d] /= D_norm;
                                }
                            }
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_directors_pts << " director triads from ASCII input file named " << directors_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }
    return;
}// readDirectorFiles

void
IBStandardInitializer::readInstrumentationFiles()
{
    std::string line_string;
    const int rank = SAMRAI_MPI::getRank();
    const int nodes = SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    int instrument_offset = 0;
    std::vector<std::string> instrument_names;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_instrument_idx[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            const std::string inst_filename = d_base_filename[ln][j] + ".inst";
            std::ifstream file_stream;
            file_stream.open(inst_filename.c_str(), std::ios::in);
            if (file_stream.is_open() && d_enable_instrumentation[ln][j])
            {
                plog << d_object_name << ":  "
                     << "processing instrumentation data from ASCII input file named " << inst_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of
                // instruments in the input file.
                int num_inst;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << inst_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_inst))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << inst_filename << std::endl);
                    }
                }

                if (num_inst <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << inst_filename << std::endl);
                }

                // The next several lines in the file indicate the names of the
                // instruments in the input file.
                for (int m = 0; m < num_inst; ++m)
                {
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << m+2 << " of file " << inst_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);

                        // trim leading whitespace
                        std::string::size_type notwhite = line_string.find_first_not_of(" \t\n");
                        line_string.erase(0,notwhite);

                        // trim trailing whitespace
                        notwhite = line_string.find_last_not_of(" \t\n");
                        line_string.erase(notwhite+1);

                        instrument_names.push_back(line_string);
                    }
                }

                // The next line in the file indicates the number of
                // instrumented IB points in the input file.
                int num_inst_pts;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << num_inst+2 << " of file " << inst_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_inst_pts))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+2 << " of file " << inst_filename << std::endl);
                    }
                }

                if (num_inst_pts <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+2 << " of file " << inst_filename << std::endl);
                }

                // Each successive line indicates the vertex number, meter
                // number, and meter node indices of each of the instrumented IB
                // points in the input file.
                std::vector<bool> encountered_instrument_idx;
                std::map<int,std::vector<bool> > encountered_node_idx;
                for (int k = 0; k < num_inst_pts; ++k)
                {
                    int n;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << num_inst+k+3 << " of file " << inst_filename << std::endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> n))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << std::endl);
                        }
                        else if ((n < 0) || (n >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << std::endl
                                       << "  vertex index " << n << " is out of range" << std::endl);
                        }

                        std::pair<int,int>& idx = d_instrument_idx[ln][j][n];

                        if (!(line_stream >> idx.first))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << std::endl);
                        }
                        else if (idx.first < 0 || idx.first >= num_inst)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << std::endl
                                       << "  meter index " << idx.first << " is out of range" << std::endl);
                        }

                        if (idx.first >= int(encountered_instrument_idx.size()))
                        {
                            encountered_instrument_idx.resize(idx.first+1,false);
                        }
                        encountered_instrument_idx[idx.first] = true;

                        if (!(line_stream >> idx.second))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << std::endl);
                        }
                        else if (idx.second < 0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << std::endl
                                       << "  meter node index is negative" << std::endl);
                        }

                        if (idx.second >= int(encountered_node_idx[idx.first].size()))
                        {
                            encountered_node_idx[idx.first].resize(idx.second+1,false);
                        }
                        encountered_node_idx[idx.first][idx.second] = true;

                        // Correct the instrument index to account for
                        // instrument indices from earlier files.
                        idx.first += instrument_offset;
                    }
                }

                // Ensure that a complete range of instrument indices were found
                // in the input file.
                for (std::vector<bool>::iterator meter_it = encountered_instrument_idx.begin();
                     meter_it != encountered_instrument_idx.end(); ++meter_it)
                {
                    const int meter_idx = std::distance(encountered_instrument_idx.begin(),meter_it);
                    if ((*meter_it) == false)
                    {
                        TBOX_ERROR(d_object_name << ":\n  "
                                   << "  Instrument index " << meter_idx << " not found in input file " << inst_filename << std::endl);
                    }

                    std::vector<bool>& meter_node_idxs = encountered_node_idx[meter_idx];
                    for (std::vector<bool>::iterator node_it = meter_node_idxs.begin();
                         node_it != meter_node_idxs.end(); ++node_it)
                    {
                        const int node_idx = std::distance(meter_node_idxs.begin(),node_it);
                        if ((*node_it) == false)
                        {
                            TBOX_ERROR(d_object_name << ":\n  "
                                       << "  Node index " << node_idx << " associated with meter index " << meter_idx << " not found in input file " << inst_filename << std::endl);
                        }
                    }
                }

                if (int(encountered_instrument_idx.size()) != num_inst)
                {
                    TBOX_ERROR(d_object_name << ":\n  "
                               << "  Not all anticipated instrument indices were found in input file " << inst_filename
                               << "  Expected to find " << num_inst << " distinct meter indices in input file" << std::endl);
                }

                // Increment the meter offset.
                instrument_offset += encountered_instrument_idx.size();

                // Close the input file.
                file_stream.close();

                plog << d_object_name << ":  "
                     << "read " << num_inst_pts << " instrumentation points from ASCII input file named " << inst_filename << std::endl
                     << "  on MPI process " << SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }
    IBInstrumentationSpec::setInstrumentNames(instrument_names);
    return;
}// readInstrumentationFiles

void
IBStandardInitializer::getPatchVertices(
    std::vector<std::pair<int,int> >& patch_vertices,
    const Pointer<Patch<NDIM> > patch,
    const int level_number,
    const bool can_be_refined) const
{
    // Loop over all of the vertices to determine the indices of those vertices
    // within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will work for
    // now.
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
        patch->getPatchGeometry();
    const double* const xLower = patch_geom->getXLower();
    const double* const xUpper = patch_geom->getXUpper();

    for (unsigned j = 0; j < d_num_vertex[level_number].size(); ++j)
    {
        for (int k = 0; k < d_num_vertex[level_number][j]; ++k)
        {
            const double* const X = &d_vertex_posn[level_number][j][k*NDIM];
            const bool patch_owns_node =
                ((  xLower[0] <= X[0])&&(X[0] < xUpper[0]))
#if (NDIM > 1)
                &&((xLower[1] <= X[1])&&(X[1] < xUpper[1]))
#if (NDIM > 2)
                &&((xLower[2] <= X[2])&&(X[2] < xUpper[2]))
#endif
#endif
                ;
            if (patch_owns_node) patch_vertices.push_back(std::make_pair(j,k));
        }
    }
    return;
}// getPatchVertices

int
IBStandardInitializer::getCanonicalLagrangianIndex(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_vertex_offset[level_number][point_index.first]+point_index.second;
}// getCanonicalLagrangianIndex

std::vector<double>
IBStandardInitializer::getVertexPosn(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    const double* posn_ptr = &d_vertex_posn[level_number][point_index.first][point_index.second*NDIM];
    return std::vector<double>(posn_ptr,posn_ptr+NDIM);
}// getVertexPosn

const IBStandardInitializer::TargetSpec&
IBStandardInitializer::getVertexTargetSpec(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_target_spec_data[level_number][point_index.first][point_index.second];
}// getVertexTargetSpec

const IBStandardInitializer::AnchorSpec&
IBStandardInitializer::getVertexAnchorSpec(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_anchor_spec_data[level_number][point_index.first][point_index.second];
}// getVertexAnchorSpec

const IBStandardInitializer::BdryMassSpec&
IBStandardInitializer::getVertexBdryMassSpec(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_bdry_mass_spec_data[level_number][point_index.first][point_index.second];
}// getVertexBdryMassSpec

const std::vector<double>&
IBStandardInitializer::getVertexDirectors(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_directors[level_number][point_index.first][point_index.second];
}// getVertexDirectors

std::pair<int,int>
IBStandardInitializer::getVertexInstrumentationIndices(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    std::map<int,std::pair<int,int> >::const_iterator it =
        d_instrument_idx[level_number][point_index.first].find(point_index.second);
    if (it != d_instrument_idx[level_number][point_index.first].end())
    {
        return it->second;
    }
    else
    {
        return std::make_pair(-1,-1);
    }
}// getVertexInstrumentationIndices

std::vector<Pointer<Streamable> >
IBStandardInitializer::initializeSpecs(
    const std::pair<int,int>& point_index,
    const int global_index_offset,
    const int level_number) const
{
    std::vector<Pointer<Streamable> > vertex_specs;

    const int j = point_index.first;
    const int mastr_idx = getCanonicalLagrangianIndex(point_index, level_number);

    // Initialize any spring specifications associated with the present vertex.
    if (d_enable_springs[level_number][j])
    {
        std::vector<int> slave_idxs, force_fcn_idxs;
        std::vector<double> stiffness, rest_length;
        std::vector<int> subdomain_idxs;
        for (std::multimap<int,Edge>::const_iterator it = d_spring_edge_map[level_number][j].lower_bound(mastr_idx);
             it != d_spring_edge_map[level_number][j].upper_bound(mastr_idx); ++it)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(mastr_idx == it->first);
#endif
            // The connectivity information.
            const Edge& e = it->second;
            if (e.first == mastr_idx)
            {
                slave_idxs.push_back(e.second+global_index_offset);
            }
            else
            {
                slave_idxs.push_back(e.first +global_index_offset);
            }

            // The material properties.
            const SpringSpec& spec_data = d_spring_spec_data[level_number][j].find(e)->second;
            stiffness     .push_back(spec_data.stiffness    );
            rest_length   .push_back(spec_data.rest_length  );
            force_fcn_idxs.push_back(spec_data.force_fcn_idx);
            subdomain_idxs.push_back(spec_data.subdomain_idx);
        }
        if (slave_idxs.size() > 0)
        {
            vertex_specs.push_back(
                new IBSpringForceSpec(
                    mastr_idx, slave_idxs, force_fcn_idxs, stiffness, rest_length, subdomain_idxs));
        }
    }

    // Initialize any beam specifications associated with the present vertex.
    if (d_enable_beams[level_number][j])
    {
        std::vector<std::pair<int,int> > beam_neighbor_idxs;
        std::vector<double> beam_bend_rigidity;
        std::vector<std::vector<double> > beam_mesh_dependent_curvature;
        std::vector<int> beam_subdomain_idxs;
        for (std::multimap<int,BeamSpec>::const_iterator it = d_beam_spec_data[level_number][j].lower_bound(mastr_idx);
             it != d_beam_spec_data[level_number][j].upper_bound(mastr_idx); ++it)
        {
            const BeamSpec& spec_data = it->second;
            beam_neighbor_idxs.push_back(spec_data.neighbor_idxs);
            beam_bend_rigidity.push_back(spec_data.bend_rigidity);
            beam_mesh_dependent_curvature.push_back(spec_data.curvature);
            beam_subdomain_idxs.push_back(spec_data.subdomain_idx);
        }
        if (!beam_neighbor_idxs.empty())
        {
            vertex_specs.push_back(
                new IBBeamForceSpec(
                    mastr_idx, beam_neighbor_idxs, beam_bend_rigidity, beam_mesh_dependent_curvature, beam_subdomain_idxs));
        }
    }

    // Initialize any rod specifications associated with the present vertex.
    if (d_enable_rods[level_number][j])
    {
        std::vector<int> rod_next_idxs;
        std::vector<std::vector<double> > rod_material_params;
        std::vector<int> rod_subdomain_idxs;
        for (std::multimap<int,Edge>::const_iterator it = d_rod_edge_map[level_number][j].lower_bound(mastr_idx);
             it != d_rod_edge_map[level_number][j].upper_bound(mastr_idx); ++it)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(mastr_idx == it->first);
#endif
            // The connectivity information.
            const Edge& e = it->second;
            if (e.first == mastr_idx)
            {
                rod_next_idxs.push_back(e.second+global_index_offset);
            }
            else
            {
                rod_next_idxs.push_back(e.first +global_index_offset);
            }

            // The material properties.
            const RodSpec& spec_data = d_rod_spec_data[level_number][j].find(e)->second;
            rod_material_params.push_back(spec_data.properties);
            rod_subdomain_idxs.push_back(spec_data.subdomain_idx);
        }
        if (!rod_next_idxs.empty())
        {
            vertex_specs.push_back(
                new IBRodForceSpec(
                    mastr_idx, rod_next_idxs, rod_material_params, rod_subdomain_idxs));
        }
    }

    // Initialize any target point specifications associated with the present
    // vertex.
    if (d_enable_target_points[level_number][j])
    {
        const TargetSpec& spec_data = getVertexTargetSpec(point_index, level_number);
        const double kappa_target = spec_data.stiffness;
        const double eta_target = spec_data.damping;
        const int subdomain_idx = spec_data.subdomain_idx;
        const std::vector<double> X_target = getVertexPosn(point_index, level_number);
        vertex_specs.push_back(
            new IBTargetPointForceSpec(
                mastr_idx, kappa_target, eta_target, X_target, subdomain_idx));
    }

    // Initialize any anchor point specifications associated with the present
    // vertex.
    if (d_enable_anchor_points[level_number][j])
    {
        const AnchorSpec& spec_data = getVertexAnchorSpec(point_index, level_number);
        const bool is_anchor_point = spec_data.is_anchor_point;
        const int subdomain_idx = spec_data.subdomain_idx;
        if (is_anchor_point)
        {
            vertex_specs.push_back(new IBAnchorPointSpec(mastr_idx, subdomain_idx));
        }
    }

    // Initialize any instrumentation specifications associated with the present
    // vertex.
    if (d_enable_instrumentation[level_number][j])
    {
        const std::pair<int,int> inst_idx = getVertexInstrumentationIndices(point_index, level_number);
        if (inst_idx.first != -1 && inst_idx.second != -1)
        {
            vertex_specs.push_back(
                new IBInstrumentationSpec(
                    mastr_idx, inst_idx.first, inst_idx.second));
        }
    }
    return vertex_specs;
}// initializeSpecs

void
IBStandardInitializer::getFromInput(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif

    // Determine whether to use "batons" to prevent multiple MPI processes from
    // reading the same file at once.
    d_use_file_batons = db->getBoolWithDefault("use_file_batons",d_use_file_batons);

    // Determine the (maximum) number of levels in the locally refined grid.
    // Note that each piece of the Lagrangian structure must be assigned to a
    // particular level of the grid.
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

    // Resize the vectors that are indexed by the level number.
    d_level_is_initialized.resize(d_max_levels,false);

    d_base_filename.resize(d_max_levels);

    d_num_vertex.resize(d_max_levels);
    d_vertex_offset.resize(d_max_levels);
    d_vertex_posn.resize(d_max_levels);

    d_enable_springs.resize(d_max_levels);
    d_spring_edge_map.resize(d_max_levels);
    d_spring_spec_data.resize(d_max_levels);
    d_using_uniform_spring_stiffness.resize(d_max_levels);
    d_uniform_spring_stiffness.resize(d_max_levels);
    d_using_uniform_spring_rest_length.resize(d_max_levels);
    d_uniform_spring_rest_length.resize(d_max_levels);
    d_using_uniform_spring_force_fcn_idx.resize(d_max_levels);
    d_uniform_spring_force_fcn_idx.resize(d_max_levels);
    d_using_uniform_spring_subdomain_idx.resize(d_max_levels);
    d_uniform_spring_subdomain_idx.resize(d_max_levels);

    d_enable_beams.resize(d_max_levels);
    d_beam_spec_data.resize(d_max_levels);
    d_using_uniform_beam_bend_rigidity.resize(d_max_levels);
    d_uniform_beam_bend_rigidity.resize(d_max_levels);
    d_using_uniform_beam_curvature.resize(d_max_levels);
    d_uniform_beam_curvature.resize(d_max_levels);
    d_using_uniform_beam_subdomain_idx.resize(d_max_levels);
    d_uniform_beam_subdomain_idx.resize(d_max_levels);

    d_enable_rods.resize(d_max_levels);
    d_rod_edge_map.resize(d_max_levels);
    d_rod_spec_data.resize(d_max_levels);
    d_using_uniform_rod_properties.resize(d_max_levels);
    d_uniform_rod_properties.resize(d_max_levels);
    d_using_uniform_rod_subdomain_idx.resize(d_max_levels);
    d_uniform_rod_subdomain_idx.resize(d_max_levels);

    d_enable_target_points.resize(d_max_levels);
    d_target_spec_data.resize(d_max_levels);
    d_using_uniform_target_stiffness.resize(d_max_levels);
    d_uniform_target_stiffness.resize(d_max_levels);
    d_using_uniform_target_damping.resize(d_max_levels);
    d_uniform_target_damping.resize(d_max_levels);
    d_using_uniform_target_subdomain_idx.resize(d_max_levels);
    d_uniform_target_subdomain_idx.resize(d_max_levels);

    d_enable_anchor_points.resize(d_max_levels);
    d_anchor_spec_data.resize(d_max_levels);
    d_using_uniform_anchor_subdomain_idx.resize(d_max_levels);
    d_uniform_anchor_subdomain_idx.resize(d_max_levels);

    d_enable_bdry_mass.resize(d_max_levels);
    d_bdry_mass_spec_data.resize(d_max_levels);
    d_using_uniform_bdry_mass.resize(d_max_levels);
    d_uniform_bdry_mass.resize(d_max_levels);
    d_using_uniform_bdry_mass_stiffness.resize(d_max_levels);
    d_uniform_bdry_mass_stiffness.resize(d_max_levels);

    d_directors.resize(d_max_levels);

    d_enable_instrumentation.resize(d_max_levels);
    d_instrument_idx.resize(d_max_levels);

    d_global_index_offset.resize(d_max_levels);

    // Determine the various input file names.
    //
    // Prefer to use the new ``structure_names'' key, but revert to the
    // level-by-level ``base_filenames'' keys if necessary.
    if (db->keyExists("structure_names"))
    {
        const int n_strcts = db->getArraySize("structure_names");
        std::vector<std::string> structure_names(n_strcts);
        db->getStringArray("structure_names", &structure_names[0], n_strcts);
        for (int n = 0; n < n_strcts; ++n)
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
                                   << "Key data `level_number' associated with structure `" << strct_name << "' is negative.");
                    }
                    else if (ln > d_max_levels)
                    {
                        TBOX_ERROR(d_object_name << ":  "
                                   << "Key data `level_number' associated with structure `" << strct_name << "' is greater than the expected maximum level number " << d_max_levels << ".");
                    }
                    d_base_filename[ln].push_back(strct_name);
                }
                else
                {
                    TBOX_ERROR(d_object_name << ":  "
                               << "Key data `level_number' not found in structure `" << strct_name << "' input.");
                }
            }
            else
            {
                TBOX_ERROR(d_object_name << ":  "
                           << "Key data `" << strct_name << "' not found in input.");
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
                const int n_files = db->getArraySize(db_key_name);
                d_base_filename[ln].resize(n_files);
                db->getStringArray(db_key_name, &d_base_filename[ln][0], n_files);
            }
            else
            {
                TBOX_WARNING(d_object_name << ":  "
                             << "Key data `" << db_key_name << "' not found in input.");
            }
        }
    }

    // Read in any shift and scale factors.
    if (db->keyExists("length_scale_factor"))
    {
        d_length_scale_factor = db->getDouble("length_scale_factor");
    }

    if (db->keyExists("posn_shift"))
    {
        db->getDoubleArray("posn_shift", &d_posn_shift[0], NDIM);
    }

    // Read in any sub-databases associated with the input file names.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();

        d_enable_springs[ln].resize(num_base_filename,true);
        d_using_uniform_spring_stiffness[ln].resize(num_base_filename,false);
        d_uniform_spring_stiffness[ln].resize(num_base_filename,-1.0);
        d_using_uniform_spring_rest_length[ln].resize(num_base_filename,false);
        d_uniform_spring_rest_length[ln].resize(num_base_filename,-1.0);
        d_using_uniform_spring_force_fcn_idx[ln].resize(num_base_filename,false);
        d_uniform_spring_force_fcn_idx[ln].resize(num_base_filename,-1);
        d_using_uniform_spring_subdomain_idx[ln].resize(num_base_filename,false);
        d_uniform_spring_subdomain_idx[ln].resize(num_base_filename,-1);

        d_enable_beams[ln].resize(num_base_filename,true);
        d_using_uniform_beam_bend_rigidity[ln].resize(num_base_filename,false);
        d_uniform_beam_bend_rigidity[ln].resize(num_base_filename,-1.0);
        d_using_uniform_beam_curvature[ln].resize(num_base_filename,false);
        d_uniform_beam_curvature[ln].resize(num_base_filename,std::vector<double>(NDIM,0.0));
        d_using_uniform_beam_subdomain_idx[ln].resize(num_base_filename,false);
        d_uniform_beam_subdomain_idx[ln].resize(num_base_filename,-1);

        d_enable_rods[ln].resize(num_base_filename,true);
        d_using_uniform_rod_properties[ln].resize(num_base_filename,false);
        d_uniform_rod_properties[ln].resize(num_base_filename,std::vector<double>(NUM_ROD_PARAMS,0.0));
        d_using_uniform_rod_subdomain_idx[ln].resize(num_base_filename,false);
        d_uniform_rod_subdomain_idx[ln].resize(num_base_filename,-1);

        d_enable_target_points[ln].resize(num_base_filename,true);
        d_using_uniform_target_stiffness[ln].resize(num_base_filename,false);
        d_uniform_target_stiffness[ln].resize(num_base_filename,-1.0);
        d_using_uniform_target_damping[ln].resize(num_base_filename,false);
        d_uniform_target_damping[ln].resize(num_base_filename,-1.0);
        d_using_uniform_target_subdomain_idx[ln].resize(num_base_filename,false);
        d_uniform_target_subdomain_idx[ln].resize(num_base_filename,-1);

        d_enable_anchor_points[ln].resize(num_base_filename,true);
        d_using_uniform_anchor_subdomain_idx[ln].resize(num_base_filename,false);
        d_uniform_anchor_subdomain_idx[ln].resize(num_base_filename,-1);

        d_enable_bdry_mass[ln].resize(num_base_filename,true);
        d_using_uniform_bdry_mass[ln].resize(num_base_filename,false);
        d_uniform_bdry_mass[ln].resize(num_base_filename,-1.0);
        d_using_uniform_bdry_mass_stiffness[ln].resize(num_base_filename,false);
        d_uniform_bdry_mass_stiffness[ln].resize(num_base_filename,-1.0);

        d_enable_instrumentation[ln].resize(num_base_filename,true);

        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string& base_filename = d_base_filename[ln][j];
            if (db->isDatabase(base_filename))
            {
                Pointer<Database> sub_db =
                    db->getDatabase(base_filename);

                // Determine whether to enable or disable any particular
                // features.
                if (sub_db->keyExists("enable_springs"))
                {
                    d_enable_springs[ln][j] = sub_db->getBool("enable_springs");
                }
                if (sub_db->keyExists("enable_beams"))
                {
                    d_enable_beams[ln][j] = sub_db->getBool("enable_beams");
                }
                if (sub_db->keyExists("enable_rods"))
                {
                    d_enable_rods[ln][j] = sub_db->getBool("enable_rods");
                }
                if (sub_db->keyExists("enable_target_points"))
                {
                    d_enable_target_points[ln][j] = sub_db->getBool("enable_target_points");
                }
                if (sub_db->keyExists("enable_anchor_points"))
                {
                    d_enable_anchor_points[ln][j] = sub_db->getBool("enable_anchor_points");
                }
                if (sub_db->keyExists("enable_bdry_mass"))
                {
                    d_enable_bdry_mass[ln][j] = sub_db->getBool("enable_bdry_mass");
                }
                if (sub_db->keyExists("enable_instrumentation"))
                {
                    d_enable_instrumentation[ln][j] = sub_db->getBool("enable_instrumentation");
                }

                // Determine whether to use uniform values for any particular
                // structure attributes.
                if (sub_db->keyExists("uniform_spring_stiffness"))
                {
                    d_using_uniform_spring_stiffness[ln][j] = true;
                    d_uniform_spring_stiffness[ln][j] = sub_db->getDouble("uniform_spring_stiffness");
                    if (d_uniform_spring_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_spring_stiffness' in database " << base_filename << std::endl
                                   << "  spring constant is negative" << std::endl);
                    }
                }
                if (sub_db->keyExists("uniform_spring_rest_length"))
                {
                    d_using_uniform_spring_rest_length[ln][j] = true;
                    d_uniform_spring_rest_length[ln][j] = sub_db->getDouble("uniform_spring_rest_length");
                    if (d_uniform_spring_rest_length[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_spring_rest_length' in database " << base_filename << std::endl
                                   << "  spring resting length is negative" << std::endl);
                    }
                }
                if (sub_db->keyExists("uniform_spring_force_fcn_idx"))
                {
                    d_using_uniform_spring_force_fcn_idx[ln][j] = true;
                    d_uniform_spring_force_fcn_idx[ln][j] = sub_db->getInteger("uniform_spring_force_fcn_idx");
                }
                if (sub_db->keyExists("uniform_spring_subdomain_idx"))
                {
                    d_using_uniform_spring_subdomain_idx[ln][j] = true;
                    d_uniform_spring_subdomain_idx[ln][j] = sub_db->getInteger("uniform_spring_subdomain_idx");
                }

                if (sub_db->keyExists("uniform_beam_bend_rigidity"))
                {
                    d_using_uniform_beam_bend_rigidity[ln][j] = true;
                    d_uniform_beam_bend_rigidity[ln][j] = sub_db->getDouble("uniform_beam_bend_rigidity");
                    if (d_uniform_beam_bend_rigidity[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_beam_bend_rigidity' in database " << base_filename << std::endl
                                   << "  beam bending rigidity is negative" << std::endl);
                    }
                }
                if (sub_db->keyExists("uniform_beam_curvature"))
                {
                    d_using_uniform_beam_curvature[ln][j] = true;
                    sub_db->getDoubleArray("uniform_beam_curvature", &d_uniform_beam_curvature[ln][j][0], NDIM);
                }
                if (sub_db->keyExists("uniform_beam_subdomain_idx"))
                {
                    d_using_uniform_beam_subdomain_idx[ln][j] = true;
                    d_uniform_beam_subdomain_idx[ln][j] = sub_db->getInteger("uniform_beam_subdomain_idx");
                }

                if (sub_db->keyExists("uniform_rod_properties"))
                {
                    d_using_uniform_rod_properties[ln][j] = true;
                    sub_db->getDoubleArray("uniform_rod_properties", &d_uniform_rod_properties[ln][j][0], NUM_ROD_PARAMS);
                }
                if (sub_db->keyExists("uniform_rod_subdomain_idx"))
                {
                    d_using_uniform_rod_subdomain_idx[ln][j] = true;
                    d_uniform_rod_subdomain_idx[ln][j] = sub_db->getInteger("uniform_rod_subdomain_idx");
                }

                if (sub_db->keyExists("uniform_target_stiffness"))
                {
                    d_using_uniform_target_stiffness[ln][j] = true;
                    d_uniform_target_stiffness[ln][j] = sub_db->getDouble("uniform_target_stiffness");
                    if (d_uniform_target_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_stiffness' in database " << base_filename << std::endl
                                   << "  target point spring constant is negative" << std::endl);
                    }
                }
                if (sub_db->keyExists("uniform_target_damping"))
                {
                    d_using_uniform_target_damping[ln][j] = true;
                    d_uniform_target_damping[ln][j] = sub_db->getDouble("uniform_target_damping");
                    if (d_uniform_target_damping[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_damping' in database " << base_filename << std::endl
                                   << "  target point spring constant is negative" << std::endl);
                    }
                }
                if (sub_db->keyExists("uniform_target_subdomain_idx"))
                {
                    d_using_uniform_target_subdomain_idx[ln][j] = true;
                    d_uniform_target_subdomain_idx[ln][j] = sub_db->getInteger("uniform_target_subdomain_idx");
                }

                if (sub_db->keyExists("uniform_anchor_subdomain_idx"))
                {
                    d_using_uniform_anchor_subdomain_idx[ln][j] = true;
                    d_uniform_anchor_subdomain_idx[ln][j] = sub_db->getInteger("uniform_anchor_subdomain_idx");
                }

                if (sub_db->keyExists("uniform_bdry_mass"))
                {
                    d_using_uniform_bdry_mass[ln][j] = true;
                    d_uniform_bdry_mass[ln][j] = sub_db->getDouble("uniform_bdry_mass");

                    if (d_uniform_bdry_mass[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_bdry_mass' in database " << base_filename << std::endl
                                   << "  boundary mass is negative" << std::endl);
                    }
                }
                if (sub_db->keyExists("uniform_bdry_mass_stiffness"))
                {
                    d_using_uniform_bdry_mass_stiffness[ln][j] = true;
                    d_uniform_bdry_mass_stiffness[ln][j] = sub_db->getDouble("uniform_bdry_mass_stiffness");

                    if (d_uniform_bdry_mass_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_bdry_mass_stiffness' in database " << base_filename << std::endl
                                   << "  boundary mass spring constant is negative" << std::endl);
                    }
                }
            }
        }
    }

    // Output the names of the input files to be read along with additional
    // debugging information.
    pout << d_object_name << ":  Reading from input files: " << std::endl;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string& base_filename = d_base_filename[ln][j];
            pout << "  base filename: " << base_filename << std::endl
                 << "  assigned to level " << ln << " of the Cartesian grid patch hierarchy" << std::endl
                 << "     required files: " << base_filename << ".vertex" << std::endl
                 << "     optional files: " << base_filename << ".spring, " << base_filename << ".beam, " << base_filename << ".rod, " << base_filename << ".target, " << base_filename << ".anchor, " << base_filename << ".mass, " << base_filename << ".director, " << base_filename << ".inst" << std::endl;
            if (!d_enable_springs[ln][j])
            {
                pout << "  NOTE: spring forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_spring_stiffness[ln][j])
                {
                    pout << "  NOTE: uniform spring stiffnesses are being employed for the structure named " << base_filename << std::endl
                         << "        any stiffness information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
                if (d_using_uniform_spring_rest_length[ln][j])
                {
                    pout << "  NOTE: uniform spring resting lengths are being employed for the structure named " << base_filename << std::endl
                         << "        any resting length information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
                if (d_using_uniform_spring_force_fcn_idx[ln][j])
                {
                    pout << "  NOTE: uniform spring force functions are being employed for the structure named " << base_filename << std::endl
                         << "        any force function index information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
                if (d_using_uniform_spring_subdomain_idx[ln][j])
                {
                    pout << "  NOTE: uniform spring subdomain indicies are being employed for the structure named " << base_filename << std::endl
                         << "        any subdomain index information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_beams[ln][j])
            {
                pout << "  NOTE: beam forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_beam_bend_rigidity[ln][j])
                {
                    pout << "  NOTE: uniform beam bending rigidities are being employed for the structure named " << base_filename << std::endl
                         << "        any stiffness information in optional file " << base_filename << ".beam will be IGNORED" << std::endl;
                }
                if (d_using_uniform_beam_curvature[ln][j])
                {
                    pout << "  NOTE: uniform beam curvatures are being employed for the structure named " << base_filename << std::endl
                         << "        any curvature information in optional file " << base_filename << ".beam will be IGNORED" << std::endl;
                }
                if (d_using_uniform_beam_subdomain_idx[ln][j])
                {
                    pout << "  NOTE: uniform beam subdomain indicies are being employed for the structure named " << base_filename << std::endl
                         << "        any subdomain index information in optional file " << base_filename << ".beam will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_rods[ln][j])
            {
                pout << "  NOTE: rod forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_rod_properties[ln][j])
                {
                    pout << "  NOTE: uniform rod material properties are being employed for the structure named " << base_filename << std::endl
                         << "        any material property information in optional file " << base_filename << ".rod will be IGNORED" << std::endl;
                }
                if (d_using_uniform_rod_subdomain_idx[ln][j])
                {
                    pout << "  NOTE: uniform rod subdomain indicies are being employed for the structure named " << base_filename << std::endl
                         << "        any subdomain index information in optional file " << base_filename << ".rod will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_target_points[ln][j])
            {
                pout << "  NOTE: target point penalty forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_target_stiffness[ln][j])
                {
                    pout << "  NOTE: uniform target point stiffnesses are being employed for the structure named " << base_filename << std::endl
                         << "        any target point stiffness information in optional file " << base_filename << ".target will be IGNORED" << std::endl;
                }
                if (d_using_uniform_target_damping[ln][j])
                {
                    pout << "  NOTE: uniform target point damping factors are being employed for the structure named " << base_filename << std::endl
                         << "        any target point damping factor information in optional file " << base_filename << ".target will be IGNORED" << std::endl;
                }
                if (d_using_uniform_target_subdomain_idx[ln][j])
                {
                    pout << "  NOTE: uniform target point subdomain indicies are being employed for the structure named " << base_filename << std::endl
                         << "        any subdomain index information in optional file " << base_filename << ".target will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_anchor_points[ln][j])
            {
                pout << "  NOTE: anchor points are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_anchor_subdomain_idx[ln][j])
                {
                    pout << "  NOTE: uniform anchor point subdomain indicies are being employed for the structure named " << base_filename << std::endl
                         << "        any subdomain index information in optional file " << base_filename << ".anchor will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_bdry_mass[ln][j])
            {
                pout << "  NOTE: massive boundary points are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_bdry_mass[ln][j])
                {
                    pout << "  NOTE: uniform boundary point masses are being employed for the structure named " << base_filename << std::endl
                         << "        any boundary point mass information in optional file " << base_filename << ".mass will be IGNORED" << std::endl;
                }
                if (d_using_uniform_bdry_mass_stiffness[ln][j])
                {
                    pout << "  NOTE: uniform massive boundary point stiffnesses are being employed for the structure named " << base_filename << std::endl
                         << "        any massive boundary point stiffness information in optional file " << base_filename << ".mass will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_instrumentation[ln][j])
            {
                pout << "  NOTE: instrumentation is DISABLED for " << base_filename << std::endl;
            }

            pout << std::endl;
        }
    }
    return;
}// getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBStandardInitializer>;

//////////////////////////////////////////////////////////////////////////////
