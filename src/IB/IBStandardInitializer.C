// Filename: IBStandardInitializer.C
// Last modified: <01.Apr.2008 17:09:37 griffith@box221.cims.nyu.edu>
// Created on 22 Nov 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)

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
#include <ibamr/IBBeamForceSpec.h>
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBTargetPointForceSpec.h>

// IBTK INCLUDES
#include <ibtk/IndexUtilities.h>
#include <ibtk/LNodeIndexData2.h>

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
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
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
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
      d_spring_stiffness(),
      d_spring_rest_length(),
      d_spring_force_fcn_idx(),
      d_using_uniform_spring_stiffness(),
      d_uniform_spring_stiffness(),
      d_using_uniform_spring_rest_length(),
      d_uniform_spring_rest_length(),
      d_using_uniform_spring_force_fcn_idx(),
      d_uniform_spring_force_fcn_idx(),
      d_enable_beams(),
      d_beam_specs(),
      d_using_uniform_beam_bend_rigidity(),
      d_uniform_beam_bend_rigidity(),
      d_enable_target_points(),
      d_target_stiffness(),
      d_using_uniform_target_stiffness(),
      d_uniform_target_stiffness(),
      d_enable_bdry_mass(),
      d_bdry_mass(),
      d_bdry_mass_stiffness(),
      d_using_uniform_bdry_mass(),
      d_uniform_bdry_mass(),
      d_using_uniform_bdry_mass_stiffness(),
      d_uniform_bdry_mass_stiffness(),
      d_enable_instrumentation(),
      d_instrument_idx(),
      d_global_index_offset()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif

    // Register the specification objects with the IBTK::StashableManager class.
    IBSpringForceSpec::registerWithStashableManager();
    IBBeamForceSpec::registerWithStashableManager();
    IBTargetPointForceSpec::registerWithStashableManager();
    IBInstrumentationSpec::registerWithStashableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Check to see if we are starting from a restart file.
    SAMRAI::tbox::RestartManager* restart_manager = SAMRAI::tbox::RestartManager::getManager();
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

        // Process the (optional) target point information.
        readTargetPointFiles();

        // Process the (optional) mass information.
        readBoundaryMassFiles();

        // Process the (optional) instrumentation information.
        readInstrumentationFiles();

        // Wait for all processes to finish.
        SAMRAI::tbox::SAMRAI_MPI::barrier();
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
    SAMRAI::tbox::Pointer<IBTK::LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!silo_writer.isNull());
#endif

    // Cache a pointer to the data writer.
    d_silo_writer = silo_writer;

    // Check to see if we are starting from a restart file.
    SAMRAI::tbox::RestartManager* restart_manager = SAMRAI::tbox::RestartManager::getManager();
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
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Loop over all patches in the specified level of the patch level and count
    // the number of local vertices.
    int local_node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of vertices whose initial locations will be within
        // the given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
    }
    return local_node_count;
}// getLocalNodeCountOnPatchLevel

int
IBStandardInitializer::initializeDataOnPatchLevel(
    const int lag_node_index_idx,
    const int global_index_offset,
    const int local_index_offset,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    (void) lag_manager;

    // Determine the extents of the physical domain.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    const double* const XLower = grid_geom->getXLower();
    const double* const XUpper = grid_geom->getXUpper();

    // Set the global index offset.  This is equal to the number of Lagrangian
    // indices that have already been initialized on the specified level.
    d_global_index_offset[level_number] = global_index_offset;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    int local_idx = -1;
    int local_node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::pdat::CellIndex<NDIM>& patch_lower = patch_box.lower();
        const SAMRAI::pdat::CellIndex<NDIM>& patch_upper = patch_box.upper();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> index_data =
            patch->getPatchData(lag_node_index_idx);

        // Initialize the vertices whose initial locations will be within the
        // given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int>& point_idx = (*it);
            const int current_global_idx = getCannonicalLagrangianIndex(
                point_idx, level_number) + global_index_offset;
            const int current_local_idx = ++local_idx + local_index_offset;

            // Get the coordinates of the present vertex.
            const std::vector<double> X = getVertexPosn(point_idx, level_number);

            // Initialize the location of the present vertex.
            double* const node_X = &(*X_data)(current_local_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                node_X[d] = X[d];

                if (SAMRAI::tbox::MathUtilities<double>::equalEps(X[d],XLower[d]))
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node intersecting lower physical boundary.\n"
                               << "  please ensure that all nodes are within the computational domain."<< std::endl);
                }
                else if (X[d] <= XLower[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node below lower physical boundary\n"
                               << "  please ensure that all nodes are within the computational domain."<< std::endl);
                }

                if (SAMRAI::tbox::MathUtilities<double>::equalEps(X[d],XUpper[d]))
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node intersecting upper physical boundary.\n"
                               << "  please ensure that all nodes are within the computational domain."<< std::endl);
                }
                else if (X[d] >= XUpper[d])
                {
                    TBOX_ERROR(d_object_name << "::initializeDataOnPatchLevel():\n"
                               << "  encountered node above upper physical boundary\n"
                               << "  please ensure that all nodes are within the computational domain."<< std::endl);
                }
            }

            // Get the index of the cell in which the present vertex is
            // initially located.
            const SAMRAI::pdat::CellIndex<NDIM> idx = IBTK::IndexUtilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);

            // Initialize the force specification object assocaited with the
            // present vertex.
            std::vector<SAMRAI::tbox::Pointer<IBTK::Stashable> > force_spec =
                initializeSpecs(
                    point_idx, global_index_offset, level_number);

            IBTK::LNodeIndexSet& node_set = (*index_data)(idx);
            node_set.push_back(
                new IBTK::LNodeIndex(current_global_idx, current_local_idx,
                                     &(*X_data)(current_local_idx), force_spec));

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
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& M_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData>& K_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    (void) lag_manager;

    // Loop over all patches in the specified level of the patch level and
    // initialize the local vertices.
    int local_idx = -1;
    int local_node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

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
            const double M = getVertexMass(point_idx, level_number);
            const double K = getVertexMassStiffness(point_idx, level_number);

            // Avoid division by zero at massless nodes.
            if (SAMRAI::tbox::MathUtilities<double>::equalEps(M,0.0))
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

void
IBStandardInitializer::tagCellsForInitialRefinement(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    // Loop over all patches in the specified level of the patch level and tag
    // cells for refinement wherever there are vertices assigned to a finer
    // level of the Cartesian grid.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::pdat::CellIndex<NDIM>& patch_lower = patch_box.lower();
        const SAMRAI::pdat::CellIndex<NDIM>& patch_upper = patch_box.upper();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

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
                const SAMRAI::pdat::CellIndex<NDIM> i = IBTK::IndexUtilities::getCellIndex(
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
    if (SAMRAI::tbox::SAMRAI_MPI::getRank() == 0)
    {
        for (unsigned j = 0; j < d_num_vertex[level_number].size(); ++j)
        {
            d_silo_writer->registerMarkerCloud(
                d_base_filename[level_number][j] + "_vertices",
                d_num_vertex[level_number][j], d_vertex_offset[level_number][j], level_number);
            if (d_spring_edge_map[level_number][j].size() > 0)
            {
                d_silo_writer->registerUnstructuredMesh(
                    d_base_filename[level_number][j] + "_mesh",
                    d_spring_edge_map[level_number][j], level_number);
            }
        }
    }
    return;
}// initializeLagSiloDataWriter

void
IBStandardInitializer::readVertexFiles()
{
    std::string line_string;
    const int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
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
            if (d_use_file_batons && rank != 0) SAMRAI::tbox::SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

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

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "processing vertex data from ASCII input filename " << vertex_filename << std::endl
                               << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

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

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "read " << d_num_vertex[ln][j] << " vertices from ASCII input filename " << vertex_filename << std::endl
                               << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI::tbox::SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI::tbox::SAMRAI_MPI::barrier();
    return;
}// readVertexFiles

void
IBStandardInitializer::readSpringFiles()
{
    std::string line_string;
    const int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_spring_edge_map[ln].resize(num_base_filename);
        d_spring_stiffness[ln].resize(num_base_filename);
        d_spring_rest_length[ln].resize(num_base_filename);
        d_spring_force_fcn_idx[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI::tbox::SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            // Ensure that the file exists.
            const std::string spring_filename = d_base_filename[ln][j] + ".spring";
            std::ifstream file_stream;
            file_stream.open(spring_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing spring data from ASCII input filename " << spring_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

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
                    int force_fcn_idx;
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
                    }

                    // Modify kappa and length according to whether spring forces are
                    // enabled, or whether uniform values are to be employed, for this
                    // particular structure.
                    if (!d_enable_springs[ln][j])
                    {
                        kappa = 0.0;
                        length = 0.0;
                    }
                    else
                    {
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
                        const Edge& other_e = (*it).second;
                        if (e.first  == other_e.first &&
                            e.second == other_e.second)
                        {
                            // This is a duplicate edge and should not be inserted into the
                            // edge map.
                            duplicate_edge = true;

                            // Ensure that the stiffness and rest length information is
                            // consistent.
                            if (!SAMRAI::tbox::MathUtilities<double>::equalEps(
                                    (*d_spring_stiffness[ln][j].find(e)).second, kappa) ||
                                !SAMRAI::tbox::MathUtilities<double>::equalEps(
                                    (*d_spring_rest_length[ln][j].find(e)).second, length) ||
                                !SAMRAI::tbox::MathUtilities<double>::equalEps(
                                    (*d_spring_force_fcn_idx[ln][j].find(e)).second, force_fcn_idx))
                            {
                                TBOX_ERROR(d_object_name << ":\n  Inconsistent duplicate edges in input file encountered on line " << k+2 << " of file " << spring_filename << std::endl
                                           << "  first vertex = " << e.first-d_vertex_offset[ln][j] << " second vertex = " << e.second-d_vertex_offset[ln][j] << std::endl
                                           << "  original spring constant = " << (*d_spring_stiffness[ln][j].find(e)).second << std::endl
                                           << "  original resting length = " << (*d_spring_rest_length[ln][j].find(e)).second << std::endl
                                           << "  original force function index = " << (*d_spring_force_fcn_idx[ln][j].find(e)).second << std::endl);
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
                        d_spring_stiffness[ln][j][e] = kappa;
                        d_spring_rest_length[ln][j][e] = length;
                        d_spring_force_fcn_idx[ln][j][e] = force_fcn_idx;
                    }
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_edges << " edges from ASCII input filename " << spring_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI::tbox::SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI::tbox::SAMRAI_MPI::barrier();
    return;
}// readSpringFiles

void
IBStandardInitializer::readBeamFiles()
{
    std::string line_string;
    const int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_beam_specs[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI::tbox::SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            const std::string beam_filename = d_base_filename[ln][j] + ".beam";
            std::ifstream file_stream;
            file_stream.open(beam_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing beam data from input filename " << beam_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

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
                    double kappa;
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

                        if (!(line_stream >> kappa))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl);
                        }
                        else if (kappa < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << std::endl
                                       << "  beam constant is negative" << std::endl);
                        }
                    }

                    // Modify kappa according to whether beam forces are
                    // enabled, or whether uniform values are to be employed,
                    // for this particular structure.
                    if (!d_enable_beams[ln][j])
                    {
                        kappa = 0.0;
                    }
                    else
                    {
                        if (d_using_uniform_beam_bend_rigidity[ln][j])
                        {
                            kappa = d_uniform_beam_bend_rigidity[ln][j];
                        }
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
                    d_beam_specs[ln][j].insert(
                        std::make_pair(
                            curr_idx, std::make_pair(
                                std::make_pair(next_idx,prev_idx),kappa)));
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_beams << " beams from input filename " << beam_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI::tbox::SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI::tbox::SAMRAI_MPI::barrier();
    return;
}// readBeamFiles

void
IBStandardInitializer::readTargetPointFiles()
{
    std::string line_string;
    const int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_target_stiffness[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI::tbox::SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            d_target_stiffness[ln][j].resize(d_num_vertex[ln][j], 0.0);

            const std::string target_point_stiffness_filename = d_base_filename[ln][j] + ".target";
            std::ifstream file_stream;
            file_stream.open(target_point_stiffness_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing target point data from input filename " << target_point_stiffness_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

                // The first line in the file indicates the number of target
                // point stiffnesses in the input file.
                int num_target_stiffness;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << target_point_stiffness_filename << std::endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_target_stiffness))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << target_point_stiffness_filename << std::endl);
                    }
                }

                if (num_target_stiffness <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << target_point_stiffness_filename << std::endl);
                }

                // Each successive line indicates the vertex number and spring
                // constant associated with any target points.
                for (int k = 0; k < num_target_stiffness; ++k)
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

                        if (!(line_stream >> d_target_stiffness[ln][j][n]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl);
                        }
                        else if (d_target_stiffness[ln][j][n] < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << std::endl
                                       << "  target point spring constant is negative" << std::endl);
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_target_stiffness << " target points from input filename " << target_point_stiffness_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;
            }

            // Modify the target point stiffness constant according to whether
            // target point penalty forces are enabled, or whether uniform
            // values are to be employed, for this particular structure.
            if (!d_enable_target_points[ln][j])
            {
                d_target_stiffness[ln][j] = std::vector<double>(
                    d_num_vertex[ln][j], 0.0);
            }
            else
            {
                if (d_using_uniform_target_stiffness[ln][j])
                {
                    d_target_stiffness[ln][j] = std::vector<double>(
                        d_num_vertex[ln][j],
                        d_uniform_target_stiffness[ln][j]);
                }
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI::tbox::SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }

    // Synchronize the processes.
    if (d_use_file_batons) SAMRAI::tbox::SAMRAI_MPI::barrier();
    return;
}// readTargetPointFiles

void
IBStandardInitializer::readBoundaryMassFiles()
{
    std::string line_string;
    const int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
    int flag = 1;
    int sz = 1;

    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        d_bdry_mass[ln].resize(num_base_filename);
        d_bdry_mass_stiffness[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            // Wait for the previous MPI process to finish reading the current file.
            if (d_use_file_batons && rank != 0) SAMRAI::tbox::SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            d_bdry_mass[ln][j].resize(d_num_vertex[ln][j], 0.0);
            d_bdry_mass_stiffness[ln][j].resize(d_num_vertex[ln][j], 0.0);

            const std::string bdry_mass_filename = d_base_filename[ln][j] + ".mass";
            std::ifstream file_stream;
            file_stream.open(bdry_mass_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing boundary mass data from input filename " << bdry_mass_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

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

                        if (!(line_stream >> d_bdry_mass[ln][j][n]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl);
                        }
                        else if (d_bdry_mass[ln][j][n] < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl
                                       << "  boundary mass is negative" << std::endl);
                        }

                        if (!(line_stream >> d_bdry_mass_stiffness[ln][j][n]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl);
                        }
                        else if (d_bdry_mass_stiffness[ln][j][n] < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << std::endl
                                       << "  boundary mass spring constant is negative" << std::endl);
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_bdry_mass_pts << " boundary mass points from input filename " << bdry_mass_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;
            }

            // Modify the boundary mass and boundary mass stiffness constant
            // according to whether boundary mass is enabled, or whether uniform
            // values are to be employed, for this particular structure.
            if (!d_enable_bdry_mass[ln][j])
            {
                d_bdry_mass[ln][j] = std::vector<double>(
                    d_num_vertex[ln][j], 0.0);
                d_bdry_mass_stiffness[ln][j] = std::vector<double>(
                    d_num_vertex[ln][j], 0.0);
            }
            else
            {
                if (d_using_uniform_bdry_mass[ln][j])
                {
                    d_bdry_mass[ln][j] = std::vector<double>(
                        d_num_vertex[ln][j],
                        d_uniform_bdry_mass[ln][j]);
                }
                if (d_using_uniform_bdry_mass_stiffness[ln][j])
                {
                    d_bdry_mass_stiffness[ln][j] = std::vector<double>(
                        d_num_vertex[ln][j],
                        d_uniform_bdry_mass_stiffness[ln][j]);
                }
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI::tbox::SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }
    return;
}// readBoundaryMassFiles

void
IBStandardInitializer::readInstrumentationFiles()
{
    std::string line_string;
    const int rank = SAMRAI::tbox::SAMRAI_MPI::getRank();
    const int nodes = SAMRAI::tbox::SAMRAI_MPI::getNodes();
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
            if (d_use_file_batons && rank != 0) SAMRAI::tbox::SAMRAI_MPI::recv(&flag, sz, rank-1, false, j);

            const std::string inst_filename = d_base_filename[ln][j] + ".inst";
            std::ifstream file_stream;
            file_stream.open(inst_filename.c_str(), std::ios::in);
            if (file_stream.is_open() && d_enable_instrumentation[ln][j])
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing instrumentation data from input filename " << inst_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;

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

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_inst_pts << " instrumentation points from input filename " << inst_filename << std::endl
                                   << "  on MPI process " << SAMRAI::tbox::SAMRAI_MPI::getRank() << std::endl;
            }

            // Free the next MPI process to start reading the current file.
            if (d_use_file_batons && rank != nodes-1) SAMRAI::tbox::SAMRAI_MPI::send(&flag, sz, rank+1, false, j);
        }
    }
    IBInstrumentationSpec::setInstrumentNames(instrument_names);
    return;
}// readInstrumentationFiles

void
IBStandardInitializer::getPatchVertices(
    std::vector<std::pair<int,int> >& patch_vertices,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const int level_number,
    const bool can_be_refined) const
{
    // Loop over all of the vertices to determine the indices of those vertices
    // within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will work for
    // now.
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
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
IBStandardInitializer::getCannonicalLagrangianIndex(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_vertex_offset[level_number][point_index.first]+point_index.second;
}// getCannonicalLagrangianIndex

std::vector<double>
IBStandardInitializer::getVertexPosn(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return std::vector<double>(
        &d_vertex_posn[level_number][point_index.first][point_index.second*NDIM     ],
        &d_vertex_posn[level_number][point_index.first][point_index.second*NDIM+NDIM]);
}// getVertexPosn

double
IBStandardInitializer::getVertexTargetStiffness(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_target_stiffness[level_number][point_index.first][point_index.second];
}// getVertexTargetStiffness

double
IBStandardInitializer::getVertexMass(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_bdry_mass[level_number][point_index.first][point_index.second];
}// getVertexMass

double
IBStandardInitializer::getVertexMassStiffness(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_bdry_mass_stiffness[level_number][point_index.first][point_index.second];
}// getVertexMassStiffness

std::pair<int,int>
IBStandardInitializer::getVertexInstrumentationIndices(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    std::map<int,std::pair<int,int> >::const_iterator it =
        d_instrument_idx[level_number][point_index.first].find(point_index.second);
    if (it != d_instrument_idx[level_number][point_index.first].end())
    {
        return (*it).second;
    }
    else
    {
        return std::make_pair(-1,-1);
    }
}// getVertexInstrumentationIndices

std::vector<SAMRAI::tbox::Pointer<IBTK::Stashable> >
IBStandardInitializer::initializeSpecs(
    const std::pair<int,int>& point_index,
    const int global_index_offset,
    const int level_number) const
{
    std::vector<SAMRAI::tbox::Pointer<IBTK::Stashable> > vertex_specs;

    const int j = point_index.first;
    const int mastr_idx = getCannonicalLagrangianIndex(point_index, level_number);

    std::vector<int> slave_idxs, force_fcn_idxs;
    std::vector<double> stiffness, rest_length;
    for (std::multimap<int,Edge>::const_iterator it = d_spring_edge_map[level_number][j].lower_bound(mastr_idx);
         it != d_spring_edge_map[level_number][j].upper_bound(mastr_idx); ++it)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(mastr_idx == (*it).first);
#endif
        // The connectivity information.
        const Edge& e = (*it).second;
        if (e.first == mastr_idx)
        {
            slave_idxs.push_back(e.second+global_index_offset);
        }
        else
        {
            slave_idxs.push_back(e.first+global_index_offset);
        }

        // The material properties.
        force_fcn_idxs.push_back((*d_spring_force_fcn_idx[level_number][j].find(e)).second);
        stiffness     .push_back((*d_spring_stiffness    [level_number][j].find(e)).second);
        rest_length   .push_back((*d_spring_rest_length  [level_number][j].find(e)).second);
    }

    if (slave_idxs.size() > 0)
    {
        vertex_specs.push_back(
            new IBSpringForceSpec(
                mastr_idx, slave_idxs, force_fcn_idxs, stiffness, rest_length));
    }

    std::vector<std::pair<int,int> > beam_neighbor_idxs;
    std::vector<double> beam_bend_rigidity;
    for (std::multimap<int,std::pair<Neighbors,double> >::const_iterator it = d_beam_specs[level_number][j].lower_bound(mastr_idx);
         it != d_beam_specs[level_number][j].upper_bound(mastr_idx); ++it)
    {
        const std::pair<int,int>& neighbor_idxs = (*it).second.first;
        const double& bend_rigidity = (*it).second.second;
        if (!SAMRAI::tbox::MathUtilities<double>::equalEps(bend_rigidity,0.0))
        {
            beam_neighbor_idxs.push_back(neighbor_idxs);
            beam_bend_rigidity.push_back(bend_rigidity);
        }
    }

    if (!beam_neighbor_idxs.empty())
    {
        vertex_specs.push_back(
            new IBBeamForceSpec(
                mastr_idx, beam_neighbor_idxs, beam_bend_rigidity));
    }

    const double kappa_target = getVertexTargetStiffness(point_index, level_number);
    const std::vector<double> X_target = getVertexPosn(point_index, level_number);

    if (!SAMRAI::tbox::MathUtilities<double>::equalEps(kappa_target,0.0))
    {
        vertex_specs.push_back(
            new IBTargetPointForceSpec(
                mastr_idx, kappa_target, X_target));
    }

    const std::pair<int,int> inst_idx = getVertexInstrumentationIndices(point_index, level_number);

    if (inst_idx.first != -1 && inst_idx.second != -1)
    {
        vertex_specs.push_back(
            new IBInstrumentationSpec(
                mastr_idx, inst_idx.first, inst_idx.second));
    }
    return vertex_specs;
}// initializeSpecs

void
IBStandardInitializer::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
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
    d_spring_stiffness.resize(d_max_levels);
    d_spring_rest_length.resize(d_max_levels);
    d_spring_force_fcn_idx.resize(d_max_levels);
    d_using_uniform_spring_stiffness.resize(d_max_levels);
    d_uniform_spring_stiffness.resize(d_max_levels);
    d_using_uniform_spring_rest_length.resize(d_max_levels);
    d_uniform_spring_rest_length.resize(d_max_levels);
    d_using_uniform_spring_force_fcn_idx.resize(d_max_levels);
    d_uniform_spring_force_fcn_idx.resize(d_max_levels);

    d_enable_beams.resize(d_max_levels);
    d_beam_specs.resize(d_max_levels);
    d_using_uniform_beam_bend_rigidity.resize(d_max_levels);
    d_uniform_beam_bend_rigidity.resize(d_max_levels);

    d_enable_target_points.resize(d_max_levels);
    d_target_stiffness.resize(d_max_levels);
    d_using_uniform_target_stiffness.resize(d_max_levels);
    d_uniform_target_stiffness.resize(d_max_levels);

    d_enable_bdry_mass.resize(d_max_levels);
    d_bdry_mass.resize(d_max_levels);
    d_bdry_mass_stiffness.resize(d_max_levels);
    d_using_uniform_bdry_mass.resize(d_max_levels);
    d_uniform_bdry_mass.resize(d_max_levels);
    d_using_uniform_bdry_mass_stiffness.resize(d_max_levels);
    d_uniform_bdry_mass_stiffness.resize(d_max_levels);

    d_enable_instrumentation.resize(d_max_levels);
    d_instrument_idx.resize(d_max_levels);

    d_global_index_offset.resize(d_max_levels);

    // Determine the various input file names.
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
                         << "Key data `" + db_key_name + "' not found in input.");
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

        d_enable_beams[ln].resize(num_base_filename,true);

        d_using_uniform_beam_bend_rigidity[ln].resize(num_base_filename,false);
        d_uniform_beam_bend_rigidity[ln].resize(num_base_filename,-1.0);

        d_enable_target_points[ln].resize(num_base_filename,true);

        d_using_uniform_target_stiffness[ln].resize(num_base_filename,false);
        d_uniform_target_stiffness[ln].resize(num_base_filename,-1.0);

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
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> sub_db =
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
                if (sub_db->keyExists("enable_target_points"))
                {
                    d_enable_target_points[ln][j] = sub_db->getBool("enable_target_points");
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
    SAMRAI::tbox::pout << d_object_name << ":  Reading from input files: " << std::endl;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = d_base_filename[ln].size();
        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string& base_filename = d_base_filename[ln][j];
            SAMRAI::tbox::pout << "  base filename: " << base_filename << std::endl
                               << "  assigned to level " << ln << " of the Cartesian grid patch hierarchy" << std::endl
                               << "     required files: " << base_filename << ".vertex" << std::endl
                               << "     optional files: " << base_filename << ".spring, " << base_filename << ".beam, " << base_filename << ".target, " << base_filename << ".mass, " << base_filename << ".inst " << std::endl;
            if (!d_enable_springs[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: spring forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_spring_stiffness[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform spring stiffnesses are being employed for the structure named " << base_filename << std::endl
                                       << "        any stiffness information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
                if (d_using_uniform_spring_rest_length[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform spring resting lengths are being employed for the structure named " << base_filename << std::endl
                                       << "        any resting length information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
                if (d_using_uniform_spring_force_fcn_idx[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform spring force functions are being employed for the structure named " << base_filename << std::endl
                                       << "        any force function index information in optional file " << base_filename << ".spring will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_beams[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: beam forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_beam_bend_rigidity[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform beam bending rigidities are being employed for the structure named " << base_filename << std::endl
                                       << "        any stiffness information in optional file " << base_filename << ".beam will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_target_points[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: target point penalty forces are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_target_stiffness[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform target point stiffnesses are being employed for the structure named " << base_filename << std::endl
                                       << "        any target point stiffness information in optional file " << base_filename << ".target will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_bdry_mass[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: massive boundary points are DISABLED for " << base_filename << std::endl;
            }
            else
            {
                if (d_using_uniform_bdry_mass[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform boundary point masses are being employed for the structure named " << base_filename << std::endl
                                       << "        any boundary point mass information in optional file " << base_filename << ".mass will be IGNORED" << std::endl;
                }
                if (d_using_uniform_bdry_mass_stiffness[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform massive boundary point stiffnesses are being employed for the structure named " << base_filename << std::endl
                                       << "        any massive boundary point stiffness information in optional file " << base_filename << ".mass will be IGNORED" << std::endl;
                }
            }

            if (!d_enable_instrumentation[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: instrumentation is DISABLED for " << base_filename << std::endl;
            }

            SAMRAI::tbox::pout << std::endl;
        }
    }
    return;
}// getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardInitializer>;

//////////////////////////////////////////////////////////////////////////////
