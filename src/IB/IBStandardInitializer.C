// Filename: IBStandardInitializer.C
// Last modified: <31.Jan.2007 01:05:13 boyce@bigboy.nyconnect.com>
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
#include <ibamr/IBStandardForceSpec.h>
#include <ibamr/LNodeIndexData.h>

// STOOLS INCLUDES
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <Index.h>
#include <tbox/MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <fstream>
#include <iostream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardInitializer::IBStandardInitializer(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_max_levels(-1),
      d_base_filename(),
      d_num_vertex(),
      d_vertex_offset(),
      d_vertex_posn(),
      d_enable_edges(),
      d_edge_map(),
      d_edge_stiffness(),
      d_edge_rest_length(),
      d_use_uniform_edge_stiffness(),
      d_uniform_edge_stiffness(),
      d_use_uniform_edge_rest_length(),
      d_uniform_edge_rest_length(),
      d_enable_target_points(),
      d_target_stiffness(),
      d_use_uniform_target_stiffness(),
      d_uniform_target_stiffness(),
      d_enable_bdry_mass(),
      d_bdry_mass(),
      d_bdry_mass_stiffness(),
      d_use_uniform_bdry_mass(),
      d_uniform_bdry_mass(),
      d_use_uniform_bdry_mass_stiffness(),
      d_uniform_bdry_mass_stiffness(),
      d_global_index_offset()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
#endif

    // Register the force specification object with the
    // StashableManager class.
    IBStandardForceSpec::registerWithStashableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Process the various input files on each MPI process.
    for (int rank = 0; rank < SAMRAI::tbox::MPI::getNodes(); ++rank)
    {
        if (rank == SAMRAI::tbox::MPI::getRank())
        {
            // Process the vertex information.
            readVertexFiles();

            // Compute the index offsets.  Note that this must be done
            // prior to processing all remaining input files.
            //
            // NOTE: A separate, independent Lagrangian numbering
            // scheme is used on each level of the locally refined
            // Cartesian grid.
            for (int ln = 0; ln < d_max_levels; ++ln)
            {
                d_vertex_offset[ln].resize(d_num_vertex[ln].size());
                d_vertex_offset[ln][0] = 0;
                for (int j = 1; j < static_cast<int>(d_num_vertex[ln].size()); ++j)
                {
                    d_vertex_offset[ln][j] = d_vertex_offset[ln][j-1]+d_num_vertex[ln][j-1];
                }
            }

            // Process the (optional) edge information.
            readEdgeFiles();

            // Process the (optional) target point information.
            readTargetPointFiles();

            // Process the (optional) mass information.
            readBoundaryMassFiles();
        }
        SAMRAI::tbox::MPI::barrier();
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
    SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!silo_writer.isNull());
#endif

    // XXXX: This code is broken if the global node offset is nonzero
    // on any of the levels of the locally refined Cartesian grid.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        if (d_global_index_offset[ln] != 0)
        {
            TBOX_ERROR("This is broken --- please submit a bug report if you encounter this error.\n");
        }
    }

    // For now, we just register the data on MPI process 0.  This will
    // fail if the structure is too large to be stored in the memory
    // allocated to a single process.
    if (SAMRAI::tbox::MPI::getRank() == 0)
    {
        for (int ln = 0; ln < d_max_levels; ++ln)
        {
            for (unsigned j = 0; j < d_num_vertex[ln].size(); ++j)
            {
                silo_writer->registerMarkerCloud(
                    d_base_filename[ln][j] + "_vertices",
                    d_num_vertex[ln][j], d_vertex_offset[ln][j], ln);

                if (d_edge_map[ln][j].size() > 0)
                {
                    silo_writer->registerUnstructuredMesh(
                        d_base_filename[ln][j] + "_mesh",
                        d_edge_map[ln][j], ln);
                }
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
    // Loop over all patches in the specified level of the patch level
    // and count the number of local vertices.
    int local_node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of vertices whose initial locations will
        // be within the given patch.
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
    SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Set the global index offset.  This is equal to the number of
    // Lagrangian indices that have already been initialized on the
    // specified level.
    d_global_index_offset[level_number] = global_index_offset;

    // Loop over all patches in the specified level of the patch level
    // and initialize the local vertices.
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

        SAMRAI::tbox::Pointer<LNodeIndexData> index_data =
            patch->getPatchData(lag_node_index_idx);

        // Initialize the vertices whose initial locations will be
        // within the given patch.
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
            const vector<double> X = getVertexPosn(point_idx, level_number);

            // Initialize the location of the present vertex.
            double* const node_X = &(*X_data)(current_local_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                node_X[d] = X[d];
            }

            // Get the index of the cell in which the present vertex
            // is initially located.
            const SAMRAI::pdat::CellIndex<NDIM> idx = STOOLS::STOOLS_Utilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);

            // Initialize the force specification object assocaited
            // with the present vertex.
            std::vector<SAMRAI::tbox::Pointer<Stashable> > force_spec =
                initializeForceSpec(
                    point_idx, global_index_offset, level_number);

            if (!index_data->isElement(idx))
            {
                index_data->appendItem(idx,LNodeIndexSet());
            }
            LNodeIndexSet* node_set = index_data->getItem(idx);
            node_set->push_back(
                new LNodeIndex(current_global_idx, current_local_idx,
                               &(*X_data)(current_local_idx), force_spec));
        }
    }
    return local_node_count;
}// initializeDataOnPatchLevel

int
IBStandardInitializer::initializeMassDataOnPatchLevel(
    const int global_index_offset,
    const int local_index_offset,
    SAMRAI::tbox::Pointer<LNodeLevelData>& M_data,
    SAMRAI::tbox::Pointer<LNodeLevelData>& K_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Loop over all patches in the specified level of the patch level
    // and initialize the local vertices.
    int local_idx = -1;
    int local_node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

        // Initialize the vertices whose initial locations will be
        // within the given patch.
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number, can_be_refined);
        local_node_count += patch_vertices.size();
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int>& point_idx = (*it);
            const int current_local_idx = ++local_idx + local_index_offset;

            // Initialize the mass and stiffness coefficient
            // corresponding to the present vertex.
            const double M = getVertexMass(point_idx, level_number);
            const double K = getVertexMassStiffness(point_idx, level_number);

            // Avoid division by zero.
            if (SAMRAI::tbox::Utilities::deq(M,0.0))
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
    // Loop over all patches in the specified level of the patch level
    // and tag cells for refinement wherever there are vertices
    // assigned to a finer level of the Cartesian grid.
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

        // Tag cells for refinement whenever there are vertices whose
        // initial locations will be within the index space of the
        // given patch, but on the next finer level of the AMR patch
        // hierarchy.
        const bool can_be_refined = level_number+2 < d_max_levels;
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number+1, can_be_refined);
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int>& point_idx = (*it);

            // Get the coordinates of the present vertex.
            const vector<double> X = getVertexPosn(point_idx, level_number+1);

            // Get the index of the cell in which the present vertex
            // is initially located.
            const SAMRAI::pdat::CellIndex<NDIM> i = STOOLS::STOOLS_Utilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);

            // Tag the cell for refinement.
            if (patch_box.contains(i)) (*tag_data)(i) = 1;
        }
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

namespace
{

inline std::string
discard_comments(
    const std::string& input_string)
{
    // Create a copy of the input string, but without any text
    // following a '!', '#', or '%' character.
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

void
IBStandardInitializer::readVertexFiles()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = static_cast<int>(d_base_filename[ln].size());
        d_num_vertex[ln].resize(num_base_filename);
        d_vertex_posn[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string vertex_filename = d_base_filename[ln][j] + ".vertex";
            std::ifstream file_stream;
            std::string line_string;
            file_stream.open(vertex_filename.c_str(), std::ios::in);
            if (!file_stream.is_open()) TBOX_ERROR(d_object_name << ":\n  Unable to open input file " << vertex_filename << endl);

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "processing vertex data from input filename " << vertex_filename << endl
                               << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

            // The first entry in the file is the number of vertices.
            if (!std::getline(file_stream, line_string))
            {
                TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << vertex_filename << endl);
            }
            else
            {
                line_string = discard_comments(line_string);
                std::istringstream line_stream(line_string);
                if (!(line_stream >> d_num_vertex[ln][j]))
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << vertex_filename << endl);
                }
            }

            // Each successive line provides the initial position of
            // each vertex in the input file.
            d_vertex_posn[ln][j].resize(d_num_vertex[ln][j]*NDIM);
            for (int k = 0; k < d_num_vertex[ln][j]; ++k)
            {
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << vertex_filename << endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (!(line_stream >> d_vertex_posn[ln][j][k*NDIM+d]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << vertex_filename << endl);
                        }
                    }
                }
            }

            // Close the input file.
            file_stream.close();

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "read " << d_num_vertex[ln][j] << " vertices from input filename " << vertex_filename << endl
                               << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;
        }
    }
    return;
}// readVertexFiles

void
IBStandardInitializer::readEdgeFiles()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = static_cast<int>(d_base_filename[ln].size());
        d_edge_map[ln].resize(num_base_filename);
        d_edge_stiffness[ln].resize(num_base_filename);
        d_edge_rest_length[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string edge_filename = d_base_filename[ln][j] + ".edge";
            std::ifstream file_stream;
            std::string line_string;
            file_stream.open(edge_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing edge data from input filename " << edge_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

                // The first line in the file indicates the number of
                // edges in the input file.
                int num_edges;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << edge_filename << endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_edges))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << edge_filename << endl);
                    }
                }

                if (num_edges <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << edge_filename << endl);
                }

                // Each successive line provides the connectivity
                // information for each edge in the structure.
                for (int k = 0; k < num_edges; ++k)
                {
                    Edge e;
                    double kappa, length;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << edge_filename << endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> e.first))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl);
                        }
                        else if ((e.first < 0) || (e.first >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl
                                       << "  vertex index " << e.first << " is out of range" << endl);
                        }

                        if (!(line_stream >> e.second))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl);
                        }
                        else if ((e.second < 0) || (e.second >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl
                                       << "  vertex index " << e.second << " is out of range" << endl);
                        }

                        if (!(line_stream >> kappa))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl);
                        }
                        else if (kappa < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl
                                       << "  spring constant is negative" << endl);
                        }

                        if (!(line_stream >> length))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl);
                        }
                        else if (length < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << edge_filename << endl
                                       << "  spring resting length is negative" << endl);
                        }
                    }

                    // Modify kappa and length according to whether
                    // edge forces are enabled, or whether uniform
                    // values are to be employed, for this particular
                    // structure.
                    if (!d_enable_edges[ln][j])
                    {
                        kappa = 0.0;
                        length = 0.0;
                    }
                    else
                    {
                        if (d_use_uniform_edge_stiffness[ln][j])
                        {
                            kappa = d_uniform_edge_stiffness[ln][j];
                        }
                        if (d_use_uniform_edge_rest_length[ln][j])
                        {
                            length = d_uniform_edge_rest_length[ln][j];
                        }
                    }

                    // Correct the edge numbers to be in the global
                    // Lagrangian indexing scheme.
                    e.first  += d_vertex_offset[ln][j];
                    e.second += d_vertex_offset[ln][j];

                    // Always place the lower index first.
                    if (e.first > e.second)
                    {
                        const int tmp = e.first;
                        e.first = e.second;
                        e.second = tmp;
                    }

                    // Check to see if the edge has already been
                    // inserted in the edge map.
                    bool duplicate_edge = false;
                    for (std::multimap<int,Edge>::const_iterator it =
                             d_edge_map[ln][j].lower_bound(e.first);
                         it != d_edge_map[ln][j].upper_bound(e.first); ++it)
                    {
                        const Edge& other_e = (*it).second;
                        if (e.first  == other_e.first &&
                            e.second == other_e.second)
                        {
                            // This is a duplicate edge and should not
                            // be inserted into the edge map.
                            duplicate_edge = true;

                            // Ensure that the stiffness and rest
                            // length information is consistent.
                            if (!SAMRAI::tbox::Utilities::deq(
                                    (*d_edge_stiffness[ln][j].find(e)).second, kappa) ||
                                !SAMRAI::tbox::Utilities::deq(
                                    (*d_edge_rest_length[ln][j].find(e)).second, length))
                            {
                                TBOX_ERROR(d_object_name << ":\n  Inconsistent duplicate edges found in file " << edge_filename <<endl);
                            }
                        }
                    }

                    if (!duplicate_edge)
                    {
                        for (std::multimap<int,Edge>::const_iterator it =
                                 d_edge_map[ln][j].lower_bound(e.second);
                             it != d_edge_map[ln][j].upper_bound(e.second); ++it)
                        {
                            const Edge& other_e = (*it).second;
                            if (e.first  == other_e.first &&
                                e.second == other_e.second)
                            {
                                TBOX_ERROR(d_object_name << ":\n  Edge map is inconsistent.  Please contact the IBAMR developers." << endl);
                            }
                        }
                    }

                    // Initialize the edge map entries corresponding
                    // to the present edge.
                    //
                    // Note that the edge is associated with both the
                    // first and the second vertex in the edge map.
                    if (!duplicate_edge)
                    {
                        d_edge_map[ln][j].insert(std::make_pair(e.first ,e));
                        d_edge_map[ln][j].insert(std::make_pair(e.second,e));
                        d_edge_stiffness[ln][j][e] = kappa;
                        d_edge_rest_length[ln][j][e] = length;
                    }
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_edges << " edges from input filename " << edge_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;
            }
        }
    }
    return;
}// readEdgeFiles

void
IBStandardInitializer::readTargetPointFiles()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = static_cast<int>(d_base_filename[ln].size());
        d_target_stiffness[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            d_target_stiffness[ln][j].resize(d_num_vertex[ln][j], 0.0);

            const std::string target_point_stiffness_filename = d_base_filename[ln][j] + ".target_stiff";
            std::ifstream file_stream;
            std::string line_string;
            file_stream.open(target_point_stiffness_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing target point data from input filename " << target_point_stiffness_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

                // The first line in the file indicates the number of
                // target point stiffnesses in the input file.
                int num_target_stiffness;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << target_point_stiffness_filename << endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_target_stiffness))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << target_point_stiffness_filename << endl);
                    }
                }

                if (num_target_stiffness <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << target_point_stiffness_filename << endl);
                }

                // Each successive line indicates the vertex number
                // and spring constant associated with any target
                // points.
                for (int k = 0; k < num_target_stiffness; ++k)
                {
                    int n;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << target_point_stiffness_filename << endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> n))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << endl);
                        }
                        else if ((n < 0) || (n >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << endl
                                       << "  vertex index " << n << " is out of range" << endl);
                        }

                        if (!(line_stream >> d_target_stiffness[ln][j][n]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << endl);
                        }
                        else if (d_target_stiffness[ln][j][n] < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << target_point_stiffness_filename << endl
                                       << "  target point spring constant is negative" << endl);
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_target_stiffness << " target points from input filename " << target_point_stiffness_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;
            }

            // Modify the target point stiffness constant according to
            // whether target point penalty forces are enabled, or
            // whether uniform values are to be employed, for this
            // particular structure.
            if (!d_enable_target_points[ln][j])
            {
                d_target_stiffness[ln][j] = std::vector<double>(
                    d_num_vertex[ln][j], 0.0);
            }
            else
            {
                if (d_use_uniform_target_stiffness[ln][j])
                {
                    d_target_stiffness[ln][j] = std::vector<double>(
                        d_num_vertex[ln][j],
                        d_uniform_target_stiffness[ln][j]);
                }
            }
        }
    }
    return;
}// readTargetPointFiles

void
IBStandardInitializer::readBoundaryMassFiles()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = static_cast<int>(d_base_filename[ln].size());
        d_bdry_mass[ln].resize(num_base_filename);
        d_bdry_mass_stiffness[ln].resize(num_base_filename);
        for (int j = 0; j < num_base_filename; ++j)
        {
            d_bdry_mass[ln][j].resize(d_num_vertex[ln][j], 0.0);
            d_bdry_mass_stiffness[ln][j].resize(d_num_vertex[ln][j], 0.0);

            const std::string bdry_mass_filename = d_base_filename[ln][j] + ".bdry_mass";
            std::ifstream file_stream;
            std::string line_string;
            file_stream.open(bdry_mass_filename.c_str(), std::ios::in);
            if (file_stream.is_open())
            {
                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing boundary mass data from input filename " << bdry_mass_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

                // The first line in the file indicates the number of
                // massive IB points in the input file.
                int num_bdry_mass_pts;
                if (!std::getline(file_stream, line_string))
                {
                    TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line 1 of file " << bdry_mass_filename << endl);
                }
                else
                {
                    line_string = discard_comments(line_string);
                    std::istringstream line_stream(line_string);
                    if (!(line_stream >> num_bdry_mass_pts))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << bdry_mass_filename << endl);
                    }
                }

                if (num_bdry_mass_pts <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line 1 of file " << bdry_mass_filename << endl);
                }

                // Each successive line indicates the vertex number
                // and spring constant associated with any massive IB
                // points.
                for (int k = 0; k < num_bdry_mass_pts; ++k)
                {
                    int n;
                    if (!std::getline(file_stream, line_string))
                    {
                        TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered before line " << k+2 << " of file " << bdry_mass_filename << endl);
                    }
                    else
                    {
                        line_string = discard_comments(line_string);
                        std::istringstream line_stream(line_string);
                        if (!(line_stream >> n))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << endl);
                        }
                        else if ((n < 0) || (n >= d_num_vertex[ln][j]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << endl
                                       << "  vertex index " << n << " is out of range" << endl);
                        }

                        if (!(line_stream >> d_bdry_mass[ln][j][n]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << endl);
                        }
                        else if (d_bdry_mass[ln][j][n] < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << endl
                                       << "  boundary mass is negative" << endl);
                        }

                        if (!(line_stream >> d_bdry_mass_stiffness[ln][j][n]))
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << endl);
                        }
                        else if (d_bdry_mass_stiffness[ln][j][n] < 0.0)
                        {
                            TBOX_ERROR(d_object_name << ":\n  Invalid entry in input file encountered on line " << k+2 << " of file " << bdry_mass_filename << endl
                                       << "  boundary mass spring constant is negative" << endl);
                        }
                    }
                }

                // Close the input file.
                file_stream.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << num_bdry_mass_pts << " boundary mass points from input filename " << bdry_mass_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;
            }

            // Modify the boundary mass and boundary mass stiffness
            // constant according to whether boundary mass is enabled,
            // or whether uniform values are to be employed, for this
            // particular structure.
            if (!d_enable_bdry_mass[ln][j])
            {
                d_bdry_mass[ln][j] = std::vector<double>(
                    d_num_vertex[ln][j], 0.0);
                d_bdry_mass_stiffness[ln][j] = std::vector<double>(
                    d_num_vertex[ln][j], 0.0);
            }
            else
            {
                if (d_use_uniform_bdry_mass[ln][j])
                {
                    d_bdry_mass[ln][j] = std::vector<double>(
                        d_num_vertex[ln][j],
                        d_uniform_bdry_mass[ln][j]);
                }
                if (d_use_uniform_bdry_mass_stiffness[ln][j])
                {
                    d_bdry_mass_stiffness[ln][j] = std::vector<double>(
                        d_num_vertex[ln][j],
                        d_uniform_bdry_mass_stiffness[ln][j]);
                }
            }
        }
    }
    return;
}// readBoundaryMassFiles

void
IBStandardInitializer::getPatchVertices(
    std::vector<std::pair<int,int> >& patch_vertices,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const int level_number,
    const bool can_be_refined) const
{
    // Loop over all of the vertices to determine the indices of those
    // vertices within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will
    // work for now.
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

std::vector<SAMRAI::tbox::Pointer<Stashable> >
IBStandardInitializer::initializeForceSpec(
    const std::pair<int,int>& point_index,
    const int global_index_offset,
    const int level_number) const
{
    std::vector<SAMRAI::tbox::Pointer<Stashable> > force_spec;

    const int j = point_index.first;
    const int lag_index = getCannonicalLagrangianIndex(point_index, level_number);

    std::vector<int> dst_idxs;
    std::vector<double> stiffness, rest_length;
    for (std::multimap<int,Edge>::const_iterator it = d_edge_map[level_number][j].lower_bound(lag_index);
         it != d_edge_map[level_number][j].upper_bound(lag_index); ++it)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(lag_index == (*it).first);
#endif
        // The connectivity information.
        const Edge& e = (*it).second;
        if (e.first == lag_index)
        {
            dst_idxs.push_back(e.second+global_index_offset);
        }
        else
        {
            dst_idxs.push_back(e.first+global_index_offset);
        }

        // The material properties.
        stiffness  .push_back((*d_edge_stiffness  [level_number][j].find(e)).second);
        rest_length.push_back((*d_edge_rest_length[level_number][j].find(e)).second);
    }

    const std::vector<double> X_target = getVertexPosn(point_index, level_number);
    const double kappa_target = getVertexTargetStiffness(point_index, level_number);

    force_spec.push_back(new IBStandardForceSpec(
                             dst_idxs, stiffness, rest_length,
                             X_target, kappa_target));

    return force_spec;
}// initializeForceSpec

void
IBStandardInitializer::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    // Determine the (maximum) number of levels in the locally refined
    // grid.  Note that each piece of the Lagrangian structure must be
    // assigned to a particular level of the grid.
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
    d_base_filename.resize(d_max_levels);
    d_num_vertex.resize(d_max_levels);
    d_vertex_offset.resize(d_max_levels);
    d_vertex_posn.resize(d_max_levels);
    d_enable_edges.resize(d_max_levels);
    d_edge_map.resize(d_max_levels);
    d_edge_stiffness.resize(d_max_levels);
    d_edge_rest_length.resize(d_max_levels);
    d_use_uniform_edge_stiffness.resize(d_max_levels);
    d_uniform_edge_stiffness.resize(d_max_levels);
    d_use_uniform_edge_rest_length.resize(d_max_levels);
    d_uniform_edge_rest_length.resize(d_max_levels);
    d_enable_target_points.resize(d_max_levels);
    d_target_stiffness.resize(d_max_levels);
    d_use_uniform_target_stiffness.resize(d_max_levels);
    d_uniform_target_stiffness.resize(d_max_levels);
    d_enable_bdry_mass.resize(d_max_levels);
    d_bdry_mass.resize(d_max_levels);
    d_bdry_mass_stiffness.resize(d_max_levels);
    d_use_uniform_bdry_mass.resize(d_max_levels);
    d_uniform_bdry_mass.resize(d_max_levels);
    d_use_uniform_bdry_mass_stiffness.resize(d_max_levels);
    d_uniform_bdry_mass_stiffness.resize(d_max_levels);
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

    // Read in any sub-databases associated with the input file names.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = static_cast<int>(d_base_filename[ln].size());

        d_enable_edges[ln].resize(num_base_filename,true);

        d_use_uniform_edge_stiffness[ln].resize(num_base_filename,false);
        d_uniform_edge_stiffness[ln].resize(num_base_filename,-1.0);

        d_use_uniform_edge_rest_length[ln].resize(num_base_filename,false);
        d_uniform_edge_rest_length[ln].resize(num_base_filename,-1.0);

        d_enable_target_points[ln].resize(num_base_filename,true);

        d_use_uniform_target_stiffness[ln].resize(num_base_filename,false);
        d_uniform_target_stiffness[ln].resize(num_base_filename,-1.0);

        d_enable_bdry_mass[ln].resize(num_base_filename,true);

        d_use_uniform_bdry_mass[ln].resize(num_base_filename,false);
        d_uniform_bdry_mass[ln].resize(num_base_filename,-1.0);

        d_use_uniform_bdry_mass_stiffness[ln].resize(num_base_filename,false);
        d_uniform_bdry_mass_stiffness[ln].resize(num_base_filename,-1.0);

        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string& base_filename = d_base_filename[ln][j];
            if (db->isDatabase(base_filename))
            {
                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> sub_db =
                    db->getDatabase(base_filename);

                // Determine whether to enable or disable any
                // particular features.
                if (sub_db->keyExists("enable_edges"))
                {
                    d_enable_edges[ln][j] = sub_db->getBool("enable_edges");
                }
                if (sub_db->keyExists("enable_target_points"))
                {
                    d_enable_target_points[ln][j] = sub_db->getBool("enable_target_points");
                }
                if (sub_db->keyExists("enable_bdry_mass"))
                {
                    d_enable_bdry_mass[ln][j] = sub_db->getBool("enable_bdry_mass");
                }

                // Determine whether to use uniform values for any
                // particular structure attributes.
                if (sub_db->keyExists("uniform_edge_stiffness"))
                {
                    d_use_uniform_edge_stiffness[ln][j] = true;
                    d_uniform_edge_stiffness[ln][j] = sub_db->getDouble("uniform_edge_stiffness");

                    if (d_uniform_edge_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_edge_stiffness' in database " << base_filename << endl
                                   << "  spring constant is negative" << endl);
                    }
                }
                if (sub_db->keyExists("uniform_edge_rest_length"))
                {
                    d_use_uniform_edge_rest_length[ln][j] = true;
                    d_uniform_edge_rest_length[ln][j] = sub_db->getDouble("uniform_edge_rest_length");

                    if (d_uniform_edge_rest_length[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_edge_rest_length' in database " << base_filename << endl
                                   << "  spring resting length is negative" << endl);
                    }
                }

                if (sub_db->keyExists("uniform_target_stiffness"))
                {
                    d_use_uniform_target_stiffness[ln][j] = true;
                    d_uniform_target_stiffness[ln][j] = sub_db->getDouble("uniform_target_stiffness");

                    if (d_uniform_target_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_target_stiffness' in database " << base_filename << endl
                                   << "  target point spring constant is negative" << endl);
                    }
                }

                if (sub_db->keyExists("uniform_bdry_mass"))
                {
                    d_use_uniform_bdry_mass[ln][j] = true;
                    d_uniform_bdry_mass[ln][j] = sub_db->getDouble("uniform_bdry_mass");

                    if (d_uniform_bdry_mass[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_bdry_mass' in database " << base_filename << endl
                                   << "  boundary mass is negative" << endl);
                    }
                }
                if (sub_db->keyExists("uniform_bdry_mass_stiffness"))
                {
                    d_use_uniform_bdry_mass_stiffness[ln][j] = true;
                    d_uniform_bdry_mass_stiffness[ln][j] = sub_db->getDouble("uniform_bdry_mass_stiffness");

                    if (d_uniform_bdry_mass_stiffness[ln][j] < 0.0)
                    {
                        TBOX_ERROR(d_object_name << ":\n  Invalid entry for key `uniform_bdry_mass_stiffness' in database " << base_filename << endl
                                   << "  boundary mass spring constant is negative" << endl);
                    }
                }
            }
        }
    }

    // Output the names of the input files to be read along with
    // additional debugging information.
    SAMRAI::tbox::pout << d_object_name << ":  Reading from input files: " << endl;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filename = static_cast<int>(d_base_filename[ln].size());
        for (int j = 0; j < num_base_filename; ++j)
        {
            const std::string& base_filename = d_base_filename[ln][j];
            SAMRAI::tbox::pout << "  base filename: " << base_filename << endl
                               << "  assigned to level " << ln << " of the Cartesian grid patch hierarchy" << endl
                               << "     required files: " << base_filename << ".vertex" << endl
                               << "     optional files: " << base_filename << ".edge, " << base_filename << ".target_stiff, " << base_filename << ".bdry_mass" << endl;
            if (!d_enable_edges[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: edge forces are DISABLED for " << base_filename << endl;
            }
            else
            {
                if (d_use_uniform_edge_stiffness[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform edge stiffnesses are being employed for " << base_filename << endl
                                       << "        any stiffness information in file " << base_filename << ".vertex will be IGNORED" << endl;
                }
                if (d_use_uniform_edge_rest_length[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform edge resting lengths are being employed for " << base_filename << endl
                                       << "        any resting length information in file " << base_filename << ".vertex will be IGNORED" << endl;
                }
            }

            if (!d_enable_target_points[ln][j])
            {
                SAMRAI::tbox::pout << "  NOTE: target point penalty forces are DISABLED for " << base_filename << endl;
            }
            else
            {
                if (d_use_uniform_edge_stiffness[ln][j])
                {
                    SAMRAI::tbox::pout << "  NOTE: uniform target point stiffnesses are being employed for " << base_filename << endl
                                       << "        any target point stiffness information in file " << base_filename << ".vertex will be IGNORED" << endl;
                }
            }

            SAMRAI::tbox::pout << endl;
        }
    }

    return;
}// getFromDatabase

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardInitializer>;

//////////////////////////////////////////////////////////////////////////////
