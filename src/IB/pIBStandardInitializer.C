// Filename: pIBStandardInitializer.C
// Last modified: <18.Jan.2007 15:39:13 boyce@bigboy.nyconnect.com>
// Created on 18 Jan 2007 by Boyce Griffith (boyce@bigboy.nyconnect.com)

#include "pIBStandardInitializer.h"

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
#include <ibamr/SpringForceSpec.h>
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

// C++ STDLIB INCLUDES
#include <cassert>
#include <fstream>
#include <iostream>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

pIBStandardInitializer::pIBStandardInitializer(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_max_levels(-1),
      d_base_filenames(),
      d_num_vertices(),
      d_vertex_offsets(),
      d_vertex_posns(),
      d_vertex_masses(),
      d_vertex_anchor_tags(),
      d_edge_map(),
      d_edge_stiffnesses(),
      d_edge_rest_lengths(),
      d_global_index_offset()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
#endif

    // Register the force specification object with the
    // StashableManager class.
    SpringForceSpec::registerWithStashableManager();

    // Determine the (maximum) number of levels in the locally refined
    // grid.  Note that each piece of the Lagrangian structure must be
    // assigned to a particular level of the grid.
    if (input_db->keyExists("max_levels"))
    {
        d_max_levels = input_db->getInteger("max_levels");
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

    d_base_filenames.resize(d_max_levels);
    d_num_vertices.resize(d_max_levels);
    d_vertex_offsets.resize(d_max_levels);
    d_vertex_posns.resize(d_max_levels);
    d_vertex_masses.resize(d_max_levels);
    d_vertex_anchor_tags.resize(d_max_levels);
    d_edge_map.resize(d_max_levels);
    d_edge_stiffnesses.resize(d_max_levels);
    d_edge_rest_lengths.resize(d_max_levels);
    d_global_index_offset.resize(d_max_levels);

    // Determine the various input file names.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        std::ostringstream db_key_name_stream;
        db_key_name_stream << "base_filenames_" << ln;
        const std::string db_key_name = db_key_name_stream.str();
        if (input_db->keyExists(db_key_name))
        {
            const int n_files = input_db->getArraySize(db_key_name);
            d_base_filenames[ln].resize(n_files);
            input_db->getStringArray(
                db_key_name, &d_base_filenames[ln][0], n_files);
        }
        else
        {
            TBOX_WARNING(d_object_name << ":  "
                         << "Key data `" + db_key_name + "' not found in input.");
        }
    }

    // Output the names of the input files to be read.
    SAMRAI::tbox::pout << d_object_name << ":  Reading from input files: " << endl;
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        for (std::vector<std::string>::const_iterator it = d_base_filenames[ln].begin();
             it != d_base_filenames[ln].end(); ++it)
        {
            SAMRAI::tbox::pout << "  base filename: " << *it << endl
                               << "  assigned to level " << ln << " of the Cartesian grid patch hierarchy" << endl
                               << "     required files: " << *it << ".vertices" << endl
                               << "     optional files: " << *it << ".edges" << endl;
        }
    }

    // Process the vertex information on each MPI process.
    for (int rank = 0; rank < SAMRAI::tbox::MPI::getNodes(); ++rank)
    {
        if (rank == SAMRAI::tbox::MPI::getRank())
        {
            readVertexFiles();
        }
        SAMRAI::tbox::MPI::barrier();
    }

    // Compute the index offsets.
    //
    // NOTE: A separate, independent Lagrangian numbering scheme is
    // used on each level of the locally refined Cartesian grid.
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        d_vertex_offsets[ln].resize(d_num_vertices[ln].size());
        d_vertex_offsets[ln][0] = 0;
        for (int j = 1; j < static_cast<int>(d_num_vertices[ln].size()); ++j)
        {
            d_vertex_offsets[ln][j] = d_vertex_offsets[ln][j-1]+d_num_vertices[ln][j-1];
        }
    }

    // Read the edge information on each MPI process.
    for (int rank = 0; rank < SAMRAI::tbox::MPI::getNodes(); ++rank)
    {
        if (rank == SAMRAI::tbox::MPI::getRank())
        {
            readEdgeFiles();
        }
        SAMRAI::tbox::MPI::barrier();
    }

    return;
}// pIBStandardInitializer

pIBStandardInitializer::~pIBStandardInitializer()
{
    // intentionally blank
    return;
}// ~pIBStandardInitializer

void
pIBStandardInitializer::registerLagSiloDataWriter(
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
            for (unsigned j = 0; j < d_num_vertices[ln].size(); ++j)
            {
                silo_writer->registerMarkerCloud(
                    d_base_filenames[ln][j] + "_vertices",
                    d_num_vertices[ln][j], d_vertex_offsets[ln][j], ln);

                if (d_edge_map[ln][j].size() > 0)
                {
                    silo_writer->registerUnstructuredMesh(
                        d_base_filenames[ln][j] + "_mesh",
                        d_edge_map[ln][j], ln);
                }
            }
        }
    }

    return;
}// registerLagSiloDataWriter

bool
pIBStandardInitializer::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    return !d_num_vertices[level_number].empty();
}// getLevelHasLagrangianData

int
pIBStandardInitializer::getLocalNodeCountOnPatchLevel(
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
pIBStandardInitializer::initializeDataOnPatchLevel(
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

void
pIBStandardInitializer::tagCellsForInitialRefinement(
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
        // given patch, but on a finer level of the AMR patch
        // hierarchy.
        for (int ln = level_number+1; ln < d_max_levels; ++ln)
        {
            const bool can_be_refined = ln+1 < d_max_levels;
            std::vector<std::pair<int,int> > patch_vertices;
            getPatchVertices(patch_vertices, patch, level_number+1, can_be_refined);
            for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
                 it != patch_vertices.end(); ++it)
            {
                const std::pair<int,int>& point_idx = (*it);

                // Get the coordinates of the present vertex.
                const vector<double> X = getVertexPosn(point_idx, ln);

                // Get the index of the cell in which the present
                // vertex is initially located.
                const SAMRAI::pdat::CellIndex<NDIM> idx = STOOLS::STOOLS_Utilities::getCellIndex(
                    X, xLower, xUpper, dx, patch_lower, patch_upper);

                // Tag the cell for refinement.
                if (patch_box.contains(idx)) (*tag_data)(idx) = 1;
            }
        }
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

namespace
{

static const char DELIM = '\n';
static const int NBUF = 256;
static char buf[NBUF];

void
ignore_rest_of_line(
    std::ifstream& is)
{
    for (char c = is.peek(); c != DELIM && !is.eof(); )
    {
        is.get(buf, NBUF, DELIM);
        c = is.peek();
    }
    return;
}// ignore_rest_of_line

}

void
pIBStandardInitializer::readVertexFiles()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filenames = static_cast<int>(d_base_filenames[ln].size());
        d_num_vertices[ln].resize(num_base_filenames);
        d_vertex_posns[ln].resize(num_base_filenames);
        d_vertex_masses[ln].resize(num_base_filenames);
        d_vertex_anchor_tags[ln].resize(num_base_filenames);
        for (int j = 0; j < num_base_filenames; ++j)
        {
            const std::string vertex_filename = d_base_filenames[ln][j] + ".vertices";
            std::ifstream is;
            is.open(vertex_filename.c_str(), std::ios::in);
            if (!is.is_open()) TBOX_ERROR(d_object_name << ":\n  Unable to open input file " << vertex_filename << endl);
            is.peek();  // try to catch premature end-of-file for empty files

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "processing vertex data from input filename " << vertex_filename << endl
                               << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

            // The first entry in the file is the number of vertices.
            if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line 1 of file " << vertex_filename << endl);
            is >> d_num_vertices[ln][j];
            ignore_rest_of_line(is);

            // Each successive line provides the initial position of
            // each vertex in the input file, its mass, and whether it
            // is anchored.
            d_vertex_posns[ln][j].resize(d_num_vertices[ln][j]*NDIM);
            d_vertex_masses[ln][j].resize(d_num_vertices[ln][j]);
            d_vertex_anchor_tags[ln][j].resize(d_num_vertices[ln][j]);
            for (int k = 0; k < d_num_vertices[ln][j]; ++k)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+2 << " of file " << vertex_filename << endl);
                    is >> d_vertex_posns[ln][j][k*NDIM+d];
                }
                if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+2 << " of file " << vertex_filename << endl);
                is >> d_vertex_masses[ln][j][k];
                if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+2 << " of file " << vertex_filename << endl);
                is >> d_vertex_anchor_tags[ln][j][k];
                ignore_rest_of_line(is);

                // We shift the value of the mass slightly to avoid
                // division by zero in IBHierarchyIntegrator.
                if (SAMRAI::tbox::Utilities::deq(d_vertex_masses[ln][j][k],0.0))
                {
                    d_vertex_masses[ln][j][k] = std::numeric_limits<double>::epsilon();
                }
            }

            // Close the input file.
            is.close();

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "read " << d_num_vertices[ln][j] << " vertices from input filename " << vertex_filename << endl
                               << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;
        }
    }
    return;
}// readVertexFiles

void
pIBStandardInitializer::readEdgeFiles()
{
    for (int ln = 0; ln < d_max_levels; ++ln)
    {
        const int num_base_filenames = static_cast<int>(d_base_filenames[ln].size());
        d_edge_map[ln].resize(num_base_filenames);
        d_edge_stiffnesses[ln].resize(num_base_filenames);
        d_edge_rest_lengths[ln].resize(num_base_filenames);
        for (int j = 0; j < num_base_filenames; ++j)
        {
            const std::string edge_filename = d_base_filenames[ln][j] + ".edges";
            std::ifstream is;
            is.open(edge_filename.c_str(), std::ios::in);
            if (is.is_open())
            {
                is.peek();  // try to catch premature end-of-file for empty files

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "processing edge data from input filename " << edge_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

                // The first entry in the file indicates whether the
                // numbering is assumed to start from 0 or from 1.
                int base_index;
                if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line 1 of file " << edge_filename << endl);
                is >> base_index;
                ignore_rest_of_line(is);

                if (base_index != 0 && base_index != 1)
                {
                    TBOX_ERROR(d_object_name << ":\n  Line 1 of filename " << edge_filename << " is invalid" << endl);
                }

                // The second line in the file indicates the number of
                // edges in the input file.
                int num_edges;
                if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line 2 of file " << edge_filename << endl);
                is >> num_edges;
                ignore_rest_of_line(is);

                if (num_edges <= 0)
                {
                    TBOX_ERROR(d_object_name << ":\n  Line 2 of filename " << edge_filename << " is invalid" << endl);
                }

                // Each successive line provides the connectivity
                // information for each edge in the structure.
                for (int k = 0; k < num_edges; ++k)
                {
                    Edge e;
                    double kappa, length;
                    if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+3 << " of file " << edge_filename << endl);
                    is >> e.first;   // first  vertex
                    if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+3 << " of file " << edge_filename << endl);
                    is >> e.second;  // second vertex
                    if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+3 << " of file " << edge_filename << endl);
                    is >> kappa;     // spring constant
                    if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+3 << " of file " << edge_filename << endl);
                    is >> length;    // rest length
                    ignore_rest_of_line(is);

                    e.first  += d_vertex_offsets[ln][j]-base_index; // correct the edge numbers to be in the
                    e.second += d_vertex_offsets[ln][j]-base_index; // global Lagrangian indexing scheme

                    // Always place the lower index first.
                    if (e.first > e.second)
                    {
                        const int tmp = e.first;
                        e.first = e.second;
                        e.second = tmp;
                    }

                    d_edge_map[ln][j].insert(std::make_pair(e.first ,e));
                    d_edge_map[ln][j].insert(std::make_pair(e.second,e));
                    d_edge_stiffnesses [ln][j][e] = kappa;
                    d_edge_rest_lengths[ln][j][e] = length;
                }

                // Close the input file.
                is.close();

                SAMRAI::tbox::plog << d_object_name << ":  "
                                   << "read " << d_num_vertices[ln][j] << " vertices from input filename " << edge_filename << endl
                                   << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;
            }
        }
    }
    return;
}// readEdgeFiles

void
pIBStandardInitializer::getPatchVertices(
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

    for (unsigned j = 0; j < d_num_vertices[level_number].size(); ++j)
    {
        for (int k = 0; k < d_num_vertices[level_number][j]; ++k)
        {
            const double* const X = &d_vertex_posns[level_number][j][k*NDIM];
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
pIBStandardInitializer::getCannonicalLagrangianIndex(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_vertex_offsets[level_number][point_index.first]+point_index.second;
}// getCannonicalLagrangianIndex

std::vector<double>
pIBStandardInitializer::getVertexPosn(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return std::vector<double>(
        &d_vertex_posns[level_number][point_index.first][point_index.second*NDIM     ],
        &d_vertex_posns[level_number][point_index.first][point_index.second*NDIM+NDIM]);
}// getVertexPosn

double
pIBStandardInitializer::getVertexMass(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return d_vertex_masses[level_number][point_index.first][point_index.second];
}// getVertexMass

bool
pIBStandardInitializer::getVertexAnchorTag(
    const std::pair<int,int>& point_index,
    const int level_number) const
{
    return (d_vertex_anchor_tags[level_number][point_index.first][point_index.second]
            == 1);
}// getVertexAnchorTag

std::vector<SAMRAI::tbox::Pointer<Stashable> >
pIBStandardInitializer::initializeForceSpec(
    const std::pair<int,int>& point_index,
    const int global_index_offset,
    const int level_number) const
{
    std::vector<SAMRAI::tbox::Pointer<Stashable> > force_spec;

    const int j = point_index.first;
    const int lag_index = getCannonicalLagrangianIndex(point_index, level_number);

    bool has_edges = false;
    std::vector<int> dst_idxs;
    std::vector<double> stiffnesses, rest_lengths;
    for (std::multimap<int,Edge>::const_iterator it = d_edge_map[level_number][j].lower_bound(lag_index);
         it != d_edge_map[level_number][j].upper_bound(lag_index); ++it)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(lag_index == (*it).first);
#endif
        has_edges = true;

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
        stiffnesses .push_back((*d_edge_stiffnesses [level_number][j].find(e)).second);
        rest_lengths.push_back((*d_edge_rest_lengths[level_number][j].find(e)).second);
    }

    // Don't bother to create a force spec if there are no edges.
    if (has_edges)
    {
        force_spec.push_back(new SpringForceSpec(dst_idxs, stiffnesses, rest_lengths));
    }

    return force_spec;
}// initializeForceSpec

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::pIBStandardInitializer>;

//////////////////////////////////////////////////////////////////////////////
