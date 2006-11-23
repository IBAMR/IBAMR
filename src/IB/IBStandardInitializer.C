// Filename: IBStandardInitializer.C
// Last modified: <23.Nov.2006 15:28:10 boyce@boyce-griffiths-powerbook-g4-15.local>
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
#include <ibamr/SpringForceSpec.h>
#include <ibamr/TargetPointForceSpec.h>
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardInitializer::IBStandardInitializer(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_num_vertices(),
      d_vertex_offsets(),
      d_vertex_posns()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
#endif

    // Register the various force specification objects with the
    // StashableManager class.
    SpringForceSpec::registerWithStashableManager();
    TargetPointForceSpec::registerWithStashableManager();

    // Get the input filenames.
    std::vector<std::string> base_filenames;
    if (input_db->keyExists("base_filenames"))
    {
        const int n_files = input_db->getArraySize("base_filenames");
        base_filenames.resize(n_files);
        input_db->getStringArray("base_filenames", &base_filenames[0], n_files);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `base_filenames' not found in input.");
    }

    // Output the names of the input files to be read.
    SAMRAI::tbox::pout << d_object_name << ":  Reading from input files: " << endl;
    for (std::vector<std::string>::const_iterator it = base_filenames.begin();
         it != base_filenames.end(); ++it)
    {
        SAMRAI::tbox::pout << "  base filename: " << *it << endl
                           << "     required files: " << *it << ".vertices" << endl
                           << "     optional files: " << *it << ".edges" << endl;
    }

    // Process the input files on each MPI process.
    for (int rank = 0; rank < SAMRAI::tbox::MPI::getNodes(); ++rank)
    {
        if (rank == SAMRAI::tbox::MPI::getRank())
        {
            readVertexFiles(base_filenames);
            //readEdgeFiles(base_filenames);
        }
        SAMRAI::tbox::MPI::barrier();
    }

    // Compute the index offsets.
    d_vertex_offsets.resize(d_num_vertices.size());
    d_vertex_offsets[0] = 0;
    for (int j = 1; j < static_cast<int>(d_num_vertices.size()); ++j)
    {
        d_vertex_offsets[j] = d_vertex_offsets[j-1]+d_num_vertices[j-1];
    }

    return;
}// IBStandardInitializer

IBStandardInitializer::~IBStandardInitializer()
{
    // intentionally blank
    return;
}// ~IBStandardInitializer

bool
IBStandardInitializer::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    // All curvilinear mesh nodes must reside on the finest level of
    // the AMR patch hierarchy.
    return !can_be_refined;
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
            const std::pair<int,int> point_idx = (*it);
            const int current_global_idx = getCannonicalLagrangianIndex(point_idx) + global_index_offset;
            const int current_local_idx = ++local_idx + local_index_offset;

            // Get the coordinates and stiffnesses of the present
            // vertex.
            const vector<double> X = getVertexPosn(point_idx);

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
            std::vector<SAMRAI::tbox::Pointer<Stashable> > force_spec;
            //force_spec.push_back(new TargetPointForceSpec(X, kappa));

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
IBStandardInitializer::tagCellsForInitialRefinement(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    // Loop over all patches in the specified level of the patch level
    // and tag cells for refinement whenever there are vertices.
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
        static const bool can_be_refined = false;
        std::vector<std::pair<int,int> > patch_vertices;
        getPatchVertices(patch_vertices, patch, level_number+1, can_be_refined);
        for (std::vector<std::pair<int,int> >::const_iterator it = patch_vertices.begin();
             it != patch_vertices.end(); ++it)
        {
            const std::pair<int,int> point_idx = (*it);

            // Get the coordinates of the present vertex.
            const vector<double> X = getVertexPosn(point_idx);

            // Get the index of the cell in which the present vertex
            // is initially located.
            const SAMRAI::pdat::CellIndex<NDIM> idx = STOOLS::STOOLS_Utilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);

            // Tag the cell for refinement.
            if (patch_box.contains(idx)) (*tag_data)(idx) = 1;
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
IBStandardInitializer::readVertexFiles(
    const std::vector<std::string>& base_filenames)
{
    const int num_base_filenames = static_cast<int>(base_filenames.size());
    d_num_vertices.resize(num_base_filenames);
    d_vertex_posns.resize(num_base_filenames);
    for (int j = 0; j < num_base_filenames; ++j)
    {
        const std::string vertex_filename = base_filenames[j] + ".vertices";
        std::ifstream is;
        is.open(vertex_filename.c_str(), std::ios::in);
        if (!is.is_open()) TBOX_ERROR(d_object_name << ":\n  Unable to open input file " << vertex_filename << endl);
        is.peek();  // try to catch premature end-of-file for empty files

        SAMRAI::tbox::plog << d_object_name << ":  "
                           << "processing vertex data from input filename " << vertex_filename << endl
                           << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

        // The first entry in the file is the number of vertices.
        if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line 1 of file " << vertex_filename << endl);
        is >> d_num_vertices[j];
        ignore_rest_of_line(is);

        // Each successive line provides the initial position of each
        // vertex in the input file.
        d_vertex_posns[j].resize(d_num_vertices[j]*NDIM);
        for (int k = 0; k < d_num_vertices[j]; ++k)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+2 << " of file " << vertex_filename << endl);
                is >> d_vertex_posns[j][k*NDIM+d];
            }
            ignore_rest_of_line(is);
        }

        SAMRAI::tbox::plog << d_object_name << ":  "
                           << "read " << d_num_vertices[j] << " vertices from input filename " << vertex_filename << endl
                           << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

        // Close the input file.
        is.close();
    }
    return;
}// readVertexFiles

#if 0
void
IBStandardInitializer::readEdgeFiles(
    const std::vector<std::string>& base_filenames)
{
    const int num_base_filenames = static_cast<int>(base_filenames.size());
    d_num_vertices.resize(num_base_filenames);
    d_vertex_posns.resize(num_base_filenames);
    for (int j = 0; j < num_base_filenames; ++j)
    {
        const std::string edge_filename = base_filenames[j] + ".vertices";
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
            if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line 1 of file " << edge_filename << endl);
            int base_index;
            is >> base_index;
            ignore_rest_of_line(is);

            if (base_index != 0 && base_index != 1)
            {
                TBOX_ERROR(d_object_name << ":\n  Line 1 of filename " << edge_filename << " is invalid" << endl);
            }

            // The second line in the file indicates the number of
            // edges in the input file.


            // Each successive line provides the connectivity
            // information for each edge in the structure.
            d_edge_posns[j].resize(d_num_vertices[j]*NDIM);
            for (int k = 0; k < d_num_vertices[j]; ++k)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    if (is.eof()) TBOX_ERROR(d_object_name << ":\n  Premature end to input file encountered on line " << k+2 << " of file " << edge_filename << endl);
                    is >> d_edge_posns[j][k*NDIM+d];
                }
                ignore_rest_of_line(is);
            }

            SAMRAI::tbox::plog << d_object_name << ":  "
                               << "read " << d_num_vertices[j] << " vertices from input filename " << edge_filename << endl
                               << "  on MPI process " << SAMRAI::tbox::MPI::getRank() << endl;

            // Close the input file.
            is.close();
        }
    }
    return;
}// readEdgeFiles
#endif

void
IBStandardInitializer::getPatchVertices(
    std::vector<std::pair<int,int> >& patch_vertices,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const int level_number,
    const bool can_be_refined) const
{
    // All curvilinear mesh nodes must reside on the finest level of
    // the AMR grid.  Consequently, if can_be_refied == true, then
    // there CANNOT be ANY points within the specified patch.
    if (can_be_refined) return;

    // Loop over all of the vertices to determine the indices of those
    // vertices within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will
    // work for now.
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch->getPatchGeometry();
    const double* const xLower = patch_geom->getXLower();
    const double* const xUpper = patch_geom->getXUpper();

    for (unsigned j = 0; j < d_num_vertices.size(); ++j)
    {
        for (int k = 0; k < d_num_vertices[j]; ++k)
        {
            const double* const X = &d_vertex_posns[j][k*NDIM];
            const bool patch_owns_node =
                ((  xLower[0] <= X[0])&&(X[0] < xUpper[0]))
#if (NDIM > 1)
                &&((xLower[1] <= X[1])&&(X[1] < xUpper[1]))
#if (NDIM > 2)
                &&((xLower[2] <= X[2])&&(X[2] < xUpper[2]))
#endif
#endif
                ;
            if (patch_owns_node) patch_vertices.push_back(std::pair<int,int>(j,k));
        }
    }

    return;
}// getPatchVertices

int
IBStandardInitializer::getCannonicalLagrangianIndex(
    const std::pair<int,int> point_index) const
{
    return d_vertex_offsets[point_index.first]+point_index.second;
}// getCannonicalLagrangianIndex

std::vector<double>
IBStandardInitializer::getVertexPosn(
    const std::pair<int,int> point_index) const
{
    return std::vector<double>(
        &d_vertex_posns[point_index.first][point_index.second*NDIM     ],
        &d_vertex_posns[point_index.first][point_index.second*NDIM+NDIM]);
}// getIBStandardPosn

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardInitializer>;

//////////////////////////////////////////////////////////////////////////////
