// Filename: TargetPointInitializer.C
// Last modified: <09.Nov.2006 01:25:06 boyce@bigboy.nyconnect.com>
// Created on 26 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)

#include "TargetPointInitializer.h"

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

namespace
{
static const int INPUT_MPI_ROOT = 0;
static const char DELIM = '\n';
static const int NBUF = 256;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

TargetPointInitializer::TargetPointInitializer(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name),
      d_num_points(0),
      d_initial_posns(),
      d_stiffnesses()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
#endif

    // Register the TargetPointForceSpec object with the
    // StashableManager class.
    TargetPointForceSpec::registerWithStashableManager();

    // Get the input filename.
    string input_filename;
    if (input_db->keyExists("input_filename"))
    {
        input_filename = input_db->getString("input_filename");
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Key data `input_filename' not found in input.");
    }

    // Process the input file.
    //
    // The following code is not very flexible or fault tolerant; it
    // could be improved dramatically.
    if (SAMRAI::tbox::MPI::getRank() == INPUT_MPI_ROOT)
    {
        char buf[NBUF];
        std::ifstream is;
        is.open(input_filename.c_str(), ios::in);

        // The first entry in the file must be the number of target
        // points.
        is >> d_num_points;

        for (char c = is.peek(); c != DELIM; )
        {
            is.get(buf, NBUF, DELIM);
            c = is.peek();
        }

        // Read in each target point location followed by the
        // associated stiffness.
        d_initial_posns.resize(d_num_points*NDIM);
        d_stiffnesses.resize(d_num_points);

        for (int k = 0; k < d_num_points; ++k)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                is >> d_initial_posns[k*NDIM+d];
            }
            is >> d_stiffnesses[k];

            for (char c = is.peek(); c != DELIM; )
            {
                is.get(buf, NBUF, DELIM);
                c = is.peek();
            }
        }
    }

    // Broadcast the input file to the remaining processors.
    if (SAMRAI::tbox::MPI::getNodes() > 1)
    {
        d_num_points = SAMRAI::tbox::MPI::bcast(d_num_points, INPUT_MPI_ROOT);

        if (SAMRAI::tbox::MPI::getRank() != INPUT_MPI_ROOT)
        {
            d_initial_posns.resize(d_num_points*NDIM);
            d_stiffnesses.resize(d_num_points);
        }

        MPI_Bcast(static_cast<void*>(&d_initial_posns[0]),
                  d_num_points*NDIM, MPI_DOUBLE, INPUT_MPI_ROOT,
                  SAMRAI::tbox::MPI::getCommunicator());
        MPI_Bcast(static_cast<void*>(&d_stiffnesses[0]),
                  d_num_points, MPI_DOUBLE, INPUT_MPI_ROOT,
                  SAMRAI::tbox::MPI::getCommunicator());

        if (SAMRAI::tbox::MPI::getRank() == INPUT_MPI_ROOT)
        {
            SAMRAI::tbox::MPI::updateOutgoingStatistics(
                2, (NDIM+1)*d_num_points);
        }
        else
        {
            SAMRAI::tbox::MPI::updateIncomingStatistics(
                2, (NDIM+1)*d_num_points);
        }
    }

    return;
}// TargetPointInitializer

TargetPointInitializer::~TargetPointInitializer()
{
    // intentionally blank
    return;
}// ~TargetPointInitializer

bool
TargetPointInitializer::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    // All target points must reside on the finest level of the AMR
    // grid.
    return !can_be_refined;
}// getLevelHasLagrangianData

int
TargetPointInitializer::getLocalNodeCountOnPatchLevel(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // Loop over all patches in the specified level of the patch level
    // and count the number of local target points.
    int node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

        // Count the number of target points whose initial locations
        // will be within the given patch.
        std::vector<int> point_indices;
        getPatchTargetPointIndices(point_indices, patch,
                                   level_number, can_be_refined);
        node_count += point_indices.size();
    }

    return node_count;
}// getLocalNodeCountOnPatchLevel

void
TargetPointInitializer::initializeDataOnPatchLevel(
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
    // and initialize the local target points.
    int local_idx = -1;
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

        // Initialize the target points whose initial locations will
        // be within the given patch.
        std::vector<int> point_indices;
        getPatchTargetPointIndices(
            point_indices, patch, level_number, can_be_refined);
        for (std::vector<int>::const_iterator it = point_indices.begin();
             it != point_indices.end(); ++it)
        {
            const int point_idx = (*it);
            const int current_global_idx = point_idx + global_index_offset;
            const int current_local_idx = ++local_idx + local_index_offset;

            // Get the coordinates and stiffnesses of the present
            // target point.
            const vector<double> X = getTargetPointPosn(point_idx);
            const double kappa = getTargetPointStiffness(point_idx);

            // Initialize the location of the present target point.
            double* const node_X = &(*X_data)(current_local_idx);
            for (int d = 0; d < NDIM; ++d)
            {
                node_X[d] = X[d];
            }

            // Get the index of the cell in which the present target
            // point is initially located.
            const SAMRAI::pdat::CellIndex<NDIM> idx = STOOLS::STOOLS_Utilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);

            // Initialize the force specification object assocaited
            // with the present target point.
            std::vector<SAMRAI::tbox::Pointer<Stashable> > force_spec;
            force_spec.push_back(new TargetPointForceSpec(X, kappa));

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
    return;
}// initializeDataOnPatchLevel

void
TargetPointInitializer::tagCellsForInitialRefinement(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    // Loop over all patches in the specified level of the patch level
    // and tag cells for refinement whenever there are target points.
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

        // Initialize the target points whose initial locations will
        // be within the given patch.
        //
        // NOTE: Here, we want to tag the locations of all target
        // points that are to be be assigned to any finer levels in
        // the Cartesian grid.
        static const bool can_be_refined = false;
        std::vector<int> point_indices;
        getPatchTargetPointIndices(
            point_indices, patch, level_number+1, can_be_refined);
        for (std::vector<int>::const_iterator it = point_indices.begin();
             it != point_indices.end(); ++it)
        {
            const int point_idx = (*it);

            // Get the coordinates of the present target point.
            const vector<double> X = getTargetPointPosn(point_idx);

            // Get the index of the cell in which the present target
            // point is initially located.
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

void
TargetPointInitializer::getPatchTargetPointIndices(
    std::vector<int>& point_indices,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
    const int level_number,
    const bool can_be_refined) const
{
    // All target points must reside on the finest level of the AMR
    // grid.  Consequently, if can_be_refied == true, then there must
    // be no points within the specified patch.
    if (can_be_refined) return;

    // Loop over all of the target points to determine the indices of
    // those points within the present patch.
    //
    // NOTE: This is clearly not the best way to do this, but it will
    // work for now.
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch->getPatchGeometry();
    const double* const xLower = patch_geom->getXLower();
    const double* const xUpper = patch_geom->getXUpper();

    for (int k = 0; k < d_num_points; ++k)
    {
        const double* const X = &d_initial_posns[k*NDIM];
        const bool patch_owns_node =
            ((  xLower[0] <= X[0])&&(X[0] < xUpper[0]))
#if (NDIM > 1)
            &&((xLower[1] <= X[1])&&(X[1] < xUpper[1]))
#if (NDIM > 2)
            &&((xLower[2] <= X[2])&&(X[2] < xUpper[2]))
#endif
#endif
            ;
        if (patch_owns_node) point_indices.push_back(k);
    }

    return;
}// getPatchTargetPointIndices

std::vector<double>
TargetPointInitializer::getTargetPointPosn(
    const int point_index) const
{
    return std::vector<double>(
        &d_initial_posns[point_index*NDIM], &d_initial_posns[(point_index+1)*NDIM]);
}// getTargetPointPosn

double
TargetPointInitializer::getTargetPointStiffness(
    const int point_index) const
{
    return d_stiffnesses[point_index];
}// getTargetPointStiffness

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::TargetPointInitializer>;

//////////////////////////////////////////////////////////////////////////////
