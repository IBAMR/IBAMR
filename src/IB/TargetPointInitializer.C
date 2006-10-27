// Filename: TargetPointInitializer.C
// Last modified: <27.Oct.2006 01:16:30 boyce@bigboy.nyconnect.com>
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

// C++ STDLIB INCLUDES
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
    // All target points must reside on the finest level of the AMR
    // grid.  Consequently, if can_be_refied == true, then the local
    // node count is 0.
    if (can_be_refined) return 0;

    // Loop over all patches in the specified level of the patch level
    // and count the number of local target points.
    int node_count = 0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();

        // Loop over all target points and increment the local node
        // counter whenever the initial position of a target point
        // lies within the present patch.
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

            if (patch_owns_node)
            {
                ++node_count;
            }
        }
    }

    return node_count;
}// getLocalNodeCountOnPatchLevel

void
TargetPointInitializer::initializeDataOnPatchLevel(
    const int lag_node_index_idx,
    SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    // All target points must reside on the finest level of the AMR
    // grid.  Consequently, if can_be_refied == true, then there is no
    // data to initialize.
    if (can_be_refined) return;

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

        // Loop over all target points and initialize all local target
        // points.
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

            if (patch_owns_node)
            {
                const SAMRAI::pdat::CellIndex<NDIM> idx = STOOLS::STOOLS_Utilities::getCellIndex(
                    X, xLower, xUpper, dx, patch_lower, patch_upper);

                const int current_lag_idx = k;
                const int current_local_idx = ++local_idx;

                std::vector<SAMRAI::tbox::Pointer<Stashable> > force_spec;
                force_spec.push_back(
                    new TargetPointForceSpec(X, d_stiffnesses[k]));

                if (!index_data->isElement(idx))
                {
                    index_data->appendItem(idx,LNodeIndexSet());
                }
                LNodeIndexSet* node_set = index_data->getItem(idx);
                node_set->push_back(
                    new LNodeIndex(current_lag_idx,
                                   current_local_idx,
                                   &(*X_data)(current_local_idx),
                                   force_spec));

                double* const node_X = &(*X_data)(current_local_idx);
                for (int d = 0; d < NDIM; ++d)
                {
                    node_X[d] = X[d];
                }
            }
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

        for (int k = 0; k < d_num_points; ++k)
        {
            const double* const X = &d_initial_posns[k*NDIM];
            const SAMRAI::pdat::CellIndex<NDIM> idx = STOOLS::STOOLS_Utilities::getCellIndex(
                X, xLower, xUpper, dx, patch_lower, patch_upper);
            if (patch_box.contains(idx))
            {
                (*tag_data)(idx) = 1;
            }
        }
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::TargetPointInitializer>;

//////////////////////////////////////////////////////////////////////////////
