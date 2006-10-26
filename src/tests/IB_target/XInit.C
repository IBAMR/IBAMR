// Filename: XInit.C
// Created on 12 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <25.Oct.2006 18:32:10 boyce@bigboy.nyconnect.com>

#include "XInit.h"

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

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

XInit::XInit(
    const string& object_name,
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geom,
    tbox::Pointer<tbox::Database> input_db)
    : d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_num_nodes(256),
      d_center(NDIM),
      d_radius(0.25),
      d_stiffness(1000.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif

    if (NDIM != 2) TBOX_ERROR("only NDIM=2 is presently implemented!\n");

    // Register the TargetPointForceSpec object with the
    // StashableManager class.
    TargetPointForceSpec::registerWithStashableManager();

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    return;
}// XInit

XInit::~XInit()
{
    // intentionally blank
    return;
}// ~XInit

bool
XInit::getLevelHasLagrangianData(
    const int level_number,
    const bool can_be_refined) const
{
    return !can_be_refined;
}// getLevelHasLagrangianData

int
XInit::getLocalNodeCountOnPatchLevel(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    if (can_be_refined) return 0;

    int node_count = 0;

    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
        const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();

        std::vector<double> X(NDIM);
        for (int l = 0; l < d_num_nodes; ++l)
        {
            getNodePosn(X,l);
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
XInit::initializeDataOnPatchLevel(
    const int lag_node_index_idx,
    tbox::Pointer<LNodeLevelData>& X_data,
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time)
{
    if (can_be_refined) return;

    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    int local_idx = -1;

    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
        const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const hier::Box<NDIM>& patch_box = patch->getBox();
        const pdat::CellIndex<NDIM>& patch_lower = patch_box.lower();
        const pdat::CellIndex<NDIM>& patch_upper = patch_box.upper();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        tbox::Pointer<LNodeIndexData> index_data =
            patch->getPatchData(lag_node_index_idx);

        std::vector<double> X(NDIM);
        pdat::CellIndex<NDIM> idx;
        for (int l = 0; l < d_num_nodes; ++l)
        {
            getNodePosn(X,l);
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
                idx = STOOLS::STOOLS_Utilities::getCellIndex(
                    X, xLower, xUpper, dx, patch_lower, patch_upper);

                const int current_lag_idx = l;
                const int current_local_idx = ++local_idx;

                const double ds = 2.0*M_PI*d_radius/static_cast<double>(d_num_nodes);
                vector< tbox::Pointer<Stashable> > force_spec;
                force_spec.push_back(
                    new TargetPointForceSpec(X, ds*d_stiffness));

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
XInit::tagCellsForInitialRefinement(
    const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index)
{
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
        const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
        const hier::Box<NDIM>& patch_box = patch->getBox();
        const pdat::CellIndex<NDIM>& patch_lower = patch_box.lower();
        const pdat::CellIndex<NDIM>& patch_upper = patch_box.upper();
        const double* const xLower = patch_geom->getXLower();
        const double* const xUpper = patch_geom->getXUpper();
        const double* const dx = patch_geom->getDx();

        tbox::Pointer< pdat::CellData<NDIM,int> > tag_data = patch->getPatchData(tag_index);

        std::vector<double> X(NDIM);
        pdat::CellIndex<NDIM> idx;
        for (int l = 0; l < d_num_nodes; ++l)
        {
            getNodePosn(X,l);

            idx = STOOLS::STOOLS_Utilities::getCellIndex(
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

void
XInit::getNodePosn(
    std::vector<double>& X,
    const int l) const
{
    const double theta =
        2.0*M_PI*static_cast<double>(l)/static_cast<double>(d_num_nodes);

    X[0] = d_center[0] + d_radius*cos(theta);
    X[1] = d_center[1] + d_radius*sin(theta);

    return;
}// getNodePosn

void
XInit::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        if (db->keyExists("center"))
        {
            d_center = db->getDoubleArray("center");
        }
        else
        {
            TBOX_ERROR("key `center' not specifed in input file");
        }

        d_num_nodes = db->getIntegerWithDefault("num_nodes", d_num_nodes);
        d_radius = db->getDoubleWithDefault("radius", d_radius);
        d_stiffness = db->getDoubleWithDefault("stiffness", d_stiffness);
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
