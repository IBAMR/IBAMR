// Filename: XInit.C
// Created on 12 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <07.Oct.2006 23:21:08 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "XInit.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

#include <ibamr/SpringForceSpec.h>
#include <ibamr/LNodeIndexData.h>

// STOOLS INCLUDES
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

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
      d_tapered(true),
      d_num_nodes(256),
      d_num_layers(6),
      d_num_stacks(128),
      d_width(3.0/64.0),
      d_alpha(0.25),
      d_beta(0.25),
      d_l(0.75),
      d_stiffness(1.0),
      d_rest_length(0.0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!grid_geom.isNull());
#endif

    // Register the SpringForceSpec object with the StashableManager
    // class.
    SpringForceSpec::registerWithStashableManager();

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

        double X[NDIM];

#if (NDIM == 3)
        for (int s = 0; s < d_num_stacks; ++s)
#endif
        {
            for (int r = 0; r < d_num_layers; ++r)
            {
                for (int l = 0; l < d_num_nodes; ++l)
                {
#if (NDIM == 2)
                    getNodePosn(X,l,r);
#endif
#if (NDIM == 3)
                    getNodePosn(X,l,r,s);
#endif
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

        double X[NDIM];
        pdat::CellIndex<NDIM> idx;

#if (NDIM == 3)
        for (int s = 0; s < d_num_stacks; ++s)
#endif
        {
            for (int r = 0; r < d_num_layers; ++r)
            {
                double r_fac;
                if (d_tapered)
                {
                    const double theta = -0.5 +
                        (static_cast<double>(r)+0.5)/static_cast<double>(d_num_layers);
                    r_fac = 1.0 + cos(2.0*M_PI*theta);
                }
                else
                {
                    r_fac = 1.0;
                }

                for (int l = 0; l < d_num_nodes; ++l)
                {
#if (NDIM == 2)
                    getNodePosn(X,l,r);
#endif
#if (NDIM == 3)
                    getNodePosn(X,l,r,s);
#endif
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

                        const int current_lag_idx = l + r*d_num_nodes
#if (NDIM == 3)
                            + s*d_num_nodes*d_num_layers
#endif
                            ;
                        const int current_local_idx = ++local_idx;

                        vector<int> dst_idxs;
                        vector<double> stiffnesses, rest_lengths;

                        double X_prev[NDIM];
                        if (l > 0)
                        {
                            dst_idxs.push_back(l-1+r*d_num_nodes
#if (NDIM == 3)
                                               + s*d_num_nodes*d_num_layers
#endif
                                               );
#if (NDIM == 2)
                            getNodePosn(X_prev,l-1,r);
#endif
#if (NDIM == 3)
                            getNodePosn(X_prev,l-1,r,s);
#endif
                        }
                        else
                        {
                            dst_idxs.push_back(d_num_nodes-1+r*d_num_nodes
#if (NDIM == 3)
                                               + s*d_num_nodes*d_num_layers
#endif
                                               );
#if (NDIM == 2)
                            getNodePosn(X_prev,d_num_nodes-1,r);
#endif
#if (NDIM == 3)
                            getNodePosn(X_prev,d_num_nodes-1,r,s);
#endif
                        }

                        double l_prev_sq = 0.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            l_prev_sq += pow(X[d]-X_prev[d],2.0);
                        }

                        // scale stiffness by 1/ds = num_nodes
                        stiffnesses.push_back(
                            r_fac*
                            d_stiffness*static_cast<double>(d_num_nodes));
                        rest_lengths.push_back(sqrt(l_prev_sq));

                        double X_next[NDIM];
                        if (l < d_num_nodes-1)
                        {
                            dst_idxs.push_back(l+1+r*d_num_nodes
#if (NDIM == 3)
                                               + s*d_num_nodes*d_num_layers
#endif
                                               );
#if (NDIM == 2)
                            getNodePosn(X_next,l+1,r);
#endif
#if (NDIM == 3)
                            getNodePosn(X_next,l+1,r,s);
#endif
                        }
                        else
                        {
                            dst_idxs.push_back(0+r*d_num_nodes
#if (NDIM == 3)
                                               + s*d_num_nodes*d_num_layers
#endif
                                               );
#if (NDIM == 2)
                            getNodePosn(X_next,0,r);
#endif
#if (NDIM == 3)
                            getNodePosn(X_next,0,r,s);
#endif
                        }

                        double l_next_sq = 0.0;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            l_next_sq += pow(X[d]-X_next[d],2.0);
                        }

                        // scale stiffness by 1/ds = num_nodes
                        stiffnesses.push_back(
                            r_fac*
                            d_stiffness*static_cast<double>(d_num_nodes));
                        rest_lengths.push_back(sqrt(l_next_sq));

#if (NDIM == 3)
                        if (s < d_num_stacks-1 )
                        {
                            dst_idxs.push_back(l+r*d_num_nodes+(s+1)*d_num_nodes*d_num_layers);

                            stiffnesses.push_back(
                                r_fac*
                                d_stiffness*static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#if 0
                            // scale stiffness by 1/ds = num_stacks
                            stiffnesses.push_back(
                                2.0*r_fac*static_cast<double>(d_num_stacks)*
                                d_stiffness*static_cast<double>(d_num_stacks)/static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#endif
                        }
                        if (s > 0 )
                        {
                            dst_idxs.push_back(l+r*d_num_nodes+(s-1)*d_num_nodes*d_num_layers);

                            stiffnesses.push_back(
                                r_fac*
                                d_stiffness*static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#if 0
                            // scale stiffness by 1/ds = num_stacks
                            stiffnesses.push_back(
                                2.0*r_fac*static_cast<double>(d_num_stacks)*
                                d_stiffness*static_cast<double>(d_num_stacks)/static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#endif
                        }

                        if (r < d_num_layers-1 )
                        {
                            dst_idxs.push_back(l+(r+1)*d_num_nodes+s*d_num_nodes*d_num_layers);

                            stiffnesses.push_back(
                                r_fac*
                                d_stiffness*static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#if 0
                            // scale stiffness by 1/dr = num_layers
                            stiffnesses.push_back(
                                2.0*8.0*r_fac*static_cast<double>(d_num_layers)*
                                d_stiffness*static_cast<double>(d_num_layers)/static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#endif
                        }
                        if (r > 0 )
                        {
                            dst_idxs.push_back(l+(r-1)*d_num_nodes+s*d_num_nodes*d_num_layers);

                            stiffnesses.push_back(
                                r_fac*
                                d_stiffness*static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#if 0
                            // scale stiffness by 1/dr = num_layers
                            stiffnesses.push_back(
                                2.0*8.0*r_fac*static_cast<double>(d_num_layers)*
                                d_stiffness*static_cast<double>(d_num_layers)/static_cast<double>(d_num_nodes));
                            rest_lengths.push_back(0.0);
#endif
                        }
#endif

                        vector< tbox::Pointer<Stashable> > force_spec;
                        force_spec.push_back(
                            new SpringForceSpec(dst_idxs,stiffnesses,rest_lengths));

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

        double X[NDIM];
        pdat::CellIndex<NDIM> idx;

#if (NDIM == 3)
        for (int s = 0; s < d_num_stacks; ++s)
#endif
        {
            for (int r = 0; r < d_num_layers; ++r)
            {
                for (int l = 0; l < d_num_nodes; ++l)
                {
#if (NDIM == 2)
                    getNodePosn(X,l,r);
#endif
#if (NDIM == 3)
                    getNodePosn(X,l,r,s);
#endif
                    idx = STOOLS::STOOLS_Utilities::getCellIndex(
                        X, xLower, xUpper, dx, patch_lower, patch_upper);

                    if (patch_box.contains(idx))
                    {
                        (*tag_data)(idx) = 1;
                    }
                }
            }
        }
    }
    return;
}// tagCellsForInitialRefinement

/////////////////////////////// PRIVATE //////////////////////////////////////

void
XInit::getNodePosn(
    double* const X,
    const int l,
    const int r
#if (NDIM == 3)
    , const int s
#endif
                   ) const
{
    const double theta =
        2.0*M_PI*static_cast<double>(l)/static_cast<double>(d_num_nodes);

    const double dr =
        d_width*((static_cast<double>(r)+0.5)/
                 static_cast<double>(d_num_layers) - 0.5);
    X[0] = 0.5 + (d_alpha+dr)*cos(theta);
    X[1] = 0.5 + (d_beta +dr)*sin(theta);
#if (NDIM == 3)
    const double ds =
        8.0*d_width*((static_cast<double>(s)+0.5)/
                     static_cast<double>(d_num_stacks) - 0.5);
    X[2] = 0.5 + ds;
#endif

#if 0
#if (NDIM == 3)
    const double phi =
        4.0*M_PI*static_cast<double>(s)/static_cast<double>(d_num_stacks);

    X[0] = 0.5 + (d_alpha + (d_l+dr)*cos(theta))*cos(phi);
    X[1] = 0.5 + (d_beta  + (d_l+dr)*cos(theta))*sin(phi);
    X[2] = 0.5 + (          (d_l+dr)*sin(theta))         ;

//    X[0] = 0.5 + (d_alpha + (4.0*d_width+dr)*cos(theta))*cos(phi);
//    X[1] = 0.5 + (d_beta  + (4.0*d_width+dr)*cos(theta))*sin(phi);
//    X[2] = 0.5 + (4.0*d_width+dr)*sin(theta);

    // X[2] = 0.5 + d_l*((static_cast<double>(s)+0.5)/static_cast<double>(d_num_stacks)-0.5);
#endif
#endif
    return;
}// getNodePosn

void
XInit::getFromInput(
    tbox::Pointer<tbox::Database> db)
{
    if (!db.isNull())
    {
        d_tapered = db->getBoolWithDefault("tapered", d_tapered);
        d_num_nodes = db->getIntegerWithDefault("num_nodes", d_num_nodes);
        d_num_layers = db->getIntegerWithDefault("num_layers", d_num_layers);
        d_num_stacks = db->getIntegerWithDefault("num_stacks", d_num_stacks);
        d_width = db->getDoubleWithDefault("width", d_width);
        d_alpha = db->getDoubleWithDefault("alpha", d_alpha);
        d_beta = db->getDoubleWithDefault("beta", d_beta);
        d_l = db->getDoubleWithDefault("l", d_l);
        d_stiffness = db->getDoubleWithDefault("stiffness",d_stiffness);
        d_rest_length = db->getDoubleWithDefault("rest_length",d_rest_length);
    }
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////
