// Filename: IBInstrumentPanel.C
// Last modified: <11.Jun.2007 19:32:07 griffith@box221.cims.nyu.edu>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

#include "IBInstrumentPanel.h"

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
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/LNodeIndexData2.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentPanel::IBInstrumentPanel(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_num_meters(0),
      d_num_meter_nodes(),
      d_X_perimeter(),
      d_X_centroid(),
      d_X_web()
{
    if (!input_db.isNull())
    {

    }
    return;
}// IBInstrumentPanel

IBInstrumentPanel::~IBInstrumentPanel()
{
    // intentionally blank
    return;
}// ~IBInstrumentPanel

void
IBInstrumentPanel::computeMeterPositions(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Loop over all local nodes to determine the locations of the meters.
    d_X_perimeter.clear();
    d_X_perimeter.resize(d_num_meters);
    for (int k = 0; k < d_num_meters; ++k)
    {
        d_X_perimeter[k].resize(NDIM*d_num_meter_nodes[k], 0.0);
    }
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (lag_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                    patch->getPatchData(lag_node_index_idx);

                for (LNodeIndexData2::Iterator it(patch_box); it; it++)
                {
                    const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                    const LNodeIndexSet& node_set = (*idx_data)(i);
                    for (LNodeIndexSet::const_iterator n = node_set.begin();
                         n != node_set.end(); ++n)
                    {
                        const LNodeIndexSet::value_type& node_idx = *n;
                        const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data =
                            node_idx->getStashData();
                        for (unsigned l = 0; l < stash_data.size(); ++l)
                        {
                            SAMRAI::tbox::Pointer<IBInstrumentationSpec> spec = stash_data[l];
                            if (!spec.isNull())
                            {
                                const int meter_index = spec->getMeterIndex();
                                const int node_index = spec->getNodeIndex();
                                copy(node_idx->getNodeLocation(),
                                     node_idx->getNodeLocation()+NDIM,
                                     &d_X_perimeter[meter_index][NDIM*node_index]);
                            }
                        }
                    }
                }
            }
        }
    }

    // Set the positions of all meter nodes on all processes.
    for (int k = 0; k < d_num_meters; ++k)
    {
        SAMRAI::tbox::MPI::sumReduction(&d_X_perimeter[k][0], NDIM*d_num_meter_nodes[k]);
    }

    // Determine the centroid of each meter.
    d_X_centroid.clear();
    d_X_centroid.resize(d_num_meters);
    for (int k = 0; k < d_num_meters; ++k)
    {

    }

    // Find the patch and index corresponding to the finest cell in the patch
    // hierarchy that contains the present position of each pressure gauge.
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++k)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
            hierarchy->getPatchLevel(ln);

        std::vector<SAMRAI::hier::Index<NDIM> > pressure_gauge_index(d_num_pressure_gauges);
        for (int k = 0; k < d_num_pressure_gauges; ++k)
        {
            const std::vector<double>& X = d_pressure_gauge_loc[k];
            pressure_gauge_index[k] = compute_pressure_gauge_index(X);
        }

        const SAMRAI::hier::BoxArray<NDIM>& box_array = level->getBoxes();
        std::vector<bool> found_pressure_gauge_box(d_num_pressure_gauges, false);
        for (int p = 0; p < box_array.size(); ++p)
        {
            const SAMRAI::hier::Box<NDIM>& box = box_array[p];
            for (int k = 0; k < d_num_pressure_gauges; ++k)
            {
                if (!found_pressure_gauge_box[k] && box.contains(pressure_gauge_index[k]))
                {
                    d_pressure_gauge_level_num[k] = ln;
                    d_pressure_gauge_patch_num[k] = p;
                    d_pressure_gauge_on_proc[k] = (rank == level->getMappingForPatch(p));
                    found_pressure_gauge_box[k] = true;
                }
            }
        }
    }
#endif
    return;
}// computeMeterPositions

void
IBInstrumentPanel::initializeHierarchyData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Determine how many flow meters/pressure gauges are present.
    int max_meter_index = 0;
    std::vector<int> max_node_index;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        if (lag_manager->levelContainsLagrangianData(ln))
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                hierarchy->getPatchLevel(ln);
            for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
                const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
                const SAMRAI::tbox::Pointer<LNodeIndexData2> idx_data =
                    patch->getPatchData(lag_node_index_idx);

                for (LNodeIndexData2::Iterator it(patch_box); it; it++)
                {
                    const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                    const LNodeIndexSet& node_set = (*idx_data)(i);
                    for (LNodeIndexSet::const_iterator n = node_set.begin();
                         n != node_set.end(); ++n)
                    {
                        const LNodeIndexSet::value_type& node_idx = *n;
                        const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data =
                            node_idx->getStashData();
                        for (unsigned l = 0; l < stash_data.size(); ++l)
                        {
                            SAMRAI::tbox::Pointer<IBInstrumentationSpec> spec = stash_data[l];
                            if (!spec.isNull())
                            {
                                int meter_index = spec->getMeterIndex();
                                max_meter_index = max(meter_index, max_meter_index);

                                int node_index = spec->getNodeIndex();
                                max_node_index.resize(max_meter_index+1,-1);
                                max_node_index[meter_index] = max(node_index, max_node_index[meter_index]);
                            }
                        }
                    }
                }
            }
        }
    }

    d_num_meters = SAMRAI::tbox::MPI::maxReduction(max_meter_index)+1;
    max_node_index.resize(d_num_meters,-1);

    d_num_meter_nodes.clear();
    d_num_meter_nodes.resize(d_num_meters,0);
    for (int k = 0; k < d_num_meters; ++k)
    {
        d_num_meter_nodes[k] = max_node_index[k]+1;
    }
    SAMRAI::tbox::MPI::maxReduction(&d_num_meter_nodes[0], d_num_meters);
    return;
}// initializeHierarchyData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBInstrumentPanel>;

//////////////////////////////////////////////////////////////////////////////
