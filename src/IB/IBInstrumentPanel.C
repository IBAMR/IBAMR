// Filename: IBPlumbingToolkit.C
// Last modified: <11.Jun.2007 16:57:54 griffith@box221.cims.nyu.edu>
// Created on 12 May 2007 by Boyce Griffith (boyce@trasnaform2.local)

#include "IBPlumbingToolkit.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBPlumbingToolkit::IBPlumbingToolkit(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_interp_type("LINEAR"),
      d_num_flow_meters(0),
      d_num_pressure_gauges(0)
{
    if (!input_db.isNull())
    {

    }
    return;
}// IBPlumbingToolkit

IBPlumbingToolkit::~IBPlumbingToolkit()
{
    // intentionally blank
    return;
}// ~IBPlumbingToolkit

void
IBPlumbingToolkit::initializeHierarchyData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // Initialize the data associated with each flow meter.
    initializeFlowMeter(hierarchy, init_data_time, initial_time, lag_manager);

    // Initialize the data associated with each pressure gauge.
    initializePressureGauges(hierarchy, init_data_time, initial_time, lag_manager);
    return;
}// initializeHierarchyData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBPlumbingToolkit::initializeFlowMeter(
    const int meter_number)
{
    return;
}// initializeFlowMeter

void
IBPlumbingToolkit::initializePressureGauges(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->
        getLNodeIndexPatchDescriptorIndex();

    // At the initial time, determine the number of pressure gauges, and the
    // number of points associated with each of the pressure gauges.
    if (!d_initialized_pressure_data)
    {
        int max_gauge_idx = 0;
        std::vector<int> max_point_idx;
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
                    const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data =
                        patch->getPatchData(lag_node_index_idx);

                    for (LNodeIndexData::Iterator it(*idx_data); it; it++)
                    {
                        if (patch_box.contains(it.getIndex()))
                        {
                            const LNodeIndexSet& node_set = *it;
                            for (LNodeIndexSet::const_iterator n = node_set.begin();
                                 n != node_set.end(); ++n)
                            {
                                const LNodeIndexSet::value_type& node_idx = *n;
                                const int& mastr_idx = node_idx->getLagrangianIndex();
                                const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data =
                                    node_idx->getStashData();
                                for (unsigned l = 0; l < stash_data.size(); ++l)
                                {
                                    SAMRAI::tbox::Pointer<IBInstrumentationSpec> spec = stash_data[l];
                                    if (!spec.isNull() && spec->getPressureGaugeIndex() != -1)
                                    {
                                        int gauge_idx = spec->getPressureGaugeIndex();
                                        max_gauge_idx = max(gauge_idx, max_gauge_idx);

                                        int point_idx = spec->getPointIndex();
                                        max_point_idx.resize(max_gauge_idx+1,-1);
                                        max_point_idx[gauge_idx] = max(point_idx, max_point_idx[gauge_idx]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        d_num_pressure_gauges = SAMRAI::tbox::MPI::maxReduction(max_gauge_idx)+1;
        max_point_idx.resize(d_num_pressure_gauges,-1);
        d_num_pressure_gauge_points.resize(d_num_pressure_gauges,0);
        for (int k = 0; k < d_num_pressure_gauges; ++k)
        {
            d_num_pressure_gauge_points[k] = max_point_idx[k]+1;
        }
        SAMRAI::tbox::MPI::maxReduction(&d_num_pressure_gauge_points[0], d_num_pressure_gauge_points);
    }

    // Loop over all local points to determine the locations of the pressure
    // gauges.
    std::vector<std::vector<double> > X(d_num_pressure_gauges);
    for (int k = 0; k < d_num_pressure_gauges; ++k)
    {
        X[k].resize(NDIM*d_num_pressure_gauge_points[k], 0.0);
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
                const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data =
                    patch->getPatchData(lag_node_index_idx);

                for (LNodeIndexData::Iterator it(*idx_data); it; it++)
                {
                    if (patch_box.contains(it.getIndex()))
                    {
                        const LNodeIndexSet& node_set = *it;
                        for (LNodeIndexSet::const_iterator n = node_set.begin();
                             n != node_set.end(); ++n)
                        {
                            const LNodeIndexSet::value_type& node_idx = *n;
                            const int& mastr_idx = node_idx->getLagrangianIndex();
                            const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data =
                                node_idx->getStashData();
                            for (unsigned l = 0; l < stash_data.size(); ++l)
                            {
                                SAMRAI::tbox::Pointer<IBPressureGaugeSpec> spec = stash_data[l];
                                if (!spec.isNull())
                                {
                                    const int gauge_idx = spec->getPressureGaugeIndex();
                                    const int point_idx = spec->getPointIndex();
                                    copy(node_index->getNodeLocation(),
                                         node_index->getNodeLocation()+NDIM,
                                         &X[gauge_idx][NDIM*point_idx]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Determine the position of each pressure gauge.
    for (int k = 0; k < d_num_pressure_gauges; ++k)
    {
        SAMRAI::tbox::MPI::sumReduction(&X[k][0], d_num_pressure_gauge_points[k]);
        compute_centroid(d_pressure_gauge_loc[k], X[k]);
    }

    // Find the patch and index corresponding to the finest cell in the patch
    // hierarchy that contains the present position of each pressure gauge.
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++k)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
            hierarchy->getPatchLevel(ln);

        std::vector<SAMRAI::hier::Index<NDIM> > pressure_gauge_idx(d_num_pressure_gauges);
        for (int k = 0; k < d_num_pressure_gauges; ++k)
        {
            const std::vector<double>& X = d_pressure_gauge_loc[k];
            pressure_gauge_idx[k] = compute_pressure_gauge_index(X);
        }

        const SAMRAI::hier::BoxArray<NDIM>& box_array = level->getBoxes();
        std::vector<bool> found_pressure_gauge_box(d_num_pressure_gauges, false);
        for (int p = 0; p < box_array.size(); ++p)
        {
            const SAMRAI::hier::Box<NDIM>& box = box_array[p];
            for (int k = 0; k < d_num_pressure_gauges; ++k)
            {
                if (!found_pressure_gauge_box[k] && box.contains(pressure_gauge_idx[k]))
                {
                    d_pressure_gauge_level_num[k] = ln;
                    d_pressure_gauge_patch_num[k] = p;
                    d_pressure_gauge_on_proc[k] = (rank == level->getMappingForPatch(p));
                    found_pressure_gauge_box[k] = true;
                }
            }
        }
    }
    return;
}// initializePressureGauges

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBPlumbingToolkit>;

//////////////////////////////////////////////////////////////////////////////
