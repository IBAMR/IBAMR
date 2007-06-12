// Filename: IBInstrumentPanel.C
// Last modified: <12.Jun.2007 18:14:01 griffith@box221.cims.nyu.edu>
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

// STOOLS INCLUDES
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec-et.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <cassert>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
inline void
init_meter_elements2d(
    blitz::Array<blitz::TinyVector<double,NDIM>,2>& X_web,
    blitz::Array<blitz::TinyVector<double,NDIM>,2>& dA_web,
    const blitz::Array<blitz::TinyVector<double,NDIM>,1>& X_perimeter,
    const blitz::TinyVector<double,NDIM>& X_centroid)
{
#if (NDIM != 2)
    assert(false);
#endif
    assert(false);  // XXXX
    return;
}// init_meter_elements2d

inline void
init_meter_elements3d(
    blitz::Array<blitz::TinyVector<double,NDIM>,2>& X_web,
    blitz::Array<blitz::TinyVector<double,NDIM>,2>& dA_web,
    const blitz::Array<blitz::TinyVector<double,NDIM>,1>& X_perimeter,
    const blitz::TinyVector<double,NDIM>& X_centroid)
{
#if (NDIM != 3)
    assert(false);
#else
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(X_web.extent(0) == X_perimeter.extent(0));
    assert(dA_web.extent(0) == X_perimeter.extent(0));
#endif
    const int num_perimeter_nodes = X_perimeter.extent(1);
    const int num_web_nodes = X_web.extent(2);
    for (int m = 0; m < num_perimeter_nodes; ++m)
    {
        const blitz::TinyVector<double,NDIM> X_perimeter0(X_perimeter(m));
        const blitz::TinyVector<double,NDIM> dX0((X_centroid-X_perimeter0)/double(num_web_nodes));

        const blitz::TinyVector<double,NDIM> X_perimeter1(X_perimeter((m+1)%num_perimeter_nodes));
        const blitz::TinyVector<double,NDIM> dX1((X_centroid-X_perimeter1)/double(num_web_nodes));

        for (int n = 0; n < num_web_nodes; ++n)
        {
            // Compute the four vertices of the quadrilateral web patch.
            //
            // Note that here the vertices are placed in "standard" (i.e.,
            // "counter-clockwise") orientation.
            const blitz::TinyVector<double,NDIM> X0(X_perimeter0+double(n  )*dX0);
            const blitz::TinyVector<double,NDIM> X1(X_perimeter1+double(n  )*dX1);
            const blitz::TinyVector<double,NDIM> X2(X_perimeter1+double(n+1)*dX1);
            const blitz::TinyVector<double,NDIM> X3(X_perimeter0+double(n+1)*dX0);
#ifdef DEBUG_CHECK_ASSERTIONS
            if (n+1 == num_web_nodes)
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    assert(SAMRAI::tbox::Utilities::deq(X2(d),X_centroid(d)));
                    assert(SAMRAI::tbox::Utilities::deq(X3(d),X_centroid(d)));
                }
            }
#endif
            // Compute the centroid of the quadrilateral web patch.
            X_web(m,n) = 0.25*(X0+X1+X2+X3);

            // Compute the area-weighted normal to the quadrilateral web patch,
            // i.e.,
            //
            //    dA = 0.5*((X2-X1) X (X0-X1) + (X0-X3) X (X2-X3))
            //
            // Note that by construction, the quadrilateral is guaranteed to lie
            // within a plane, and that the vectors are chosen so that the
            // resulting normal vectors have the same orientation.
            dA_web(m,n) = 0.5*(
                cross(blitz::TinyVector<double,NDIM>(X2-X1), blitz::TinyVector<double,NDIM>(X0-X1)) +
                cross(blitz::TinyVector<double,NDIM>(X0-X3), blitz::TinyVector<double,NDIM>(X2-X3)));
        }
    }
#endif
    return;
}// init_meter_elements3d
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentPanel::IBInstrumentPanel(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_num_meters(0),
      d_num_perimeter_nodes(),
      d_X_centroid(),
      d_X_perimeter(),
      d_X_web(),
      d_dA_web(),
      d_web_patch_map(),
      d_meter_centroid_map()
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
IBInstrumentPanel::initializeHierarchyIndependentData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const lag_manager)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Determine how many flow meters/pressure gauges are present in the local
    // data.
    int max_meter_index = -1;
    std::vector<int> max_node_index;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
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
                                const int m = spec->getMeterIndex();
                                max_meter_index = max(m, max_meter_index);

                                const int n = spec->getNodeIndex();
                                max_node_index.resize(max_meter_index+1,-1);
                                max_node_index[m] = max(n, max_node_index[m]);
                            }
                        }
                    }
                }
            }
        }
    }

    // Communicate local data to all processes.
    d_num_meters = SAMRAI::tbox::MPI::maxReduction(max_meter_index)+1;
    max_node_index.resize(d_num_meters,-1);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_meters >= 0);
#endif
    d_num_perimeter_nodes.clear();
    d_num_perimeter_nodes.resize(d_num_meters,-1);
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_num_perimeter_nodes[m] = max_node_index[m]+1;
    }
    SAMRAI::tbox::MPI::maxReduction(&d_num_perimeter_nodes[0], d_num_meters);
#ifdef DEBUG_CHECK_ASSERTIONS
    for (int m = 0; m < d_num_meters; ++m)
    {
        assert(d_num_perimeter_nodes[m] > 0);
    }
#endif

    // Resize arrays.
    d_X_centroid.resize(d_num_meters);
    d_X_perimeter.resize(d_num_meters);
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_X_perimeter[m].resize(d_num_perimeter_nodes[m]);
    }
    d_X_web.resize(d_num_meters);
    d_dA_web.resize(d_num_meters);
    return;
}// initializeHierarchyIndependentData

void
IBInstrumentPanel::initializeHierarchyDependentData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const lag_manager)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Loop over all local nodes to determine the positions of the local
    // perimeter nodes.
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            d_X_perimeter[m](n) = 0.0;
        }
    }
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
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
                                const int m = spec->getMeterIndex();
                                const int n = spec->getNodeIndex();
                                copy(node_idx->getNodeLocation(),
                                     node_idx->getNodeLocation()+NDIM,
                                     d_X_perimeter[m](n).data());
                            }
                        }
                    }
                }
            }
        }
    }

    // Set the positions of all perimeter nodes on all processes.
    std::vector<double> X_perimeter_flattened;
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            X_perimeter_flattened.insert(
                X_perimeter_flattened.end(),
                d_X_perimeter[m](n).data(),d_X_perimeter[m](n).data()+NDIM);
        }
    }

    SAMRAI::tbox::MPI::sumReduction(&X_perimeter_flattened[0],X_perimeter_flattened.size());

    for (int m = 0, k = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n, ++k)
        {
            copy(&X_perimeter_flattened[NDIM*k],(&X_perimeter_flattened[NDIM*k])+NDIM,
                 d_X_perimeter[m](n).data());
        }
    }

    // Determine the centroid of each perimeter.
    for (int m = 0; m < d_num_meters; ++m)
    {
        d_X_centroid[m] = 0.0;
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            d_X_centroid[m] += d_X_perimeter[m](n);
        }
        d_X_centroid[m] /= double(d_num_perimeter_nodes[m]);
    }

    // Determine the maximum distance from perimeter nodes to centroids.
    std::vector<double> r_max(d_num_meters,0.0);
    for (int m = 0; m < d_num_meters; ++m)
    {
        for (int n = 0; n < d_num_perimeter_nodes[m]; ++n)
        {
            const blitz::TinyVector<double,NDIM> r(d_X_perimeter[m](n)-d_X_centroid[m]);
            r_max[m] = std::max(r_max[m],sqrt(dot(r,r)));
        }
    }

    // Determine the finest grid spacing.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom =
        hierarchy->getGridGeometry();
    const double* const domainXLower = grid_geom->getXLower();
    const double* const domainXUpper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    assert(grid_geom->getDomainIsSingleBox());
    const SAMRAI::hier::Box<NDIM> domain_box = grid_geom->getPhysicalDomain()(0);

    const SAMRAI::hier::IntVector<NDIM>& ratio_to_level_zero =
        hierarchy->getPatchLevel(finest_ln)->getRatio();
    std::vector<double> dx_finest(NDIM,0.0);
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d]/double(ratio_to_level_zero(d));
    }
    const double h_finest = *(std::min_element(dx_finest.begin(),dx_finest.end()));

    // Build the meter web patch centroids and area elements.
    for (int m = 0; m < d_num_meters; ++m)
    {
        int num_web_nodes = 2*static_cast<int>(ceil(r_max[m]/h_finest));
        d_X_web[m].resize(d_num_perimeter_nodes[m],num_web_nodes);
        d_dA_web[m].resize(d_num_perimeter_nodes[m],num_web_nodes);
#if (NDIM == 2)
        init_meter_elements2d(d_X_web[m],d_dA_web[m],d_X_perimeter[m],d_X_centroid[m]);
#endif
#if (NDIM == 3)
        init_meter_elements3d(d_X_web[m],d_dA_web[m],d_X_perimeter[m],d_X_centroid[m]);
#endif
    }

    // Setup the mappings from cell indices to the web patch data.
    //
    // NOTE: Each web patch is assigned to precisely one cell in precisely one
    // level.  In particular, each web patch is assigned to whichever cell is
    // the finest cell that contains the region of physical space in which the
    // web patch centroid is located.
    d_web_patch_map.clear();
    d_web_patch_map.resize(finest_ln+1);
    d_meter_centroid_map.clear();
    d_meter_centroid_map.resize(finest_ln+1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        const SAMRAI::hier::IntVector<NDIM>& ratio = level->getRatio();
        const SAMRAI::hier::Box<NDIM> domain_box_level = SAMRAI::hier::Box<NDIM>::refine(
            domain_box, ratio);
        const SAMRAI::hier::Index<NDIM>& domain_box_level_lower = domain_box_level.lower();
        const SAMRAI::hier::Index<NDIM>& domain_box_level_upper = domain_box_level.upper();
        double dx[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d]/double(ratio(d));
        }

        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > finer_level =
            (ln < finest_ln
             ? hierarchy->getPatchLevel(ln+1)
             : SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL));
        const SAMRAI::hier::IntVector<NDIM>& finer_ratio =
            (ln < finest_ln
             ? finer_level->getRatio()
             : SAMRAI::hier::IntVector<NDIM>(1));
        const SAMRAI::hier::Box<NDIM> finer_domain_box_level = SAMRAI::hier::Box<NDIM>::refine(
            domain_box, finer_ratio);
        const SAMRAI::hier::Index<NDIM>& finer_domain_box_level_lower = finer_domain_box_level.lower();
        const SAMRAI::hier::Index<NDIM>& finer_domain_box_level_upper = finer_domain_box_level.upper();
        double finer_dx[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            finer_dx[d] = dx_coarsest[d]/double(finer_ratio(d));
        }

        for (int l = 0; l < d_num_meters; ++l)
        {
            // Setup the web patch mapping.
            for (int m = 0; m < d_X_web[l].extent(0); ++m)
            {
                for (int n = 0; n < d_X_web[l].extent(1); ++n)
                {
                    const double* const X = d_X_web[l](m,n).data();
                    const SAMRAI::hier::Index<NDIM> i =
                        STOOLS::STOOLS_Utilities::getCellIndex(
                            X, domainXLower, domainXUpper, dx,
                            domain_box_level_lower, domain_box_level_upper);
                    const SAMRAI::hier::Index<NDIM> finer_i =
                        (ln < finest_ln
                         ? STOOLS::STOOLS_Utilities::getCellIndex(
                             X, domainXLower, domainXUpper, finer_dx,
                             finer_domain_box_level_lower, finer_domain_box_level_upper)
                         : SAMRAI::hier::Index<NDIM>(-1));

                    if (level->getBoxes().contains(i) &&
                        (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
                    {
                        WebPatch p;
                        p.meter_num = l;
                        p.X = &d_X_web[l](m,n);
                        p.dA = &d_dA_web[l](m,n);
                        d_web_patch_map[ln].insert(std::make_pair(i,p));
                    }
                }
            }

            // Setup the web centroid mapping.
            {
                const double* const X = d_X_centroid[l].data();
                const SAMRAI::hier::Index<NDIM> i =
                    STOOLS::STOOLS_Utilities::getCellIndex(
                        X, domainXLower, domainXUpper, dx,
                        domain_box_level_lower, domain_box_level_upper);
                const SAMRAI::hier::Index<NDIM> finer_i =
                    (ln < finest_ln
                     ? STOOLS::STOOLS_Utilities::getCellIndex(
                         X, domainXLower, domainXUpper, finer_dx,
                         finer_domain_box_level_lower, finer_domain_box_level_upper)
                     : SAMRAI::hier::Index<NDIM>(-1));

                if (level->getBoxes().contains(i) &&
                    (ln == finest_ln || !finer_level->getBoxes().contains(finer_i)))
                {
                    MeterCentroid c;
                    c.meter_num = l;
                    c.X = &d_X_centroid[l];
                    d_meter_centroid_map[ln].insert(std::make_pair(i,c));
                }
            }
        }
    }
    return;
}// initializeHierarchyDependentData

void
IBInstrumentPanel::readMeterData(
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_var,
    const int U_data_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_var,
    const int P_data_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const lag_manager)
{
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();

    // Compute the local contributions to the flux of U through the flow meter,
    // and the value of P at the centroid of the meter.
    std::vector<double> U_dA(d_num_meters,0.0);
    std::vector<double> P_centroid(d_num_meters,0.0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
            hierarchy->getPatchLevel(ln);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();

            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_data =
                patch->getPatchData(U_data_idx);
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_data =
                patch->getPatchData(P_data_idx);

            for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& i = b();

                std::pair<WebPatchMap::const_iterator,WebPatchMap::const_iterator> r1 =
                    d_web_patch_map[ln].equal_range(i);
                for (WebPatchMap::const_iterator it = r1.first; it != r1.second; ++it)
                {
                    const int& meter_num = (*it).second.meter_num;
                    const blitz::TinyVector<double,NDIM>& X = *((*it).second.X);
                    const blitz::TinyVector<double,NDIM>& dA = *((*it).second.dA);

                    // XXXX: Should perform linear interpolation here!!!
                    (void) X;
                    const blitz::TinyVector<double,NDIM> U(
                        (*U_data)(i,0),(*U_data)(i,1)
#if (NDIM == 3)
                        ,(*U_data)(i,2)
#endif
                                                           );
                    U_dA[meter_num] += dot(U,dA);
                }

                std::pair<MeterCentroidMap::const_iterator,MeterCentroidMap::const_iterator> r2 =
                    d_meter_centroid_map[ln].equal_range(i);
                for (MeterCentroidMap::const_iterator it = r2.first; it != r2.second; ++it)
                {
                    const int& meter_num = (*it).second.meter_num;
                    const blitz::TinyVector<double,NDIM>& X = *((*it).second.X);

                    // XXXX: Should perform linear interpolation here!!!
                    (void) X;
                    P_centroid[meter_num] = (*P_data)(i);
                }
            }
        }
    }

    // Synchronize the values across all processes.
    SAMRAI::tbox::MPI::sumReduction(&U_dA[0],d_num_meters);
    SAMRAI::tbox::MPI::sumReduction(&P_centroid[0],d_num_meters);
    return;
}// readMeterData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBInstrumentPanel>;

//////////////////////////////////////////////////////////////////////////////
