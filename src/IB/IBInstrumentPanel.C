// Filename: IBInstrumentPanel.C
// Last modified: <12.Jun.2007 02:01:20 boyce@bigboy.nyconnect.com>
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

namespace
{
inline void
init_meter_elements2d(
    std::vector<std::vector<double> >& X_web,
    std::vector<std::vector<double> >& dA_web,
    const int interp_level,
    const int num_perimeter_nodes,
    const std::vector<double>& X_perimeter,
    const std::vector<double>& X_centroid)
{
#if (NDIM != 2)
    assert(false);
#endif
    assert(false);  // XXXX
    return;
}// init_meter_elements2d

inline void
cross_product(
    std::vector<double>& z,
    const std::vector<double>& x,
    const std::vector<double>& y)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(x.size() == NDIM);
    assert(y.size() == NDIM);
    assert(z.size() == NDIM);
#endif
    z[0] = x[1]*y[2]-x[2]*y[1];
    z[1] = x[2]*y[0]-x[0]*y[2];
    z[2] = x[0]*y[1]-x[1]*y[0];
    return;
}// cross_product

inline void
init_meter_elements3d(
    std::vector<std::vector<double> >& X_web,
    std::vector<std::vector<double> >& dA_web,
    const int interp_level,
    const int num_perimeter_nodes,
    const std::vector<double>& X_perimeter,
    const std::vector<double>& X_centroid)
{
#if (NDIM != 3)
    assert(false);
#else
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(static_cast<int>( X_web.size()) == num_perimeter_nodes);
    assert(static_cast<int>(dA_web.size()) == num_perimeter_nodes);
    assert(static_cast<int>(X_perimeter.size()) == NDIM*num_perimeter_nodes);
    assert(static_cast<int>( X_centroid.size()) == NDIM);
#endif
    for (int n = 0; n < num_perimeter_nodes; ++n)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(static_cast<int>( X_web[n].size()) == NDIM*interp_level);
        assert(static_cast<int>(dA_web[n].size()) == NDIM*interp_level);
#endif
        const std::vector<double> X_perimeter0(&X_perimeter[NDIM*n],(&X_perimeter[NDIM*n])+NDIM);
        const std::vector<double> X_perimeter1(&X_perimeter[NDIM*((n+1)%num_perimeter_nodes)],(&X_perimeter[NDIM*((n+1)%num_perimeter_nodes)])+NDIM);
        std::vector<double> dX0(NDIM);
        std::vector<double> dX1(NDIM);
        for (int d = 0; d < NDIM; ++d)
        {
            dX0[d] = (X_centroid[d]-X_perimeter0[d])/double(interp_level+1);
            dX1[d] = (X_centroid[d]-X_perimeter1[d])/double(interp_level+1);
        }
        std::vector<double> X0(NDIM), X1(NDIM), X2(NDIM), X3(NDIM);
        std::vector<double> a(NDIM), b(NDIM), c(NDIM);
        for (int l = 0; l < interp_level; ++l)
        {
            // Compute the four vertices of the quadrilateral.
            for (int d = 0; d < NDIM; ++d)
            {
                X0[d] = X_perimeter0[d]+double(l  )*dX0[d];
                X1[d] = X_perimeter1[d]+double(l  )*dX1[d];
                X2[d] = X_perimeter1[d]+double(l+1)*dX1[d];
                X3[d] = X_perimeter0[d]+double(l+1)*dX0[d];
            }

            // Compute the centroid of the quatiladeral.
            for (int d = 0; d < NDIM; ++d)
            {
                X_web[n][NDIM*l+d] = 0.25*(X0[d]+X1[d]+X2[d]+X3[d]);
            }

            // Compute the area weighted normal to the quadrilateral, i.e.,
            //
            //    dA = 0.5*((X2-X1) X (X0-X1) + (X0-X3) X (X2-X3))
            //
            // Note that the quadrilateral is guaranteed to lie within a plane,
            // and that the cross-products are setup so that the normal vectors
            // have the same orientation.
            for (int d = 0; d < NDIM; ++d)
            {
                a[d] = X2[d]-X1[d];
                b[d] = X0[d]-X1[d];
            }
            cross_product(c,a,b);
            for (int d = 0; d < NDIM; ++d)
            {
                dA_web[n][NDIM*l+d] = 0.5*c[d];
            }

            for (int d = 0; d < NDIM; ++d)
            {
                a[d] = X0[d]-X3[d];
                b[d] = X2[d]-X3[d];
            }
            cross_product(c,a,b);
            for (int d = 0; d < NDIM; ++d)
            {
                dA_web[n][NDIM*l+d] += 0.5*c[d];
            }
        }
    }
#endif
    return;
}// init_meter_elements3d
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentPanel::IBInstrumentPanel(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_interp_level(10),
      d_num_meters(0),
      d_num_perimeter_nodes(),
      d_X_perimeter(),
      d_X_centroid(),
      d_X_web(),
      d_dA_web()
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
IBInstrumentPanel::readMeters(
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > U_var,
    const int U_data_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_var,
    const int P_data_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const lag_manager)
{
    // Compute the flux of U through the flow meter, and the value of P at the
    // centroid of the meter.
    std::vector<double> U_dA(d_num_meters,0.0);
    std::vector<double> P_centroid(d_num_meters,0.0);
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
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
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,bool> > is_refined_data =
                patch->getPatchData(d_is_refined_data_idx);

            for (SAMRAI::hier::Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& i = b();
                if (!(*is_refined_data)(i))
                {
                    std::pair<WebPatchMap::const_iterator,WebPatchMap::const_iterator> r;

                    r = d_web_patch_map[ln].equal_range(i);
                    for (WebPatchMap::const_iterator it = r.first; it != r.second; ++it)
                    {
                        const int& meter_num = (*it).second.meter_num;
                        const double* const X = (*it).second.X;
                        const double* const dA = (*it).second.dA;

                        (void) X;

                        // XXXX: Should perform linear interpolation here!!!
                        for (int d = 0; d < NDIM; ++d)
                        {
                            U_dA[meter_num] += (*U_data)(i,d)*dA[d];
                        }
                    }

                    r = d_centroid_map[ln].equal_range(i);
                    for (WebPatchMap::const_iterator it = r.first; it != r.second; ++it)
                    {
                        const int& meter_num = (*it).second.meter_num;
                        const double* const X = (*it).second.X;
                        const double* const dA = (*it).second.dA;

                        (void) X;
                        (void) dA;

                        // XXXX: Should perform linear interpolation here!!!
                        P_centroid[meter_num] = (*P_data)(i);
                    }
                }
            }
        }
    }

    SAMRAI::tbox::MPI::sumReduction(&U_dA[0],d_num_meters);
    SAMRAI::tbox::MPI::sumReduction(&P_centroid[0],d_num_meters);

    return;
}// readMeters

void
IBInstrumentPanel::initializeMeterData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const lag_manager)
{
    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager-> getLNodeIndexPatchDescriptorIndex();

    // Loop over all local nodes to determine the positions of the local
    // perimeter nodes.
    d_X_perimeter.clear();
    d_X_perimeter.resize(d_num_meters);
    for (int k = 0; k < d_num_meters; ++k)
    {
        d_X_perimeter[k].resize(NDIM*d_num_perimeter_nodes[k], 0.0);
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

    // Set the positions of all perimeter nodes on all processes.
    for (int k = 0; k < d_num_meters; ++k)
    {
        SAMRAI::tbox::MPI::sumReduction(&d_X_perimeter[k][0], NDIM*d_num_perimeter_nodes[k]);
    }

    // Determine the centroid of each perimeter.
    d_X_centroid.clear();
    d_X_centroid.resize(d_num_meters);
    for (int k = 0; k < d_num_meters; ++k)
    {
        d_X_centroid[k].resize(NDIM, 0.0);
        for (int n = 0; n < d_num_perimeter_nodes[k]; ++n)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                d_X_centroid[k][d] += d_X_perimeter[k][NDIM*d_num_perimeter_nodes[k]+d];
            }
        }
        for (int d = 0; d < NDIM; ++d)
        {
            d_X_centroid[k][d] /= double(d_num_perimeter_nodes[k]);
        }
    }

    // Build the meter web centroids and area elements.
    d_X_web.clear();
    d_X_web.resize(d_num_meters);
    d_dA_web.clear();
    d_dA_web.resize(d_num_meters);
    for (int k = 0; k < d_num_meters; ++k)
    {
        d_X_web[k].resize(d_num_perimeter_nodes[k],std::vector<double>(NDIM*d_interp_level,0.0));
        d_dA_web[k].resize(d_num_perimeter_nodes[k],std::vector<double>(NDIM*d_interp_level,0.0));
#if (NDIM == 2)
        init_meter_elements2d(d_X_web[k],d_dA_web[k],d_interp_level,d_num_perimeter_nodes[k],d_X_perimeter[k],d_X_centroid[k]);
#endif
#if (NDIM == 3)
        init_meter_elements3d(d_X_web[k],d_dA_web[k],d_interp_level,d_num_perimeter_nodes[k],d_X_perimeter[k],d_X_centroid[k]);
#endif
    }
    return;
}// initializeMeterData

void
IBInstrumentPanel::initializeHierarchyData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
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

    d_num_perimeter_nodes.clear();
    d_num_perimeter_nodes.resize(d_num_meters,0);
    for (int k = 0; k < d_num_meters; ++k)
    {
        d_num_perimeter_nodes[k] = max_node_index[k]+1;
    }
    SAMRAI::tbox::MPI::maxReduction(&d_num_perimeter_nodes[0], d_num_meters);
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
