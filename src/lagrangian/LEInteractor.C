//
// LEInteractor.C
//
// Created on 14 Jul 2004
//         by Boyce Griffith (boyce@trasnaform.speakeasy.net).
//
// Last modified: <27.May.2006 17:27:55 boyce@boyce.cims.nyu.edu>
//

#include "LEInteractor.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// STL INCLUDES
//
#include <algorithm>
#include <functional>
#include <vector>

// SAMRAI-tools INCLUDES
//
#include "LNodeIndex.h"
#include "LNodeIndexSet.h"
#include "SAMRAI_tools_Utilities.h"

// SAMRAI INCLUDES
//
#include "CellIndex.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

// FORTRAN ROUTINES
//
#include <math.h>
#include <signal.h>

// Link between C/C++ and Fortran files
//       name in                         name in
//      C/C++ code                     Fortran code
//      ----------                     ------------
#define FORT_PW_CUBIC_INTERP           lagpwcubicinterp_
#define FORT_PW_CUBIC_SPREAD           lagpwcubicspread_

#define FORT_IB_4_INTERP               lagib4interp_
#define FORT_IB_4_SPREAD               lagib4spread_

#define FORT_IB_6_INTERP               lagib6interp_
#define FORT_IB_6_SPREAD               lagib6spread_

#define FORT_WIDE_IB_4_INTERP          lagwib4interp_
#define FORT_WIDE_IB_4_SPREAD          lagwib4spread_

// Function interfaces
extern "C" {
    void FORT_PW_CUBIC_INTERP(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
        );

    void FORT_PW_CUBIC_SPREAD(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
        );



    void FORT_IB_4_INTERP(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
        );

    void FORT_IB_4_SPREAD(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
        );



    void FORT_IB_6_INTERP(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
        );

    void FORT_IB_6_SPREAD(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
        );



    void FORT_WIDE_IB_4_INTERP(
        const double* , const double* , const double* , const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        const double* ,
        const int* , const double* , const int& ,
        const double* , double*
        );

    void FORT_WIDE_IB_4_SPREAD(
        const double* , const double* , const double* , const int& ,
        const int* , const double* , const int& ,
        const double* , const double* ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
#endif
        double*
        );
}

/////////////////////////////// INLINE ///////////////////////////////////////

//#ifdef DEBUG_NO_INLINE
//#include "LEInteractor.I"
//#endif

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
    // Timers.
    static tbox::Pointer<tbox::Timer> t_interpolate;
    static tbox::Pointer<tbox::Timer> t_interpolate_f77;
    static tbox::Pointer<tbox::Timer> t_spread;
    static tbox::Pointer<tbox::Timer> t_spread_f77;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void LEInteractor::initializeTimers()
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_interpolate = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::LEInteractor::interpolate()");
        t_interpolate_f77 = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::LEInteractor::interpolate()[fortran]");
        t_spread = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::LEInteractor::spread()");
        t_spread_f77 = tbox::TimerManager::getManager()->
            getTimer("SAMRAI-tools::LEInteractor::spread()[fortran]");
        timers_need_init = false;
    }    
    return;
}// initializeTimers

namespace
{
    struct GetLocalPETScIndex
        : unary_function<tbox::Pointer<LNodeIndex>,int>
    {
        int operator()(
            const tbox::Pointer<LNodeIndex>& index) const
            {
                return(index->getLocalPETScIndex());
            }
    };        
}

int LEInteractor::getStencilSize(
    const string& weighting_fcn)
{
    if (weighting_fcn == "PIECEWISE_CUBIC") return(4);
    if (weighting_fcn == "IB_4") return(4);
    if (weighting_fcn == "IB_6") return(6);
    if (weighting_fcn == "WIDE_IB_4") return(8);

    TBOX_ERROR("LEInteractor::getStencilSize()\n"
               << "  Unknown weighting function "
               << weighting_fcn << endl);

    return(-1);
}// getStencilSize

void LEInteractor::interpolate(
    tbox::Pointer<LNodeLevelData>& Q_data,
    const tbox::Pointer<LNodeLevelData>& X_data,
    const tbox::Pointer<LNodeIndexData>& idx_data,
    const tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
    const tbox::Pointer<hier::Patch<NDIM> >& patch,
    const hier::Box<NDIM>& box,
    const string& interp_fcn,
    const bool enforce_periodic_bcs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_data.isNull());
    assert(!q_data.isNull());
    assert(!X_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    assert(Q_data->getDepth() == q_data->getDepth());
    assert(X_data->getDepth() == NDIM);
#endif

    interpolate(&(*Q_data)(0), Q_data->getDepth(),
                &(*X_data)(0), X_data->getDepth(),
                idx_data,
                q_data,
                patch, box, interp_fcn, enforce_periodic_bcs);

    return;
}// interpolate

void LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const tbox::Pointer<LNodeIndexData>& idx_data,
    const tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
    const tbox::Pointer<hier::Patch<NDIM> >& patch,
    const hier::Box<NDIM>& box,
    const string& interp_fcn,
    const bool enforce_periodic_bcs)
{
    t_interpolate->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    assert(Q_depth == q_data->getDepth());
    assert(X_depth == NDIM);
#endif

    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& ilower = patch_box.lower();
    const hier::Index<NDIM>& iupper = patch_box.upper();
    
    const hier::Box<NDIM>& idx_ghost_box = idx_data->getGhostBox();
    const hier::IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    
    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    const int depth = q_data->getDepth();
    
    // Generate a list of local indices which lie in the specified
    // box.
    vector<int> local_indices;

    if (box == patch_box)
    {
        local_indices = idx_data->getInteriorLocalIndices();
    }
    else if (box == idx_ghost_box)
    {
        local_indices = idx_data->getInteriorLocalIndices();
        local_indices.insert(local_indices.end(),
                             idx_data->getGhostLocalIndices().begin(),
                             idx_data->getGhostLocalIndices().end());
    }
    else
    {
        local_indices.reserve(idx_data->getInteriorLocalIndices().size() +
                              idx_data->getGhostLocalIndices().size());
        
        for (LNodeIndexData::Iterator it(*idx_data); it; it++)
        {
            if (box.contains(it.getIndex()))
            {
                const LNodeIndexSet& node_set = *it;
                const LNodeIndexSet::size_type& num_ids = node_set.size();
                
                local_indices.resize(local_indices.size()+num_ids);
                transform(node_set.begin(), node_set.end(),
                          local_indices.end()-num_ids,
                          GetLocalPETScIndex());
            }
        }
    }

    // Generate the periodic offsets when necessary.
    vector<double> periodic_offsets(NDIM*local_indices.size(),0.0);

    if (enforce_periodic_bcs)
    {
        if (box == patch_box)
        {
            // intentionally blank
        }
        else if (box == idx_ghost_box)
        {
            int k = idx_data->getInteriorLocalIndices().size();
            for (LNodeIndexData::Iterator it(*idx_data); it; it++)
            {
                if (!patch_box.contains(it.getIndex()))
                {
                    const LNodeIndexSet& node_set = *it;
                    const hier::IntVector<NDIM>& offset =
                        node_set.getPeriodicOffset();

                    for (LNodeIndexSet::size_type n = 0;
                         n < node_set.size(); ++n, ++k)
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            periodic_offsets[NDIM*k+d] = static_cast<double>(offset(d))*dx[d];
                        }
                    }
                }
            }
        }
        else
        {
            int k = 0;
            for (LNodeIndexData::Iterator it(*idx_data); it; it++)
            {
                if (box.contains(it.getIndex()))
                {
                    const LNodeIndexSet& node_set = *it;
                    const hier::IntVector<NDIM>& offset =
                        node_set.getPeriodicOffset();

                    for (LNodeIndexSet::size_type n = 0;
                         n < node_set.size(); ++n, ++k)
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            periodic_offsets[NDIM*k+d] = static_cast<double>(offset(d))*dx[d];
                        }
                    }
                }
            }
        }
    }
    
    // Interpolate.
    t_interpolate_f77->start();

    const int local_indices_size = local_indices.size();
    if (local_indices_size > 0)
    {
        if (interp_fcn == "PIECEWISE_CUBIC")
        {
            FORT_PW_CUBIC_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else if (interp_fcn == "IB_4")
        {
            FORT_IB_4_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else if (interp_fcn == "IB_6")
        {
            FORT_IB_6_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else if (interp_fcn == "WIDE_IB_4")
        {
            FORT_WIDE_IB_4_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else
        {
            TBOX_ERROR("LEInteractor::interpolate()\n" <<
                       "  Unknown interpolation weighting function "
                       << interp_fcn << endl);
        }
    }
    t_interpolate_f77->stop();
    
    t_interpolate->stop();
    return;
}// interpolate

void LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const int num_vals,
    const tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
    const tbox::Pointer<hier::Patch<NDIM> >& patch,
    const hier::Box<NDIM>& box,
    const string& interp_fcn)
{
    t_interpolate->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
    assert(!patch.isNull());
    assert(Q_depth == q_data->getDepth());
    assert(X_depth == NDIM);
    assert(num_vals > 0);
#endif

    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& ilower = patch_box.lower();
    const hier::Index<NDIM>& iupper = patch_box.upper();
    
    const hier::IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    
    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    const int depth = q_data->getDepth();
    
    // Generate a list of local indices which lie in the specified
    // box.
    vector<int> local_indices;
    for (int l = 0; l < num_vals; ++l)
    {
        const pdat::CellIndex<NDIM> idx =
            SAMRAI_tools_Utilities::getCellIndex(
                &X_data[l*NDIM],xLower,dx,ilower);
        if (box.contains(idx)) local_indices.push_back(l);
    }
    
    // This routine does not have enough data to generate periodic
    // offsets.
    vector<double> periodic_offsets(NDIM*local_indices.size(),0.0);
    
    // Interpolate.
    t_interpolate_f77->start();

    const int local_indices_size = local_indices.size();
    if (local_indices_size > 0)
    {
        if (interp_fcn == "PIECEWISE_CUBIC")
        {
            FORT_PW_CUBIC_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else if (interp_fcn == "IB_4")
        {
            FORT_IB_4_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else if (interp_fcn == "IB_6")
        {
            FORT_IB_6_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else if (interp_fcn == "WIDE_IB_4")
        {
            FORT_WIDE_IB_4_INTERP(
                dx,xLower,xUpper,depth,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer(),
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data);
        }
        else
        {
            TBOX_ERROR("LEInteractor::interpolate()\n" <<
                       "  Unknown interpolation weighting function "
                       << interp_fcn << endl);
        }
    }
    t_interpolate_f77->stop();
    
    t_interpolate->stop();
    return;
}// interpolate

void LEInteractor::spread(
    tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
    const tbox::Pointer<LNodeLevelData>& Q_data,
    const tbox::Pointer<LNodeLevelData>& X_data,
    const tbox::Pointer<LNodeIndexData>& idx_data,
    const tbox::Pointer<hier::Patch<NDIM> >& patch,
    const hier::Box<NDIM>& box,
    const string& spread_fcn,
    const bool enforce_periodic_bcs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_data.isNull());
    assert(!q_data.isNull());
    assert(!X_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    assert(Q_data->getDepth() == q_data->getDepth());
    assert(X_data->getDepth() == NDIM);
#endif

    spread(q_data,
           &(*Q_data)(0), Q_data->getDepth(),
           &(*X_data)(0), X_data->getDepth(),
           idx_data,
           patch, box, spread_fcn, enforce_periodic_bcs);
    
    return;
}// spread

void LEInteractor::spread(
    tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const tbox::Pointer<LNodeIndexData>& idx_data,
    const tbox::Pointer<hier::Patch<NDIM> >& patch,
    const hier::Box<NDIM>& box,
    const string& spread_fcn,
    const bool enforce_periodic_bcs)
{
    t_spread->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    assert(Q_depth == q_data->getDepth());
    assert(X_depth == NDIM);
#endif

    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& ilower = patch_box.lower();
    const hier::Index<NDIM>& iupper = patch_box.upper();
    
    const hier::Box<NDIM>& idx_ghost_box = idx_data->getGhostBox();
    const hier::IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    
    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    const int depth = q_data->getDepth();
    
    // Generate a list of local indices which lie in the specified
    // box.
    vector<int> local_indices;

    if (box == patch_box)
    {
        local_indices = idx_data->getInteriorLocalIndices();
    }
    else if (box == idx_ghost_box)
    {
        local_indices = idx_data->getInteriorLocalIndices();
        local_indices.insert(local_indices.end(),
                             idx_data->getGhostLocalIndices().begin(),
                             idx_data->getGhostLocalIndices().end());
    }
    else
    {
        local_indices.reserve(idx_data->getInteriorLocalIndices().size() +
                              idx_data->getGhostLocalIndices().size());
        
        for (LNodeIndexData::Iterator it(*idx_data); it; it++)
        {
            if (box.contains(it.getIndex()))
            {
                const LNodeIndexSet& node_set = *it;
                const LNodeIndexSet::size_type& num_ids = node_set.size();
                
                local_indices.resize(local_indices.size()+num_ids);
                transform(node_set.begin(), node_set.end(),
                          local_indices.end()-num_ids,
                          GetLocalPETScIndex());
            }
        }
    }

    // Generate the periodic offsets when necessary.
    vector<double> periodic_offsets(NDIM*local_indices.size(),0.0);

    if (enforce_periodic_bcs)
    {
        if (box == patch_box)
        {
            // intentionally blank
        }
        else if (box == idx_ghost_box)
        {
            int k = idx_data->getInteriorLocalIndices().size();
            for (LNodeIndexData::Iterator it(*idx_data); it; it++)
            {
                if (!patch_box.contains(it.getIndex()))
                {
                    const LNodeIndexSet& node_set = *it;
                    const hier::IntVector<NDIM>& offset =
                        node_set.getPeriodicOffset();

                    for (LNodeIndexSet::size_type n = 0;
                         n < node_set.size(); ++n, ++k)
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            periodic_offsets[NDIM*k+d] = static_cast<double>(offset(d))*dx[d];
                        }
                    }
                }
            }
        }
        else
        {
            int k = 0;
            for (LNodeIndexData::Iterator it(*idx_data); it; it++)
            {
                if (box.contains(it.getIndex()))
                {
                    const LNodeIndexSet& node_set = *it;
                    const hier::IntVector<NDIM>& offset =
                        node_set.getPeriodicOffset();

                    for (LNodeIndexSet::size_type n = 0;
                         n < node_set.size(); ++n, ++k)
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            periodic_offsets[NDIM*k+d] = static_cast<double>(offset(d))*dx[d];
                        }
                    }
                }
            }
        }
    }
    
    // Spread.
    t_spread_f77->start();

    const int local_indices_size = local_indices.size();
    if (local_indices_size > 0)
    {
        if (spread_fcn == "PIECEWISE_CUBIC")
        {
            FORT_PW_CUBIC_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else if (spread_fcn == "IB_4")
        {
            FORT_IB_4_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else if (spread_fcn == "IB_6")
        {
            FORT_IB_6_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data, Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else if (spread_fcn == "WIDE_IB_4")
        {
            FORT_WIDE_IB_4_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread()\n" <<
                       "  Unknown spreading weighting function "
                       << spread_fcn << endl);
        }
    }
    t_spread_f77->stop();
    
    t_spread->stop();
    return;
}// spread

void LEInteractor::spread(
    tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const int num_vals,
    const tbox::Pointer<hier::Patch<NDIM> >& patch,
    const hier::Box<NDIM>& box,
    const string& spread_fcn)
{
    t_spread->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
    assert(!patch.isNull());
    assert(Q_depth == q_data->getDepth());
    assert(X_depth == NDIM);
    assert(num_vals > 0);
#endif

    const hier::Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& ilower = patch_box.lower();
    const hier::Index<NDIM>& iupper = patch_box.upper();
    
    const hier::IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    
    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    const int depth = q_data->getDepth();
    
    // Generate a list of local indices which lie in the specified
    // box.
    vector<int> local_indices;
    for (int l = 0; l < num_vals; ++l)
    {
        const pdat::CellIndex<NDIM> idx =
            SAMRAI_tools_Utilities::getCellIndex(
                &X_data[l*NDIM],xLower,dx,ilower);
        if (box.contains(idx)) local_indices.push_back(l);
    }

    // This routine does not have enough data to generate periodic
    // offsets.
    vector<double> periodic_offsets(NDIM*local_indices.size(),0.0);
    
    // Spread.
    t_spread_f77->start();

    const int local_indices_size = local_indices.size();
    if (local_indices_size > 0)
    {
        if (spread_fcn == "PIECEWISE_CUBIC")
        {
            FORT_PW_CUBIC_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else if (spread_fcn == "IB_4")
        {
            FORT_IB_4_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else if (spread_fcn == "IB_6")
        {
            FORT_IB_6_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else if (spread_fcn == "WIDE_IB_4")
        {
            FORT_WIDE_IB_4_SPREAD(
                dx,xLower,xUpper,depth,
                &local_indices[0], &periodic_offsets[0], local_indices_size,
                X_data,Q_data,
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                q_gcw(0),q_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                q_gcw(0),q_gcw(1),q_gcw(2),
#endif
                q_data->getPointer());
        }
        else
        {
            TBOX_ERROR("LEInteractor::spread()\n" <<
                       "  Unknown spreading weighting function "
                       << spread_fcn << endl);
        }
    }
    t_spread_f77->stop();
    
    t_spread->stop();
    return;
}// spread

//////////////////////////////////////////////////////////////////////////////
