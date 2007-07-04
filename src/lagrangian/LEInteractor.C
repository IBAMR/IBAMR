// Filename: LEInteractor.C
// Last modified: <04.Jul.2007 13:38:00 boyce@bigboy.nyconnect.com>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "LEInteractor.h"

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
#include <ibamr/LNodeIndex.h>
#include <ibamr/LNodeIndexSet.h>

// STOOLS INCLUDES
#include <stools/STOOLS_Utilities.h>

// SAMRAI INCLUDES
#include <CellIndex.h>
#include <CartesianPatchGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define LAGRANGIAN_PWCUBIC_INTERP_F77 F77_FUNC_(lagrangian_pwcubic_interp2d, LAGRANGIAN_PWCUBIC_INTERP2D)
#define LAGRANGIAN_PWCUBIC_SPREAD_F77 F77_FUNC_(lagrangian_pwcubic_spread2d, LAGRANGIAN_PWCUBIC_SPREAD2D)

#define LAGRANGIAN_IB4_INTERP_F77 F77_FUNC_(lagrangian_ib4_interp2d, LAGRANGIAN_IB4_INTERP2D)
#define LAGRANGIAN_IB4_SPREAD_F77 F77_FUNC_(lagrangian_ib4_spread2d, LAGRANGIAN_IB4_SPREAD2D)

#define LAGRANGIAN_WIB4_INTERP_F77 F77_FUNC_(lagrangian_wib4_interp2d, LAGRANGIAN_WIB4_INTERP2D)
#define LAGRANGIAN_WIB4_SPREAD_F77 F77_FUNC_(lagrangian_wib4_spread2d, LAGRANGIAN_WIB4_SPREAD2D)

#define LAGRANGIAN_IB6_INTERP_F77 F77_FUNC_(lagrangian_ib6_interp2d, LAGRANGIAN_IB6_INTERP2D)
#define LAGRANGIAN_IB6_SPREAD_F77 F77_FUNC_(lagrangian_ib6_spread2d, LAGRANGIAN_IB6_SPREAD2D)

#endif

#if (NDIM == 3)
#define LAGRANGIAN_PWCUBIC_INTERP_F77 F77_FUNC_(lagrangian_pwcubic_interp3d, LAGRANGIAN_PWCUBIC_INTERP3D)
#define LAGRANGIAN_PWCUBIC_SPREAD_F77 F77_FUNC_(lagrangian_pwcubic_spread3d, LAGRANGIAN_PWCUBIC_SPREAD3D)

#define LAGRANGIAN_IB4_INTERP_F77 F77_FUNC_(lagrangian_ib4_interp3d, LAGRANGIAN_IB4_INTERP3D)
#define LAGRANGIAN_IB4_SPREAD_F77 F77_FUNC_(lagrangian_ib4_spread3d, LAGRANGIAN_IB4_SPREAD3D)

#define LAGRANGIAN_WIB4_INTERP_F77 F77_FUNC_(lagrangian_wib4_interp3d, LAGRANGIAN_WIB4_INTERP3D)
#define LAGRANGIAN_WIB4_SPREAD_F77 F77_FUNC_(lagrangian_wib4_spread3d, LAGRANGIAN_WIB4_SPREAD3D)

#define LAGRANGIAN_IB6_INTERP_F77 F77_FUNC_(lagrangian_ib6_interp3d, LAGRANGIAN_IB6_INTERP3D)
#define LAGRANGIAN_IB6_SPREAD_F77 F77_FUNC_(lagrangian_ib6_spread3d, LAGRANGIAN_IB6_SPREAD3D)

#endif

extern "C"
{
    void
    LAGRANGIAN_PWCUBIC_INTERP_F77(
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

    void
    LAGRANGIAN_PWCUBIC_SPREAD_F77(
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



    void
    LAGRANGIAN_IB4_INTERP_F77(
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

    void
    LAGRANGIAN_IB4_SPREAD_F77(
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



    void
    LAGRANGIAN_IB6_INTERP_F77(
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

    void
    LAGRANGIAN_IB6_SPREAD_F77(
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



    void
    LAGRANGIAN_WIB4_INTERP_F77(
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

    void
    LAGRANGIAN_WIB4_SPREAD_F77(
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_interpolate;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_interpolate_f77;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_spread;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_spread_f77;

struct GetLocalPETScIndex
    : std::unary_function<SAMRAI::tbox::Pointer<LNodeIndex>,int>
{
    inline int
    operator()(
        const SAMRAI::tbox::Pointer<LNodeIndex>& index) const
        {
            return index->getLocalPETScIndex();
        }
};
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
LEInteractor::initializeTimers()
{
    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_interpolate = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LEInteractor::interpolate()");
        t_interpolate_f77 = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LEInteractor::interpolate()[fortran]");
        t_spread = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LEInteractor::spread()");
        t_spread_f77 = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::LEInteractor::spread()[fortran]");
        timers_need_init = false;
    }
    return;
}// initializeTimers

int
LEInteractor::getStencilSize(
    const std::string& weighting_fcn)
{
    if (weighting_fcn == "PIECEWISE_CUBIC") return 4;
    if (weighting_fcn == "IB_4") return 4;
    if (weighting_fcn == "IB_6") return 6;
    if (weighting_fcn == "WIDE_IB_4") return 8;

    TBOX_ERROR("LEInteractor::getStencilSize()\n"
               << "  Unknown weighting function "
               << weighting_fcn << endl);

    return -1;
}// getStencilSize

void
LEInteractor::interpolate(
    SAMRAI::tbox::Pointer<LNodeLevelData>& Q_data,
    const SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
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
                patch, box, periodic_shift, interp_fcn);

    return;
}// interpolate

void
LEInteractor::interpolate(
    double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& interp_fcn)
{
    t_interpolate->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    assert(Q_depth == q_data->getDepth());
    assert(X_depth == NDIM);
#endif

    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
    const SAMRAI::hier::Index<NDIM>& ilower = patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch_box.upper();

    const SAMRAI::hier::IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    const int depth = q_data->getDepth();

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    getLocalIndices(local_indices, box, idx_data);

    // Generate periodic offsets.
    std::vector<double> periodic_offsets(NDIM*local_indices.size(),0.0);
    getPeriodicOffsets(periodic_offsets, box, patch, periodic_shift, idx_data);

    // Interpolate.
    t_interpolate_f77->start();

    const int local_indices_size = local_indices.size();
    if (local_indices_size > 0)
    {
        if (interp_fcn == "PIECEWISE_CUBIC")
        {
            LAGRANGIAN_PWCUBIC_INTERP_F77(
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
            LAGRANGIAN_IB4_INTERP_F77(
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
            LAGRANGIAN_IB6_INTERP_F77(
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
            LAGRANGIAN_WIB4_INTERP_F77(
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

void
LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const SAMRAI::tbox::Pointer<LNodeLevelData>& Q_data,
    const SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
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
           patch, box, periodic_shift, spread_fcn);
    return;
}// spread

void
LEInteractor::spread(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
    const double* const Q_data,
    const int Q_depth,
    const double* const X_data,
    const int X_depth,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
    t_spread->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    assert(Q_depth == q_data->getDepth());
    assert(X_depth == NDIM);
#endif
    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
    const SAMRAI::hier::Index<NDIM>& ilower = patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch_box.upper();

    const SAMRAI::hier::IntVector<NDIM>& q_gcw = q_data->getGhostCellWidth();
    const int depth = q_data->getDepth();

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    // Generate a list of local indices which lie in the specified box.
    std::vector<int> local_indices;
    getLocalIndices(local_indices, box, idx_data);

    // Generate periodic offsets.
    std::vector<double> periodic_offsets(NDIM*local_indices.size(),0.0);
    getPeriodicOffsets(periodic_offsets, box, patch, periodic_shift, idx_data);

    // Spread.
    t_spread_f77->start();

    const int local_indices_size = local_indices.size();
    if (local_indices_size > 0)
    {
        if (spread_fcn == "PIECEWISE_CUBIC")
        {
            LAGRANGIAN_PWCUBIC_SPREAD_F77(
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
            LAGRANGIAN_IB4_SPREAD_F77(
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
            LAGRANGIAN_IB6_SPREAD_F77(
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
            LAGRANGIAN_WIB4_SPREAD_F77(
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

void
LEInteractor::spreadReflectedForces(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data,
    const SAMRAI::tbox::Pointer<LNodeLevelData>& F_data,
    const SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const std::string& spread_fcn)
{
    t_spread->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!f_data.isNull());
    assert(!idx_data.isNull());
    assert(!patch.isNull());
    const int F_depth = F_data->getDepth();
    const int X_depth = X_data->getDepth();
    assert(F_depth == f_data->getDepth());
    assert(F_depth == NDIM);
    assert(X_depth == NDIM);
#endif
    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
    const SAMRAI::hier::Index<NDIM>& ilower = patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch_box.upper();

    const SAMRAI::hier::IntVector<NDIM>& f_gcw = f_data->getGhostCellWidth();
    const int depth = f_data->getDepth();

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const xLower = pgeom->getXLower();
    const double* const xUpper = pgeom->getXUpper();
    const double* const dx = pgeom->getDx();

    // Generate reflected position and force vectors.
    //
    // XXXX: THIS DOES NOT ACCOUNT FOR THE CASE THAT POINTS ARE NEAR EDGES OR
    // CORNERS!!!!
    const int stencil_size = getStencilSize(spread_fcn);
    const int ghosts = static_cast<int>(floor(0.5*double(stencil_size))+1);
    const SAMRAI::hier::Box<NDIM>& stencil_box = SAMRAI::hier::Box<NDIM>::grow(patch_box,ghosts);

    const double* const X_arr = &(*X_data)(0);
    const double* const F_arr = &(*F_data)(0);
    std::vector<double> reflected_X_data, reflected_F_data;
    std::vector<int> reflected_local_indices;
    int node_counter = 0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (int side = 0; side <= 1; ++side)
        {
            const double* const X_bdry = (side == 0 ? xLower : xUpper);
            if (pgeom->getTouchesRegularBoundary(axis,side))
            {
                for (LNodeIndexData2::Iterator it(box); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = *it;
                    SAMRAI::hier::Index<NDIM> reflected_i = i;
                    reflected_i(axis) = (side == 0
                                         ? ilower(axis) - (ilower(axis) - i(axis)    )
                                         : iupper(axis) - (iupper(axis) - i(axis) - 1));
                    if (stencil_box.contains(reflected_i))
                    {
                        const LNodeIndexSet& node_set = (*idx_data)(i);
                        for (LNodeIndexSet::const_iterator n = node_set.begin();
                             n != node_set.end(); ++n)
                        {
                            const LNodeIndexSet::value_type& node_idx = *n;
                            const int& petsc_idx = node_idx->getLocalPETScIndex();
                            const double* const X = &X_arr[NDIM*petsc_idx];
                            const double* const F = &F_arr[NDIM*petsc_idx];

                            // Reflect X and F across the appropriate boundary.
                            std::vector<double> X_reflect(X,X+NDIM);
                            std::vector<double> F_reflect(F,F+NDIM);
                            X_reflect[axis] = X_bdry[axis] - (X[axis] - X_bdry[axis]);
                            F_reflect[axis] = -F[axis];

                            // Insert the data into the vectors.
                            reflected_X_data.insert(reflected_X_data.end(),X_reflect.begin(),X_reflect.end());
                            reflected_F_data.insert(reflected_F_data.end(),F_reflect.begin(),F_reflect.end());
                            reflected_local_indices.push_back(node_counter);
                            ++node_counter;
                        }
                    }
                }
            }
        }
    }

    // Set the perodic offsets to equal zero.
    std::vector<double> periodic_offsets(NDIM*reflected_local_indices.size(),0.0);

    // Spread.
    t_spread_f77->start();

    const int reflected_local_indices_size = reflected_local_indices.size();
    if (reflected_local_indices_size > 0)
    {
        if (spread_fcn == "PIECEWISE_CUBIC")
        {
            LAGRANGIAN_PWCUBIC_SPREAD_F77(
                dx,xLower,xUpper,depth,
                &reflected_local_indices[0], &periodic_offsets[0], reflected_local_indices_size,
                &reflected_X_data[0], &reflected_F_data[0],
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                f_gcw(0),f_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                f_gcw(0),f_gcw(1),f_gcw(2),
#endif
                f_data->getPointer());
        }
        else if (spread_fcn == "IB_4")
        {
            LAGRANGIAN_IB4_SPREAD_F77(
                dx,xLower,xUpper,depth,
                &reflected_local_indices[0], &periodic_offsets[0], reflected_local_indices_size,
                &reflected_X_data[0], &reflected_F_data[0],
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                f_gcw(0),f_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                f_gcw(0),f_gcw(1),f_gcw(2),
#endif
                f_data->getPointer());
        }
        else if (spread_fcn == "IB_6")
        {
            LAGRANGIAN_IB6_SPREAD_F77(
                dx,xLower,xUpper,depth,
                &reflected_local_indices[0], &periodic_offsets[0], reflected_local_indices_size,
                &reflected_X_data[0], &reflected_F_data[0],
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                f_gcw(0),f_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                f_gcw(0),f_gcw(1),f_gcw(2),
#endif
                f_data->getPointer());
        }
        else if (spread_fcn == "WIDE_IB_4")
        {
            LAGRANGIAN_WIB4_SPREAD_F77(
                dx,xLower,xUpper,depth,
                &reflected_local_indices[0], &periodic_offsets[0], reflected_local_indices_size,
                &reflected_X_data[0], &reflected_F_data[0],
#if (NDIM == 2)
                ilower(0),iupper(0),ilower(1),iupper(1),
                f_gcw(0),f_gcw(1),
#endif
#if (NDIM == 3)
                ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                f_gcw(0),f_gcw(1),f_gcw(2),
#endif
                f_data->getPointer());
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
}// spreadReflectedForces

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
LEInteractor::getLocalIndices(
    std::vector<int>& local_indices,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data)
{
    const SAMRAI::hier::Box<NDIM>& patch_box = idx_data->getBox();
    const SAMRAI::hier::Box<NDIM>& ghost_box = idx_data->getGhostBox();

    if (box == patch_box)
    {
        local_indices = idx_data->getInteriorLocalIndices();
    }
    else if (box == ghost_box)
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
        for (LNodeIndexData2::Iterator it(box); it; it++)
        {
            const SAMRAI::hier::Index<NDIM>& i = *it;
            const LNodeIndexSet& node_set = (*idx_data)(i);
            const LNodeIndexSet::size_type& num_ids = node_set.size();
            if (num_ids > 0)
            {
                local_indices.resize(local_indices.size()+num_ids);
                std::transform(node_set.begin(), node_set.end(),
                               local_indices.end()-num_ids,
                               GetLocalPETScIndex());
            }
        }
    }
    return;
}// getLocalIndices

void
LEInteractor::getPeriodicOffsets(
    std::vector<double>& periodic_offsets,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
    const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
    const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data)
{
    const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
    const SAMRAI::hier::Index<NDIM>& ilower = patch_box.lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch_box.upper();

    const SAMRAI::hier::Box<NDIM>& ghost_box = idx_data->getGhostBox();

    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();

    if (pgeom->getTouchesPeriodicBoundary() && box != patch_box)
    {
        if (box == ghost_box)
        {
            periodic_offsets.resize(NDIM*idx_data->getInteriorLocalIndices().size());
            periodic_offsets.reserve(NDIM*(idx_data->getInteriorLocalIndices().size()+
                                           idx_data->getGhostLocalIndices().size()));
            SAMRAI::hier::BoxList<NDIM> ghost_boxes;
            ghost_boxes.removeIntersections(ghost_box,patch_box);
            for (SAMRAI::hier::BoxList<NDIM>::Iterator b(ghost_boxes); b; b++)
            {
                for (LNodeIndexData2::Iterator it(*b); it; it++)
                {
                    const SAMRAI::hier::Index<NDIM>& i = *it;
                    const LNodeIndexSet& node_set = (*idx_data)(i);
                    SAMRAI::hier::IntVector<NDIM> offset = 0;
                    static const int lower = 0;
                    static const int upper = 1;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if      (pgeom->getTouchesPeriodicBoundary(d,lower) && i(d) < ilower(d))
                        {
                            offset(d) = -periodic_shift(d);  // X is ABOVE the top of the patch --- need to shift DOWN
                        }
                        else if (pgeom->getTouchesPeriodicBoundary(d,upper) && i(d) > iupper(d))
                        {
                            offset(d) = +periodic_shift(d);  // X is BELOW the bottom of the patch --- need to shift UP
                        }
                    }
                    for (LNodeIndexSet::size_type n = 0; n < node_set.size(); ++n)
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            periodic_offsets.push_back(double(offset(d))*dx[d]);
                        }
                    }
                }
            }
        }
        else
        {
            periodic_offsets.resize(0);
            for (LNodeIndexData2::Iterator it(box); it; it++)
            {
                const SAMRAI::hier::Index<NDIM>& i = *it;
                const LNodeIndexSet& node_set = (*idx_data)(i);
                SAMRAI::hier::IntVector<NDIM> offset = 0;
                static const int lower = 0;
                static const int upper = 1;
                for (int d = 0; d < NDIM; ++d)
                {
                    if      (pgeom->getTouchesPeriodicBoundary(d,lower) && i(d) < ilower(d))
                    {
                        offset(d) = -periodic_shift(d);  // X is ABOVE the top of the patch --- need to shift DOWN
                    }
                    else if (pgeom->getTouchesPeriodicBoundary(d,upper) && i(d) > iupper(d))
                    {
                        offset(d) = +periodic_shift(d);  // X is BELOW the bottom of the patch --- need to shift UP
                    }
                }
                for (LNodeIndexSet::size_type n = 0; n < node_set.size(); ++n)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        periodic_offsets.push_back(double(offset(d))*dx[d]);
                    }
                }
            }
        }
    }
    return;
}// getPeriodicOffsets

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
