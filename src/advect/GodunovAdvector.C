// Filename: GodunovAdvector.C
// Last modified: <12.Feb.2008 21:19:24 griffith@box221.cims.nyu.edu>
// Created on 14 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "GodunovAdvector.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <ArrayData.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 1)
#define ADVECT_DERIVATIVE_F77 F77_FUNC_(advect_derivative1d, ADVECT_DERIVATIVE1D)
#define ADVECT_FLUX_F77 F77_FUNC_(advect_flux1d, ADVECT_FLUX1D)
#define ADVECT_STABLEDT_F77 F77_FUNC_(advect_stabledt1d, ADVECT_STABLEDT1D)
#define GODUNOV_PREDICT_F77 F77_FUNC_(godunov_predict1d, GODUNOV_PREDICT1D)
#define GODUNOV_PREDICTWITHSOURCE_F77 F77_FUNC_(godunov_predictwithsource1d, GODUNOV_PREDICTWITHSOURCE1D)
#endif

#if (NDIM == 2)
#define ADVECT_DERIVATIVE_F77 F77_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define ADVECT_FLUX_F77 F77_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define ADVECT_STABLEDT_F77 F77_FUNC_(advect_stabledt2d, ADVECT_STABLEDT2D)
#define GODUNOV_INCOMPRESSIBILITYFIX_F77 F77_FUNC_(godunov_incompressibilityfix2d, GODUNOV_INCOMPRESSIBILITYFIX2D)
#define GODUNOV_PREDICT_F77 F77_FUNC_(godunov_predict2d, GODUNOV_PREDICT2D)
#define GODUNOV_PREDICTWITHSOURCE_F77 F77_FUNC_(godunov_predictwithsource2d, GODUNOV_PREDICTWITHSOURCE2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_F77 F77_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define ADVECT_FLUX_F77 F77_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define ADVECT_STABLEDT_F77 F77_FUNC_(advect_stabledt3d, ADVECT_STABLEDT3D)
#define GODUNOV_INCOMPRESSIBILITYFIX_F77 F77_FUNC_(godunov_incompressibilityfix3d, GODUNOV_INCOMPRESSIBILITYFIX3D)
#define GODUNOV_PREDICT_F77 F77_FUNC_(godunov_predict3d, GODUNOV_PREDICT3D)
#define GODUNOV_PREDICTWITHSOURCE_F77 F77_FUNC_(godunov_predictwithsource3d, GODUNOV_PREDICTWITHSOURCE3D)
#endif

extern "C"
{
    void
    ADVECT_DERIVATIVE_F77(
        const double*,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const double* ,
        const double* ,
        const int& ,
        double*
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        const int& , const int& ,
        double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        const int& , const int& , const int& ,
        double*
#endif
                               );

    void
    ADVECT_FLUX_F77(
        const double& ,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const int& ,
        const double* ,
        const double* ,
        double*
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        const double* , const double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        const double* , const double* , const double* ,
        double* , double* , double*
#endif
                         );

    void
    ADVECT_STABLEDT_F77(
        const double*,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const double* ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double&
                             );

#if ((NDIM == 2) || (NDIM == 3))
    void
    GODUNOV_INCOMPRESSIBILITYFIX_F77(
        const int& ,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double*
#endif
                                          );
#endif

    void
    GODUNOV_PREDICT_F77(
        const double* , const double& ,
#if (NDIM == 3)
        const unsigned int& ,
#endif
        const int& ,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const double* ,
        const int& ,
        const int& ,
        const double* ,
        double*
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , double* , double* ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double* ,
        double* , double* , double*
#endif
                             );

    void
    GODUNOV_PREDICTWITHSOURCE_F77(
        const double* , const double& ,
#if (NDIM == 3)
        const unsigned int& ,
#endif
        const int& ,
#if (NDIM == 1)
        const int& , const int& ,
        const int& ,
        const int& ,
        const double* ,
        const double* ,
        const int& ,
        const int& ,
        const double* ,
        double*
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , double* ,
        const double* , double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , double* , double* ,
        const double* , double* , double* ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double* ,
        double* , double* , double*
#endif
                                       );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_stable_dt_on_patch;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_advective_derivative;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_flux;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict_value;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict_value_with_source_term;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict_normal_velocity;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict_normal_velocity_with_source_term;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_enforce_incompressibility;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_predict_with_source_term;

// Number of ghosts cells used for each variable quantity.
static const int FACEG = 1;

// Version of GodunovAdvector restart file data
static const int GODUNOV_ADVECTOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

GodunovAdvector::GodunovAdvector(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    const bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_limiter_type(FOURTH_ORDER)
#if (NDIM==3)
    , d_using_full_ctu(true)
#endif
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from given input/restart databases.
    bool is_from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (is_from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, is_from_restart);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_stable_dt_on_patch = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::computeStableDtOnPatch()");
        t_compute_advective_derivative = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::computeAdvectiveDerivative()");
        t_compute_flux = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::computeFlux()");
        t_predict_value = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::predictValue()");
        t_predict_value_with_source_term = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::predictValueWithSourceTerm()");
        t_predict_normal_velocity = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::predictNormalVelocity()");
        t_predict_normal_velocity_with_source_term = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::predictNormalVelocityWithSourceTerm()");
        t_enforce_incompressibility = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::enforceIncompressibility()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::putToDatabase()");
        t_predict = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::predict()");
        t_predict_with_source_term = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::GodunovAdvector::predictWithSourceTerm()");
        timers_need_init = false;
    }
    return;
}// GodunovAdvector

GodunovAdvector::~GodunovAdvector()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    return;
}// ~GodunovAdvector

const std::string&
GodunovAdvector::getName() const
{
    return d_object_name;
}// getName

double
GodunovAdvector::computeStableDtOnPatch(
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::hier::Patch<NDIM>& patch) const
{
    t_compute_stable_dt_on_patch->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());
#endif
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::IntVector<NDIM>& u_ghost_cells = u_ADV.getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();

#if (NDIM == 1)
    ADVECT_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),
        u_ghost_cells(0),
        u_ADV.getPointer(0),
        stable_dt);
#endif
#if (NDIM == 2)
    ADVECT_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        u_ghost_cells(0),u_ghost_cells(1),
        u_ADV.getPointer(0),u_ADV.getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    ADVECT_STABLEDT_F77(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        u_ghost_cells(0),u_ghost_cells(1),u_ghost_cells(2),
        u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
        stable_dt);
#endif

    t_compute_stable_dt_on_patch->stop();
    return stable_dt;
}// computeStableDtOnPatch

void
GodunovAdvector::computeAdvectiveDerivative(
    SAMRAI::pdat::CellData<NDIM,double>& N,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::FaceData<NDIM,double>& q_half,
    const SAMRAI::hier::Patch<NDIM>& patch) const
{
    t_compute_advective_derivative->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(N.getDepth() == q_half.getDepth());
    TBOX_ASSERT(N.getBox()   == patch.getBox());

    TBOX_ASSERT(q_half.getBox() == patch.getBox());
#endif
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::IntVector<NDIM>&  u_ADV_ghost_cells = u_ADV .getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>&      N_ghost_cells = N     .getGhostCellWidth();

    for (int depth = 0; depth < q_half.getDepth() ; ++depth)
    {
#if (NDIM == 1)
        ADVECT_DERIVATIVE_F77(
            dx,
            ilower(0),iupper(0),
            u_ADV_ghost_cells(0),
            q_half_ghost_cells(0),
            u_ADV.getPointer(0),
            q_half.getPointer(0,depth),
            N_ghost_cells(0),
            N.getPointer(depth));
#endif
#if (NDIM == 2)
        ADVECT_DERIVATIVE_F77(
            dx,
            ilower(0),iupper(0),ilower(1),iupper(1),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),
            q_half_ghost_cells(0),q_half_ghost_cells(1),
            u_ADV.getPointer(0),u_ADV.getPointer(1),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),
            N_ghost_cells(0),N_ghost_cells(1),
            N.getPointer(depth));
#endif
#if (NDIM == 3)
        ADVECT_DERIVATIVE_F77(
            dx,
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),u_ADV_ghost_cells(2),
            q_half_ghost_cells(0),q_half_ghost_cells(1),q_half_ghost_cells(2),
            u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),q_half.getPointer(2,depth),
            N_ghost_cells(0),N_ghost_cells(1),N_ghost_cells(2),
            N.getPointer(depth));
#endif
    }

    t_compute_advective_derivative->stop();
    return;
}// computeAdvectiveDerivative

void
GodunovAdvector::computeFlux(
    SAMRAI::pdat::FaceData<NDIM,double>& flux,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::FaceData<NDIM,double>& q_half,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_compute_flux->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(flux.getDepth() == q_half.getDepth());
    TBOX_ASSERT(flux.getBox()   == patch.getBox());

    TBOX_ASSERT(q_half.getBox() == patch.getBox());
#endif
    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::IntVector<NDIM>&   flux_ghost_cells = flux  .getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>&  u_ADV_ghost_cells = u_ADV .getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();

    for (int depth = 0; depth < q_half.getDepth(); ++depth)
    {
#if (NDIM == 1)
        ADVECT_FLUX_F77(
            dt,
            ilower(0),iupper(0),
            u_ADV_ghost_cells(0),
            q_half_ghost_cells(0),
            flux_ghost_cells(0),
            u_ADV.getPointer(0),
            q_half.getPointer(0,depth),
            flux.getPointer(0,depth));
#endif
#if (NDIM == 2)
        ADVECT_FLUX_F77(
            dt,
            ilower(0),iupper(0),ilower(1),iupper(1),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),
            q_half_ghost_cells(0),q_half_ghost_cells(1),
            flux_ghost_cells(0),flux_ghost_cells(1),
            u_ADV.getPointer(0),u_ADV.getPointer(1),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),
            flux.getPointer(0,depth),flux.getPointer(1,depth));
#endif
#if (NDIM == 3)
        ADVECT_FLUX_F77(
            dt,
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),u_ADV_ghost_cells(2),
            q_half_ghost_cells(0),q_half_ghost_cells(1),q_half_ghost_cells(2),
            flux_ghost_cells(0),flux_ghost_cells(1),flux_ghost_cells(2),
            u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),q_half.getPointer(2,depth),
            flux.getPointer(0,depth),flux.getPointer(1,depth),flux.getPointer(2,depth));
#endif
    }

    t_compute_flux->stop();
    return;
}// computeFlux

void
GodunovAdvector::predictValue(
    SAMRAI::pdat::FaceData<NDIM,double>& q_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::CellData<NDIM,double>& Q,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_predict_value->start();

    predict(q_half, u_ADV, Q, patch, dt);

    t_predict_value->stop();
    return;
}// predictValue

void
GodunovAdvector::predictValueWithSourceTerm(
    SAMRAI::pdat::FaceData<NDIM,double>& q_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::CellData<NDIM,double>& Q,
    const SAMRAI::pdat::CellData<NDIM,double>& F,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_predict_value_with_source_term->start();

    predictWithSourceTerm(q_half, u_ADV, Q, F, patch, dt);

    t_predict_value_with_source_term->stop();
    return;
}// predictValueWithSourceTerm

void
GodunovAdvector::predictNormalVelocity(
    SAMRAI::pdat::FaceData<NDIM,double>& v_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::CellData<NDIM,double>& V,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_predict_normal_velocity->start();

    SAMRAI::pdat::FaceData<NDIM,double> v_half_tmp(v_half.getBox(), NDIM, SAMRAI::hier::IntVector<NDIM>(FACEG));

    predict(v_half_tmp, u_ADV, V, patch, dt);

    for (int axis = 0; axis < NDIM; ++axis)
    {
        SAMRAI::pdat::ArrayData<NDIM,double>& v_half_arr = v_half.getArrayData(axis);
        const SAMRAI::pdat::ArrayData<NDIM,double>& v_half_tmp_arr =
            v_half_tmp.getArrayData(axis);
        const SAMRAI::hier::Box<NDIM> box = (v_half_arr.getBox())*(v_half_tmp_arr.getBox());
        v_half_arr.copyDepth(0, v_half_tmp_arr, axis, box);
    }

    t_predict_normal_velocity->stop();
    return;
}// predictNormalVelocity

void
GodunovAdvector::predictNormalVelocityWithSourceTerm(
    SAMRAI::pdat::FaceData<NDIM,double>& v_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::CellData<NDIM,double>& V,
    const SAMRAI::pdat::CellData<NDIM,double>& F,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_predict_normal_velocity_with_source_term->start();

    SAMRAI::pdat::FaceData<NDIM,double> v_half_tmp(v_half.getBox(), NDIM, SAMRAI::hier::IntVector<NDIM>(FACEG));

    predictWithSourceTerm(v_half_tmp, u_ADV, V, F, patch, dt);

    for (int axis = 0; axis < NDIM; ++axis)
    {
        SAMRAI::pdat::ArrayData<NDIM,double>& v_half_arr = v_half.getArrayData(axis);
        const SAMRAI::pdat::ArrayData<NDIM,double>& v_half_tmp_arr =
            v_half_tmp.getArrayData(axis);
        const SAMRAI::hier::Box<NDIM> box = (v_half_arr.getBox())*(v_half_tmp_arr.getBox());
        v_half_arr.copyDepth(0, v_half_tmp_arr, axis, box);
    }

    t_predict_normal_velocity_with_source_term->stop();
    return;
}// predictNormalVelocityWithSourceTerm

void
GodunovAdvector::enforceIncompressibility(
    SAMRAI::pdat::FaceData<NDIM,double>& v_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::FaceData<NDIM,double>& grad_phi,
    const SAMRAI::hier::Patch<NDIM>& patch) const
{
    t_enforce_incompressibility->start();

#if (NDIM != 1)

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(grad_phi.getDepth() == 1);
    TBOX_ASSERT(grad_phi.getBox()   == patch.getBox());

    TBOX_ASSERT(v_half.getBox()   == patch.getBox());
    TBOX_ASSERT(v_half.getDepth() == NDIM);
#endif
    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::IntVector<NDIM>& grad_phi_ghost_cells = grad_phi.getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>&   v_half_ghost_cells = v_half  .getGhostCellWidth();

    for (int depth = 0; depth < NDIM; ++depth)
    {
#if (NDIM == 2)
        GODUNOV_INCOMPRESSIBILITYFIX_F77(
            depth,
            ilower(0),iupper(0),ilower(1),iupper(1),
            grad_phi_ghost_cells(0),grad_phi_ghost_cells(1),
            v_half_ghost_cells(0),v_half_ghost_cells(1),
            grad_phi.getPointer(0),grad_phi.getPointer(1),
            v_half.getPointer(0,depth),v_half.getPointer(1,depth));
#endif
#if (NDIM == 3)
        GODUNOV_INCOMPRESSIBILITYFIX_F77(
            depth,
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            grad_phi_ghost_cells(0),grad_phi_ghost_cells(1),grad_phi_ghost_cells(2),
            v_half_ghost_cells(0),v_half_ghost_cells(1),v_half_ghost_cells(2),
            grad_phi.getPointer(0),grad_phi.getPointer(1),grad_phi.getPointer(2),
            v_half.getPointer(0,depth),v_half.getPointer(1,depth),v_half.getPointer(2,depth));
#endif
    }

#endif

    t_enforce_incompressibility->stop();
    return;
}// enforceIncompressibility

void
GodunovAdvector::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("GODUNOV_ADVECTOR_VERSION",GODUNOV_ADVECTOR_VERSION);
    db->putInteger("d_limiter_type", d_limiter_type);
#if (NDIM == 3)
    db->putBool("d_using_full_ctu", d_using_full_ctu);
#endif

    t_put_to_database->stop();
    return;
}// putToDatabase

void
GodunovAdvector::printClassData(
    std::ostream &os) const
{
    os << "\nGodunovAdvector::printClassData..." << std::endl;
    os << "GodunovAdvector: this = " << const_cast<GodunovAdvector*>(this) << "\n";
    os << "d_object_name = " << d_object_name << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart << std::endl;
    os << "d_limiter_type = "<< d_limiter_type << "\n";
#if (NDIM == 3)
    os << "d_using_full_ctu = " << d_using_full_ctu << std::endl;
#endif

    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
GodunovAdvector::predict(
    SAMRAI::pdat::FaceData<NDIM,double>& q_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::CellData<NDIM,double>& Q,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_predict->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(q_half.getDepth() == Q.getDepth());
    TBOX_ASSERT(q_half.getBox()   == patch.getBox());

    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(Q.getBox()        == patch.getBox());
#endif
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom =
        patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::IntVector<NDIM>&  u_ADV_ghost_cells = u_ADV .getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>&      Q_ghost_cells = Q     .getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();

#if (NDIM>1)
    SAMRAI::pdat::CellData<NDIM,double> Q_temp1(patch.getBox(),1,Q_ghost_cells);
    SAMRAI::pdat::FaceData<NDIM,double> q_half_temp(patch.getBox(),1,q_half_ghost_cells);
#endif

#if (NDIM>2)
    SAMRAI::pdat::CellData<NDIM,double> Q_temp2(patch.getBox(),1,Q_ghost_cells);
#endif

    for (int depth = 0; depth < Q.getDepth(); ++depth)
    {
#if (NDIM == 1)
        GODUNOV_PREDICT_F77(
            dx,dt,
            d_limiter_type,
            ilower(0),iupper(0),
            Q_ghost_cells(0),
            Q.getPointer(depth),
            u_ADV_ghost_cells(0),
            q_half_ghost_cells(0),
            u_ADV.getPointer(0),
            q_half.getPointer(0,depth));
#endif
#if (NDIM == 2)
        GODUNOV_PREDICT_F77(
            dx,dt,
            d_limiter_type,
            ilower(0),iupper(0),ilower(1),iupper(1),
            Q_ghost_cells(0),Q_ghost_cells(1),
            Q.getPointer(depth),Q_temp1.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),
            q_half_ghost_cells(0),q_half_ghost_cells(1),
            u_ADV.getPointer(0),u_ADV.getPointer(1),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth));
#endif
#if (NDIM == 3)
        GODUNOV_PREDICT_F77(
            dx,dt,
            static_cast<unsigned int>(d_using_full_ctu),
            d_limiter_type,
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            Q_ghost_cells(0),Q_ghost_cells(1),Q_ghost_cells(2),
            Q.getPointer(depth),Q_temp1.getPointer(0),Q_temp2.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),u_ADV_ghost_cells(2),
            q_half_ghost_cells(0),q_half_ghost_cells(1),q_half_ghost_cells(2),
            u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),q_half_temp.getPointer(2),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),q_half.getPointer(2,depth));
#endif
    }

    t_predict->stop();
    return;
}// predict

void
GodunovAdvector::predictWithSourceTerm(
    SAMRAI::pdat::FaceData<NDIM,double>& q_half,
    const SAMRAI::pdat::FaceData<NDIM,double>& u_ADV,
    const SAMRAI::pdat::CellData<NDIM,double>& Q,
    const SAMRAI::pdat::CellData<NDIM,double>& F,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const double dt) const
{
    t_predict_with_source_term->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(q_half.getDepth() == Q.getDepth());
    TBOX_ASSERT(q_half.getDepth() == F.getDepth());
    TBOX_ASSERT(q_half.getBox()   == patch.getBox());

    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(Q.getBox() == patch.getBox());

    TBOX_ASSERT(F.getBox() == patch.getBox());
#endif
    const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const SAMRAI::hier::Index<NDIM>& ilower = patch.getBox().lower();
    const SAMRAI::hier::Index<NDIM>& iupper = patch.getBox().upper();

    const SAMRAI::hier::IntVector<NDIM>& u_ADV_ghost_cells = u_ADV.getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>& Q_ghost_cells = Q.getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>& F_ghost_cells = F.getGhostCellWidth();
    const SAMRAI::hier::IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();

#if (NDIM>1)
    SAMRAI::pdat::CellData<NDIM,double> Q_temp1(patch.getBox(),1,Q_ghost_cells);
    SAMRAI::pdat::CellData<NDIM,double> F_temp1(patch.getBox(),1,F_ghost_cells);
    SAMRAI::pdat::FaceData<NDIM,double> q_half_temp(patch.getBox(),1,q_half_ghost_cells);
#endif

#if (NDIM>2)
    SAMRAI::pdat::CellData<NDIM,double> Q_temp2(patch.getBox(),1,Q_ghost_cells);
    SAMRAI::pdat::CellData<NDIM,double> F_temp2(patch.getBox(),1,F_ghost_cells);
#endif

    for (int depth = 0; depth < Q.getDepth(); ++depth)
    {
#if (NDIM == 1)
        GODUNOV_PREDICTWITHSOURCE_F77(
            dx,dt,
            d_limiter_type,
            ilower(0),iupper(0),
            Q_ghost_cells(0),
            F_ghost_cells(0),
            Q.getPointer(depth),
            F.getPointer(depth),
            u_ADV_ghost_cells(0),
            q_half_ghost_cells(0),
            u_ADV.getPointer(0),
            q_half.getPointer(0,depth));
#endif
#if (NDIM == 2)
        GODUNOV_PREDICTWITHSOURCE_F77(
            dx,dt,
            d_limiter_type,
            ilower(0),iupper(0),ilower(1),iupper(1),
            Q_ghost_cells(0),Q_ghost_cells(1),
            F_ghost_cells(0),F_ghost_cells(1),
            Q.getPointer(depth),Q_temp1.getPointer(0),
            F.getPointer(depth),F_temp1.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),
            q_half_ghost_cells(0),q_half_ghost_cells(1),
            u_ADV.getPointer(0),u_ADV.getPointer(1),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth));
#endif
#if (NDIM == 3)
        GODUNOV_PREDICTWITHSOURCE_F77(
            dx,dt,
            static_cast<unsigned int>(d_using_full_ctu),
            d_limiter_type,
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            Q_ghost_cells(0),Q_ghost_cells(1),Q_ghost_cells(2),
            F_ghost_cells(0),F_ghost_cells(1),F_ghost_cells(2),
            Q.getPointer(depth),Q_temp1.getPointer(0),Q_temp2.getPointer(0),
            F.getPointer(depth),F_temp1.getPointer(0),F_temp2.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),u_ADV_ghost_cells(2),
            q_half_ghost_cells(0),q_half_ghost_cells(1),q_half_ghost_cells(2),
            u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),q_half_temp.getPointer(2),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),q_half.getPointer(2,depth));
#endif
    }

    t_predict_with_source_term->stop();
    return;
}// predictWithSourceTerm

void
GodunovAdvector::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    (void) is_from_restart;
    d_limiter_type = db->getIntegerWithDefault("limiter_type", d_limiter_type);
#if (NDIM == 3)
    d_using_full_ctu = db->getBoolWithDefault("using_full_ctu", d_using_full_ctu);
#endif
    return;
}// getFromInput

void
GodunovAdvector::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> root_db =
        SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;

    if (root_db->isDatabase(d_object_name))
    {
        db = root_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::getFromRestart():\n"
                   << "  Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("GODUNOV_ADVECTOR_VERSION");
    if (ver != GODUNOV_ADVECTOR_VERSION)
    {
        TBOX_ERROR(d_object_name << "::getFromRestart():\n"
                   << "  Restart file version different than class version.");
    }

    d_limiter_type = db->getInteger("d_limiter_type");
#if (NDIM == 3)
    d_using_full_ctu = db->getBool("d_using_full_ctu");
#endif

    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::GodunovAdvector>;

//////////////////////////////////////////////////////////////////////////////
