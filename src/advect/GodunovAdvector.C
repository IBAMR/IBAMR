// Filename: GodunovAdvector.C
// Created on 14 Feb 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

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

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <ArrayData.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <Index.h>
#include <IntVector.h>
#include <tbox/RestartManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative2d, ADVECT_DERIVATIVE2D)
#define ADVECT_FLUX_FC FC_FUNC_(advect_flux2d, ADVECT_FLUX2D)
#define ADVECT_STABLEDT_FC FC_FUNC_(advect_stabledt2d, ADVECT_STABLEDT2D)
#define GODUNOV_INCOMPRESSIBILITY_FIX_FC FC_FUNC_(godunov_incompressibility_fix2d, GODUNOV_INCOMPRESSIBILITY_FIX2D)
#define GODUNOV_PREDICT_FC FC_FUNC_(godunov_predict2d, GODUNOV_PREDICT2D)
#define GODUNOV_PREDICT_WITH_SOURCE_FC FC_FUNC_(godunov_predict_with_source2d, GODUNOV_PREDICT_WITH_SOURCE2D)
#endif

#if (NDIM == 3)
#define ADVECT_DERIVATIVE_FC FC_FUNC_(advect_derivative3d, ADVECT_DERIVATIVE3D)
#define ADVECT_FLUX_FC FC_FUNC_(advect_flux3d, ADVECT_FLUX3D)
#define ADVECT_STABLEDT_FC FC_FUNC_(advect_stabledt3d, ADVECT_STABLEDT3D)
#define GODUNOV_INCOMPRESSIBILITY_FIX_FC FC_FUNC_(godunov_incompressibility_fix3d, GODUNOV_INCOMPRESSIBILITY_FIX3D)
#define GODUNOV_PREDICT_FC FC_FUNC_(godunov_predict3d, GODUNOV_PREDICT3D)
#define GODUNOV_PREDICT_WITH_SOURCE_FC FC_FUNC_(godunov_predict_with_source3d, GODUNOV_PREDICT_WITH_SOURCE3D)
#endif

extern "C"
{
    void
    ADVECT_DERIVATIVE_FC(
        const double*,
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
    ADVECT_FLUX_FC(
        const double& ,
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
    ADVECT_STABLEDT_FC(
        const double*,
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
    GODUNOV_INCOMPRESSIBILITY_FIX_FC(
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
    GODUNOV_PREDICT_FC(
        const double* , const double& ,
#if (NDIM == 3)
        const unsigned int& ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , double* , double* , double* , double* ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
        double* , double* ,
        double* , double*
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , double* , double* , double* , double* , double* ,
        const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
        double* , double* , double* ,
        double* , double* , double*
#endif
                       );

    void
    GODUNOV_PREDICT_WITH_SOURCE_FC(
        const double* , const double& ,
#if (NDIM == 3)
        const unsigned int& ,
#endif
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const int& , const int& ,
        const double* , double* , double* , double* , double* ,
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
        const double* , double* , double* , double* , double* , double* ,
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
// Number of ghosts cells used for each variable quantity.
static const int FACEG = 1;

// Version of GodunovAdvector restart file data
static const int GODUNOV_ADVECTOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

GodunovAdvector::GodunovAdvector(
    const std::string& object_name,
    Pointer<Database> input_db,
    const bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart)
#if (NDIM == 3)
    , d_using_full_ctu(true)
#endif
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from given input/restart databases.
    bool is_from_restart = RestartManager::getManager()->isFromRestart();
    if (is_from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, is_from_restart);
    return;
}// GodunovAdvector

GodunovAdvector::~GodunovAdvector()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
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
    const FaceData<NDIM,double>& u_ADV,
    const Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());
#endif
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
        patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const IntVector<NDIM>& u_ghost_cells = u_ADV.getGhostCellWidth();

    double stable_dt = std::numeric_limits<double>::max();

#if (NDIM == 2)
    ADVECT_STABLEDT_FC(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),
        u_ghost_cells(0),u_ghost_cells(1),
        u_ADV.getPointer(0),u_ADV.getPointer(1),
        stable_dt);
#endif
#if (NDIM == 3)
    ADVECT_STABLEDT_FC(
        dx,
        ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
        u_ghost_cells(0),u_ghost_cells(1),u_ghost_cells(2),
        u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
        stable_dt);
#endif
    return stable_dt;
}// computeStableDtOnPatch

void
GodunovAdvector::computeAdvectiveDerivative(
    CellData<NDIM,double>& N,
    const FaceData<NDIM,double>& u_ADV,
    const FaceData<NDIM,double>& q_half,
    const Patch<NDIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(N.getDepth() == q_half.getDepth());
    TBOX_ASSERT(N.getBox()   == patch.getBox());

    TBOX_ASSERT(q_half.getBox() == patch.getBox());
#endif
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
        patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const IntVector<NDIM>&  u_ADV_ghost_cells = u_ADV .getGhostCellWidth();
    const IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();
    const IntVector<NDIM>&      N_ghost_cells = N     .getGhostCellWidth();

    for (int depth = 0; depth < q_half.getDepth() ; ++depth)
    {
#if (NDIM == 2)
        ADVECT_DERIVATIVE_FC(
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
        ADVECT_DERIVATIVE_FC(
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
    return;
}// computeAdvectiveDerivative

void
GodunovAdvector::computeFlux(
    FaceData<NDIM,double>& flux,
    const FaceData<NDIM,double>& u_ADV,
    const FaceData<NDIM,double>& q_half,
    const Patch<NDIM>& patch,
    const double dt) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(flux.getDepth() == q_half.getDepth());
    TBOX_ASSERT(flux.getBox()   == patch.getBox());

    TBOX_ASSERT(q_half.getBox() == patch.getBox());
#endif
    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const IntVector<NDIM>&   flux_ghost_cells = flux  .getGhostCellWidth();
    const IntVector<NDIM>&  u_ADV_ghost_cells = u_ADV .getGhostCellWidth();
    const IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();

    for (int depth = 0; depth < q_half.getDepth(); ++depth)
    {
#if (NDIM == 2)
        ADVECT_FLUX_FC(
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
        ADVECT_FLUX_FC(
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
    return;
}// computeFlux

void
GodunovAdvector::predictValue(
    FaceData<NDIM,double>& q_half,
    const FaceData<NDIM,double>& u_ADV,
    const CellData<NDIM,double>& Q,
    const Patch<NDIM>& patch,
    const double dt) const
{
    predict(q_half, u_ADV, Q, patch, dt);
    return;
}// predictValue

void
GodunovAdvector::predictValueWithSourceTerm(
    FaceData<NDIM,double>& q_half,
    const FaceData<NDIM,double>& u_ADV,
    const CellData<NDIM,double>& Q,
    const CellData<NDIM,double>& F,
    const Patch<NDIM>& patch,
    const double dt) const
{
    predictWithSourceTerm(q_half, u_ADV, Q, F, patch, dt);
    return;
}// predictValueWithSourceTerm

void
GodunovAdvector::predictNormalVelocity(
    FaceData<NDIM,double>& v_half,
    const FaceData<NDIM,double>& u_ADV,
    const CellData<NDIM,double>& V,
    const Patch<NDIM>& patch,
    const double dt) const
{
    FaceData<NDIM,double> v_half_tmp(v_half.getBox(), NDIM, IntVector<NDIM>(FACEG));

    predict(v_half_tmp, u_ADV, V, patch, dt);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        ArrayData<NDIM,double>& v_half_arr = v_half.getArrayData(axis);
        const ArrayData<NDIM,double>& v_half_tmp_arr =
            v_half_tmp.getArrayData(axis);
        const Box<NDIM> box = (v_half_arr.getBox())*(v_half_tmp_arr.getBox());
        v_half_arr.copyDepth(0, v_half_tmp_arr, axis, box);
    }
    return;
}// predictNormalVelocity

void
GodunovAdvector::predictNormalVelocityWithSourceTerm(
    FaceData<NDIM,double>& v_half,
    const FaceData<NDIM,double>& u_ADV,
    const CellData<NDIM,double>& V,
    const CellData<NDIM,double>& F,
    const Patch<NDIM>& patch,
    const double dt) const
{
    FaceData<NDIM,double> v_half_tmp(v_half.getBox(), NDIM, IntVector<NDIM>(FACEG));

    predictWithSourceTerm(v_half_tmp, u_ADV, V, F, patch, dt);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        ArrayData<NDIM,double>& v_half_arr = v_half.getArrayData(axis);
        const ArrayData<NDIM,double>& v_half_tmp_arr =
            v_half_tmp.getArrayData(axis);
        const Box<NDIM> box = (v_half_arr.getBox())*(v_half_tmp_arr.getBox());
        v_half_arr.copyDepth(0, v_half_tmp_arr, axis, box);
    }
    return;
}// predictNormalVelocityWithSourceTerm

void
GodunovAdvector::enforceIncompressibility(
    FaceData<NDIM,double>& v_half,
    const FaceData<NDIM,double>& u_ADV,
    const FaceData<NDIM,double>& grad_phi,
    const Patch<NDIM>& patch) const
{
#if (NDIM != 1)

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(grad_phi.getDepth() == 1);
    TBOX_ASSERT(grad_phi.getBox()   == patch.getBox());

    TBOX_ASSERT(v_half.getBox()   == patch.getBox());
    TBOX_ASSERT(v_half.getDepth() == NDIM);
#else
    NULL_USE(u_ADV);
#endif
    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const IntVector<NDIM>& grad_phi_ghost_cells = grad_phi.getGhostCellWidth();
    const IntVector<NDIM>&   v_half_ghost_cells = v_half  .getGhostCellWidth();

    for (unsigned int depth = 0; depth < NDIM; ++depth)
    {
#if (NDIM == 2)
        GODUNOV_INCOMPRESSIBILITY_FIX_FC(
            depth,
            ilower(0),iupper(0),ilower(1),iupper(1),
            grad_phi_ghost_cells(0),grad_phi_ghost_cells(1),
            v_half_ghost_cells(0),v_half_ghost_cells(1),
            grad_phi.getPointer(0),grad_phi.getPointer(1),
            v_half.getPointer(0,depth),v_half.getPointer(1,depth));
#endif
#if (NDIM == 3)
        GODUNOV_INCOMPRESSIBILITY_FIX_FC(
            depth,
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            grad_phi_ghost_cells(0),grad_phi_ghost_cells(1),grad_phi_ghost_cells(2),
            v_half_ghost_cells(0),v_half_ghost_cells(1),v_half_ghost_cells(2),
            grad_phi.getPointer(0),grad_phi.getPointer(1),grad_phi.getPointer(2),
            v_half.getPointer(0,depth),v_half.getPointer(1,depth),v_half.getPointer(2,depth));
#endif
    }

#endif
    return;
}// enforceIncompressibility

void
GodunovAdvector::putToDatabase(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("GODUNOV_ADVECTOR_VERSION",GODUNOV_ADVECTOR_VERSION);
#if (NDIM == 3)
    db->putBool("d_using_full_ctu", d_using_full_ctu);
#endif
    return;
}// putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
GodunovAdvector::predict(
    FaceData<NDIM,double>& q_half,
    const FaceData<NDIM,double>& u_ADV,
    const CellData<NDIM,double>& Q,
    const Patch<NDIM>& patch,
    const double dt) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(q_half.getDepth() == Q.getDepth());
    TBOX_ASSERT(q_half.getBox()   == patch.getBox());

    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(Q.getBox()        == patch.getBox());
#endif
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
        patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const IntVector<NDIM>&  u_ADV_ghost_cells = u_ADV .getGhostCellWidth();
    const IntVector<NDIM>&      Q_ghost_cells = Q     .getGhostCellWidth();
    const IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();

    CellData<NDIM,double> dQ(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> Q_L(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> Q_R(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> Q_temp1(patch.getBox(),1,Q_ghost_cells);
    FaceData<NDIM,double> q_half_temp(patch.getBox(),1,q_half_ghost_cells);
#if (NDIM > 2)
    CellData<NDIM,double> Q_temp2(patch.getBox(),1,Q_ghost_cells);
#endif

    for (int depth = 0; depth < Q.getDepth(); ++depth)
    {
#if (NDIM == 2)
        GODUNOV_PREDICT_FC(
            dx,dt,
            ilower(0),iupper(0),ilower(1),iupper(1),
            Q_ghost_cells(0),Q_ghost_cells(1),
            Q.getPointer(depth),Q_temp1.getPointer(0),
            dQ.getPointer(0),Q_L.getPointer(0),Q_R.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),
            q_half_ghost_cells(0),q_half_ghost_cells(1),
            u_ADV.getPointer(0),u_ADV.getPointer(1),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth));
#endif
#if (NDIM == 3)
        GODUNOV_PREDICT_FC(
            dx,dt,
            static_cast<unsigned int>(d_using_full_ctu),
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            Q_ghost_cells(0),Q_ghost_cells(1),Q_ghost_cells(2),
            Q.getPointer(depth),Q_temp1.getPointer(0),Q_temp2.getPointer(0),
            dQ.getPointer(0),Q_L.getPointer(0),Q_R.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),u_ADV_ghost_cells(2),
            q_half_ghost_cells(0),q_half_ghost_cells(1),q_half_ghost_cells(2),
            u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),q_half_temp.getPointer(2),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),q_half.getPointer(2,depth));
#endif
    }
    return;
}// predict

void
GodunovAdvector::predictWithSourceTerm(
    FaceData<NDIM,double>& q_half,
    const FaceData<NDIM,double>& u_ADV,
    const CellData<NDIM,double>& Q,
    const CellData<NDIM,double>& F,
    const Patch<NDIM>& patch,
    const double dt) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(q_half.getDepth() == Q.getDepth());
    TBOX_ASSERT(q_half.getDepth() == F.getDepth());
    TBOX_ASSERT(q_half.getBox()   == patch.getBox());

    TBOX_ASSERT(u_ADV.getDepth() == 1);
    TBOX_ASSERT(u_ADV.getBox()   == patch.getBox());

    TBOX_ASSERT(Q.getBox() == patch.getBox());

    TBOX_ASSERT(F.getBox() == patch.getBox());
#endif
    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();

    const Index<NDIM>& ilower = patch.getBox().lower();
    const Index<NDIM>& iupper = patch.getBox().upper();

    const IntVector<NDIM>& u_ADV_ghost_cells = u_ADV.getGhostCellWidth();
    const IntVector<NDIM>& Q_ghost_cells = Q.getGhostCellWidth();
    const IntVector<NDIM>& F_ghost_cells = F.getGhostCellWidth();
    const IntVector<NDIM>& q_half_ghost_cells = q_half.getGhostCellWidth();

    CellData<NDIM,double> dQ(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> Q_L(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> Q_R(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> Q_temp1(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> F_temp1(patch.getBox(),1,F_ghost_cells);
    FaceData<NDIM,double> q_half_temp(patch.getBox(),1,q_half_ghost_cells);
#if (NDIM > 2)
    CellData<NDIM,double> Q_temp2(patch.getBox(),1,Q_ghost_cells);
    CellData<NDIM,double> F_temp2(patch.getBox(),1,F_ghost_cells);
#endif

    for (int depth = 0; depth < Q.getDepth(); ++depth)
    {
#if (NDIM == 2)
        GODUNOV_PREDICT_WITH_SOURCE_FC(
            dx,dt,
            ilower(0),iupper(0),ilower(1),iupper(1),
            Q_ghost_cells(0),Q_ghost_cells(1),
            F_ghost_cells(0),F_ghost_cells(1),
            Q.getPointer(depth),Q_temp1.getPointer(0),
            dQ.getPointer(0),Q_L.getPointer(0),Q_R.getPointer(0),
            F.getPointer(depth),F_temp1.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),
            q_half_ghost_cells(0),q_half_ghost_cells(1),
            u_ADV.getPointer(0),u_ADV.getPointer(1),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth));
#endif
#if (NDIM == 3)
        GODUNOV_PREDICT_WITH_SOURCE_FC(
            dx,dt,
            static_cast<unsigned int>(d_using_full_ctu),
            ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
            Q_ghost_cells(0),Q_ghost_cells(1),Q_ghost_cells(2),
            F_ghost_cells(0),F_ghost_cells(1),F_ghost_cells(2),
            Q.getPointer(depth),Q_temp1.getPointer(0),Q_temp2.getPointer(0),
            dQ.getPointer(0),Q_L.getPointer(0),Q_R.getPointer(0),
            F.getPointer(depth),F_temp1.getPointer(0),F_temp2.getPointer(0),
            u_ADV_ghost_cells(0),u_ADV_ghost_cells(1),u_ADV_ghost_cells(2),
            q_half_ghost_cells(0),q_half_ghost_cells(1),q_half_ghost_cells(2),
            u_ADV.getPointer(0),u_ADV.getPointer(1),u_ADV.getPointer(2),
            q_half_temp.getPointer(0),q_half_temp.getPointer(1),q_half_temp.getPointer(2),
            q_half.getPointer(0,depth),q_half.getPointer(1,depth),q_half.getPointer(2,depth));
#endif
    }
    return;
}// predictWithSourceTerm

void
GodunovAdvector::getFromInput(
    Pointer<Database> db,
    bool /*is_from_restart*/)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
#if (NDIM == 3)
    if (db->keyExists("using_full_ctu")) d_using_full_ctu = db->getBool("using_full_ctu");
#else
    NULL_USE(db);
#endif
    return;
}// getFromInput

void
GodunovAdvector::getFromRestart()
{
    Pointer<Database> root_db =
        RestartManager::getManager()->getRootDatabase();

    Pointer<Database> db;

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
#if (NDIM == 3)
    d_using_full_ctu = db->getBool("d_using_full_ctu");
#endif
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
