// Filename: INSProjectionBcCoef.C
// Created on 22 Feb 2007 by Boyce Griffith
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

#include "INSProjectionBcCoef.h"

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

// IBTK INCLUDES
#include <ibtk/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CellData.h>
#include <FaceData.h>
#include <SideData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_FC FC_FUNC_(navier_stokes_homogeneous_projection_bc_coefs2d,NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS2D)
#define NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC FC_FUNC_(navier_stokes_fc_inhomogeneous_projection_bc_coefs2d,NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS2D)
#define NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC FC_FUNC_(navier_stokes_sc_inhomogeneous_projection_bc_coefs2d,NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS2D)
#endif
#if (NDIM == 3)
#define NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_FC FC_FUNC_(navier_stokes_homogeneous_projection_bc_coefs3d,NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS3D)
#define NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC FC_FUNC_(navier_stokes_fc_inhomogeneous_projection_bc_coefs3d,NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS3D)
#define NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC FC_FUNC_(navier_stokes_sc_inhomogeneous_projection_bc_coefs3d,NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS3D)
#endif

// Function interfaces
extern "C"
{
    void
    NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_FC(
        double* acoef, double* bcoef,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1
#if (NDIM == 3)
        ,const int& blower2, const int& bupper2
#endif
                                                     );

    void
    NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC(
        const double* u0, const double* u1,
#if (NDIM == 3)
        const double* u2,
#endif
        const int& u_gcw,
        const double* P, const int& P_gcw,
        const double* acoef, const double* bcoef, double* gcoef, const double* P_bdry,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
#if (NDIM == 3)
        const int& blower2, const int& bupper2,
#endif
        const int& location_index,
        const int& using_pressure_increment,
        const double& rho,
        const double& dt
                                                          );

    void
    NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC(
        const double* u0, const double* u1,
#if (NDIM == 3)
        const double* u2,
#endif
        const int& u_gcw,
        const double* P, const int& P_gcw,
        const double* acoef, const double* bcoef, double* gcoef, const double* P_bdry,
        const int& ilower0, const int& iupper0,
        const int& ilower1, const int& iupper1,
#if (NDIM == 3)
        const int& ilower2, const int& iupper2,
#endif
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1,
#if (NDIM == 3)
        const int& blower2, const int& bupper2,
#endif
        const int& location_index,
        const int& using_pressure_increment,
        const double& rho,
        const double& dt
                                                          );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSProjectionBcCoef::INSProjectionBcCoef(
    const int P_idx,
    RobinBcCoefStrategy<NDIM>* const P_bc_coef,
    const ProjectionMethodType& projection_type,
    const int u_idx,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_P_idx(-1),
      d_P_bc_coef(NULL),
      d_projection_type(),
      d_u_idx(-1),
      d_u_bc_coefs(static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_homogeneous_bc(false),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN())
{
    setCurrentPressurePatchDataIndex(P_idx);
    setPressurePhysicalBcCoef(P_bc_coef);
    setProjectionType(projection_type);
    setIntermediateVelocityPatchDataIndex(u_idx);
    setVelocityPhysicalBcCoefs(u_bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSProjectionBcCoef

INSProjectionBcCoef::~INSProjectionBcCoef()
{
    // intentionally blank
    return;
}// ~INSProjectionBcCoef

void
INSProjectionBcCoef::setProblemCoefs(
    const double rho,
    const double dt)
{
    d_rho = rho;
    d_dt = dt;
    return;
}// setProblemCoefs

void
INSProjectionBcCoef::setCurrentPressurePatchDataIndex(
    const int P_idx)
{
    d_P_idx = P_idx;
    return;
}// setCurrentPressurePatchDataIndex

void
INSProjectionBcCoef::setProjectionType(
    const ProjectionMethodType& projection_type)
{
    d_projection_type = projection_type;
    return;
}// setProjectionType

void
INSProjectionBcCoef::setPressurePhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const P_bc_coef)
{
    d_P_bc_coef = P_bc_coef;
    return;
}// setPressurePhysicalBcCoef

void
INSProjectionBcCoef::setIntermediateVelocityPatchDataIndex(
    const int u_idx)
{
    d_u_idx = u_idx;
    return;
}// setIntermediateVelocityPatchDataIndex

void
INSProjectionBcCoef::setVelocityPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& u_bc_coefs)
{
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

void
INSProjectionBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    // intentionally blank
    return;
}// setTargetPatchDataIndex

void
INSProjectionBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSProjectionBcCoef::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_u_bc_coefs[d] != NULL);
    }
#endif
    const unsigned int location_index   = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index/2;
//  const bool is_lower        = location_index%2 == 0;
    const Box<NDIM>& patch_box = patch.getBox();
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();

    // Set the normal velocity bc coefs.
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // Set the corresponding projection Poisson problem homogeneous Robin
    // coefficients.
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
#endif
    NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_FC(
        acoef_data->getPointer(), bcoef_data->getPointer(),
        bc_coef_box.lower(0), bc_coef_box.upper(0),
        bc_coef_box.lower(1), bc_coef_box.upper(1)
#if (NDIM == 3)
        ,bc_coef_box.lower(2), bc_coef_box.upper(2)
#endif
                                                     );

    if (d_homogeneous_bc && !gcoef_data.isNull()) gcoef_data->fillAll(0.0);

    // Do not further modify the boundary condition coefficients unless we are
    // setting inhomogeneous boundary conditions.
    if (d_homogeneous_bc) return;

    // Loop over the boundary box and reset the inhomogeneous coefficients.
    Pointer<FaceData<NDIM,double> > u_fc_data = patch.getPatchData(d_u_idx);
    Pointer<SideData<NDIM,double> > u_sc_data = patch.getPatchData(d_u_idx);
    Pointer<CellData<NDIM,double> > P_data = patch.getPatchData(d_P_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_fc_data.isNull() || !u_sc_data.isNull());
    TBOX_ASSERT(u_fc_data.isNull() || u_fc_data->getGhostCellWidth().max() == u_fc_data->getGhostCellWidth().min());
    TBOX_ASSERT(u_sc_data.isNull() || u_sc_data->getGhostCellWidth().max() == u_sc_data->getGhostCellWidth().min());
    TBOX_ASSERT(!P_data.isNull());
    TBOX_ASSERT(P_data->getGhostCellWidth().max() == P_data->getGhostCellWidth().min());
    TBOX_ASSERT(!gcoef_data.isNull());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif

    ArrayData<NDIM,double> acoef_data_P(bc_coef_box, 1);
    ArrayData<NDIM,double> bcoef_data_P(bc_coef_box, 1);
    ArrayData<NDIM,double> gcoef_data_P(bc_coef_box, 1);

    Pointer<ArrayData<NDIM,double> > acoef_data_P_ptr(&acoef_data_P, false);
    Pointer<ArrayData<NDIM,double> > bcoef_data_P_ptr(&bcoef_data_P, false);
    Pointer<ArrayData<NDIM,double> > gcoef_data_P_ptr(&gcoef_data_P, false);

    if (d_P_bc_coef != NULL)
    {
        d_P_bc_coef->setBcCoefs(
            acoef_data_P_ptr, bcoef_data_P_ptr, gcoef_data_P_ptr, variable, patch, bdry_box, fill_time);
    }

    if (!u_fc_data.isNull())
    {
        const int u_ghosts = (u_fc_data->getGhostCellWidth()).max();
        const int P_ghosts = (P_data->getGhostCellWidth()).max();
        const int using_pressure_increment = (d_projection_type == PRESSURE_INCREMENT ? 1 : 0);
        NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC(
            u_fc_data->getPointer(0), u_fc_data->getPointer(1),
#if (NDIM == 3)
            u_fc_data->getPointer(2),
#endif
            u_ghosts,
            P_data->getPointer(), P_ghosts,
            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(), gcoef_data_P.getPointer(),
            patch_box.lower(0), patch_box.upper(0),
            patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
            patch_box.lower(2), patch_box.upper(2),
#endif
            bc_coef_box.lower(0), bc_coef_box.upper(0),
            bc_coef_box.lower(1), bc_coef_box.upper(1),
#if (NDIM == 3)
            bc_coef_box.lower(2), bc_coef_box.upper(2),
#endif
            location_index, using_pressure_increment,
            d_rho, d_dt);
    }
    else if (!u_sc_data.isNull())
    {
        const int u_ghosts = (u_sc_data->getGhostCellWidth()).max();
        const int P_ghosts = (P_data->getGhostCellWidth()).max();
        const int using_pressure_increment = (d_projection_type == PRESSURE_INCREMENT ? 1 : 0);
        NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_FC(
            u_sc_data->getPointer(0), u_sc_data->getPointer(1),
#if (NDIM == 3)
            u_sc_data->getPointer(2),
#endif
            u_ghosts,
            P_data->getPointer(), P_ghosts,
            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(), gcoef_data_P.getPointer(),
            patch_box.lower(0), patch_box.upper(0),
            patch_box.lower(1), patch_box.upper(1),
#if (NDIM == 3)
            patch_box.lower(2), patch_box.upper(2),
#endif
            bc_coef_box.lower(0), bc_coef_box.upper(0),
            bc_coef_box.lower(1), bc_coef_box.upper(1),
#if (NDIM == 3)
            bc_coef_box.lower(2), bc_coef_box.upper(2),
#endif
            location_index, using_pressure_increment,
            d_rho, d_dt);
    }
    else
    {
        TBOX_ERROR("this statement should not be reached!\n");
    }
    return;
}// setBcCoefs

IntVector<NDIM>
INSProjectionBcCoef::numberOfExtensionsFillable() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_u_bc_coefs[d] != NULL);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_u_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
