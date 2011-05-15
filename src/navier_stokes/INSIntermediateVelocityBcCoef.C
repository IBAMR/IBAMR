// Filename: INSIntermediateVelocityBcCoef.C
// Created on 30 Aug 2007 by Boyce Griffith
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

#include "INSIntermediateVelocityBcCoef.h"

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
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_OPEN_BC_COEFS_FC FC_FUNC_(navier_stokes_open_bc_coefs2d,NAVIER_STOKES_OPEN_BC_COEFS2D)
#define NAVIER_STOKES_TANGENTIAL_BC_COEFS_FC FC_FUNC_(navier_stokes_tangential_bc_coefs2d,NAVIER_STOKES_TANGENTIAL_BC_COEFS2D)
#endif
#if (NDIM == 3)
#define NAVIER_STOKES_OPEN_BC_COEFS_FC FC_FUNC_(navier_stokes_open_bc_coefs3d,NAVIER_STOKES_OPEN_BC_COEFS3D)
#define NAVIER_STOKES_TANGENTIAL_BC_COEFS_FC FC_FUNC_(navier_stokes_tangential_bc_coefs3d,NAVIER_STOKES_TANGENTIAL_BC_COEFS3D)
#endif

// Function interfaces
extern "C"
{
    void
    NAVIER_STOKES_OPEN_BC_COEFS_FC(
        const double* U, const int& U_gcw,
        const double* acoef, const double* bcoef, double* gcoef,
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
        const int& comp_idx,
        const double* dx);

    void
    NAVIER_STOKES_TANGENTIAL_BC_COEFS_FC(
        const double* Phi, const int& Phi_gcw,
        const double* acoef, const double* bcoef, double* gcoef,
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
        const int& comp_idx,
        const double& rho, const double* dx, const double& dt);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSIntermediateVelocityBcCoef::INSIntermediateVelocityBcCoef(
    const int comp_idx,
    const int Phi_idx,
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx),
      d_target_idx(-1),
      d_Phi_idx(-1),
      d_u_bc_coefs(static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_homogeneous_bc(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_using_intermediate_velocity_bc_coefs(false),
      d_velocity_correction(false)
{
    useTrueVelocityBcCoefs();
    setPhiPatchDataIndex(Phi_idx);
    setVelocityPhysicalBcCoefs(u_bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSIntermediateVelocityBcCoef

INSIntermediateVelocityBcCoef::~INSIntermediateVelocityBcCoef()
{
    // intentionally blank
    return;
}// ~INSIntermediateVelocityBcCoef

void
INSIntermediateVelocityBcCoef::useTrueVelocityBcCoefs()
{
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_rho = std::numeric_limits<double>::quiet_NaN();
    d_using_intermediate_velocity_bc_coefs = false;
    d_velocity_correction = false;
    return;
}// useTrueVelocityBcCoef

void
INSIntermediateVelocityBcCoef::useIntermediateVelocityBcCoefs(
    const double current_time,
    const double new_time,
    const double rho,
    const bool velocity_correction)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_rho = rho;
    d_using_intermediate_velocity_bc_coefs = true;
    d_velocity_correction = velocity_correction;
    return;
}// useIntermediateVelocityBcCoef

void
INSIntermediateVelocityBcCoef::setPhiPatchDataIndex(
    const int Phi_idx)
{
    d_Phi_idx = Phi_idx;
    return;
}// setPhiPatchDataIndex

void
INSIntermediateVelocityBcCoef::setVelocityPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& u_bc_coefs)
{
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

void
INSIntermediateVelocityBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSIntermediateVelocityBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSIntermediateVelocityBcCoef::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_u_bc_coefs[d] != NULL);
    }
#endif
    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;
//  const bool is_lower        = location_index%2 == 0;
    const Box<NDIM>& patch_box = patch.getBox();
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();

    // Set the "true" velocity bc coefs.
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data beyond this point.
    if (gcoef_data.isNull()) return;

    // Patch box information.
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Enforce homogeneous boundary conditions.
    if (d_homogeneous_bc || d_velocity_correction) gcoef_data->fillAll(0.0);

    // At "open" boundaries, modify the normal velocity boundary conditions to
    // enforce div u = 0, and modify the tangential velocity boundary conditions
    // to enforce zero stress.  This is done by specifying a normal flux F at
    // the boundary.
    //
    // Note that this flux F may be non-zero even in the case that we are
    // employing homogeneous boundary conditions.
    Pointer<CellData<NDIM,double> > U_data =
        patch.checkAllocated(d_target_idx)
        ? patch.getPatchData(d_target_idx)
        : Pointer<PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!U_data.isNull());
    TBOX_ASSERT(U_data->getGhostCellWidth().max() == U_data->getGhostCellWidth().min());
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    const int U_ghosts = (U_data->getGhostCellWidth()).max();
    NAVIER_STOKES_OPEN_BC_COEFS_FC(
        U_data->getPointer(), U_ghosts,
        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
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
        location_index, d_comp_idx,
        dx);

    // Do not further modify the boundary condition coefficients unless we are
    // setting inhomogeneous boundary conditions for the tangential components
    // of U^{*}.
    if ((d_comp_idx == bdry_normal_axis) || d_homogeneous_bc || !d_using_intermediate_velocity_bc_coefs || !(fill_time > d_current_time))
    {
        return;
    }

    // Modify the inhomogeneous coefficients for the tangential components of
    // the intermediate velocity.
    Pointer<CellData<NDIM,double> > Phi_data =
        patch.checkAllocated(d_Phi_idx)
        ? patch.getPatchData(d_Phi_idx)
        : Pointer<PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Phi_data.isNull());
    TBOX_ASSERT(Phi_data->getGhostCellWidth().max() == Phi_data->getGhostCellWidth().min());
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    const int Phi_ghosts = (Phi_data->getGhostCellWidth()).max();
    const double dt = d_new_time - d_current_time;
    NAVIER_STOKES_TANGENTIAL_BC_COEFS_FC(
        Phi_data->getPointer(), Phi_ghosts,
        acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(),
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
        location_index, d_comp_idx,
        d_rho, dx, dt);
    return;
}// setBcCoefs

IntVector<NDIM>
INSIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_u_bc_coefs[d] != NULL);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (int d = 0; d < NDIM; ++d)
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
