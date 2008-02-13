// Filename: INSIntermediateVelocityBcCoef.C
// Last modified: <12.Feb.2008 21:22:13 griffith@box221.cims.nyu.edu>
// Created on 30 Aug 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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

// STOOLS INCLUDES
#include <stools/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define NAVIER_STOKES_OPEN_BC_COEFS_F77 F77_FUNC_(navier_stokes_open_bc_coefs2d,NAVIER_STOKES_OPEN_BC_COEFS2D)
#define NAVIER_STOKES_TANGENTIAL_BC_COEFS_F77 F77_FUNC_(navier_stokes_tangential_bc_coefs2d,NAVIER_STOKES_TANGENTIAL_BC_COEFS2D)
#endif
#if (NDIM == 3)
#define NAVIER_STOKES_OPEN_BC_COEFS_F77 F77_FUNC_(navier_stokes_open_bc_coefs3d,NAVIER_STOKES_OPEN_BC_COEFS3D)
#define NAVIER_STOKES_TANGENTIAL_BC_COEFS_F77 F77_FUNC_(navier_stokes_tangential_bc_coefs3d,NAVIER_STOKES_TANGENTIAL_BC_COEFS3D)
#endif

// Function interfaces
extern "C"
{
    void
    NAVIER_STOKES_OPEN_BC_COEFS_F77(
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
    NAVIER_STOKES_TANGENTIAL_BC_COEFS_F77(
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
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx),
      d_target_idx(-1),
      d_Phi_idx(-1),
      d_u_bc_coefs(NDIM,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL)),
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
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSIntermediateVelocityBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    TBOX_ERROR("INSIntermediateVelocityBcCoef::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << std::endl);
#else
    setBcCoefs_private(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif
    return;
}// setBcCoefs

void
INSIntermediateVelocityBcCoef::setBcCoefs(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > bcoef_data =
        (acoef_data.isNull()
         ? NULL
         : new SAMRAI::pdat::ArrayData<NDIM,double>(acoef_data->getBox(), acoef_data->getDepth()));
    setBcCoefs_private(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    TBOX_ERROR("INSIntermediateVelocityBcCoef::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << std::endl);
#endif
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
INSIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_bc_coefs.size() == NDIM);
    for (unsigned l = 0; l < d_u_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(d_u_bc_coefs[l] != NULL);
    }
#endif
    SAMRAI::hier::IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (int d = 0; d < NDIM; ++d)
    {
        ret_val = SAMRAI::hier::IntVector<NDIM>::min(
            ret_val, d_u_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSIntermediateVelocityBcCoef::setBcCoefs_private(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_bc_coefs.size() == NDIM);
    for (unsigned l = 0; l < d_u_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(d_u_bc_coefs[l] != NULL);
    }
#endif
    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;
//  const bool is_lower        = location_index%2 == 0;
    const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Box<NDIM>& bc_coef_box = acoef_data->getBox();

    // Set the "true" velocity bc coefs.
#if USING_OLD_ROBIN_BC_INTERFACE
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data beyond this point.
    if (gcoef_data.isNull()) return;

    // Patch box information.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_data =
        patch.checkAllocated(d_target_idx)
        ? patch.getPatchData(d_target_idx)
        : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
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
    NAVIER_STOKES_OPEN_BC_COEFS_F77(
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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Phi_data =
        patch.checkAllocated(d_Phi_idx)
        ? patch.getPatchData(d_Phi_idx)
        : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
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
    NAVIER_STOKES_TANGENTIAL_BC_COEFS_F77(
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
}// setBcCoefs_private

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
