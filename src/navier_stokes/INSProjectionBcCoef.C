// Filename: INSProjectionBcCoef.C
// Last modified: <22.Jul.2008 16:59:16 griffith@box230.cims.nyu.edu>
// Created on 22 Feb 2007 by Boyce Griffith (boyce@trasnaform2.local)

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
#define NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_F77 F77_FUNC_(navier_stokes_homogeneous_projection_bc_coefs2d,NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS2D)
#define NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77 F77_FUNC_(navier_stokes_fc_inhomogeneous_projection_bc_coefs2d,NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS2D)
#define NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77 F77_FUNC_(navier_stokes_sc_inhomogeneous_projection_bc_coefs2d,NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS2D)
#endif
#if (NDIM == 3)
#define NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_F77 F77_FUNC_(navier_stokes_homogeneous_projection_bc_coefs3d,NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS3D)
#define NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77 F77_FUNC_(navier_stokes_fc_inhomogeneous_projection_bc_coefs3d,NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS3D)
#define NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77 F77_FUNC_(navier_stokes_sc_inhomogeneous_projection_bc_coefs3d,NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS3D)
#endif

// Function interfaces
extern "C"
{
    void
    NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_F77(
        double* acoef, double* bcoef,
        const int& blower0, const int& bupper0,
        const int& blower1, const int& bupper1
#if (NDIM == 3)
        ,const int& blower2, const int& bupper2
#endif
                                                      );

    void
    NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77(
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
    NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77(
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
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const P_bc_coef,
    const std::string& projection_type,
    const int u_idx,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_P_idx(-1),
      d_P_bc_coef(NULL),
      d_projection_type(),
      d_u_idx(-1),
      d_u_bc_coefs(NDIM,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL)),
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
    const std::string& projection_type)
{
    if (projection_type != "pressure_increment" && projection_type != "pressure_update")
    {
        TBOX_ERROR("INSProjectionBcCoef::setProjectionType():\n"
                   << "  invalid velocity projection type: " << projection_type << "\n"
                   << "  valid choices are: ``pressure_increment'' or ``pressure_update''" << std::endl);
    }
    d_projection_type = projection_type;
    return;
}// setProjectionType

void
INSProjectionBcCoef::setPressurePhysicalBcCoef(
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const P_bc_coef)
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
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSProjectionBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
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
    NAVIER_STOKES_HOMOGENEOUS_PROJECTION_BC_COEFS_F77(
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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_fc_data = patch.getPatchData(d_u_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_sc_data = patch.getPatchData(d_u_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_data = patch.getPatchData(d_P_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_fc_data.isNull() || !u_sc_data.isNull());
    TBOX_ASSERT(u_fc_data.isNull() || u_fc_data->getGhostCellWidth().max() == u_fc_data->getGhostCellWidth().min());
    TBOX_ASSERT(u_sc_data.isNull() || u_sc_data->getGhostCellWidth().max() == u_sc_data->getGhostCellWidth().min());
    TBOX_ASSERT(!P_data.isNull());
    TBOX_ASSERT(P_data->getGhostCellWidth().max() == P_data->getGhostCellWidth().min());
    TBOX_ASSERT(!gcoef_data.isNull());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif

    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > acoef_data_P =
        new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > bcoef_data_P =
        new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > gcoef_data_P =
        new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);

    if (d_P_bc_coef != NULL)
    {
        d_P_bc_coef->setBcCoefs(
            acoef_data_P, bcoef_data_P, gcoef_data_P, variable, patch, bdry_box, fill_time);
    }

    if (!u_fc_data.isNull())
    {
        const int u_ghosts = (u_fc_data->getGhostCellWidth()).max();
        const int P_ghosts = (P_data->getGhostCellWidth()).max();
        const int using_pressure_increment = (d_projection_type == "pressure_increment" ? 1 : 0);
        NAVIER_STOKES_FC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77(
            u_fc_data->getPointer(0), u_fc_data->getPointer(1),
#if (NDIM == 3)
            u_fc_data->getPointer(2),
#endif
            u_ghosts,
            P_data->getPointer(), P_ghosts,
            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(), gcoef_data_P->getPointer(),
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
        const int using_pressure_increment = (d_projection_type == "pressure_increment" ? 1 : 0);
        NAVIER_STOKES_SC_INHOMOGENEOUS_PROJECTION_BC_COEFS_F77(
            u_sc_data->getPointer(0), u_sc_data->getPointer(1),
#if (NDIM == 3)
            u_sc_data->getPointer(2),
#endif
            u_ghosts,
            P_data->getPointer(), P_ghosts,
            acoef_data->getPointer(), bcoef_data->getPointer(), gcoef_data->getPointer(), gcoef_data_P->getPointer(),
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

SAMRAI::hier::IntVector<NDIM>
INSProjectionBcCoef::numberOfExtensionsFillable() const
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
