// Filename: IntermediateVelocityRobinBcCoefs.C
// Last modified: <09.Mar.2007 22:43:47 griffith@box221.cims.nyu.edu>
// Created on 30 Sep 2006 by Boyce Griffith (boyce@boyce-griffiths-powerbook-g4-15.local)

#include "IntermediateVelocityRobinBcCoefs.h"

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
#include <ArrayDataBasicOps.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IntermediateVelocityRobinBcCoefs::IntermediateVelocityRobinBcCoefs(
    int velocity_depth,
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef)
    : d_velocity_depth(velocity_depth),
      d_bc_coef(bc_coef),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_P_idx(-1),
      d_Phi_idx(-1)
{
    // intentionally blank
    return;
}// IntermediateVelocityRobinBcCoefs

IntermediateVelocityRobinBcCoefs::~IntermediateVelocityRobinBcCoefs()
{
    // intentionally blank
    return;
}// ~IntermediateVelocityRobinBcCoefs

void
IntermediateVelocityRobinBcCoefs::setCurrentTime(
    const double current_time)
{
    d_current_time = current_time;
    return;
}// setCurrentTime

void
IntermediateVelocityRobinBcCoefs::setNewTime(
    const double new_time)
{
    d_new_time = new_time;
    return;
}// setNewTime

void
IntermediateVelocityRobinBcCoefs::setRho(
    const double rho)
{
    d_rho = rho;
    return;
}// setRho

void
IntermediateVelocityRobinBcCoefs::setPressureIndex(
    const int P_idx)
{
    d_P_idx = P_idx;
    return;
}// setPressureIndex

void
IntermediateVelocityRobinBcCoefs::setPhiIndex(
    const int Phi_idx)
{
    d_Phi_idx = Phi_idx;
    return;
}// setPhiIndex

void
IntermediateVelocityRobinBcCoefs::setBcCoefs(
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
     const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
     const SAMRAI::hier::Patch<NDIM>& patch,
     const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
     double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    TBOX_ERROR("IntermediateVelocityRobinBcCoefs::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#else
    d_bc_coef->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif
    if (d_P_idx != -1 && d_Phi_idx != -1 &&
        !SAMRAI::tbox::Utilities::deq(fill_time, d_current_time))
    {
        correctBcCoefs(acoef_data, bcoef_data, gcoef_data, patch, bdry_box);
    }
    return;
}// setBcCoefs

void
IntermediateVelocityRobinBcCoefs::setBcCoefs(
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
     const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
     const SAMRAI::hier::Patch<NDIM>& patch,
     const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
     double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    d_bc_coef->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    TBOX_ERROR("IntermediateVelocityRobinBcCoefs::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#endif
    const SAMRAI::hier::Box<NDIM>& bc_coef_box = acoef_data->getBox();
    SAMRAI::math::ArrayDataBasicOps<NDIM,double> array_ops;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > bcoef_data =
        new SAMRAI::pdat::ArrayData<NDIM,double>(bc_coef_box, 1);
    array_ops.scale(*bcoef_data, -1.0, *acoef_data, bc_coef_box);
    array_ops.addScalar(*bcoef_data, *bcoef_data, 1.0, bc_coef_box);

    if (d_P_idx != -1 && d_Phi_idx != -1 &&
        !SAMRAI::tbox::Utilities::deq(fill_time, d_current_time))
    {
        correctBcCoefs(acoef_data, bcoef_data, gcoef_data, patch, bdry_box);
    }
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
IntermediateVelocityRobinBcCoefs::numberOfExtensionsFillable() const
{
    return d_bc_coef->numberOfExtensionsFillable();
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IntermediateVelocityRobinBcCoefs::correctBcCoefs(
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
     const SAMRAI::hier::Patch<NDIM>& patch,
     const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box) const
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_data =
        patch.getPatchData(d_P_idx);
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Phi_data =
        patch.getPatchData(d_Phi_idx);

    // IMPORTANT NOTE: This is a clumsy and kludgey mechanism to
    // detect when homogeneous boundary conditions are being employed.
    if (P_data.isNull()) return;

    // NOTE: At this point, presumably inhomogeneous boundary
    // conditions are being employed.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double dt = d_new_time - d_current_time;

    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2;
    const bool bdry_upper_side = (location_index % 2) != 0;

    const SAMRAI::hier::BoundaryBox<NDIM> trimmed_bdry_box =
        STOOLS::PhysicalBoundaryUtilities::trimBoundaryCodim1Box(bdry_box, patch);
    const SAMRAI::hier::Box<NDIM> bc_coef_box =
        STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);

    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i_s_bdry = b();
        const double& a = (*acoef_data)(i_s_bdry,0);
        const double& b = (*bcoef_data)(i_s_bdry,0);
        double& g = (*gcoef_data)(i_s_bdry,0);
        assert(b == 0.0);

        // i_s_bdry: side index located on physical boundary
        //
        // i_c_intr: cell index located adjacent to physical boundary
        // in the patch interior
        //
        // i_c_bdry: cell index located adjacent to physical boundary
        // in the patch exterior
        SAMRAI::pdat::CellIndex<NDIM> i_c_intr = i_s_bdry;
        SAMRAI::pdat::CellIndex<NDIM> i_c_bdry = i_c_intr;
        if (bdry_upper_side)
        {
            i_c_intr(bdry_normal_axis) -= 1;
            i_c_bdry(bdry_normal_axis) = i_c_intr(bdry_normal_axis)+1;
        }
        else
        {
            i_c_bdry(bdry_normal_axis) = i_c_intr(bdry_normal_axis)-1;
        }

        if (d_velocity_depth == bdry_normal_axis)
        {
            // Correct boundary conditions for the normal velocity.
            const double dP_dn   = ((*  P_data)(i_c_bdry)-(*  P_data)(i_c_intr))/dx[bdry_normal_axis];
            const double dPhi_dn = ((*Phi_data)(i_c_bdry)-(*Phi_data)(i_c_intr))/dx[bdry_normal_axis];

            if (bdry_upper_side)
            {
                g -= a*(dt/d_rho)*dP_dn;
            }
            else
            {
                g += a*(dt/d_rho)*dP_dn;
            }
        }
        else
        {
            // Correct boundary conditions for the tangential velocity.
            double dP_dt   = 0.0;
            double dPhi_dt = 0.0;
            i_c_intr(d_velocity_depth) += 1;
            i_c_bdry(d_velocity_depth) += 1;
            dP_dt   += 0.5*((*  P_data)(i_c_bdry)+(*  P_data)(i_c_intr))/dx[d_velocity_depth];
            dPhi_dt += 0.5*((*Phi_data)(i_c_bdry)+(*Phi_data)(i_c_intr))/dx[d_velocity_depth];
            i_c_intr(d_velocity_depth) -= 2;
            i_c_bdry(d_velocity_depth) -= 2;
            dP_dt   -= 0.5*((*  P_data)(i_c_bdry)+(*  P_data)(i_c_intr))/dx[d_velocity_depth];
            dPhi_dt -= 0.5*((*Phi_data)(i_c_bdry)+(*Phi_data)(i_c_intr))/dx[d_velocity_depth];
            i_c_intr(d_velocity_depth) += 1;
            i_c_bdry(d_velocity_depth) += 1;

            g -= a*(dt/d_rho)*dP_dt;
        }
    }

    return;
}// correctBcCoefs

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
