// Filename: INSStaggeredPressureBcCoef.C
// Last modified: <03.Oct.2008 17:49:53 griffith@box230.cims.nyu.edu>
// Created on 23 Jul 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredPressureBcCoef.h"

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
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <SideData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredPressureBcCoef::INSStaggeredPressureBcCoef(
    const INSCoefs& problem_coefs,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_problem_coefs(problem_coefs),
      d_u_current_idx(-1),
      d_u_new_idx(-1),
      d_u_bc_coefs(NDIM,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_target_idx(-1),
      d_homogeneous_bc(false)
{
    setVelocityPhysicalBcCoefs(u_bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredPressureBcCoef

INSStaggeredPressureBcCoef::~INSStaggeredPressureBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredPressureBcCoef

void
INSStaggeredPressureBcCoef::setVelocityCurrentPatchDataIndex(
    const int u_current_idx)
{
    d_u_current_idx = u_current_idx;
    return;
}// setVelocityCurrentPatchDataIndex

void
INSStaggeredPressureBcCoef::setVelocityNewPatchDataIndex(
    const int u_new_idx)
{
    d_u_new_idx = u_new_idx;
    return;
}// setVelocityNewPatchDataIndex

void
INSStaggeredPressureBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSStaggeredPressureBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

void
INSStaggeredPressureBcCoef::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

void
INSStaggeredPressureBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSStaggeredPressureBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredPressureBcCoef::setBcCoefs(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    const double half_time = 0.5*(d_current_time+d_new_time);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_bc_coefs.size() == NDIM);
    for (unsigned l = 0; l < d_u_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(d_u_bc_coefs[l] != NULL);
    }
    TBOX_ASSERT(SAMRAI::tbox::MathUtilities<double>::equalEps(fill_time,d_new_time) ||
                SAMRAI::tbox::MathUtilities<double>::equalEps(fill_time,half_time));
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
    TBOX_ASSERT(!gcoef_data.isNull());
#endif
    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;
    const bool is_lower        = location_index%2 == 0;
    const SAMRAI::hier::Box<NDIM>& bc_coef_box = acoef_data->getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif

    // Set the unmodified velocity bc coefs.
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, half_time);

    // Modify the velocity boundary conditions to correspond to pressure
    // boundary conditions.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_current_data =
        patch.checkAllocated(d_u_current_idx)
        ? patch.getPatchData(d_u_current_idx)
        : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_current_data.isNull());
    TBOX_ASSERT(u_current_data->getGhostCellWidth().max() == u_current_data->getGhostCellWidth().min());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_new_data =
        patch.checkAllocated(d_u_new_idx)
        ? patch.getPatchData(d_u_new_idx)
        : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_new_data.isNull());
    TBOX_ASSERT(u_new_data->getGhostCellWidth().max() == u_new_data->getGhostCellWidth().min());
#endif
    const SAMRAI::hier::Box<NDIM> ghost_box = u_current_data->getGhostBox() * u_new_data->getGhostBox();
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double mu = d_problem_coefs.getMu();
    for (SAMRAI::hier::Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& i = it();
        double& alpha = (*acoef_data)(i,0);
        double& beta  = (*bcoef_data)(i,0);
        double& gamma = (*gcoef_data)(i,0);

        const bool velocity_bc = SAMRAI::tbox::MathUtilities<double>::equalEps(alpha,1.0);
        const bool traction_bc = SAMRAI::tbox::MathUtilities<double>::equalEps(beta ,1.0);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (velocity_bc)
        {
            // Set the boundary condition coefficients to correspond to
            // homogeneous Neumann boundary conditions on the pressure.
            alpha = 0.0;
            beta  = 1.0;
            gamma = 0.0;
        }
        else if (traction_bc)
        {
            // Compute (d/dn)u_norm at the boundary by extrapolating the
            // divergence free condition to the boundary.
            SAMRAI::hier::Index<NDIM> i_intr0 = i;
            SAMRAI::hier::Index<NDIM> i_intr1 = i;

            if (is_lower)
            {
                i_intr0(bdry_normal_axis) += 0;
                i_intr1(bdry_normal_axis) += 1;
            }
            else
            {
                i_intr0(bdry_normal_axis) -= 1;
                i_intr1(bdry_normal_axis) -= 2;
            }

            for (int d = 0; d < NDIM; ++d)
            {
                if (d != bdry_normal_axis)
                {
                    i_intr0(d) = std::max(i_intr0(d),ghost_box.lower()(d));
                    i_intr0(d) = std::min(i_intr0(d),ghost_box.upper()(d));

                    i_intr1(d) = std::max(i_intr1(d),ghost_box.lower()(d));
                    i_intr1(d) = std::min(i_intr1(d),ghost_box.upper()(d));
                }
            }

            double du_norm_current_dn = 0.0;
            double du_norm_new_dn     = 0.0;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                if (axis != bdry_normal_axis)
                {
                    const SAMRAI::pdat::SideIndex<NDIM> i_s_intr0_upper(i_intr0, axis, SAMRAI::pdat::SideIndex<NDIM>::Upper);
                    const SAMRAI::pdat::SideIndex<NDIM> i_s_intr1_upper(i_intr1, axis, SAMRAI::pdat::SideIndex<NDIM>::Upper);
                    const double u_tan_current_upper = 1.5*(*u_current_data)(i_s_intr0_upper)-0.5*(*u_current_data)(i_s_intr1_upper);
                    const double u_tan_new_upper     = 1.5*(*u_new_data    )(i_s_intr0_upper)-0.5*(*u_new_data    )(i_s_intr1_upper);

                    const SAMRAI::pdat::SideIndex<NDIM> i_s_intr0_lower(i_intr0, axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                    const SAMRAI::pdat::SideIndex<NDIM> i_s_intr1_lower(i_intr1, axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                    const double u_tan_current_lower = 1.5*(*u_current_data)(i_s_intr0_lower)-0.5*(*u_current_data)(i_s_intr1_lower);
                    const double u_tan_new_lower     = 1.5*(*u_new_data    )(i_s_intr0_lower)-0.5*(*u_new_data    )(i_s_intr1_lower);

                    du_norm_current_dn += (is_lower ? +1.0 : -1.0)*(u_tan_current_upper-u_tan_current_lower)/dx[axis];
                    du_norm_new_dn     += (is_lower ? +1.0 : -1.0)*(u_tan_new_upper    -u_tan_new_lower    )/dx[axis];
                }
            }

            // Set the boundary condition coefficients to correspond to either
            // homogeneous or inhomogeneous Dirichlet boundary conditions on the
            // pressure.
            alpha = 1.0;
            beta  = 0.0;
// XXXX     gamma = (is_lower ? -1.0 : +1.0)*(mu*du_norm_new_dn + (d_homogeneous_bc ? 0.0 : mu*du_norm_current_dn - gamma));
            gamma = (is_lower ? -1.0 : +1.0)*((d_homogeneous_bc ? 0.0 : - gamma));
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
INSStaggeredPressureBcCoef::numberOfExtensionsFillable() const
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
