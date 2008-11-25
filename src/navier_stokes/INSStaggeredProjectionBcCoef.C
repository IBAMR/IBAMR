// Filename: INSStaggeredProjectionBcCoef.C
// Last modified: <24.Nov.2008 15:39:20 griffith@box230.cims.nyu.edu>
// Created on 23 Jul 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredProjectionBcCoef.h"

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
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredProjectionBcCoef::INSStaggeredProjectionBcCoef(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_u_bc_coefs(NDIM,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_target_idx(-1),
      d_homogeneous_bc(false)
{
    setVelocityPhysicalBcCoefs(u_bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredProjectionBcCoef

INSStaggeredProjectionBcCoef::~INSStaggeredProjectionBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredProjectionBcCoef

void
INSStaggeredProjectionBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSStaggeredProjectionBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

void
INSStaggeredProjectionBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSStaggeredProjectionBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredProjectionBcCoef::setBcCoefs(
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
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
#endif
    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;
//  const bool is_lower        = location_index%2 == 0;
    const SAMRAI::hier::Box<NDIM>& bc_coef_box = acoef_data->getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    if (!gcoef_data.isNull()) TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif

    // Set the unmodified velocity bc coefs.
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // Modify the velocity boundary conditions to correspond to pressure
    // boundary conditions.
    const bool set_acoef_vals = !acoef_data.isNull();
    const bool set_bcoef_vals = !bcoef_data.isNull();
    const bool set_gcoef_vals = !gcoef_data.isNull();
    for (SAMRAI::hier::Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& i = it();
        double dummy_val;
        double& alpha = set_acoef_vals ? (*acoef_data)(i,0) : dummy_val;
        double& beta  = set_bcoef_vals ? (*bcoef_data)(i,0) : dummy_val;
        double& gamma = set_gcoef_vals ? (*gcoef_data)(i,0) : dummy_val;

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
            // Set the boundary condition coefficients to correspond to
            // homogeneous Dirichlet boundary conditions on the pressure.
            alpha = 1.0;
            beta  = 0.0;
            gamma = 0.0;
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
INSStaggeredProjectionBcCoef::numberOfExtensionsFillable() const
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
        ret_val = SAMRAI::hier::IntVector<NDIM>::min(ret_val, d_u_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
