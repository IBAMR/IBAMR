// Filename: INSStaggeredVelocityBcCoef.C
// Last modified: <22.Jul.2008 19:30:13 griffith@box230.cims.nyu.edu>
// Created on 22 Jul 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredVelocityBcCoef.h"

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

INSStaggeredVelocityBcCoef::INSStaggeredVelocityBcCoef(
    const int comp_idx,
    const int P_idx,
    const double mu,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx),
      d_target_idx(-1),
      d_P_idx(-1),
      d_mu(mu),
      d_u_bc_coefs(NDIM,static_cast<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_homogeneous_bc(false)
{
    setPressurePatchDataIndex(P_idx);
    setVelocityPhysicalBcCoefs(u_bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredVelocityBcCoef

INSStaggeredVelocityBcCoef::~INSStaggeredVelocityBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredVelocityBcCoef

void
INSStaggeredVelocityBcCoef::setPressurePatchDataIndex(
    const int P_idx)
{
    d_P_idx = P_idx;
    return;
}// setPressurePatchDataIndex

void
INSStaggeredVelocityBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSStaggeredVelocityBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

void
INSStaggeredVelocityBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSStaggeredVelocityBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredVelocityBcCoef::setBcCoefs(
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
    // Set the unmodified velocity bc coefs.
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data beyond this point.
    if (gcoef_data.isNull()) return;

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc) gcoef_data->fillAll(0.0);

    // Modify Neumann boundary conditions to correspond to traction (stress)
    // boundary conditions.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_data =
        patch.checkAllocated(d_target_idx)
        ? patch.getPatchData(d_target_idx)
        : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_data.isNull());
    TBOX_ASSERT(u_data->getGhostCellWidth().max() == u_data->getGhostCellWidth().min());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_data =
        patch.checkAllocated(d_P_idx)
        ? patch.getPatchData(d_P_idx)
        : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!P_data.isNull());
    TBOX_ASSERT(P_data->getGhostCellWidth().max() == P_data->getGhostCellWidth().min());
    TBOX_ASSERT(!acoef_data.isNull());
    TBOX_ASSERT(!bcoef_data.isNull());
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
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (SAMRAI::hier::Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const SAMRAI::hier::Index<NDIM>& i = it();
        const double& alpha = (*acoef_data)(i,0);
        const double& beta  = (*bcoef_data)(i,0);
        double& g = (*gcoef_data)(i,0);
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(SAMRAI::tbox::MathUtilities<double>::equalEps(alpha+beta,1.0));
        TBOX_ASSERT(SAMRAI::tbox::MathUtilities<double>::equalEps(alpha,1.0) || SAMRAI::tbox::MathUtilities<double>::equalEps(beta,1.0));
#endif
        const bool traction_bc = SAMRAI::tbox::MathUtilities<double>::equalEps(beta ,1.0);
        if (traction_bc)
        {
            if (d_comp_idx == bdry_normal_axis)
            {
                // Compute the pressure at the boundary using normal
                // extrapolation.
                SAMRAI::hier::Index<NDIM> i_bdry = i;
                SAMRAI::hier::Index<NDIM> i_intr = i;
                if (is_lower)
                {
                    i_bdry(bdry_normal_axis) += 0;
                    i_intr(bdry_normal_axis) += 1;
                }
                else
                {
                    i_bdry(bdry_normal_axis) -= 1;
                    i_intr(bdry_normal_axis) -= 2;
                }
                const double P = 1.5*(*P_data)(i_bdry) - 0.5*(*P_data)(i_intr);

                // Correct the boundary condition value.
                g = 0.5*(g+P)/d_mu;
            }
            else
            {
                // Compute the tangential derivative of the normal component of
                // the velocity at the boundary.
                SAMRAI::hier::Index<NDIM> i_upper = i;
                SAMRAI::hier::Index<NDIM> i_lower = i;
                i_upper(d_comp_idx) += 1;
                i_lower(d_comp_idx) -= 0;
                const SAMRAI::pdat::SideIndex<NDIM> i_s_upper(i_upper, bdry_normal_axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                const SAMRAI::pdat::SideIndex<NDIM> i_s_lower(i_lower, bdry_normal_axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);
                const double du_norm_dtan = ((*u_data)(i_s_upper)-(*u_data)(i_s_lower))/dx[d_comp_idx];

                // Correct the boundary condition value.
                g = g/d_mu-du_norm_dtan;
            }
        }
    }
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
INSStaggeredVelocityBcCoef::numberOfExtensionsFillable() const
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
