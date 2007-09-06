// Filename: INSIntermediateVelocityBcCoef.C
// Last modified: <06.Sep.2007 01:53:01 griffith@box221.cims.nyu.edu>
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
#include <cassert>
#include <limits>

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
      d_homo_patch_data_idxs(),
      d_inhomo_patch_data_idxs(),
      d_u_bc_coefs(NDIM,NULL),
      d_homogeneous_bc(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_using_intermediate_velocity_bc_coefs(false)
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
    return;
}// useTrueVelocityBcCoef

void
INSIntermediateVelocityBcCoef::useIntermediateVelocityBcCoefs(
    const double current_time,
    const double new_time,
    const double rho)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_rho = rho;
    d_using_intermediate_velocity_bc_coefs = true;
    return;
}// useIntermediateVelocityBcCoef

void
INSIntermediateVelocityBcCoef::setPhiPatchDataIndex(
    const int Phi_idx)
{
    d_Phi_idx = Phi_idx;
    d_inhomo_patch_data_idxs.clear();
    d_inhomo_patch_data_idxs.insert(d_Phi_idx);
    return;
}// setPhiPatchDataIndex

void
INSIntermediateVelocityBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSIntermediateVelocityBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condiiton objects must be provided." << endl);
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

const std::set<int>&
INSIntermediateVelocityBcCoef::getHomogeneousBcFillDataIndices() const
{
    return d_homo_patch_data_idxs;
}// getHomogeneousBcFillDataIndices

const std::set<int>&
INSIntermediateVelocityBcCoef::getInhomogeneousBcFillDataIndices() const
{
    return d_inhomo_patch_data_idxs;
}// getInhomogeneousBcFillDataIndices

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
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
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
         ? SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >(NULL)
         : new SAMRAI::pdat::ArrayData<NDIM,double>(acoef_data->getBox(), acoef_data->getDepth()));
    setBcCoefs_private(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    TBOX_ERROR("INSIntermediateVelocityBcCoef::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#endif
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
INSIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_u_bc_coefs.size() == NDIM);
    for (unsigned l = 0; l < d_u_bc_coefs.size(); ++l)
    {
        assert(d_u_bc_coefs[l] != NULL);
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
    assert(d_u_bc_coefs.size() == NDIM);
    for (unsigned l = 0; l < d_u_bc_coefs.size(); ++l)
    {
        assert(d_u_bc_coefs[l] != NULL);
    }
#endif
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis =  location_index / 2      ;
    const bool bdry_lower_side = (location_index % 2) == 0;

    // Set the "true" velocity bc coefs.
#if USING_OLD_ROBIN_BC_INTERFACE
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif

    // Set the homogeneous boundary conditions.
    if (d_homogeneous_bc && !gcoef_data.isNull()) gcoef_data->fillAll(0.0);

    // Modify the normal velocity boundary condition to enforce div u = 0 at the
    // boundary to O(dx^2) at "open" boundaries.
    if (bdry_normal_axis == d_comp_idx && !gcoef_data.isNull())
    {
        const SAMRAI::hier::Box<NDIM>& patch_box = patch.getBox();
        const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
        const SAMRAI::hier::Index<NDIM>& patch_upper = patch_box.upper();

        const SAMRAI::hier::Box<NDIM> bc_coef_box =
            STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > U_data =
            patch.checkAllocated(d_target_idx)
            ? patch.getPatchData(d_target_idx)
            : SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(NULL);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!U_data.isNull());
        assert(!acoef_data.isNull() || !bcoef_data.isNull());
#endif
        for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const SAMRAI::hier::Index<NDIM>& i = b();
            if ((!acoef_data.isNull() && SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),0.0)) ||
                (!bcoef_data.isNull() && SAMRAI::tbox::Utilities::deq((*bcoef_data)(i,0),1.0)))
            {
                // Check to see if the boundary point is near a edge or corner
                // in the computational domain.
                bool touches_edge_or_corner = false;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    if (axis != bdry_normal_axis)
                    {
                        touches_edge_or_corner = touches_edge_or_corner ||
                            (pgeom->getTouchesRegularBoundary(axis,0) && i(axis) <= patch_lower(axis)) ||
                            (pgeom->getTouchesRegularBoundary(axis,1) && i(axis) >= patch_upper(axis));
                    }
                }

                // Compute the divergence of U at the boundary, (div U)_b, via
                // linear extrapolation using the extrapolated ghost cell value.
                double div_U_b = 0.0;
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    SAMRAI::hier::Index<NDIM> i_lower(i), i_upper(i);
                    if (axis == bdry_normal_axis)
                    {
                        if (bdry_lower_side)
                        {
                            i_lower(bdry_normal_axis) -= 1;
                            i_upper(bdry_normal_axis) += 1;
                        }
                        else
                        {
                            i_lower(bdry_normal_axis) -= 2;
                            i_upper(bdry_normal_axis) += 0;
                        }
                    }
                    else
                    {
                        i_lower(axis) -= 1;
                        i_upper(axis) += 1;
                        if (!bdry_lower_side)
                        {
                            i_lower(bdry_normal_axis) -= 1;
                            i_upper(bdry_normal_axis) -= 1;
                        }
                    }

                    div_U_b += ((*U_data)(i_upper,axis) - (*U_data)(i_lower,axis))/(double(i_upper(axis)-i_lower(axis))*dx[axis]);
#if 0
                    if (bdry_lower_side)
                    {
                        i_lower(axis) += 1;
                        i_upper(axis) += 1;
                    }
                    else
                    {
                        i_lower(axis) -= 1;
                        i_upper(axis) -= 1;
                    }

                    div_U_b -= 0.25*((*U_data)(i_upper,axis) - (*U_data)(i_lower,axis))/(double(i_upper(axis)-i_lower(axis))*dx[axis]);
#endif
                }

                // Rescale div_U_b and subtract the extrapolated ghost cell
                // value from the approximation to (div U)_b.
                double F = 2.0*dx[bdry_normal_axis]*div_U_b;

                SAMRAI::hier::Index<NDIM> i_intr(i), i_bdry(i);
                if (bdry_lower_side)
                {
                    i_bdry(bdry_normal_axis) -= 1;
                    F += (*U_data)(i_bdry,bdry_normal_axis);
                }
                else
                {
                    i_intr(bdry_normal_axis) -= 1;
                    F -= (*U_data)(i_bdry,bdry_normal_axis);
                }

                // Determine the ghost value and the appropriate value for the
                // boundary condition coefficient.
                const double& U_intr = (*U_data)(i_intr,bdry_normal_axis);
                const double& U_bdry = (bdry_lower_side ? +1.0 : -1.0)*F;
                if (!acoef_data.isNull()) (*acoef_data)(i,0) = 1.0;
                if (!bcoef_data.isNull()) (*bcoef_data)(i,0) = 0.0;
                (*gcoef_data)(i,0) = touches_edge_or_corner ? 0.0 : (U_bdry + U_intr)/2.0;
            }
        }
    }

    // XXXX: Get rid of Phi????

#if 0
    // Loop over the boundary box and reset the inhomogeneous coefficients for
    // the tangential components of the intermediate velocity.
    const double dt = d_new_time - d_current_time;
    if (!d_homogeneous_bc && !gcoef_data.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!acoef_data.isNull());
#endif
        const SAMRAI::hier::Box<NDIM> bc_coef_box =
            STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > Phi_data = patch.getPatchData(d_Phi_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!Phi_data.isNull());
#endif
        const SAMRAI::hier::Box<NDIM> ghost_box = Phi_data->getGhostBox();
        const SAMRAI::hier::Index<NDIM> ghost_upper = ghost_box.upper();
        const SAMRAI::hier::Index<NDIM> ghost_lower = ghost_box.lower();
        for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const SAMRAI::hier::Index<NDIM>& i = b();
            if (bdry_normal_axis != d_comp_idx)
            {
                if (SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),1.0))
                {
                    SAMRAI::hier::Index<NDIM> i_up0(i), i_up1(i), i_down0(i), i_down1(i);
                    i_up0  (d_comp_idx) += 1;
                    i_up1  (d_comp_idx) += 1;

                    i_down0(d_comp_idx) -= 1;
                    i_down1(d_comp_idx) -= 1;

                    i_up1  (bdry_normal_axis) -= 1;
                    i_down1(bdry_normal_axis) -= 1;

                    if (i_up0(d_comp_idx) > ghost_upper(d_comp_idx))
                    {
                        i_up0  (d_comp_idx) -= 1;
                        i_up1  (d_comp_idx) -= 1;
                        i_down0(d_comp_idx) -= 1;
                        i_down1(d_comp_idx) -= 1;
                    }
                    else if (i_down0(d_comp_idx) < ghost_lower(d_comp_idx))
                    {
                        i_up0  (d_comp_idx) += 1;
                        i_up1  (d_comp_idx) += 1;
                        i_down0(d_comp_idx) += 1;
                        i_down1(d_comp_idx) += 1;
                    }

                    const double grad_Phi = 0.5*(
                        0.5*((*Phi_data)(i_up0)-(*Phi_data)(i_down0))/dx[d_comp_idx] +
                        0.5*((*Phi_data)(i_up1)-(*Phi_data)(i_down1))/dx[d_comp_idx]);

                    (*gcoef_data)(i,0) += (dt/d_rho)*grad_Phi;
                }
                else
                {
                    assert(false);
                }
            }
        }
    }
#endif

    return;
}// setBcCoefs_private

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
