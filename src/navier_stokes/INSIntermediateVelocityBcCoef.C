// Filename: INSIntermediateVelocityBcCoef.C
// Last modified: <12.Oct.2007 02:53:29 griffith@box221.cims.nyu.edu>
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
    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;
    const bool is_lower        = location_index%2 == 0;

    const SAMRAI::hier::Box<NDIM> bc_coef_box =
        STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);

    const bool fill_acoef_data = !acoef_data.isNull();
    const bool fill_bcoef_data = !bcoef_data.isNull();
    const bool fill_gcoef_data = !gcoef_data.isNull();

    // Set the "true" velocity bc coefs.
#if USING_OLD_ROBIN_BC_INTERFACE
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data.
    if (!fill_gcoef_data)
    {
        return;
    }

    // Patch box information.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    // Enforce homogeneous boundary conditions.
    if (d_homogeneous_bc || d_velocity_correction) gcoef_data->fillAll(0.0);

#if 0  // XXXX
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
    assert(!U_data.isNull());
    assert(fill_acoef_data || fill_bcoef_data);
#endif
    const SAMRAI::hier::Box<NDIM>& U_ghost_box = U_data->getGhostBox();
    const SAMRAI::hier::Index<NDIM>& U_ghost_lower = U_ghost_box.lower();
    const SAMRAI::hier::Index<NDIM>& U_ghost_upper = U_ghost_box.upper();
    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i = b();
        if ((fill_acoef_data && SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),0.0)) ||
            (fill_bcoef_data && SAMRAI::tbox::Utilities::deq((*bcoef_data)(i,0),1.0)))
        {
            // Compute the flux to satisfy either div U = 0 or t * sigma * n = 0
            // at the boundary.
            double F = 0.0;
            for (int axis = 0; axis < NDIM; ++axis)
            {
                if (axis != bdry_normal_axis)
                {
                    // Determine the component of the velocity that we are
                    // differencing.
                    const int comp_idx = (d_comp_idx == bdry_normal_axis ? axis : bdry_normal_axis);

                    // Setup the basic difference stencil, located inside the
                    // computational domain.
                    SAMRAI::hier::Index<NDIM> i_lower(i), i_upper(i);
                    i_lower(axis) -= 1;
                    i_upper(axis) += 1;

                    if (is_lower)
                    {
                        // intentionally blank
                    }
                    else
                    {
                        i_lower(bdry_normal_axis) -= 1;
                        i_upper(bdry_normal_axis) -= 1;
                    }

                    // Ensure the difference stencil lies within the ghost cell
                    // region of the patch data.
                    if      (i_lower(axis) <= U_ghost_lower(axis))
                    {
                        i_lower(axis) = U_ghost_lower(axis);
                        i_upper(axis) = U_ghost_lower(axis)+1;
                    }
                    else if (i_upper(axis) >= U_ghost_upper(axis))
                    {
                        i_lower(axis) = U_ghost_upper(axis)-1;
                        i_upper(axis) = U_ghost_upper(axis);
                    }

                    if ((d_comp_idx == bdry_normal_axis) || axis == d_comp_idx)
                    {
                        const double& U_upper = (*U_data)(i_upper,comp_idx);
                        const double& U_lower = (*U_data)(i_lower,comp_idx);
                        F += 0.5*(U_upper-U_lower)/(double(i_upper(axis)-i_lower(axis))*dx[axis]);
                    }

                    // Shift the difference stencil outside the computational
                    // domain.
                    if (is_lower)
                    {
                        i_lower(bdry_normal_axis) -= 1;
                        i_upper(bdry_normal_axis) -= 1;
                    }
                    else
                    {
                        i_lower(bdry_normal_axis) += 1;
                        i_upper(bdry_normal_axis) += 1;
                    }

                    if ((d_comp_idx == bdry_normal_axis) || axis == d_comp_idx)
                    {
                        const double& U_upper = (*U_data)(i_upper,comp_idx);
                        const double& U_lower = (*U_data)(i_lower,comp_idx);
                        F += 0.5*(U_upper-U_lower)/(double(i_upper(axis)-i_lower(axis))*dx[axis]);
                    }
                }
            }

            // Set the normal flux.
            (*gcoef_data)(i,0) = (is_lower ? +1.0 : -1.0)*F;
        }
    }
#endif

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
    assert(!Phi_data.isNull());
    assert(fill_acoef_data || fill_bcoef_data);
#endif
    const SAMRAI::hier::Box<NDIM>& Phi_ghost_box = Phi_data->getGhostBox();
    const SAMRAI::hier::Index<NDIM>& Phi_ghost_lower = Phi_ghost_box.lower();
    const SAMRAI::hier::Index<NDIM>& Phi_ghost_upper = Phi_ghost_box.upper();
    const double dt = d_new_time - d_current_time;
    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i = b();
        if ((fill_acoef_data && SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),1.0)) ||
            (fill_bcoef_data && SAMRAI::tbox::Utilities::deq((*bcoef_data)(i,0),0.0)))
        {
            // Set the boundary condition coefficients to satisfy U*t = g at the
            // boundary.
            const int axis = d_comp_idx;

            // Setup the basic difference stencil, located inside the
            // computational domain.
            SAMRAI::hier::Index<NDIM> i_lower(i), i_upper(i);
            i_lower(axis) -= 1;
            i_upper(axis) += 1;

            if (!is_lower)
            {
                i_lower(bdry_normal_axis) -= 1;
                i_upper(bdry_normal_axis) -= 1;
            }

            // Ensure the difference stencil lies within the ghost cell
            // region of the patch data.
            if      (i_lower(axis) <= Phi_ghost_lower(axis))
            {
                i_lower(axis) = Phi_ghost_lower(axis);
                i_upper(axis) = Phi_ghost_lower(axis)+1;
            }
            else if (i_upper(axis) >= Phi_ghost_upper(axis))
            {
                i_lower(axis) = Phi_ghost_upper(axis)-1;
                i_upper(axis) = Phi_ghost_upper(axis);
            }

            const double grad_Phi = ((*Phi_data)(i_upper)-(*Phi_data)(i_lower))/(double(i_upper(axis)-i_lower(axis))*dx[axis]);

            // Approximately enforce U*t = g.
            (*gcoef_data)(i,0) += (dt/d_rho)*grad_Phi;
        }
        else
        {
            // Set the boundary condition coefficeints to satisfy zero
            // tangential stress at the boundary.
            const int axis = d_comp_idx;

            // Setup the basic difference stencil.
            SAMRAI::hier::Index<NDIM> i_intr_lower(i), i_bdry_lower(i), i_intr_upper(i), i_bdry_upper(i);

            i_intr_lower(axis) -= 1;
            i_bdry_lower(axis) -= 1;

            i_intr_upper(axis) += 1;
            i_bdry_upper(axis) += 1;

            if (is_lower)
            {
                i_bdry_lower(bdry_normal_axis) -= 1;
                i_bdry_upper(bdry_normal_axis) -= 1;
            }
            else
            {
                i_intr_lower(bdry_normal_axis) -= 1;
                i_intr_upper(bdry_normal_axis) -= 1;
            }

            // Ensure the difference stencil lies within the ghost cell
            // region of the patch data.
            if      (i_intr_lower(axis) <= Phi_ghost_lower(axis) ||
                     i_bdry_lower(axis) <= Phi_ghost_lower(axis))
            {
                i_intr_lower(axis) = Phi_ghost_lower(axis);
                i_bdry_lower(axis) = Phi_ghost_lower(axis);

                i_intr_upper(axis) = Phi_ghost_lower(axis)+1;
                i_bdry_upper(axis) = Phi_ghost_lower(axis)+1;
            }
            else if (i_intr_upper(axis) >= Phi_ghost_upper(axis) ||
                     i_bdry_upper(axis) >= Phi_ghost_upper(axis))
            {
                i_intr_lower(axis) = Phi_ghost_upper(axis)-1;
                i_bdry_lower(axis) = Phi_ghost_upper(axis)-1;

                i_intr_upper(axis) = Phi_ghost_upper(axis);
                i_bdry_upper(axis) = Phi_ghost_upper(axis);
            }

            const double t_dot_grad_grad_Phi_dot_n =
                (((*Phi_data)(i_bdry_upper)-(*Phi_data)(i_intr_upper))/(double(i_bdry_upper(bdry_normal_axis)-i_intr_upper(bdry_normal_axis))*dx[bdry_normal_axis]) -
                 ((*Phi_data)(i_bdry_lower)-(*Phi_data)(i_intr_lower))/(double(i_bdry_lower(bdry_normal_axis)-i_intr_lower(bdry_normal_axis))*dx[bdry_normal_axis]))/
                (2.0*double(i_intr_upper(axis)-i_intr_lower(axis))*dx[axis]);

            // Approximately enforce the tangential stress boundary condition.
            (*gcoef_data)(i,0) += 2.0*(dt/d_rho)*t_dot_grad_grad_Phi_dot_n;
        }
    }
    return;
}// setBcCoefs_private

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
