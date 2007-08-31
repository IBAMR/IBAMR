// Filename: INSIntermediateVelocityBcCoef.C
// Last modified: <30.Aug.2007 23:52:57 griffith@box221.cims.nyu.edu>
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
    : d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_intermediate_velocity_fix(false),
      d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_comp_idx(comp_idx),
      d_Phi_idx(-1),
      d_homo_patch_data_idxs(),
      d_inhomo_patch_data_idxs(),
      d_u_bc_coefs(),
      d_homogeneous_bc(false)
{
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
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < u_bc_coefs.size(); ++l)
    {
        assert(u_bc_coefs[l] != NULL);
    }
#endif
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

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
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();

    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;

    // Set the normal velocity bc coefs.
#if USING_OLD_ROBIN_BC_INTERFACE
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif

    if (!d_intermediate_velocity_fix || SAMRAI::tbox::Utilities::deq(fill_time, d_current_time)) return;

    // Loop over the boundary box and reset the inhomogeneous coefficients.
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

                    (*gcoef_data)(i,0) += (d_dt/d_rho)*grad_Phi;
                }
                else
                {
                    assert(false);
                }
            }
        }
    }

    return;
}// setBcCoefs_private

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
