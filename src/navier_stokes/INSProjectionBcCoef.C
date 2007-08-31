// Filename: INSProjectionBcCoef.C
// Last modified: <30.Aug.2007 20:19:42 griffith@box221.cims.nyu.edu>
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

// STOOLS INCLUDES
#include <stools/PhysicalBoundaryUtilities.h>

// SAMRAI INCLUDES
#include <FaceData.h>
#include <SideData.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSProjectionBcCoef::INSProjectionBcCoef(
    const int u_idx,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_rho(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_u_idx(-1),
      d_homo_patch_data_idxs(),
      d_inhomo_patch_data_idxs(),
      d_u_bc_coefs(),
      d_homogeneous_bc(false)
{
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
INSProjectionBcCoef::setIntermediateVelocityPatchDataIndex(
    const int u_idx)
{
    d_u_idx = u_idx;
    d_inhomo_patch_data_idxs.clear();
    d_inhomo_patch_data_idxs.insert(d_u_idx);
    return;
}// setIntermediateVelocityPatchDataIndex

void
INSProjectionBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSProjectionBcCoef::setVelocityPhysicalBcCoefs():\n"
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
INSProjectionBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

const std::set<int>&
INSProjectionBcCoef::getHomogeneousBcFillDataIndices() const
{
    return d_homo_patch_data_idxs;
}// getHomogeneousBcFillDataIndices

const std::set<int>&
INSProjectionBcCoef::getInhomogeneousBcFillDataIndices() const
{
    return d_inhomo_patch_data_idxs;
}// getInhomogeneousBcFillDataIndices

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
#if USING_OLD_ROBIN_BC_INTERFACE
    TBOX_ERROR("INSProjectionBcCoef::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#else
    setBcCoefs_private(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif
    return;
}// setBcCoefs

void
INSProjectionBcCoef::setBcCoefs(
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
    TBOX_ERROR("INSProjectionBcCoef::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#endif
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
INSProjectionBcCoef::numberOfExtensionsFillable() const
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
INSProjectionBcCoef::setBcCoefs_private(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
    const SAMRAI::hier::Patch<NDIM>& patch,
    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    const int location_index   = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index/2;
    const bool is_lower        = location_index%2 == 0;

    // Set the normal velocity bc coefs.
#if USING_OLD_ROBIN_BC_INTERFACE
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif

    // Loop over the boundary box and reset the homogeneous coefficients.
    const SAMRAI::hier::Box<NDIM> bc_coef_box =
        STOOLS::PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_box);

    const bool fill_acoef_data = !acoef_data.isNull();
    const bool fill_bcoef_data = !bcoef_data.isNull();
    const bool fill_gcoef_data = !gcoef_data.isNull();

    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i = b();

        if (fill_acoef_data) (*acoef_data)(i,0) = 1.0-(*acoef_data)(i,0);
        if (fill_bcoef_data) (*bcoef_data)(i,0) = 1.0-(*bcoef_data)(i,0);
        if (fill_gcoef_data && d_homogeneous_bc) (*gcoef_data)(i,0) = 0.0;
    }

    // Loop over the boundary box and reset the inhomogeneous coefficients.
    if (!d_homogeneous_bc && !gcoef_data.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!acoef_data.isNull());
#endif
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_fc_data = patch.getPatchData(d_u_idx);
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > u_sc_data = patch.getPatchData(d_u_idx);

        if (!u_fc_data.isNull())
        {
            for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& i = b();
                const SAMRAI::pdat::FaceIndex<NDIM> i_f(
                    i, bdry_normal_axis, SAMRAI::pdat::FaceIndex<NDIM>::Lower);

                if (SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),0.0))
                {
                    (*gcoef_data)(i,0) =
                        (is_lower ? -1.0 : +1.0)*(d_rho/d_dt)*((*u_fc_data)(i_f) - (*gcoef_data)(i,0));
                }
                else
                {
                    assert(false);
                }
            }
        }
        else if (!u_sc_data.isNull())
        {
            for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
            {
                const SAMRAI::hier::Index<NDIM>& i = b();
                const SAMRAI::pdat::SideIndex<NDIM> i_s(
                    i, bdry_normal_axis, SAMRAI::pdat::SideIndex<NDIM>::Lower);

                if (SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),0.0))
                {
                    (*gcoef_data)(i,0) =
                        (is_lower ? -1.0 : +1.0)*(d_rho/d_dt)*((*u_sc_data)(i_s) - (*gcoef_data)(i,0));
                }
                else
                {
                    assert(false);
                }
            }
        }
        else
        {
            TBOX_ERROR("INSProjectionBcCoef::setBcCoefs():\n"
                       << "  intermediate velocity data required to set inhomogenous boundary conditions is not available" << std::endl);
        }
    }

    return;
}// setBcCoefs_private

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
