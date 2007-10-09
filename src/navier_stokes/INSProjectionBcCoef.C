// Filename: INSProjectionBcCoef.C
// Last modified: <09.Oct.2007 17:44:39 griffith@box221.cims.nyu.edu>
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
#include <CellData.h>
#include <FaceData.h>
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
      d_u_bc_coefs(NDIM,NULL),
      d_homo_patch_data_idxs(),
      d_inhomo_patch_data_idxs(),
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
    d_inhomo_patch_data_idxs.clear();
    if (d_u_idx != -1) d_inhomo_patch_data_idxs.insert(d_u_idx);
    d_inhomo_patch_data_idxs.insert(d_P_idx);
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
                   << "  valid choices are: ``pressure_increment'' or ``pressure_update''" << endl);
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
    d_inhomo_patch_data_idxs.clear();
    d_inhomo_patch_data_idxs.insert(d_u_idx);
    if (d_P_idx != -1) d_inhomo_patch_data_idxs.insert(d_P_idx);
    return;
}// setIntermediateVelocityPatchDataIndex

void
INSProjectionBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSProjectionBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << endl);
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
INSProjectionBcCoef::setBcCoefs_private(
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

    // Set the normal velocity bc coefs.
#if USING_OLD_ROBIN_BC_INTERFACE
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    d_u_bc_coefs[bdry_normal_axis]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif

    // Loop over the boundary box and reset the homogeneous coefficients.
    for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
    {
        const SAMRAI::hier::Index<NDIM>& i = b();

        if (fill_acoef_data) (*acoef_data)(i,0) = 1.0-(*acoef_data)(i,0);
        if (fill_bcoef_data) (*bcoef_data)(i,0) = 1.0-(*bcoef_data)(i,0);
        if (fill_gcoef_data && d_homogeneous_bc) (*gcoef_data)(i,0) = 0.0;
    }

    // Loop over the boundary box and reset the inhomogeneous coefficients.
    if (!d_homogeneous_bc && fill_gcoef_data)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM,double> > u_data = patch.getPatchData(d_u_idx);
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > P_data = patch.getPatchData(d_P_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!u_data.isNull());
        assert(!P_data.isNull());
        assert(fill_acoef_data || fill_bcoef_data);
#endif
        SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > pcoef_data =
            (acoef_data.isNull()
             ? new SAMRAI::pdat::ArrayData<NDIM,double>(bcoef_data->getBox(), bcoef_data->getDepth())
             : new SAMRAI::pdat::ArrayData<NDIM,double>(acoef_data->getBox(), acoef_data->getDepth()));

        if (d_P_bc_coef != NULL)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> > null_data(NULL);
#if USING_OLD_ROBIN_BC_INTERFACE
            d_P_bc_coef->setBcCoefs(
                null_data, pcoef_data, variable, patch, bdry_box, fill_time);
#else
            d_P_bc_coef->setBcCoefs(
                null_data, null_data, pcoef_data, variable, patch, bdry_box, fill_time);
#endif
        }
        else
        {
            pcoef_data->fillAll(0.0);
        }

        const bool using_pressure_increment = (d_projection_type == "pressure_increment");
        if (using_pressure_increment)
        {
            SAMRAI::tbox::plog << "P_idx = " << d_P_idx << "\n";
            P_data->print(P_data->getGhostBox());
        }

        for (SAMRAI::hier::Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const SAMRAI::hier::Index<NDIM>& i = b();
            const SAMRAI::pdat::FaceIndex<NDIM> i_f(
                i, bdry_normal_axis, SAMRAI::pdat::FaceIndex<NDIM>::Lower);

            if ((fill_acoef_data && SAMRAI::tbox::Utilities::deq((*acoef_data)(i,0),1.0)) ||
                (fill_bcoef_data && SAMRAI::tbox::Utilities::deq((*bcoef_data)(i,0),0.0)))
            {
                SAMRAI::hier::Index<NDIM> i_intr0(i), i_intr1(i);
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

                const double P_bdry =
                    (using_pressure_increment
                     ? 1.5*(*P_data)(i_intr0) - 0.5*(*P_data)(i_intr1)
                     : 0.0);
                (*gcoef_data)(i,0) = (*pcoef_data)(i,0) - P_bdry;
            }
            else
            {
                (*gcoef_data)(i,0) = (is_lower ? -1.0 : +1.0)*(d_rho/d_dt)*((*u_data)(i_f) - (*gcoef_data)(i,0));
            }
        }
    }
    return;
}// setBcCoefs_private

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
