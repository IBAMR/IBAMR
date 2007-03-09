// Filename: IntermediateVelocityRobinBcCoefs.C
// Last modified: <09.Mar.2007 18:44:14 griffith@box221.cims.nyu.edu>
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

// SAMRAI INCLUDES
#include <ArrayDataBasicOps.h>
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
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef)
    : d_velocity_depth(velocity_depth),
      d_bc_coef(d_bc_coef)
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
    correctBcCoefs(acoef_data, bcoef_data, gcoef_data, patch, bdry_box);
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

    correctBcCoefs(acoef_data, bcoef_data, gcoef_data, patch, bdry_box);
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

    return;
}// correctBcCoefs

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
