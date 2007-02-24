// Filename: NormalVelocityBcCoefAdapter.C
// Last modified: <23.Feb.2007 22:10:12 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 30 Sep 2006 by Boyce Griffith (boyce@boyce-griffiths-powerbook-g4-15.local)

#include "NormalVelocityBcCoefAdapter.h"

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
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

NormalVelocityBcCoefAdapter::NormalVelocityBcCoefAdapter(
    const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    : d_bc_coefs()
{
    setRobinBcCoefStrategies(bc_coefs);
    return;
}// NormalVelocityBcCoefAdapter

NormalVelocityBcCoefAdapter::~NormalVelocityBcCoefAdapter()
{
    // intentionally blank
    return;
}// ~NormalVelocityBcCoefAdapter

void
NormalVelocityBcCoefAdapter::setRobinBcCoefStrategies(
    const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    if (bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("NormalVelocityBcCoefAdapter::setRobinBcCoefStrategies():\n"
                   << "  precisely NDIM boundary condiiton objects must be provided." << endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < bc_coefs.size(); ++l)
    {
        assert(bc_coefs[l] != NULL);
    }
#endif
    d_bc_coefs = bc_coefs;
    return;
}// setRobinBcCoefStrategies

void
NormalVelocityBcCoefAdapter::setBcCoefs(
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& bcoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
     const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
     const SAMRAI::hier::Patch<NDIM>& patch,
     const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
     double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    TBOX_ERROR("NormalVelocityBcCoefAdapter::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#else
    const int depth = bdry_box.getLocationIndex()/2;
    d_bc_coefs[depth]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#endif
    return;
}// setBcCoefs

void
NormalVelocityBcCoefAdapter::setBcCoefs(
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& acoef_data,
     SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM,double> >& gcoef_data,
     const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
     const SAMRAI::hier::Patch<NDIM>& patch,
     const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
     double fill_time) const
{
#if USING_OLD_ROBIN_BC_INTERFACE
    const int depth = bdry_box.getLocationIndex()/2;
    d_bc_coefs[depth]->setBcCoefs(
        acoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
#else
    TBOX_ERROR("NormalVelocityBcCoefAdapter::setBcCoefs():\n"
               << "  using incorrect SAMRAI::solv::RobinBcCoefStrategy interface." << endl);
#endif
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
NormalVelocityBcCoefAdapter::numberOfExtensionsFillable() const
{
    SAMRAI::hier::IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (int d = 0; d < NDIM; ++d)
    {
        ret_val = SAMRAI::hier::IntVector<NDIM>::min(
            ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
