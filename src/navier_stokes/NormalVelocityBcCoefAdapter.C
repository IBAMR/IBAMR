// Filename: NormalVelocityBcCoefAdapter.C
// Last modified: <13.Feb.2008 13:32:07 griffith@box221.cims.nyu.edu>
// Created on 22 Feb 2007 by Boyce Griffith (boyce@trasnaform2.local)

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
    setPhysicalBcCoefs(bc_coefs);
    return;
}// NormalVelocityBcCoefAdapter

NormalVelocityBcCoefAdapter::~NormalVelocityBcCoefAdapter()
{
    // intentionally blank
    return;
}// ~NormalVelocityBcCoefAdapter

void
NormalVelocityBcCoefAdapter::setPhysicalBcCoefs(
    const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    if (bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("NormalVelocityBcCoefAdapter::setPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condiiton objects must be provided." << std::endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(bc_coefs[l] != NULL);
    }
#endif
    d_bc_coefs = bc_coefs;
    return;
}// setPhysicalBcCoefs

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
    const int depth = bdry_box.getLocationIndex()/2;
    d_bc_coefs[depth]->setBcCoefs(
        acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
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
