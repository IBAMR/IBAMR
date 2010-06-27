// Filename: INSStaggeredIntermediateVelocityBcCoef.C
// Last modified: <27.Jun.2010 15:27:08 griffith@griffith-macbook-pro.local>
// Created on 24 Jul 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredIntermediateVelocityBcCoef.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredIntermediateVelocityBcCoef::INSStaggeredIntermediateVelocityBcCoef(
    const int comp_idx,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
    const bool homogeneous_bc)
    : d_comp_idx(comp_idx),
      d_u_bc_coefs(NDIM,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_target_idx(-1),
      d_homogeneous_bc(false)
{
    setVelocityPhysicalBcCoefs(u_bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
}// INSStaggeredIntermediateVelocityBcCoef

INSStaggeredIntermediateVelocityBcCoef::~INSStaggeredIntermediateVelocityBcCoef()
{
    // intentionally blank
    return;
}// ~INSStaggeredIntermediateVelocityBcCoef

void
INSStaggeredIntermediateVelocityBcCoef::setVelocityPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR("INSStaggeredIntermediateVelocityBcCoef::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condition objects must be provided." << std::endl);
    }
    d_u_bc_coefs = u_bc_coefs;
    return;
}// setVelocityPhysicalBcCoefs

void
INSStaggeredIntermediateVelocityBcCoef::setTargetPatchDataIndex(
    const int target_idx)
{
    d_target_idx = target_idx;
    return;
}// setTargetPatchDataIndex

void
INSStaggeredIntermediateVelocityBcCoef::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredIntermediateVelocityBcCoef::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
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
    d_u_bc_coefs[d_comp_idx]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // Ensure homogeneous boundary conditions are prescribed.
    if (!gcoef_data.isNull()) gcoef_data->fillAll(0.0);
    return;
}// setBcCoefs

IntVector<NDIM>
INSStaggeredIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_bc_coefs.size() == NDIM);
    for (unsigned l = 0; l < d_u_bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(d_u_bc_coefs[l] != NULL);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_u_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
