// Filename: Stashable.C
// Last modified: <24.Jun.2007 21:25:17 griffith@box221.cims.nyu.edu>
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "Stashable.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

Stashable::Stashable()
{
    // intentionally blank
    return;
}// Stashable

Stashable::~Stashable()
{
    // intentionally blank
    return;
}// ~Stashable

void
Stashable::registerPeriodicShift(
    const SAMRAI::hier::IntVector<NDIM>& offset,
    const std::vector<double>& displacement)
{
    // intentionally blank
    return;
}// registerPeriodicShift

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::Stashable>;

//////////////////////////////////////////////////////////////////////////////
