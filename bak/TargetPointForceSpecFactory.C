// Filename: TargetPointForceSpecFactory.C
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <25.Oct.2006 18:30:48 boyce@bigboy.nyconnect.com>

#include "TargetPointForceSpecFactory.h"

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
#include <ibamr/TargetPointForceSpec.h>
#include <ibamr/StashableManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int TargetPointForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

TargetPointForceSpecFactory::TargetPointForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// TargetPointForceSpecFactory

TargetPointForceSpecFactory::~TargetPointForceSpecFactory()
{
    // intentionally blank
    return;
}// ~TargetPointForceSpecFactory

int
TargetPointForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
TargetPointForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<Stashable>
TargetPointForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    std::vector<double> X(NDIM);
    stream.unpack(&X[0], NDIM);
    double kappa;
    stream.unpack(&kappa, 1);

    return new TargetPointForceSpec(X, kappa);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::TargetPointForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
