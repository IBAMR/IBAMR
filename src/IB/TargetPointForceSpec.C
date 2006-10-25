// Filename: TargetPointForceSpec.C
// Created on 23 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)
// Last modified: <24.Oct.2006 14:42:30 boyce@bigboy.nyconnect.com>

#include "TargetPointForceSpec.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#define included_IBAMR_config
#include <IBAMR_config.h>
#endif

#ifndef included_SAMRAI_config
#define included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

// IBAMR INCLUDES
#include <ibamr/TargetPointForceSpecFactory.h>
#include <ibamr/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/MPI.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool TargetPointForceSpec::s_registered_factory = false;
int  TargetPointForceSpec::s_stashable_id = -1;

void
TargetPointForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes
    // actually register the stashable factory with the stashable
    // manager, and to ensure that all processes employ the same
    // stashable id for the TargetPointForceSpec object.
    SAMRAI::tbox::MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new TargetPointForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

bool
TargetPointForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

TargetPointForceSpec::TargetPointForceSpec(
    const std::vector<double>& X,
    const double kappa)
    : d_X(X),
      d_kappa(kappa)
{
#ifdef DEBUG_CHECK_ASSERTSIONS
    assert(d_X.size() == NDIM);
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("TargetPointForceSpec::TargetPointForceSpec():\n"
                   << "  must call TargetPointForceSpec::registerWithStashableManager() before\n"
                   << "  creating any TargetPointForceSpec objects.\n");
    }
    return;
}// TargetPointForceSpec

TargetPointForceSpec::~TargetPointForceSpec()
{
    // intentionally blank
    return;
}// ~TargetPointForceSpec

const std::vector<double>&
TargetPointForceSpec::getPosition() const
{
    return d_X;
}// getPosition

const double&
TargetPointForceSpec::getStiffness() const
{
    return d_kappa;
}// getStiffness

int
TargetPointForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

size_t
TargetPointForceSpec::getDataStreamSize() const
{
    return (NDIM+1)*SAMRAI::tbox::AbstractStream::sizeofDouble();
}// getDataStreamSize

void
TargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_X[0], NDIM);
    stream.pack(&d_kappa, 1);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::TargetPointForceSpec>;

//////////////////////////////////////////////////////////////////////////////
