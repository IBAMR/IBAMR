// Filename: TargetPointForceSpec.C
// Created on 23 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)
// Last modified: <18.Jan.2007 16:40:23 boyce@bigboy.nyconnect.com>

#include "TargetPointForceSpec.h"

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
    const std::vector<double>& X_target,
    const double kappa_target)
    : d_X_target(X_target),
      d_kappa_target(kappa_target)
{
#ifdef DEBUG_CHECK_ASSERTSIONS
    assert(d_X_target.size() == NDIM);
    assert(d_kappa_target >= 0.0);
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("TargetPointForceSpec::TargetPointForceSpec():\n"
                   << "  must call TargetPointForceSpec::registerWithStashableManager() before\n"
                   << "  creating any TargetPointForceSpec objects.\n");
    }
    return;
}// TargetPointForceSpec

TargetPointForceSpec::TargetPointForceSpec(
    const double* const X_target,
    const double kappa_target)
    : d_X_target(X_target,X_target+NDIM),
      d_kappa_target(kappa_target)
{
#ifdef DEBUG_CHECK_ASSERTSIONS
    assert(d_X_target.size() == NDIM);
    assert(d_kappa_target >= 0.0);
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
    return d_X_target;
}// getPosition

const double&
TargetPointForceSpec::getStiffness() const
{
    return d_kappa_target;
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
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_X_target.size() == NDIM);
    assert(d_kappa_target >= 0.0);
#endif
    stream.pack(&d_X_target[0],NDIM);
    stream.pack(&d_kappa_target,1);
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
