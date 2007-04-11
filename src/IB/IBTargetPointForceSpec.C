// Filename: IBTargetPointForceSpec.C
// Last modified: <11.Apr.2007 04:11:10 boyce@trasnaform2.local>
// Created on 21 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBTargetPointForceSpec.h"

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
#include <ibamr/IBTargetPointForceSpecFactory.h>
#include <ibamr/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/MPI.h>
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool IBTargetPointForceSpec::s_registered_factory = false;
int  IBTargetPointForceSpec::s_stashable_id = -1;

void
IBTargetPointForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the
    // IBTargetPointForceSpec object.
    SAMRAI::tbox::MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new IBTargetPointForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

bool
IBTargetPointForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBTargetPointForceSpec::IBTargetPointForceSpec(
    const int master_idx,
    const double& kappa_target,
    const std::vector<double>& X_target)
    : d_master_idx(master_idx),
      d_kappa_target(kappa_target),
      d_X_target(X_target)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(NDIM == static_cast<int>(d_X_target.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBTargetPointForceSpec::IBTargetPointForceSpec():\n"
                   << "  must call IBTargetPointForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBTargetPointForceSpec objects.\n");
    }
    return;
}// IBTargetPointForceSpec

IBTargetPointForceSpec::~IBTargetPointForceSpec()
{
    // intentionally blank
    return;
}// ~IBTargetPointForceSpec

const int&
IBTargetPointForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

int&
IBTargetPointForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

const double&
IBTargetPointForceSpec::getStiffness() const
{
    return d_kappa_target;
}// getStiffness

double&
IBTargetPointForceSpec::getStiffness()
{
    return d_kappa_target;
}// getStiffness

const std::vector<double>&
IBTargetPointForceSpec::getTargetPointPosition() const
{
    return d_X_target;
}// getTargetPointPosition

std::vector<double>&
IBTargetPointForceSpec::getTargetPointPosition()
{
    return d_X_target;
}// getTargetPointPosition

int
IBTargetPointForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

size_t
IBTargetPointForceSpec::getDataStreamSize() const
{
    return ((1     )*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (1+NDIM)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

void
IBTargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(NDIM == static_cast<int>(d_X_target.size()));
#endif
    stream.pack(&d_master_idx,1);
    stream.pack(&d_kappa_target,1);
    stream.pack(&d_X_target[0],NDIM);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBTargetPointForceSpec>;

//////////////////////////////////////////////////////////////////////////////
