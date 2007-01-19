// Filename: IBStandardForceSpec.C
// Last modified: <18.Jan.2007 16:44:56 boyce@bigboy.nyconnect.com>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBStandardForceSpec.h"

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
#include <ibamr/IBStandardForceSpecFactory.h>
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

bool IBStandardForceSpec::s_registered_factory = false;
int  IBStandardForceSpec::s_stashable_id = -1;

void
IBStandardForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes
    // actually register the stashable factory with the stashable
    // manager, and to ensure that all processes employ the same
    // stashable id for the IBStandardForceSpec object.
    SAMRAI::tbox::MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new IBStandardForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

bool
IBStandardForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceSpec::IBStandardForceSpec(
    const std::vector<int>& dst_idxs,
    const std::vector<double>& stiffnesses,
    const std::vector<double>& rest_lengths,
    const std::vector<double>& X_target,
    const double kappa_target)
    : d_num_links(static_cast<int>(dst_idxs.size())),
      d_dst_idxs(dst_idxs),
      d_stiffnesses(stiffnesses),
      d_rest_lengths(rest_lengths),
      d_X_target(X_target.empty() ||
                 SAMRAI::tbox::Utilities::deq(kappa_target,0.0)
                 ? std::vector<double>(NDIM,0.0) : X_target),
      d_kappa_target(kappa_target)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_links == static_cast<int>(d_dst_idxs.size()));
    assert(d_num_links == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_links == static_cast<int>(d_rest_lengths.size()));
    assert(d_X_target.size() == NDIM);
    assert(d_kappa_target >= 0.0);
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBStandardForceSpec::IBStandardForceSpec():\n"
                   << "  must call IBStandardForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBStandardForceSpec objects.\n");
    }
    return;
}// IBStandardForceSpec

IBStandardForceSpec::~IBStandardForceSpec()
{
    // intentionally blank
    return;
}// ~IBStandardForceSpec

int
IBStandardForceSpec::getNumberOfLinks() const
{
    return d_num_links;
}// getNumberOfLinks

const std::vector<int>&
IBStandardForceSpec::getDestinationNodeIndices() const
{
    return d_dst_idxs;
}// getDestinationNodeIndices

vector<int>&
IBStandardForceSpec::getDestinationNodeIndices()
{
    return d_dst_idxs;
}// getDestinationNodeIndices

const std::vector<double>&
IBStandardForceSpec::getStiffnesses() const
{
    return d_stiffnesses;
}// getStiffnesses

vector<double>&
IBStandardForceSpec::getStiffnesses()
{
    return d_stiffnesses;
}// getStiffnesses

const std::vector<double>&
IBStandardForceSpec::getRestingLengths() const
{
    return d_rest_lengths;
}// getRestingLengths

vector<double>&
IBStandardForceSpec::getRestingLengths()
{
    return d_rest_lengths;
}// getRestingLengths

const std::vector<double>&
IBStandardForceSpec::getTargetPosition() const
{
    return d_X_target;
}// getTargetPosition

const double&
IBStandardForceSpec::getTargetStiffness() const
{
    return d_kappa_target;
}// getTargetStiffness

int
IBStandardForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

size_t
IBStandardForceSpec::getDataStreamSize() const
{
    return ((1+d_num_links)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (2*d_num_links)*SAMRAI::tbox::AbstractStream::sizeofDouble() +
            (NDIM+1)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

void
IBStandardForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_links == static_cast<int>(d_dst_idxs.size()));
    assert(d_num_links == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_links == static_cast<int>(d_rest_lengths.size()));
    assert(d_X_target.size() == NDIM);
    assert(d_kappa_target >= 0.0);
#endif
    stream.pack(&d_num_links,1);
    stream.pack(&d_dst_idxs[0],d_num_links);
    stream.pack(&d_stiffnesses[0],d_num_links);
    stream.pack(&d_rest_lengths[0],d_num_links);
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
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardForceSpec>;

//////////////////////////////////////////////////////////////////////////////
