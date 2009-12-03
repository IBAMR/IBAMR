// Filename: SpringForceSpec.C
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <25.Oct.2006 18:30:20 boyce@bigboy.nyconnect.com>

#include "SpringForceSpec.h"

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
#include <ibamr/SpringForceSpecFactory.h>
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

bool SpringForceSpec::s_registered_factory = false;
int  SpringForceSpec::s_stashable_id = -1;

void
SpringForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes
    // actually register the stashable factory with the stashable
    // manager, and to ensure that all processes employ the same
    // stashable id for the SpringForceSpec object.
    SAMRAI::tbox::MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new SpringForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

bool
SpringForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

SpringForceSpec::SpringForceSpec(
    const std::vector<int>& dst_idxs,
    const std::vector<double>& stiffnesses,
    const std::vector<double>& rest_lengths)
    : d_num_links(static_cast<int>(dst_idxs.size())),
      d_dst_idxs(dst_idxs),
      d_stiffnesses(stiffnesses),
      d_rest_lengths(rest_lengths)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_links == static_cast<int>(d_dst_idxs.size()));
    assert(d_num_links == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_links == static_cast<int>(d_rest_lengths.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("SpringForceSpec::SpringForceSpec():\n"
                   << "  must call SpringForceSpec::registerWithStashableManager() before\n"
                   << "  creating any SpringForceSpec objects.\n");
    }
    return;
}// SpringForceSpec

SpringForceSpec::~SpringForceSpec()
{
    // intentionally blank
    return;
}// ~SpringForceSpec

int
SpringForceSpec::getNumberOfLinks() const
{
    return d_num_links;
}// getNumberOfLinks

const std::vector<int>&
SpringForceSpec::getDestinationNodeIndices() const
{
    return d_dst_idxs;
}// getDestinationNodeIndices

vector<int>&
SpringForceSpec::getDestinationNodeIndices()
{
    return d_dst_idxs;
}// getDestinationNodeIndices

const std::vector<double>&
SpringForceSpec::getStiffnesses() const
{
    return d_stiffnesses;
}// getStiffnesses

vector<double>&
SpringForceSpec::getStiffnesses()
{
    return d_stiffnesses;
}// getStiffnesses

const std::vector<double>&
SpringForceSpec::getRestingLengths() const
{
    return d_rest_lengths;
}// getRestingLengths

vector<double>&
SpringForceSpec::getRestingLengths()
{
    return d_rest_lengths;
}// getRestingLengths

int
SpringForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

size_t
SpringForceSpec::getDataStreamSize() const
{
    return ((1+d_num_links)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (2*d_num_links)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

void
SpringForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_links == static_cast<int>(d_dst_idxs.size()));
    assert(d_num_links == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_links == static_cast<int>(d_rest_lengths.size()));
#endif
    stream.pack(&d_num_links,1);
    stream.pack(&d_dst_idxs[0],d_num_links);
    stream.pack(&d_stiffnesses[0],d_num_links);
    stream.pack(&d_rest_lengths[0],d_num_links);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::SpringForceSpec>;

//////////////////////////////////////////////////////////////////////////////
