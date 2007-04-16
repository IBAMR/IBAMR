// Filename: IBSpringForceSpec.C
// Last modified: <16.Apr.2007 05:38:43 boyce@bigboy.nyconnect.com>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBSpringForceSpec.h"

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
#include <ibamr/IBSpringForceSpecFactory.h>
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

bool IBSpringForceSpec::s_registered_factory = false;
int  IBSpringForceSpec::s_stashable_id = -1;

void
IBSpringForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the IBSpringForceSpec
    // object.
    SAMRAI::tbox::MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new IBSpringForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

bool
IBSpringForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceSpec::IBSpringForceSpec(
    const int master_idx,
    const std::vector<int>& slave_idxs,
    const std::vector<int>& force_fcn_idxs,
    const std::vector<double>& stiffnesses,
    const std::vector<double>& rest_lengths)
    : d_num_springs(static_cast<int>(slave_idxs.size())),
      d_master_idx(master_idx),
      d_slave_idxs(slave_idxs),
      d_force_fcn_idxs(force_fcn_idxs),
      d_stiffnesses(stiffnesses),
      d_rest_lengths(rest_lengths)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_springs == static_cast<int>(d_slave_idxs.size()));
    assert(d_num_springs == static_cast<int>(d_force_fcn_idxs.size()));
    assert(d_num_springs == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_springs == static_cast<int>(d_rest_lengths.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBSpringForceSpec::IBSpringForceSpec():\n"
                   << "  must call IBSpringForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBSpringForceSpec objects.\n");
    }
    return;
}// IBSpringForceSpec

IBSpringForceSpec::~IBSpringForceSpec()
{
    // intentionally blank
    return;
}// ~IBSpringForceSpec

unsigned
IBSpringForceSpec::getNumberOfSprings() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_springs == static_cast<int>(d_slave_idxs.size()));
    assert(d_num_springs == static_cast<int>(d_force_fcn_idxs.size()));
    assert(d_num_springs == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_springs == static_cast<int>(d_rest_lengths.size()));
#endif
    return d_num_springs;
}// getNumberOfSprings

const int&
IBSpringForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

int&
IBSpringForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

const std::vector<int>&
IBSpringForceSpec::getSlaveNodeIndices() const
{
    return d_slave_idxs;
}// getSlaveNodeIndices

std::vector<int>&
IBSpringForceSpec::getSlaveNodeIndices()
{
    return d_slave_idxs;
}// getSlaveNodeIndices

const std::vector<int>&
IBSpringForceSpec::getForceFunctionIndices() const
{
    return d_force_fcn_idxs;
}// getForceFunctionIndices

std::vector<int>&
IBSpringForceSpec::getForceFunctionIndices()
{
    return d_force_fcn_idxs;
}// getForceFunctionIndices

const std::vector<double>&
IBSpringForceSpec::getStiffnesses() const
{
    return d_stiffnesses;
}// getStiffnesses

std::vector<double>&
IBSpringForceSpec::getStiffnesses()
{
    return d_stiffnesses;
}// getStiffnesses

const std::vector<double>&
IBSpringForceSpec::getRestingLengths() const
{
    return d_rest_lengths;
}// getRestingLengths

std::vector<double>&
IBSpringForceSpec::getRestingLengths()
{
    return d_rest_lengths;
}// getRestingLengths

int
IBSpringForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

size_t
IBSpringForceSpec::getDataStreamSize() const
{
    return ((2+2*d_num_springs)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (  2*d_num_springs)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

void
IBSpringForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_springs == static_cast<int>(d_slave_idxs.size()));
    assert(d_num_springs == static_cast<int>(d_force_fcn_idxs.size()));
    assert(d_num_springs == static_cast<int>(d_stiffnesses.size()));
    assert(d_num_springs == static_cast<int>(d_rest_lengths.size()));
#endif
    stream.pack(&d_num_springs,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&d_slave_idxs[0],d_num_springs);
    stream.pack(&d_force_fcn_idxs[0],d_num_springs);
    stream.pack(&d_stiffnesses[0],d_num_springs);
    stream.pack(&d_rest_lengths[0],d_num_springs);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBSpringForceSpec>;

//////////////////////////////////////////////////////////////////////////////
