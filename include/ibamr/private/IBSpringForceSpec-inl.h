// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IBSpringForceSpec_inl
#define included_IBAMR_IBSpringForceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBSpringForceSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBSpringForceSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBSpringForceSpec::IBSpringForceSpec(const unsigned int num_springs)
    : d_slave_idxs(num_springs), d_force_fcn_idxs(num_springs), d_parameters(num_springs)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBSpringForceSpec::IBSpringForceSpec():\n"
                   << "  must call IBSpringForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBSpringForceSpec objects.\n");
    }
#endif
    return;
} // IBSpringForceSpec

inline IBSpringForceSpec::IBSpringForceSpec(const int master_idx,
                                            const std::vector<int>& slave_idxs,
                                            const std::vector<int>& force_fcn_idxs,
                                            const std::vector<std::vector<double> >& parameters)
    : d_master_idx(master_idx), d_slave_idxs(slave_idxs), d_force_fcn_idxs(force_fcn_idxs), d_parameters(parameters)
{
#if !defined(NDEBUG)
    const size_t num_springs = d_slave_idxs.size();
    TBOX_ASSERT(num_springs == d_force_fcn_idxs.size());
    TBOX_ASSERT(num_springs == d_parameters.size());
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBSpringForceSpec::IBSpringForceSpec():\n"
                   << "  must call IBSpringForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBSpringForceSpec objects.\n");
    }
#endif
    return;
} // IBSpringForceSpec

inline unsigned int
IBSpringForceSpec::getNumberOfSprings() const
{
    const auto num_springs = static_cast<unsigned int>(d_slave_idxs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(num_springs == d_force_fcn_idxs.size());
    TBOX_ASSERT(num_springs == d_parameters.size());
#endif
    return num_springs;
} // getNumberOfSprings

inline const int&
IBSpringForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
} // getMasterNodeIndex

inline int&
IBSpringForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
} // getMasterNodeIndex

inline const std::vector<int>&
IBSpringForceSpec::getSlaveNodeIndices() const
{
    return d_slave_idxs;
} // getSlaveNodeIndices

inline std::vector<int>&
IBSpringForceSpec::getSlaveNodeIndices()
{
    return d_slave_idxs;
} // getSlaveNodeIndices

inline const std::vector<int>&
IBSpringForceSpec::getForceFunctionIndices() const
{
    return d_force_fcn_idxs;
} // getForceFunctionIndices

inline std::vector<int>&
IBSpringForceSpec::getForceFunctionIndices()
{
    return d_force_fcn_idxs;
} // getForceFunctionIndices

inline const std::vector<std::vector<double> >&
IBSpringForceSpec::getParameters() const
{
    return d_parameters;
} // getParameters

inline std::vector<std::vector<double> >&
IBSpringForceSpec::getParameters()
{
    return d_parameters;
} // getParameters

inline int
IBSpringForceSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBSpringForceSpec::getDataStreamSize() const
{
    const size_t num_springs = d_slave_idxs.size();
#if !defined(NDEBUG)
    TBOX_ASSERT(num_springs == d_force_fcn_idxs.size());
    TBOX_ASSERT(num_springs == d_parameters.size());
#endif
    size_t size = (2 + 2 * num_springs) * SAMRAI::tbox::AbstractStream::sizeofInt();
    for (unsigned int k = 0; k < num_springs; ++k)
    {
        size += SAMRAI::tbox::AbstractStream::sizeofInt() +
                d_parameters[k].size() * SAMRAI::tbox::AbstractStream::sizeofDouble();
    }
    return size;
} // getDataStreamSize

inline void
IBSpringForceSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    const auto num_springs = static_cast<unsigned int>(d_slave_idxs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(num_springs == d_force_fcn_idxs.size());
    TBOX_ASSERT(num_springs == d_parameters.size());
#endif
    stream << static_cast<int>(num_springs);
    stream.pack(&d_master_idx, 1);
    stream.pack(&d_slave_idxs[0], num_springs);
    stream.pack(&d_force_fcn_idxs[0], num_springs);
    for (unsigned int k = 0; k < num_springs; ++k)
    {
        const int num_parameters = static_cast<int>(d_parameters[k].size());
        stream << num_parameters;
        stream.pack(&d_parameters[k][0], num_parameters);
    }
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBSpringForceSpec_inl
