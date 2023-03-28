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

#ifndef included_IBAMR_IBRodForceSpec_inl
#define included_IBAMR_IBRodForceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBRodForceSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

#include <array>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBRodForceSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBRodForceSpec::IBRodForceSpec(const unsigned int num_rods) : d_next_idxs(num_rods), d_material_params(num_rods)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBRodForceSpec::IBRodForceSpec():\n"
                   << "  must call IBRodForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBRodForceSpec objects.\n");
    }
#endif
    return;
} // IBRodForceSpec

inline IBRodForceSpec::IBRodForceSpec(
    const int master_idx,
    const std::vector<int>& next_idxs,
    const std::vector<std::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& material_params)
    : d_master_idx(master_idx), d_next_idxs(next_idxs), d_material_params(material_params)
{
#if !defined(NDEBUG)
    const size_t num_rods = d_next_idxs.size();
    TBOX_ASSERT(num_rods == d_material_params.size());
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBRodForceSpec::IBRodForceSpec():\n"
                   << "  must call IBRodForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBRodForceSpec objects.\n");
    }
#endif
    return;
} // IBRodForceSpec

inline IBRodForceSpec::~IBRodForceSpec()
{
    // intentionally blank
    return;
} // ~IBRodForceSpec

inline unsigned int
IBRodForceSpec::getNumberOfRods() const
{
    const auto num_rods = static_cast<unsigned int>(d_next_idxs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(num_rods == d_material_params.size());
#endif
    return num_rods;
} // getNumberOfRods

inline const int&
IBRodForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
} // getMasterNodeIndex

inline int&
IBRodForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
} // getMasterNodeIndex

inline const std::vector<int>&
IBRodForceSpec::getNextNodeIndices() const
{
    return d_next_idxs;
} // getNextNodeIndices

inline std::vector<int>&
IBRodForceSpec::getNextNodeIndices()
{
    return d_next_idxs;
} // getNextNodeIndices

inline const std::vector<std::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >&
IBRodForceSpec::getMaterialParams() const
{
    return d_material_params;
} // getMaterialParams

inline std::vector<std::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >&
IBRodForceSpec::getMaterialParams()
{
    return d_material_params;
} // getMaterialParams

inline int
IBRodForceSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBRodForceSpec::getDataStreamSize() const
{
    const size_t num_rods = d_next_idxs.size();
#if !defined(NDEBUG)
    TBOX_ASSERT(num_rods == d_material_params.size());
#endif
    return ((2 + num_rods) * SAMRAI::tbox::AbstractStream::sizeofInt() +
            (NUM_MATERIAL_PARAMS * num_rods) * SAMRAI::tbox::AbstractStream::sizeofDouble());
} // getDataStreamSize

inline void
IBRodForceSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    const auto num_rods = static_cast<unsigned int>(d_next_idxs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(num_rods == d_material_params.size());
#endif
    stream << static_cast<int>(num_rods);
    stream.pack(&d_master_idx, 1);
    stream.pack(&d_next_idxs[0], num_rods);
    for (unsigned int n = 0; n < num_rods; ++n)
    {
        stream.pack(d_material_params[n].data(), NUM_MATERIAL_PARAMS);
    }
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBRodForceSpec_inl
