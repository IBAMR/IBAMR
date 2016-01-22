// Filename: IBRodForceSpec-inl.h
// Created on 23 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBRodForceSpec_inl
#define included_IBRodForceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "boost/array.hpp"
#include "ibamr/IBRodForceSpec.h"
#include "ibtk/StreamableManager.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

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

inline IBRodForceSpec::IBRodForceSpec(const unsigned int num_rods)
    : d_master_idx(-1), d_next_idxs(num_rods), d_material_params(num_rods)
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
    const std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >& material_params)
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
    const unsigned int num_rods = static_cast<unsigned int>(d_next_idxs.size());
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

inline const std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >&
IBRodForceSpec::getMaterialParams() const
{
    return d_material_params;
} // getMaterialParams

inline std::vector<boost::array<double, IBRodForceSpec::NUM_MATERIAL_PARAMS> >&
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
    const unsigned int num_rods = static_cast<unsigned int>(d_next_idxs.size());
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

#endif //#ifndef included_IBRodForceSpec_inl
