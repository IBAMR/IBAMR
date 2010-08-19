// Filename: IBRodForceSpec.I
// Created on 23 Jun 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBRodForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBRodForceSpec::IBRodForceSpec(
    const int master_idx,
    const std::vector<int>& next_idxs,
    const std::vector<std::vector<double> >& material_params)
    : d_master_idx(master_idx),
      d_next_idxs(next_idxs),
      d_material_params(material_params)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    const int num_rods = d_next_idxs.size();
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBRodForceSpec::IBRodForceSpec():\n"
                   << "  must call IBRodForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBRodForceSpec objects.\n");
    }
    return;
}// IBRodForceSpec

inline
IBRodForceSpec::~IBRodForceSpec()
{
    // intentionally blank
    return;
}// ~IBRodForceSpec

inline unsigned
IBRodForceSpec::getNumberOfRods() const
{
    const int num_rods = d_next_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    return num_rods;
}// getNumberOfRods

inline const int&
IBRodForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBRodForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const std::vector<int>&
IBRodForceSpec::getNextNodeIndices() const
{
    return d_next_idxs;
}// getNextNodeIndices

inline std::vector<int>&
IBRodForceSpec::getNextNodeIndices()
{
    return d_next_idxs;
}// getNextNodeIndices

inline const std::vector<std::vector<double> >&
IBRodForceSpec::getMaterialParams() const
{
    return d_material_params;
}// getMaterialParams

inline std::vector<std::vector<double> >&
IBRodForceSpec::getMaterialParams()
{
    return d_material_params;
}// getMaterialParams

inline int
IBRodForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBRodForceSpec::getDataStreamSize() const
{
    const int num_rods = d_next_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    return (( 2+num_rods)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (10*num_rods)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBRodForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    const int num_rods = d_next_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_rods == int(d_material_params.size()));
#endif
    stream.pack(&num_rods,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&d_next_idxs[0],num_rods);
    for (int n = 0; n < num_rods; ++n)
    {
        stream.pack(&d_material_params[n][0],10);
    }
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
