// Filename: IBTargetPointForceSpec.C
// Created on 11 Jun 2007 by Boyce Griffith
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
IBTargetPointForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBTargetPointForceSpec::IBTargetPointForceSpec(
    const int master_idx,
    const double& kappa_target,
    const double& eta_target,
    const std::vector<double>& X_target)
    : d_master_idx(master_idx),
      d_kappa_target(kappa_target),
      d_eta_target(eta_target),
      d_X_target(X_target)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == int(d_X_target.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBTargetPointForceSpec::IBTargetPointForceSpec():\n"
                   << "  must call IBTargetPointForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBTargetPointForceSpec objects.\n");
    }
    return;
}// IBTargetPointForceSpec

inline
IBTargetPointForceSpec::~IBTargetPointForceSpec()
{
    // intentionally blank
    return;
}// ~IBTargetPointForceSpec

inline const int&
IBTargetPointForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBTargetPointForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const double&
IBTargetPointForceSpec::getStiffness() const
{
    return d_kappa_target;
}// getStiffness

inline double&
IBTargetPointForceSpec::getStiffness()
{
    return d_kappa_target;
}// getStiffness

inline const double&
IBTargetPointForceSpec::getDamping() const
{
    return d_eta_target;
}// getDamping

inline double&
IBTargetPointForceSpec::getDamping()
{
    return d_eta_target;
}// getDamping

inline const std::vector<double>&
IBTargetPointForceSpec::getTargetPointPosition() const
{
    return d_X_target;
}// getTargetPointPosition

inline std::vector<double>&
IBTargetPointForceSpec::getTargetPointPosition()
{
    return d_X_target;
}// getTargetPointPosition

inline int
IBTargetPointForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBTargetPointForceSpec::getDataStreamSize() const
{
    return ((1     )*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (2+NDIM)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBTargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == int(d_X_target.size()));
#endif
    stream.pack(&d_master_idx,1);
    stream.pack(&d_kappa_target,1);
    stream.pack(&d_eta_target,1);
    stream.pack(&d_X_target[0],NDIM);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
