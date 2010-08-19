// Filename: IBMovingTargetPointForceSpec.C
// Created on 14 Aug 2008 by Boyce Griffith
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
IBMovingTargetPointForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBMovingTargetPointForceSpec::IBMovingTargetPointForceSpec(
    const int master_idx,
    const double& kappa_target,
    const double& eta_target,
    const int spec_fcn_idx,
    const std::vector<double>& periodic_shift)
    : d_master_idx(master_idx),
      d_kappa_target(kappa_target),
      d_eta_target(eta_target),
      d_spec_fcn_idx(spec_fcn_idx),
      d_periodic_shift(periodic_shift)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == int(d_periodic_shift.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBMovingTargetPointForceSpec::IBMovingTargetPointForceSpec():\n"
                   << "  must call IBMovingTargetPointForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBMovingTargetPointForceSpec objects.\n");
    }
    return;
}// IBMovingTargetPointForceSpec

inline
IBMovingTargetPointForceSpec::~IBMovingTargetPointForceSpec()
{
    // intentionally blank
    return;
}// ~IBMovingTargetPointForceSpec

inline const int&
IBMovingTargetPointForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBMovingTargetPointForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const double&
IBMovingTargetPointForceSpec::getStiffness() const
{
    return d_kappa_target;
}// getStiffness

inline double&
IBMovingTargetPointForceSpec::getStiffness()
{
    return d_kappa_target;
}// getStiffness

inline const double&
IBMovingTargetPointForceSpec::getDamping() const
{
    return d_eta_target;
}// getDamping

inline double&
IBMovingTargetPointForceSpec::getDamping()
{
    return d_eta_target;
}// getDamping

inline const int&
IBMovingTargetPointForceSpec::getPositionAndVelocityFunctionIndex() const
{
    return d_spec_fcn_idx;
}// getPositionAndVelocityFunctionIndex

inline int&
IBMovingTargetPointForceSpec::getPositionAndVelocityFunctionIndex()
{
    return d_spec_fcn_idx;
}// getPositionAndVelocityFunctionIndex

inline const std::vector<double>&
IBMovingTargetPointForceSpec::getPeriodicShift() const
{
    return d_periodic_shift;
}// getPeriodicShift

inline std::vector<double>&
IBMovingTargetPointForceSpec::getPeriodicShift()
{
    return d_periodic_shift;
}// getPeriodicShift

inline int
IBMovingTargetPointForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBMovingTargetPointForceSpec::getDataStreamSize() const
{
    return ((2     )*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (2+NDIM)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBMovingTargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(NDIM == int(d_periodic_shift.size()));
#endif
    stream.pack(&d_master_idx,1);
    stream.pack(&d_kappa_target,1);
    stream.pack(&d_eta_target,1);
    stream.pack(&d_spec_fcn_idx,1);
    stream.pack(&d_periodic_shift[0],NDIM);
    return;
}// packStream

inline void
IBMovingTargetPointForceSpec::registerPeriodicShift(
    const SAMRAI::hier::IntVector<NDIM>& offset,
    const std::vector<double>& displacement)
{
    for (int d = 0; d < NDIM; ++d)
    {
        d_periodic_shift[d] += displacement[d];
    }
    return;
}// registerPeriodicShift

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
