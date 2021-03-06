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

#ifndef included_IBAMR_IBTargetPointForceSpec_inl
#define included_IBAMR_IBTargetPointForceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBTargetPointForceSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBTargetPointForceSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBTargetPointForceSpec::IBTargetPointForceSpec(const int master_idx,
                                                      const double kappa_target,
                                                      const double eta_target,
                                                      const IBTK::Point& X_target)
    : d_master_idx(master_idx), d_kappa_target(kappa_target), d_eta_target(eta_target), d_X_target(X_target)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBTargetPointForceSpec::IBTargetPointForceSpec():\n"
                   << "  must call IBTargetPointForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBTargetPointForceSpec objects.\n");
    }
#endif
    return;
} // IBTargetPointForceSpec

inline const int&
IBTargetPointForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
} // getMasterNodeIndex

inline int&
IBTargetPointForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
} // getMasterNodeIndex

inline const double&
IBTargetPointForceSpec::getStiffness() const
{
    return d_kappa_target;
} // getStiffness

inline double&
IBTargetPointForceSpec::getStiffness()
{
    return d_kappa_target;
} // getStiffness

inline const double&
IBTargetPointForceSpec::getDamping() const
{
    return d_eta_target;
} // getDamping

inline double&
IBTargetPointForceSpec::getDamping()
{
    return d_eta_target;
} // getDamping

inline const IBTK::Point&
IBTargetPointForceSpec::getTargetPointPosition() const
{
    return d_X_target;
} // getTargetPointPosition

inline IBTK::Point&
IBTargetPointForceSpec::getTargetPointPosition()
{
    return d_X_target;
} // getTargetPointPosition

inline int
IBTargetPointForceSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBTargetPointForceSpec::getDataStreamSize() const
{
    return ((1) * SAMRAI::tbox::AbstractStream::sizeofInt() +
            (2 + NDIM) * SAMRAI::tbox::AbstractStream::sizeofDouble());
} // getDataStreamSize

inline void
IBTargetPointForceSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx, 1);
    stream.pack(&d_kappa_target, 1);
    stream.pack(&d_eta_target, 1);
    stream.pack(d_X_target.data(), NDIM);
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBTargetPointForceSpec_inl
