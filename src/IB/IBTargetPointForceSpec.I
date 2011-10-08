// Filename: IBTargetPointForceSpec.C
// Created on 11 Jun 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBTargetPointForceSpec::getIsRegisteredWithStreamableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
IBTargetPointForceSpec::IBTargetPointForceSpec(
    const int master_idx,
    const double kappa_target,
    const double eta_target,
    const blitz::TinyVector<double,NDIM>& X_target
#if ENABLE_SUBDOMAIN_INDICES
    ,const int subdomain_idx
#endif
                                               )
    : d_master_idx(master_idx),
      d_kappa_target(kappa_target),
      d_eta_target(eta_target),
      d_X_target(X_target)
#if ENABLE_SUBDOMAIN_INDICES
    , d_subdomain_idx(subdomain_idx)
#endif
{
#ifdef DEBUG_CHECK_ASSERTIONS
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBTargetPointForceSpec::IBTargetPointForceSpec():\n"
                   << "  must call IBTargetPointForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBTargetPointForceSpec objects.\n");
    }
#endif
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

inline const blitz::TinyVector<double,NDIM>&
IBTargetPointForceSpec::getTargetPointPosition() const
{
    return d_X_target;
}// getTargetPointPosition

inline blitz::TinyVector<double,NDIM>&
IBTargetPointForceSpec::getTargetPointPosition()
{
    return d_X_target;
}// getTargetPointPosition

#if ENABLE_SUBDOMAIN_INDICES
inline const int&
IBTargetPointForceSpec::getSubdomainIndex() const
{
    return d_subdomain_idx;
}// getSubdomainIndex

inline int&
IBTargetPointForceSpec::getSubdomainIndex()
{
    return d_subdomain_idx;
}// getSubdomainIndex
#endif

inline int
IBTargetPointForceSpec::getStreamableClassID() const
{
    return s_class_id;
}// getStreamableClassID

inline size_t
IBTargetPointForceSpec::getDataStreamSize() const
{
    return ((1     )*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (2+NDIM)*SAMRAI::tbox::AbstractStream::sizeofDouble())
#if ENABLE_SUBDOMAIN_INDICES
        + SAMRAI::tbox::AbstractStream::sizeofInt()
#endif
        ;
}// getDataStreamSize

inline void
IBTargetPointForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    stream.pack(&d_master_idx,1);
    stream.pack(&d_kappa_target,1);
    stream.pack(&d_eta_target,1);
    stream.pack(d_X_target.data(),NDIM);
#if ENABLE_SUBDOMAIN_INDICES
    stream.pack(&d_subdomain_idx,1);
#endif
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
