// Filename: LNode.I
// Created on 05 May 2011 by Boyce Griffith
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

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/StreamableManager.h>
#include <ibtk/ibtk_utilities.h>

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
LNode::LNode(
    const int lagrangian_nidx,
    const int global_petsc_nidx,
    const int local_petsc_nidx,
    const double* const X_ptr,
    const SAMRAI::hier::IntVector<NDIM>& periodic_offset,
    const std::vector<double>& periodic_displacement,
    const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data)
    : LNodeIndex(lagrangian_nidx, global_petsc_nidx, local_petsc_nidx,
                 X_ptr, periodic_offset, periodic_displacement, node_data)
{
    // intentionally blank
    return;
}// LNode

inline
LNode::LNode(
    const LNode& from)
    : LNodeIndex(from)
{
    // intentionally blank
    return;
}// LNode

inline
LNode::~LNode()
{
    // intentionally blank
    return;
}// ~LNode

inline LNode&
LNode::operator=(
    const LNode& that)
{
    if (this != &that)
    {
        LNodeIndex::operator=(that);
    }
    return *this;
}// operator=

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

inline bool
operator<(
    const LNode& lhs,
    const LNode& rhs)
{
    return lhs.getGlobalPETScIndex() < rhs.getGlobalPETScIndex();
}// operator<

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
