// Filename: VecCellDataFactory-inl.h
// Created on 09 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#ifndef included_VecCellDataFactory_inl_h
#define included_VecCellDataFactory_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "MultiblockCellDataTranslator.h"
#include "ibtk/VecCellDataFactory.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class TYPE>
inline int
VecCellDataFactory<TYPE>::getDefaultDepth() const
{
    return d_depth;
}// getDefaultDepth

template<class TYPE>
inline void
VecCellDataFactory<TYPE>::setDefaultDepth(
    const int depth)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth > 0);
#endif
    d_depth = depth;
    return;
}// setDefaultDepth

template<class TYPE>
inline SAMRAI::hier::MultiblockDataTranslator<NDIM>*
VecCellDataFactory<TYPE>::getMultiblockDataTranslator()
{
    if (!d_mb_trans)
    {
        d_mb_trans = new SAMRAI::pdat::MultiblockCellDataTranslator<NDIM,TYPE>();
    }
    return d_mb_trans;
}// getMultiblockDataTranslator

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VecCellDataFactory_inl_h
