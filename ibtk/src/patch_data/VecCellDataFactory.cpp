// Filename: VecCellDataFactory.C
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>

#include "BoxGeometry.h"
#include "CellGeometry.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchData.h"
#include "SAMRAI_config.h"
#include "VecCellDataFactory.h"
#include "ibtk/VecCellData.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/ArenaManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template<class TYPE>
VecCellDataFactory<TYPE>::VecCellDataFactory(
    int depth,
    const IntVector<NDIM>& ghosts)
    : PatchDataFactory<NDIM>(ghosts),
      d_depth(depth),
      d_mb_trans(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(depth > 0);
    TBOX_ASSERT(ghosts.min() >= 0);
#endif
    return;
}// VecCellDataFactory

template<class TYPE>
VecCellDataFactory<TYPE>::~VecCellDataFactory()
{
    if (d_mb_trans)
    {
        delete d_mb_trans;
    }
    return;
}// ~VecCellDataFactory

template<class TYPE>
Pointer<PatchDataFactory<NDIM> >
VecCellDataFactory<TYPE>::cloneFactory(
    const IntVector<NDIM>& ghosts)
{
    return new VecCellDataFactory<TYPE>(d_depth, ghosts);
}// cloneFactory

template<class TYPE>
Pointer<PatchData<NDIM> >
VecCellDataFactory<TYPE>::allocate(
    const Box<NDIM>& box,
    Pointer<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }

    PatchData<NDIM>* patchdata = new (pool) VecCellData<TYPE>(box, d_depth, d_ghosts);
    return Pointer<PatchData<NDIM> >(patchdata, pool);
}// allocate

template<class TYPE>
Pointer<PatchData<NDIM> >
VecCellDataFactory<TYPE>::allocate(
    const Patch<NDIM>& patch,
    Pointer<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
}// allocate

template<class TYPE>
Pointer<BoxGeometry<NDIM> >
VecCellDataFactory<TYPE>::getBoxGeometry(
    const Box<NDIM>& box) const
{
    BoxGeometry<NDIM>* boxgeometry = new CellGeometry<NDIM>(box, d_ghosts);
    return Pointer<BoxGeometry<NDIM> >(boxgeometry);
}

template<class TYPE>
size_t
VecCellDataFactory<TYPE>::getSizeOfMemory(
    const Box<NDIM>& box) const
{
    const size_t obj = Arena::align(sizeof(VecCellData<TYPE>));
    const size_t data = VecCellData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);
    return obj+data;
}// getSizeOfMemory

template<class TYPE>
bool
VecCellDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
    return true;
}// fineBoundaryRepresentsVariable

template<class TYPE>
bool
VecCellDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
    return false;
}// dataLivesOnPatchBorder

template<class TYPE>
bool
VecCellDataFactory<TYPE>::validCopyTo(
    const Pointer<PatchDataFactory<NDIM> >& dst_pdf) const
{
    /*
     * The only valid option is VecCellDataFactory.
     */
    Pointer<VecCellDataFactory<TYPE> > cdf = dst_pdf;
    if (cdf)
    {
        return true;
    }
    else
    {
        return false;
    }
}// validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::VecCellDataFactory<double>;
template class Pointer<IBTK::VecCellDataFactory<double> >;

//////////////////////////////////////////////////////////////////////////////
