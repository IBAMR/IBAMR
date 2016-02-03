// Filename: LIndexSetDataFactory.cpp
// Created on 13 May 2011 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LIndexSetDataFactory.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/ArenaManager.h"

namespace IBTK
{
class LNode;
class LNodeIndex;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LIndexSetDataFactory<T>::LIndexSetDataFactory(const IntVector<NDIM>& ghosts)
    : LSetDataFactory<T>(ghosts)
{
    // intentionally blank
    return;
} // LIndexSetDataFactory

template <class T>
LIndexSetDataFactory<T>::~LIndexSetDataFactory()
{
    // intentionally blank
    return;
} // ~LIndexSetDataFactory

template <class T>
Pointer<PatchDataFactory<NDIM> >
LIndexSetDataFactory<T>::cloneFactory(const IntVector<NDIM>& ghosts)
{
    return new LIndexSetDataFactory<T>(ghosts);
} // cloneFactory

template <class T>
Pointer<PatchData<NDIM> >
LIndexSetDataFactory<T>::allocate(const Box<NDIM>& box, Pointer<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    PatchData<NDIM>* pd = new (pool) LIndexSetData<T>(box, LSetDataFactory<T>::getGhostCellWidth());
    return Pointer<PatchData<NDIM> >(pd, pool);
} // allocate

template <class T>
Pointer<PatchData<NDIM> >
LIndexSetDataFactory<T>::allocate(const Patch<NDIM>& patch, Pointer<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LIndexSetDataFactory<T>::getSizeOfMemory(const Box<NDIM>& /*box*/) const
{
    return Arena::align(sizeof(LIndexSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LIndexSetDataFactory<T>::validCopyTo(const Pointer<PatchDataFactory<NDIM> >& dst_pdf) const
{
    const Pointer<LIndexSetDataFactory<T> > lnidf = dst_pdf;
    return lnidf;
} // validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LIndexSetDataFactory<IBTK::LNode>;
template class IBTK::LIndexSetDataFactory<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
