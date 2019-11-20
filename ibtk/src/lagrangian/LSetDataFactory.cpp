// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LSetData.h"
#include "ibtk/LSetDataFactory.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "Patch.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "tbox/ArenaManager.h"

namespace IBTK
{
template <class T>
class LSet;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetDataFactory<T>::LSetDataFactory(IntVector<NDIM> ghosts)
    : IndexDataFactory<NDIM, LSet<T>, CellGeometry<NDIM> >(std::move(ghosts))
{
    // intentionally blank
    return;
} // LSetDataFactory

template <class T>
Pointer<PatchDataFactory<NDIM> >
LSetDataFactory<T>::cloneFactory(const IntVector<NDIM>& ghosts)
{
    return new LSetDataFactory<T>(ghosts);
} // cloneFactory

template <class T>
Pointer<PatchData<NDIM> >
LSetDataFactory<T>::allocate(const Box<NDIM>& box, Pointer<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    PatchData<NDIM>* pd =
        new (pool) LSetData<T>(box, IndexDataFactory<NDIM, LSet<T>, CellGeometry<NDIM> >::getGhostCellWidth());
    return Pointer<PatchData<NDIM> >(pd, pool);
} // allocate

template <class T>
Pointer<PatchData<NDIM> >
LSetDataFactory<T>::allocate(const Patch<NDIM>& patch, Pointer<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LSetDataFactory<T>::getSizeOfMemory(const Box<NDIM>& /*box*/) const
{
    return Arena::align(sizeof(LSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LSetDataFactory<T>::validCopyTo(const Pointer<PatchDataFactory<NDIM> >& dst_pdf) const
{
    Pointer<LSetDataFactory<T> > lnidf = dst_pdf;
    return lnidf;
} // validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LSetDataFactory<IBTK::LMarker>;
template class IBTK::LSetDataFactory<IBTK::LNode>;
template class IBTK::LSetDataFactory<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
