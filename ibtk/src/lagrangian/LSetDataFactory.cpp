// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSetDataFactory.h"

#include "Box.h"
#include "CellGeometry.h"
#include "IndexDataFactory.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "tbox/Arena.h"
#include "tbox/ArenaManager.h"
#include "tbox/Pointer.h"

#include <algorithm>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetDataFactory<T>::LSetDataFactory(IntVectorNd ghosts) : IndexDataFactoryNd<LSet<T>, CellGeometryNd>(std::move(ghosts))
{
    // intentionally blank
    return;
} // LSetDataFactory

template <class T>
Pointer<PatchDataFactoryNd>
LSetDataFactory<T>::cloneFactory(const IntVectorNd& ghosts)
{
    return new LSetDataFactory<T>(ghosts);
} // cloneFactory

template <class T>
Pointer<PatchDataNd>
LSetDataFactory<T>::allocate(const BoxNd& box, Pointer<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    PatchDataNd* pd = new (pool) LSetData<T>(box, IndexDataFactoryNd<LSet<T>, CellGeometryNd>::getGhostCellWidth());
    return Pointer<PatchDataNd>(pd, pool);
} // allocate

template <class T>
Pointer<PatchDataNd>
LSetDataFactory<T>::allocate(const PatchNd& patch, Pointer<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LSetDataFactory<T>::getSizeOfMemory(const BoxNd& /*box*/) const
{
    return Arena::align(sizeof(LSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LSetDataFactory<T>::validCopyTo(const Pointer<PatchDataFactoryNd>& dst_pdf) const
{
    Pointer<LSetDataFactory<T> > lnidf = dst_pdf;
    return lnidf;
} // validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LSetDataFactory<IBTK::LNode>;
template class IBTK::LSetDataFactory<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
