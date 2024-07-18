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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/LIndexSetData.h"
#include "ibtk/LIndexSetDataFactory.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSetDataFactory.h"

#include "Box.h"
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
LIndexSetDataFactory<T>::LIndexSetDataFactory(IntVectorNd ghosts) : LSetDataFactory<T>(std::move(ghosts))
{
    // intentionally blank
    return;
} // LIndexSetDataFactory

template <class T>
SAMRAIPointer<PatchDataFactoryNd>
LIndexSetDataFactory<T>::cloneFactory(const IntVectorNd& ghosts)
{
    return new LIndexSetDataFactory<T>(ghosts);
} // cloneFactory

template <class T>
SAMRAIPointer<PatchDataNd>
LIndexSetDataFactory<T>::allocate(const BoxNd& box, SAMRAIPointer<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    PatchDataNd* pd = new (pool) LIndexSetData<T>(box, LSetDataFactory<T>::getGhostCellWidth());
    return SAMRAIPointer<PatchDataNd>(pd, pool);
} // allocate

template <class T>
SAMRAIPointer<PatchDataNd>
LIndexSetDataFactory<T>::allocate(const PatchNd& patch, SAMRAIPointer<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LIndexSetDataFactory<T>::getSizeOfMemory(const BoxNd& /*box*/) const
{
    return Arena::align(sizeof(LIndexSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LIndexSetDataFactory<T>::validCopyTo(const SAMRAIPointer<PatchDataFactoryNd>& dst_pdf) const
{
    const SAMRAIPointer<LIndexSetDataFactory<T> > lnidf = dst_pdf;
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
