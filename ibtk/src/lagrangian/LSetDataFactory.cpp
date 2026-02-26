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

#include <ibtk/LNode.h>
#include <ibtk/LNodeIndex.h>
#include <ibtk/LSet.h>
#include <ibtk/LSetData.h>
#include <ibtk/LSetDataFactory.h>
#include <ibtk/samrai_compatibility_names.h>

#include <tbox/Arena.h>
#include <tbox/ArenaManager.h>

#include <SAMRAIBox.h>
#include <SAMRAICellGeometry.h>
#include <SAMRAIIndexDataFactory.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchData.h>
#include <SAMRAIPatchDataFactory.h>
#include <SAMRAIPointer.h>

#include <algorithm>
#include <utility>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetDataFactory<T>::LSetDataFactory(SAMRAIIntVector ghosts)
    : SAMRAIIndexDataFactory<LSet<T>, SAMRAICellGeometry>(std::move(ghosts))
{
    // intentionally blank
    return;
} // LSetDataFactory

template <class T>
SAMRAIPointer<SAMRAIPatchDataFactory>
LSetDataFactory<T>::cloneFactory(const SAMRAIIntVector& ghosts)
{
    return new LSetDataFactory<T>(ghosts);
} // cloneFactory

template <class T>
SAMRAIPointer<SAMRAIPatchData>
LSetDataFactory<T>::allocate(const SAMRAIBox& box, SAMRAIPointer<Arena> pool) const
{
    if (!pool)
    {
        pool = ArenaManager::getManager()->getStandardAllocator();
    }
    SAMRAIPatchData* pd =
        new (pool) LSetData<T>(box, SAMRAIIndexDataFactory<LSet<T>, SAMRAICellGeometry>::getGhostCellWidth());
    return SAMRAIPointer<SAMRAIPatchData>(pd, pool);
} // allocate

template <class T>
SAMRAIPointer<SAMRAIPatchData>
LSetDataFactory<T>::allocate(const SAMRAIPatch& patch, SAMRAIPointer<Arena> pool) const
{
    return allocate(patch.getBox(), pool);
} // allocate

template <class T>
size_t
LSetDataFactory<T>::getSizeOfMemory(const SAMRAIBox& /*box*/) const
{
    return Arena::align(sizeof(LSetData<T>));
} // getSizeOfMemory

template <class T>
bool
LSetDataFactory<T>::validCopyTo(const SAMRAIPointer<SAMRAIPatchDataFactory>& dst_pdf) const
{
    SAMRAIPointer<LSetDataFactory<T>> lnidf = dst_pdf;
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
