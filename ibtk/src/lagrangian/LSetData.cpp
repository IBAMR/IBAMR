// Filename: LSetData.cpp
// Created on 04 Jun 2007 by Boyce Griffith
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

#include "Box.h"
#include "IndexData.h"
#include "IndexDataFactory.h"
#include "IndexVariable.h"
#include "IntVector.h"
#include "ibtk/LSetData.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetData<T>::LSetData(const Box<NDIM>& box, const IntVector<NDIM>& ghosts)
    : IndexData<NDIM, LSet<T>, CellGeometry<NDIM> >(box, ghosts)
{
    // intentionally blank
    return;
} // LSetData

template <class T>
LSetData<T>::~LSetData()
{
    // intentionally blank
    return;
} // ~LSetData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include "ibtk/LMarkerSet.h"
template class SAMRAI::pdat::IndexData<NDIM, IBTK::LMarkerSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexDataFactory<NDIM, IBTK::LMarkerSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexDataNode<NDIM, IBTK::LMarkerSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexIterator<NDIM, IBTK::LMarkerSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexVariable<NDIM, IBTK::LMarkerSet, CellGeometry<NDIM> >;
template class IBTK::LSetData<IBTK::LMarker>;

#include "ibtk/LNodeSet.h"
template class SAMRAI::pdat::IndexData<NDIM, IBTK::LNodeSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexDataFactory<NDIM, IBTK::LNodeSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexDataNode<NDIM, IBTK::LNodeSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexIterator<NDIM, IBTK::LNodeSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexVariable<NDIM, IBTK::LNodeSet, CellGeometry<NDIM> >;
template class IBTK::LSetData<IBTK::LNode>;

#include "ibtk/LNodeIndexSet.h"
template class SAMRAI::pdat::IndexData<NDIM, IBTK::LNodeIndexSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexDataFactory<NDIM, IBTK::LNodeIndexSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexDataNode<NDIM, IBTK::LNodeIndexSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexIterator<NDIM, IBTK::LNodeIndexSet, CellGeometry<NDIM> >;
template class SAMRAI::pdat::IndexVariable<NDIM, IBTK::LNodeIndexSet, CellGeometry<NDIM> >;
template class IBTK::LSetData<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
