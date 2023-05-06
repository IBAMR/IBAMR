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

#include "ibtk/LNode.h"
#include "ibtk/LNodeIndex.h"
#include "ibtk/LSetData.h"

#include "Box.h"
#include "CellGeometry.h"
#include "IndexData.h"
#include "IndexDataFactory.h"
#include "IndexVariable.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

#include <algorithm>
#include <ostream>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetData<T>::LSetData(Box<NDIM> box, IntVector<NDIM> ghosts)
    : IndexData<NDIM, LSet<T>, CellGeometry<NDIM> >(std::move(box), std::move(ghosts))
{
    // intentionally blank
    return;
} // LSetData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////


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
