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
#include "ibtk/LSetDataFactory.h"
#include "ibtk/LSetVariable.h"

#include "Variable.h"

#include <string>
#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class T>
LSetVariable<T>::LSetVariable(std::string name)
    : Variable<NDIM>(std::move(name), new LSetDataFactory<T>(IntVector<NDIM>(0)))
{
    // intentionally blank
    return;
} // LSetVariable

template <class T>
bool
LSetVariable<T>::dataLivesOnPatchBorder() const
{
    return false;
} // dataLivesOnPatchBorder

template <class T>
bool
LSetVariable<T>::fineBoundaryRepresentsVariable() const
{
    return true;
} // fineBoundaryRepresentsVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

template class IBTK::LSetVariable<IBTK::LNode>;
template class IBTK::LSetVariable<IBTK::LNodeIndex>;

//////////////////////////////////////////////////////////////////////////////
