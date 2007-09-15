// Filename: IBMarker.C
// Last modified: <12.Sep.2007 23:20:12 griffith@box221.cims.nyu.edu>
// Created on 12 Sep 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBMarker.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <IndexData.C>
#include <IndexDataFactory.C>
#include <IndexVariable.C>
#include <tbox/Array.C>
#include <tbox/List.C>
#include <tbox/Pointer.C>

template class SAMRAI::tbox::Pointer<IBAMR::IBMarker>;

#if (NDIM == 1)
template class SAMRAI::pdat::IndexData<1,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexDataFactory<1,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexDataNode<1,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexIterator<1,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexVariable<1,IBAMR::IBMarker>;
template class SAMRAI::tbox::Array<IBAMR::IBMarker>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<1,IBAMR::IBMarker> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<1,IBAMR::IBMarker> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<1,IBAMR::IBMarker> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<1,IBAMR::IBMarker> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<1,IBAMR::IBMarker> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<1,IBAMR::IBMarker> >;
#endif

#if (NDIM == 2)
template class SAMRAI::pdat::IndexData<2,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexDataFactory<2,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexDataNode<2,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexIterator<2,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexVariable<2,IBAMR::IBMarker>;
template class SAMRAI::tbox::Array<IBAMR::IBMarker>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<2,IBAMR::IBMarker> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<2,IBAMR::IBMarker> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<2,IBAMR::IBMarker> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<2,IBAMR::IBMarker> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<2,IBAMR::IBMarker> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<2,IBAMR::IBMarker> >;
#endif

#if (NDIM == 3)
template class SAMRAI::pdat::IndexData<3,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexDataFactory<3,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexDataNode<3,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexIterator<3,IBAMR::IBMarker>;
template class SAMRAI::pdat::IndexVariable<3,IBAMR::IBMarker>;
template class SAMRAI::tbox::Array<IBAMR::IBMarker>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<3,IBAMR::IBMarker> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<3,IBAMR::IBMarker> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<3,IBAMR::IBMarker> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<3,IBAMR::IBMarker> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<3,IBAMR::IBMarker> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<3,IBAMR::IBMarker> >;
#endif

//////////////////////////////////////////////////////////////////////////////
