// Filename: LNodeIndexSet.C
// Last modified: <04.Jun.2007 16:27:34 griffith@box221.cims.nyu.edu>
// Created on 29 Feb 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "LNodeIndexSet.h"

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

#include <ArrayData.C>
#include <ibamr/ArrayData_specialized_LNodeIndexSet.I>
#include <CellData.C>
#include <CellDataFactory.C>
#include <CellVariable.C>
#include <IndexData.C>
#include <IndexDataFactory.C>
#include <IndexVariable.C>
#include <tbox/Array.C>
#include <tbox/List.C>
#include <tbox/Pointer.C>

template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexSet>;

#if (NDIM == 1)
template class SAMRAI::pdat::ArrayData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellDataFactory<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellVariable<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexData<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<1,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<1,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<1,IBAMR::LNodeIndexSet> >;
#endif

#if (NDIM == 2)
template class SAMRAI::pdat::ArrayData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellDataFactory<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellVariable<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexData<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<2,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<2,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<2,IBAMR::LNodeIndexSet> >;
#endif

#if (NDIM == 3)
template class SAMRAI::pdat::ArrayData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellDataFactory<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::CellVariable<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexData<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataFactory<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexIterator<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::pdat::IndexVariable<3,IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<IBAMR::LNodeIndexSet>;
template class SAMRAI::tbox::Array<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::List<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListIterator<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::ListNode<SAMRAI::pdat::IndexDataNode<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexData<3,IBAMR::LNodeIndexSet> >;
template class SAMRAI::tbox::Pointer<SAMRAI::pdat::IndexVariable<3,IBAMR::LNodeIndexSet> >;
#endif

//////////////////////////////////////////////////////////////////////////////
