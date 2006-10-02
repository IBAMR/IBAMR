//
// LNodeIndexData.C
//
// Created on 01 Mar 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <07.Mar.2005 13:38:39 boyce@trasnaform.cims.nyu.edu>
//

#include "LNodeIndexData.h"

#include "LNodeIndexVariable.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

/////////////////////////////// INLINE ///////////////////////////////////////

#ifdef DEBUG_NO_INLINE
#include "LNodeIndexData.I"
#endif

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexData::LNodeIndexData(
    const hier::Box<NDIM>& box,
    const hier::IntVector<NDIM>& ghosts)
        : pdat::IndexData<NDIM,LNodeIndexSet>(box,ghosts)
{
    // intentionally blank
    return;
}// LNodeIndexData

LNodeIndexData::~LNodeIndexData()
{
    // intentionally blank
    return;
}// ~LNodeIndexData

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the LNodeIndexData class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<LNodeIndexData>;

#endif

//////////////////////////////////////////////////////////////////////////////
