// Filename: LNodeIndexData.C
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <24.Oct.2006 14:16:56 boyce@bigboy.nyconnect.com>

#include "LNodeIndexData.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#define included_IBAMR_config
#include <IBAMR_config.h>
#endif

#ifndef included_SAMRAI_config
#define included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

// IBAMR INCLUDES
#include <ibamr/LNodeIndexVariable.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexData::LNodeIndexData(
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& ghosts)
    : SAMRAI::pdat::IndexData<NDIM,LNodeIndexSet>(box,ghosts),
      d_interior_local_indices(),
      d_ghost_local_indices()
{
    // intentionally blank
    return;
}// LNodeIndexData

LNodeIndexData::~LNodeIndexData()
{
    // intentionally blank
    return;
}// ~LNodeIndexData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexData>;

//////////////////////////////////////////////////////////////////////////////
