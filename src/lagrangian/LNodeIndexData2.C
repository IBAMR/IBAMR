// Filename: LNodeIndexData2.C
// Last modified: <19.Jun.2007 16:49:15 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "LNodeIndexData2.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// SAMRAI INCLUDES
#include <CellOverlap.h>

// ZLIB INCLUDES
#include <zlib.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int DEFAULT_DEPTH = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexData2::LNodeIndexData2(
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::hier::IntVector<NDIM>& ghosts,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool)
    : SAMRAI::pdat::CellData<NDIM,LNodeIndexSet>(box,DEFAULT_DEPTH,ghosts,pool),
      d_interior_local_indices(),
      d_ghost_local_indices()
{
    // intentionally blank
    return;
}// LNodeIndexData2

LNodeIndexData2::~LNodeIndexData2()
{
    // intentionally blank
    return;
}// ~LNodeIndexData2

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexData2>;

//////////////////////////////////////////////////////////////////////////////
