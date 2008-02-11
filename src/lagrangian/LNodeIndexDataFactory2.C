// Filename: LNodeIndexDataFactory2.C
// Last modified: <07.Feb.2008 00:30:26 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "LNodeIndexDataFactory2.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/LNodeIndexData2.h>

// SAMRAI INCLUDES
#include <tbox/ArenaManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int DEFAULT_DEPTH = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexDataFactory2::LNodeIndexDataFactory2(
    const SAMRAI::hier::IntVector<NDIM>& ghosts)
    : SAMRAI::pdat::CellDataFactory<NDIM,LNodeIndexSet>(DEFAULT_DEPTH,ghosts)
{
    // intentionally blank
    return;
}// LNodeIndexDataFactory2

LNodeIndexDataFactory2::~LNodeIndexDataFactory2()
{
    // intentionally blank
    return;
}// ~LNodeIndexDataFactory2

SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >
LNodeIndexDataFactory2::cloneFactory(
    const SAMRAI::hier::IntVector<NDIM>& ghosts)
{
    return new LNodeIndexDataFactory2(ghosts);
}// cloneFactory

SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
LNodeIndexDataFactory2::allocate(
    const SAMRAI::hier::Box<NDIM>& box,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool) const
{
    if (pool.isNull())
    {
        pool = SAMRAI::tbox::ArenaManager::getManager()->getStandardAllocator();
    }
    SAMRAI::hier::PatchData<NDIM>* pd = new (pool) LNodeIndexData2(box,getGhostCellWidth(),pool);
    return SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(pd, pool);
}// allocate

SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
LNodeIndexDataFactory2::allocate(
    const SAMRAI::hier::Patch<NDIM>& patch,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool) const
{
    return allocate(patch.getBox(), pool);
}// allocate

size_t
LNodeIndexDataFactory2::getSizeOfMemory(
    const SAMRAI::hier::Box<NDIM>& box) const
{
    (void) box;
    return SAMRAI::tbox::Arena::align(sizeof(LNodeIndexData2));
}// getSizeOfMemory

bool
LNodeIndexDataFactory2::validCopyTo(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >& dst_pdf) const
{
    SAMRAI::tbox::Pointer<LNodeIndexDataFactory2> lnidf = dst_pdf;
    return !lnidf.isNull();
}// validCopyTo

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexDataFactory2>;

#include <MultiblockCellDataTranslator.C>
template class SAMRAI::pdat::MultiblockCellDataTranslator<NDIM,IBAMR::LNodeIndexSet>;

//////////////////////////////////////////////////////////////////////////////
