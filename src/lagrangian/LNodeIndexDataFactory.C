// Filename: LNodeIndexDataFactory.C
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <24.Oct.2006 14:17:28 boyce@bigboy.nyconnect.com>

#include "LNodeIndexDataFactory.h"

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
#include <ibamr/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <tbox/ArenaManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexDataFactory::LNodeIndexDataFactory(
    const SAMRAI::hier::IntVector<NDIM>& ghosts)
    : SAMRAI::pdat::IndexDataFactory<NDIM,LNodeIndexSet>(ghosts)
{
    // intentionally blank
    return;
}// LNodeIndexDataFactory

LNodeIndexDataFactory::~LNodeIndexDataFactory()
{
    // intentionally blank
    return;
}// ~LNodeIndexDataFactory

SAMRAI::tbox::Pointer<SAMRAI::hier::PatchDataFactory<NDIM> >
LNodeIndexDataFactory::cloneFactory()
{
    return new LNodeIndexDataFactory(getDefaultGhostCellWidth());
}// cloneFactory

SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >
LNodeIndexDataFactory::allocate(
    const SAMRAI::hier::Box<NDIM>& box,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Arena> pool) const
{
    if (pool.isNull())
    {
        pool = SAMRAI::tbox::ArenaManager::getManager()->getStandardAllocator();
    }
    SAMRAI::hier::PatchData<NDIM>* pd = new (pool) LNodeIndexData(box,getDefaultGhostCellWidth());
    return SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> >(pd, pool);
}// allocate

size_t
LNodeIndexDataFactory::getSizeOfMemory(
    const SAMRAI::hier::Box<NDIM>& box) const
{
    (void) box;
    return SAMRAI::tbox::Arena::align(sizeof(LNodeIndexData));
}// getSizeOfMemory

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexDataFactory>;

//////////////////////////////////////////////////////////////////////////////
