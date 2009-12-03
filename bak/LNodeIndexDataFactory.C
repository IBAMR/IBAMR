// Filename: LNodeIndexDataFactory.C
// Last modified: <17.Apr.2007 18:20:08 griffith@box221.cims.nyu.edu>
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "LNodeIndexDataFactory.h"

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
