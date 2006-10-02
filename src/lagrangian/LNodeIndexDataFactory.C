//
// LNodeIndexDataFactory.C
//
// Created on 01 Mar 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <14.Jun.2005 22:31:28 boyce@bigboy.verizon.net>
//

#include "LNodeIndexDataFactory.h"

// SAMRAI-tools INCLUDES
//
#include "LNodeIndexData.h"

// SAMRAI INCLUDES
//
#include "tbox/ArenaManager.h"

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexDataFactory::LNodeIndexDataFactory(
    const hier::IntVector<NDIM>& ghosts)
    : pdat::IndexDataFactory<NDIM,LNodeIndexSet>(ghosts)
{
    // intentionally blank
    return;
}// LNodeIndexDataFactory

LNodeIndexDataFactory::~LNodeIndexDataFactory()
{
    // intentionally blank
    return;
}// ~LNodeIndexDataFactory

tbox::Pointer<hier::PatchDataFactory<NDIM> > LNodeIndexDataFactory::cloneFactory()
{
    return(new LNodeIndexDataFactory(getDefaultGhostCellWidth()));
}// cloneFactory

tbox::Pointer<hier::PatchData<NDIM> > LNodeIndexDataFactory::allocate(
    const hier::Box<NDIM>& box,
    tbox::Pointer<tbox::Arena> pool) const
{
    if (pool.isNull())
    {
        pool = tbox::ArenaManager::getManager()->getStandardAllocator();
    }
    hier::PatchData<NDIM> *pd = new (pool) LNodeIndexData(box,getDefaultGhostCellWidth());
    return(tbox::Pointer<hier::PatchData<NDIM> >(pd, pool));
}// allocate

size_t LNodeIndexDataFactory::getSizeOfMemory(
    const hier::Box<NDIM>& box) const
{
    (void) box;
    return(tbox::Arena::align(sizeof(LNodeIndexData)));
}// getSizeOfMemory

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the LNodeIndexDataFactory
/// class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<LNodeIndexDataFactory>;

#endif

//////////////////////////////////////////////////////////////////////////////
