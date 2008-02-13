// Filename: LNodeIndexVariable2.C
// Last modified: <12.Feb.2008 21:20:49 griffith@box221.cims.nyu.edu>
// Created on 04 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "LNodeIndexVariable2.h"

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
#include <ibamr/LNodeIndexDataFactory2.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexVariable2::LNodeIndexVariable2(
    const std::string& name)
    : SAMRAI::hier::Variable<NDIM>(name, new LNodeIndexDataFactory2(SAMRAI::hier::IntVector<NDIM>(0)))
{
    // intentionally blank
    return;
}// LNodeIndexVariable2

LNodeIndexVariable2::~LNodeIndexVariable2()
{
    // intentionally blank
    return;
}// ~LNodeIndexVariable2

bool
LNodeIndexVariable2::dataLivesOnPatchBorder() const
{
    return false;
}// dataLivesOnPatchBorder

bool
LNodeIndexVariable2::fineBoundaryRepresentsVariable() const
{
    return true;
}// fineBoundaryRepresentsVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexVariable2>;

//////////////////////////////////////////////////////////////////////////////
