// Filename: LNodeIndexVariable.C
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <25.Oct.2006 18:28:03 boyce@bigboy.nyconnect.com>

#include "LNodeIndexVariable.h"

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
#include <ibamr/LNodeIndexDataFactory.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexVariable::LNodeIndexVariable(
    const std::string& name)
    : SAMRAI::hier::Variable<NDIM>(name, new LNodeIndexDataFactory(SAMRAI::hier::IntVector<NDIM>(0)))
{
    // intentionally blank
    return;
}// LNodeIndexVariable

LNodeIndexVariable::~LNodeIndexVariable()
{
    // intentionally blank
    return;
}// ~LNodeIndexVariable

bool
LNodeIndexVariable::dataLivesOnPatchBorder() const
{
    return false;
}// dataLivesOnPatchBorder

bool
LNodeIndexVariable::fineBoundaryRepresentsVariable() const
{
    return true;
}// fineBoundaryRepresentsVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::LNodeIndexVariable>;

//////////////////////////////////////////////////////////////////////////////
