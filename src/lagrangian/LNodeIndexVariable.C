// Filename: LNodeIndexVariable.C
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <02.Oct.2006 13:16:26 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "LNodeIndexVariable.h"

// IBAMR INCLUDES
#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#endif

#include <ibamr/LNodeIndexDataFactory.h>

// SAMRAI INCLUDES
#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

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
