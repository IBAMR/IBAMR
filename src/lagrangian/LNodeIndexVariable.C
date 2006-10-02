//
// LNodeIndexVariable.C
//
// Created on 01 Mar 2004
//         by Boyce Griffith (boyce@bigboy.speakeasy.net).
//
// Last modified: <07.Mar.2005 13:38:30 boyce@trasnaform.cims.nyu.edu>
//

#include "LNodeIndexVariable.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif

// SAMRAI-tools INCLUDES
//
#include "LNodeIndexDataFactory.h"

/////////////////////////////// PUBLIC ///////////////////////////////////////

LNodeIndexVariable::LNodeIndexVariable(
    const string& name)
    : hier::Variable<NDIM>(name,new LNodeIndexDataFactory(hier::IntVector<NDIM>(0)))
{
    // intentionally blank
    return;
}// LNodeIndexVariable

LNodeIndexVariable::~LNodeIndexVariable()
{
    // intentionally blank
    return;
}// ~LNodeIndexVariable

bool LNodeIndexVariable::dataLivesOnPatchBorder() const
{
    return(false);
}// dataLivesOnPatchBorder

bool LNodeIndexVariable::fineBoundaryRepresentsVariable() const
{
    return(true);
}// fineBoundaryRepresentsVariable

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#ifndef LACKS_EXPLICIT_TEMPLATE_INSTANTIATION

#include "tbox/Pointer.C"

//////////////////////////////////////////////////////////////////////
///
/// These declarations are required to use the LNodeIndexVariable
/// class.
///
//////////////////////////////////////////////////////////////////////

template class tbox::Pointer<LNodeIndexVariable>;

#endif

//////////////////////////////////////////////////////////////////////////////
