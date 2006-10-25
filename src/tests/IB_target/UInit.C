// Filename: UInit.C
// Last modified: <24.Oct.2006 21:38:24 boyce@bigboy.nyconnect.com>
// Created on 24 Oct 2006 by Boyce Griffith (boyce@bigboy.nyconnect.com)

#include "UInit.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#define included_IBAMR_config
#include <IBAMR_config.h>
#endif

#ifndef included_SAMRAI_config
#define included_SAMRAI_config
#include <SAMRAI_config.h>
#endif

// SAMRAI INCLUDES
#include <CellData.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

UInit::UInit(
    const string& object_name)
    : SetDataStrategy(object_name),
      d_object_name(object_name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
#endif
    // intentionally blank
    return;
}// UInit

UInit::~UInit()
{
    // intentionally blank
    return;
}// ~UInit

void
UInit::setDataOnPatch(
    const int data_idx,
    tbox::Pointer<hier::Variable<NDIM> > var,
    hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    tbox::Pointer< pdat::CellData<NDIM,double> > U_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!U_data.isNull());
#endif
    U_data->fillAll(0.0, patch.getBox());
    U_data->fill(1.0, patch.getBox(), 0);
    return;
}// setDataOnPatch

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
