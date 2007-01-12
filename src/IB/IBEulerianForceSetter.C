// Filename: IBEulerianForceSetter.C
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <11.Jan.2007 16:33:51 griffith@box221.cims.nyu.edu>

#include "IBEulerianForceSetter.h"

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
#include <CellData.h>
#include <PatchCellDataOpsReal.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBEulerianForceSetter::IBEulerianForceSetter(
    const string& object_name,
    const int F_idx)
    : SetDataStrategy(object_name),
      d_F_idx(F_idx)
{
    // intentionally blank
    return;
}// IBEulerianForceSetter

IBEulerianForceSetter::~IBEulerianForceSetter()
{
    // intentionally blank
    return;
}// ~IBEulerianForceSetter

bool
IBEulerianForceSetter::isTimeDependent() const
{
    return(true);
}// isTimeDependent

void
IBEulerianForceSetter::setDataOnPatch(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!f_data.isNull());
#endif
    if (initial_time)
    {
        f_data->fillAll(0.0);
    }
    else
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > my_f_data = patch.getPatchData(d_F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!my_f_data.isNull());
#endif
        f_data->copy(*my_f_data);
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBEulerianForceSetter>;

//////////////////////////////////////////////////////////////////////////////
