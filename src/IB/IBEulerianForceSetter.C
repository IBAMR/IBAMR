// Filename: IBEulerianForceSetter.C
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <30.Mar.2007 16:09:49 griffith@box221.cims.nyu.edu>

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
    const int F_current_idx,
    const int F_new_idx,
    const int F_half_idx)
    : STOOLS::SetDataStrategy(object_name),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_F_current_idx(F_current_idx),
      d_F_new_idx(F_new_idx),
      d_F_half_idx(F_half_idx)
{
    // intentionally blank
    return;
}// IBEulerianForceSetter

IBEulerianForceSetter::~IBEulerianForceSetter()
{
    // intentionally blank
    return;
}// ~IBEulerianForceSetter

void
IBEulerianForceSetter::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

bool
IBEulerianForceSetter::isTimeDependent() const
{
    return true;
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
    else if (SAMRAI::tbox::Utilities::deq(data_time, d_current_time))
    {
        if (d_F_current_idx != -1)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_current_data = patch.getPatchData(d_F_current_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!f_current_data.isNull());
#endif
            f_data->copy(*f_current_data);
        }
    }
    else if (SAMRAI::tbox::Utilities::deq(data_time, d_new_time))
    {
        if (d_F_new_idx != -1)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_new_data = patch.getPatchData(d_F_new_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!f_new_data.isNull());
#endif
            f_data->copy(*f_new_data);
        }
    }
    else if (SAMRAI::tbox::Utilities::deq(data_time, 0.5*(d_current_time+d_new_time)))
    {
        if (d_F_half_idx != -1)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_half_data = patch.getPatchData(d_F_half_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!f_half_data.isNull());
#endif
            f_data->copy(*f_half_data);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                   << "  data time " << data_time << " is not the current, new, or half time." << endl);
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
