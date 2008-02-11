// Filename: IBEulerianSourceSetter.C
// Last modified: <04.Feb.2008 21:46:23 griffith@box221.cims.nyu.edu>
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)

#include "IBEulerianSourceSetter.h"

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
#include <tbox/MathUtilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBEulerianSourceSetter::IBEulerianSourceSetter(
    const std::string& object_name,
    const int Q_current_idx,
    const int Q_new_idx,
    const int Q_half_idx)
    : STOOLS::SetDataStrategy(object_name),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_Q_current_idx(Q_current_idx),
      d_Q_new_idx(Q_new_idx),
      d_Q_half_idx(Q_half_idx)
{
    // intentionally blank
    return;
}// IBEulerianSourceSetter

IBEulerianSourceSetter::~IBEulerianSourceSetter()
{
    // intentionally blank
    return;
}// ~IBEulerianSourceSetter

void
IBEulerianSourceSetter::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

bool
IBEulerianSourceSetter::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
IBEulerianSourceSetter::setDataOnPatch(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!q_data.isNull());
#endif
    if (initial_time)
    {
        q_data->fillAll(0.0);
    }
    else if (SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        if (d_Q_current_idx != -1)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_current_data = patch.getPatchData(d_Q_current_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!q_current_data.isNull());
#endif
            q_data->copy(*q_current_data);
        }
    }
    else if (SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        if (d_Q_new_idx != -1)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_new_data = patch.getPatchData(d_Q_new_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!q_new_data.isNull());
#endif
            q_data->copy(*q_new_data);
        }
    }
    else if (SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, 0.5*(d_current_time+d_new_time)))
    {
        if (d_Q_half_idx != -1)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_half_data = patch.getPatchData(d_Q_half_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!q_half_data.isNull());
#endif
            q_data->copy(*q_half_data);
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                   << "  data time " << data_time << " is not the current, new, or half time." << std::endl);
    }
    return;
}// setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBEulerianSourceSetter>;

//////////////////////////////////////////////////////////////////////////////
