// Filename: IBEulerianSourceSetter.C
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)
// Last modified: <24.Oct.2006 14:40:46 boyce@bigboy.nyconnect.com>

#include "IBEulerianSourceSetter.h"

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

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBEulerianSourceSetter::IBEulerianSourceSetter(
    const string& object_name,
    const int Q_idx)
    : SetDataStrategy(object_name),
      d_Q_idx(Q_idx)
{
    // intentionally blank
    return;
}// IBEulerianSourceSetter

IBEulerianSourceSetter::~IBEulerianSourceSetter()
{
    // intentionally blank
    return;
}// ~IBEulerianSourceSetter

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
    else
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > my_q_data = patch.getPatchData(d_Q_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!my_q_data.isNull());
#endif
        q_data->copy(*my_q_data);
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
