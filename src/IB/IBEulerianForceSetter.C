// Filename: IBEulerianForceSetter.C
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <24.Oct.2006 14:40:34 boyce@bigboy.nyconnect.com>

#include "IBEulerianForceSetter.h"

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
      d_F_idx(F_idx),
      d_body_force_set(NULL)
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
IBEulerianForceSetter::registerBodyForce(
    SAMRAI::tbox::Pointer<SetDataStrategy> body_force_set)
{
    d_body_force_set = body_force_set;
    return;
}// registerBodyForce

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
    else if (d_body_force_set.isNull())
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > my_f_data = patch.getPatchData(d_F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!my_f_data.isNull());
#endif
        f_data->copy(*my_f_data);
    }
    else
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > my_f_data = patch.getPatchData(d_F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(!my_f_data.isNull());
#endif
        d_body_force_set->setDataOnPatch(
            data_idx, var, patch, data_time, initial_time);

        SAMRAI::math::PatchCellDataOpsReal<NDIM,double> patch_ops;
        patch_ops.add(f_data, f_data, my_f_data, patch.getBox());
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
