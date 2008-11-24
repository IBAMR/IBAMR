// Filename: IBEulerianForceSetter.C
// Last modified: <20.Nov.2008 20:02:53 griffith@box230.cims.nyu.edu>
// Created on 28 Sep 2004 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

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
#include <PatchCellDataBasicOps.h>
#include <PatchSideDataBasicOps.h>
#include <SideData.h>
#include <tbox/MathUtilities.h>

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

IBEulerianForceSetter::IBEulerianForceSetter(
    const std::string& object_name,
    const int F_current_idx,
    const int F_new_idx,
    const int F_half_idx)
    : IBTK::SetDataStrategy(object_name),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_F_current_idx(F_current_idx),
      d_F_new_idx(F_new_idx),
      d_F_half_idx(F_half_idx),
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
IBEulerianForceSetter::registerBodyForceSpecification(
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_set)
{
    d_body_force_set = F_set;
    return;
}// registerBodyForceSpecification

void
IBEulerianForceSetter::setDataOnPatch(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
    SAMRAI::hier::Patch<NDIM>& patch,
    const double data_time,
    const bool initial_time)
{
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > f_data = patch.getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_data.isNull());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_cc_data = f_data;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_sc_data = f_data;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_cc_data.isNull() || !f_sc_data.isNull());
#endif
    if (!d_body_force_set.isNull())
    {
        d_body_force_set->setDataOnPatch(data_idx, var, patch, data_time, initial_time);
    }
    else
    {
        if (!f_cc_data.isNull()) f_cc_data->fillAll(0.0);
        if (!f_sc_data.isNull()) f_sc_data->fillAll(0.0);
    }

    if (initial_time) return;

    int ib_data_idx = -1;
    if (SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, d_current_time))
    {
        ib_data_idx = d_F_current_idx;
    }
    else if (SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, d_new_time))
    {
        ib_data_idx = d_F_new_idx;
    }
    else if (SAMRAI::tbox::MathUtilities<double>::equalEps(data_time, 0.5*(d_current_time+d_new_time)))
    {
        ib_data_idx = d_F_half_idx;
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch():\n"
                   << "  data time " << data_time << " is not the current, new, or half time." << std::endl);
    }

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > f_ib_data = patch.getPatchData(ib_data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_ib_data.isNull());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_ib_cc_data = f_ib_data;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > f_ib_sc_data = f_ib_data;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!f_ib_cc_data.isNull() || !f_ib_sc_data.isNull());
    TBOX_ASSERT((!f_ib_cc_data.isNull() && !f_cc_data.isNull()) || (!f_ib_sc_data.isNull() && !f_sc_data.isNull()));
#endif
    if (!f_cc_data.isNull())
    {
        SAMRAI::math::PatchCellDataBasicOps<NDIM,double> patch_ops;
        patch_ops.add(f_cc_data, f_cc_data, f_ib_cc_data, patch.getBox());
    }
    if (!f_sc_data.isNull())
    {
        SAMRAI::math::PatchSideDataBasicOps<NDIM,double> patch_ops;
        patch_ops.add(f_sc_data, f_sc_data, f_ib_sc_data, patch.getBox());
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
