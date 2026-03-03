// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/samrai_compatibility_names.h>

#include "CartGridBodyForce.h"

#include <SAMRAIPatch.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIVariable.h>

////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAICellData.h>
#include <SAMRAICellVariable.h>
#include <SAMRAISideData.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIUtilities.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAI_config.h>

#include <ibamr/namespaces.h>

namespace IBTK
{
CartGridBodyForce::CartGridBodyForce(const int body_force_idx) : d_body_force_idx(body_force_idx)
{
    // this is intentionally left blank
    return;
} // CartGridBodyForce

bool
CartGridBodyForce::isTimeDependent() const
{
    return true;

} // isTimeDependent

void
CartGridBodyForce::setDataOnPatch(const int data_idx,
                                  Pointer<SAMRAIVariable> var,
                                  Pointer<SAMRAIPatch> patch,
                                  const double /*data_time*/,
                                  const bool /*initial_time*/,
                                  Pointer<SAMRAIPatchLevel> /*patch_level*/)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(!var.isNull());
#endif

    Pointer<SAMRAISideVariable<double>> copy_to_sc_var = var;
    Pointer<SAMRAICellVariable<double>> copy_to_cc_var = var;

    if (!copy_to_sc_var.isNull())
    {
        Pointer<SAMRAISideData<double>> copy_to_sc_data = patch->getPatchData(data_idx);
        const Pointer<SAMRAISideData<double>> copy_from_sc_data = patch->getPatchData(d_body_force_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!copy_to_sc_data.isNull());
        TBOX_ASSERT(!copy_from_sc_data.isNull());
#endif

        copy_to_sc_data->copy(*copy_from_sc_data);
    }
    else if (!copy_to_cc_var.isNull())
    {
        Pointer<SAMRAICellData<double>> copy_to_cc_data = patch->getPatchData(data_idx);
        const Pointer<SAMRAICellData<double>> copy_from_cc_data = patch->getPatchData(d_body_force_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!copy_to_cc_data.isNull());
        TBOX_ASSERT(!copy_from_cc_data.isNull());
#endif

        copy_to_cc_data->copy(*copy_from_cc_data);
    }
    else
    {
        TBOX_ERROR("CartGridBodyForce::setDataOnPatch() "
                   << "UNKNOWN DATA TYPE ENCOUNTERED" << std::endl);
    }

    return;

} // setDataonPatch

} // namespace IBTK
