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

// SAMRAI INCLUDES
#include "ibtk/samrai_compatibility_names.h"

#include "CartGridBodyForce.h"
#include "SAMRAICellData.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideData.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIUtilities.h"
#include "SAMRAIVariable.h"
#include "SAMRAIVariableDatabase.h"

////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

#include "ibamr/namespaces.h"

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
                                  SAMRAIPointer<SAMRAIVariable> var,
                                  SAMRAIPointer<SAMRAIPatch> patch,
                                  const double /*data_time*/,
                                  const bool /*initial_time*/,
                                  SAMRAIPointer<SAMRAIPatchLevel> /*patch_level*/)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!patch.isNull());
    TBOX_ASSERT(!var.isNull());
#endif

    SAMRAIPointer<SAMRAISideVariable<double>> copy_to_sc_var = var;
    SAMRAIPointer<SAMRAICellVariable<double>> copy_to_cc_var = var;

    if (!copy_to_sc_var.isNull())
    {
        SAMRAIPointer<SAMRAISideData<double>> copy_to_sc_data = patch->getPatchData(data_idx);
        const SAMRAIPointer<SAMRAISideData<double>> copy_from_sc_data = patch->getPatchData(d_body_force_idx);

#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(!copy_to_sc_data.isNull());
        TBOX_ASSERT(!copy_from_sc_data.isNull());
#endif

        copy_to_sc_data->copy(*copy_from_sc_data);
    }
    else if (!copy_to_cc_var.isNull())
    {
        SAMRAIPointer<SAMRAICellData<double>> copy_to_cc_data = patch->getPatchData(data_idx);
        const SAMRAIPointer<SAMRAICellData<double>> copy_from_cc_data = patch->getPatchData(d_body_force_idx);

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
