// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

#include "GravityForcing.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

// SAMRAI INCLUDES
#include <SAMRAIBox.h>
#include <SAMRAIHierarchyDataOpsManager.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideData.h>
#include <SAMRAISideGeometry.h>
#include <SAMRAISideIndex.h>
#include <SAMRAIVariable.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

GravityForcing::GravityForcing(const std::string& object_name,
                               SAMRAIPointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                               std::vector<double> grav_const)
    : d_object_name(object_name), d_ins_hierarchy_integrator(ins_hierarchy_integrator), d_grav_const(grav_const)
{
    // intentionally blank
    return;
} // GravityForcing

GravityForcing::~GravityForcing()
{
    // intentionally blank
    return;
} // ~GravityForcing

bool
GravityForcing::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
GravityForcing::setDataOnPatchHierarchy(const int data_idx,
                                        SAMRAIPointer<SAMRAIVariable> /*var*/,
                                        SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
                                        const double /*data_time*/,
                                        const bool /*initial_time*/,
                                        const int coarsest_ln_in,
                                        const int finest_ln_in)
{
    const int coarsest_ln = (coarsest_ln_in == -1 ? 0 : coarsest_ln_in);
    const int finest_ln = (finest_ln_in == -1 ? hierarchy->getFinestLevelNumber() : finest_ln_in);

    // Get interpolated density variable
    const int rho_ins_idx = d_ins_hierarchy_integrator->getLinearOperatorRhoPatchDataIndex();

#if !defined(NDEBUG)
    TBOX_ASSERT(rho_ins_idx >= 0);
#endif
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(ln);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
            const SAMRAIBox& box = patch->getBox();
            SAMRAIPointer<SAMRAISideData<double>> f_data = patch->getPatchData(data_idx);
            const SAMRAIPointer<SAMRAISideData<double>> rho_data = patch->getPatchData(rho_ins_idx);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SAMRAIBox::Iterator it(SAMRAISideGeometry::toSideBox(box, axis)); it; it++)
                {
                    SAMRAISideIndex s_i(it(), axis, SAMRAISideIndex::Lower);
                    (*f_data)(s_i) = ((*rho_data)(s_i)) * d_grav_const[axis];
                }
            }
        }
    }
    return;
} // setDataOnPatchHierarchy

void
GravityForcing::setDataOnPatch(const int data_idx,
                               SAMRAIPointer<SAMRAIVariable> /*var*/,
                               SAMRAIPointer<SAMRAIPatch> patch,
                               const double /*data_time*/,
                               const bool initial_time,
                               SAMRAIPointer<SAMRAIPatchLevel> /*patch_level*/)
{
    if (initial_time)
    {
        SAMRAIPointer<SAMRAISideData<double>> f_data = patch->getPatchData(data_idx);
        f_data->fillAll(0.0);
    }
    // Intentionally left blank

} // setDataOnPatch

//////////////////////////////////////////////////////////////////////////////
