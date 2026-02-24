// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// APPLICATION INCLUDES
#include <ibtk/samrai_compatibility_names.h>
// SAMRAI INCLUDES
#include <ibtk/HierarchyMathOps.h>

#include "SetFluidProperties.h"

#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideData.h>
#include <SAMRAISideGeometry.h>
#include <SAMRAISideIndex.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIVariable.h>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidDensityCallbackFunction(int rho_idx,
                                    SAMRAIPointer<SAMRAIVariable> rho_var,
                                    SAMRAIPointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int cycle_num,
                                    const double time,
                                    const double current_time,
                                    const double new_time,
                                    void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidDensityCallbackFunction

// Various options to setting side-centered densities
#define SMOOTH_SC_RHO 1
#define DESJARDINS_SC_RHO 0

void
callSetFluidViscosityCallbackFunction(int mu_idx,
                                      SAMRAIPointer<SAMRAIVariable> mu_var,
                                      SAMRAIPointer<IBTK::HierarchyMathOps> hier_math_ops,
                                      const int cycle_num,
                                      const double time,
                                      const double current_time,
                                      const double new_time,
                                      void* ctx)
{
    // Set the density from the level set information
    static SetFluidProperties* ptr_SetFluidProperties = static_cast<SetFluidProperties*>(ctx);
    ptr_SetFluidProperties->setViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidViscosityCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidProperties::SetFluidProperties(const std::string& object_name, const double rho, const double mu)
    : d_object_name(object_name), d_rho(rho), d_mu(mu)
{
    // intentionally left blank
    return;
} // SetFluidProperties

SetFluidProperties::~SetFluidProperties()
{
    // intentionally left blank
    return;

} //~SetFluidProperties

void
SetFluidProperties::setDensityPatchData(int rho_idx,
                                        SAMRAIPointer<SAMRAIVariable> rho_var,
                                        SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                        const int /*cycle_num*/,
                                        const double /*time*/,
                                        const double /*current_time*/,
                                        const double /*new_time*/)
{
    SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    SAMRAIPointer<SAMRAICellVariable<double>> rho_cc_var = rho_var;
    if (rho_cc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
                const SAMRAIBox& patch_box = patch->getBox();
                SAMRAIPointer<SAMRAICellData<double>> rho_data = patch->getPatchData(rho_idx);

                for (SAMRAIBox::Iterator it(patch_box); it; it++)
                {
                    SAMRAICellIndex ci(it());
                    (*rho_data)(ci) = d_rho;
                }
            }
        }
    }

    SAMRAIPointer<SAMRAISideVariable<double>> rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();

                const SAMRAIBox& patch_box = patch->getBox();
                SAMRAIPointer<SAMRAISideData<double>> rho_data = patch->getPatchData(rho_idx);

                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAIBox::Iterator it(SAMRAISideGeometry::toSideBox(patch_box, axis)); it; it++)
                    {
                        SAMRAISideIndex si(it(), axis, 0);

                        (*rho_data)(si) = d_rho;
                    }
                }
            }
        }
    }

    return;
} // setDensityPatchData

void
SetFluidProperties::setViscosityPatchData(int mu_idx,
                                          SAMRAIPointer<SAMRAIVariable> /*mu_var*/,
                                          SAMRAIPointer<HierarchyMathOps> hier_math_ops,
                                          const int /*cycle_num*/,
                                          const double /*time*/,
                                          const double /*current_time*/,
                                          const double /*new_time*/)
{
    // Set the constant viscosity
    SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
            SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();

            const SAMRAIBox& patch_box = patch->getBox();
            SAMRAIPointer<SAMRAICellData<double>> mu_data = patch->getPatchData(mu_idx);

            for (SAMRAIBox::Iterator it(patch_box); it; it++)
            {
                SAMRAICellIndex ci(it());

                (*mu_data)(ci) = d_mu;
            }
        }
    }

    return;
} // setViscosityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
