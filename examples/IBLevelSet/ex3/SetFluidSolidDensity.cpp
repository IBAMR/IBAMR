// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/HierarchyMathOps.h>

#include "SetFluidSolidDensity.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidSolidDensityCallbackFunction(int rho_idx,
                                         Pointer<Variable<NDIM> > rho_var,
                                         Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         const int cycle_num,
                                         const double time,
                                         const double current_time,
                                         const double new_time,
                                         void* ctx)
{
    // Set the density from the level set information
    static SetFluidSolidDensity* ptr_SetFluidSolidDensity = static_cast<SetFluidSolidDensity*>(ctx);
    ptr_SetFluidSolidDensity->setDensityPatchData(
        rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidSolidDensityCallBackFunction

// Various options to setting side-centered densities

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidSolidDensity::SetFluidSolidDensity(const std::string& object_name, const double rho_fluid)
    : d_object_name(object_name), d_rho_fluid(rho_fluid)
{
    // intentionally left blank
    return;
} // SetFluidGasSolidDensity

void
SetFluidSolidDensity::setDensityPatchData(int rho_idx,
                                          Pointer<Variable<NDIM> > rho_var,
                                          Pointer<HierarchyMathOps> hier_math_ops,
                                          const int /*cycle_num*/,
                                          const double /*time*/,
                                          const double /*current_time*/,
                                          const double /*new_time*/)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    Pointer<SideVariable<NDIM, double> > rho_sc_var = rho_var;
    if (rho_sc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);
                rho_data->fillAll(d_rho_fluid);
            }
        }
    }
    else
    {
        // Erroring out if any other centered is used for rho
        TBOX_ERROR("This statement should not have been reached");
    }

    return;
} // setDensityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
