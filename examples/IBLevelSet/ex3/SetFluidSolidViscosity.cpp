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

#include "SetFluidSolidViscosity.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidSolidViscosityCallbackFunction(int mu_idx,
                                           Pointer<Variable<NDIM> > mu_var,
                                           Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           const int cycle_num,
                                           const double time,
                                           const double current_time,
                                           const double new_time,
                                           void* ctx)
{
    // Set the density from the level set information
    static SetFluidSolidViscosity* ptr_SetFluidSolidViscosity = static_cast<SetFluidSolidViscosity*>(ctx);
    ptr_SetFluidSolidViscosity->setViscosityPatchData(
        mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidSolidViscosityCallBackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidSolidViscosity::SetFluidSolidViscosity(const std::string& object_name, const double mu_fluid)
    : d_object_name(object_name), d_mu_fluid(mu_fluid)
{
    // intentionally left blank
    return;
} // SetFluidSolidViscosity

void
SetFluidSolidViscosity::setViscosityPatchData(int mu_idx,
                                              Pointer<Variable<NDIM> > mu_var,
                                              Pointer<HierarchyMathOps> hier_math_ops,
                                              const int /*cycle_num*/,
                                              const double /*time*/,
                                              const double /*current_time*/,
                                              const double /*new_time*/)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    Pointer<CellVariable<NDIM, double> > mu_cc_var = mu_var;
    if (mu_cc_var)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);
                mu_data->fillAll(d_mu_fluid);
            }
        }
    }
    else
    {
        // Erroring out if any other centered is used for mu
        TBOX_ERROR("This statement should not have been reached");
    }

    return;
} // setViscosityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
