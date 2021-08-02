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

// APPLICATION INCLUDES
#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyMathOps.h>

#include "SetFluidProperties.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetFluidDensityCallbackFunction(int rho_idx,
                                    Pointer<Variable<NDIM> > rho_var,
                                    Pointer<IBTK::HierarchyMathOps> hier_math_ops,
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

void
callSetFluidViscosityCallbackFunction(int mu_idx,
                                      Pointer<Variable<NDIM> > mu_var,
                                      Pointer<IBTK::HierarchyMathOps> hier_math_ops,
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

SetFluidProperties::SetFluidProperties(const std::string& object_name,
                                       Pointer<IBTK::CartGridFunction> rho_fcn,
                                       Pointer<IBTK::CartGridFunction> mu_fcn)
    : d_object_name(object_name), d_rho_fcn(rho_fcn), d_mu_fcn(mu_fcn)
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
                                        Pointer<Variable<NDIM> > rho_var,
                                        Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        const int /*cycle_num*/,
                                        const double time,
                                        const double /*current_time*/,
                                        const double /*new_time*/)
{
    // Set the density based on the function specified in the input file
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    d_rho_fcn->setDataOnPatchHierarchy(rho_idx, rho_var, patch_hierarchy, time);

    return;
} // setDensityPatchData

void
SetFluidProperties::setViscosityPatchData(int mu_idx,
                                          Pointer<Variable<NDIM> > mu_var,
                                          Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                          const int /*cycle_num*/,
                                          const double time,
                                          const double /*current_time*/,
                                          const double /*new_time*/)
{
    // Set the viscosity based on the function specified in the input file
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    d_mu_fcn->setDataOnPatchHierarchy(mu_idx, mu_var, patch_hierarchy, time);

    return;
} // setViscosityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
