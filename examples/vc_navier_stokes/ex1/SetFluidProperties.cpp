// Filename: SetFluidProperties.cpp
// Created on March 6, 2018 by Nishant Nangia
//
// Copyright (c) 2018, Nishant Nangia and Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// APPLICATION INCLUDES
#include <ibamr/app_namespaces.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/HierarchyMathOps.h>

#include "SetFluidProperties.h"

#include <CartesianGridGeometry.h>

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
