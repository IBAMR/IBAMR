// Filename: SetLSProperties.cpp
// Created on Dec 17, 2017 by Nishant Nangia
//
// Copyright (c) 2002-2017, Nishant Nangia and Amneet Bhalla
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

#include <ibtk/HierarchyMathOps.h>

#include "SetLSProperties.h"

#include <CartesianGridGeometry.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetSolidLSCallbackFunction(int ls_solid_idx,
                               Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int integrator_step,
                               const double current_time,
                               const bool initial_time,
                               const bool regrid_time,
                               void* ctx)
{
    // Set the density from the level set information
    static SetLSProperties* ptr_SetLSProperties = static_cast<SetLSProperties*>(ctx);
    ptr_SetLSProperties->setLSSolidPatchData(
        ls_solid_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;

} // callSetSolidLSCallbackFunction

void
callSetGasLSCallbackFunction(int ls_gas_idx,
                             Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int integrator_step,
                             const double current_time,
                             const bool initial_time,
                             const bool regrid_time,
                             void* ctx)
{
    // Set the density from the level set information
    static SetLSProperties* ptr_SetLSProperties = static_cast<SetLSProperties*>(ctx);
    ptr_SetLSProperties->setLSGasPatchData(
        ls_gas_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;

} // callSetGasLSCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetLSProperties::SetLSProperties(const std::string& object_name,
                                 Pointer<LSInitStrategy> ls_solid_ops,
                                 Pointer<LSInitStrategy> ls_gas_ops)
    : d_object_name(object_name), d_ls_solid_ops(ls_solid_ops), d_ls_gas_ops(ls_gas_ops)
{
    // intentionally left blank
    return;
} // SetLSProperties

SetLSProperties::~SetLSProperties()
{
    // intentionally left blank
    return;

} //~SetLSProperties

void
SetLSProperties::setLSSolidPatchData(int ls_solid_idx,
                                     SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                     const int integrator_step,
                                     const double current_time,
                                     const bool initial_time,
                                     const bool regrid_time)
{
    // If at the regrid time, force reinitialization
    d_ls_solid_ops->setReinitializeLSData(regrid_time);
    d_ls_solid_ops->initializeLSData(ls_solid_idx, hier_math_ops, integrator_step, current_time, initial_time);

    return;
} // setLSSolidPatchData

void
SetLSProperties::setLSGasPatchData(int ls_gas_idx,
                                   SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                   const int integrator_step,
                                   const double current_time,
                                   const bool initial_time,
                                   const bool regrid_time)
{
    // If at the regrid time, force reinitialization
    d_ls_gas_ops->setReinitializeLSData(regrid_time);
    d_ls_gas_ops->initializeLSData(ls_gas_idx, hier_math_ops, integrator_step, current_time, initial_time);

    return;
} // setLSSolidPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
