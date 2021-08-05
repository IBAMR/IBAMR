// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#include "SetLSProperties.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

/////////////////////////////// STATIC ///////////////////////////////////////

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

SetLSProperties::SetLSProperties(const std::string& object_name, Pointer<LSInitStrategy> ls_gas_ops)
    : d_object_name(object_name), d_ls_gas_ops(ls_gas_ops)
{
    // intentionally left blank
    return;
} // SetLSProperties

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
