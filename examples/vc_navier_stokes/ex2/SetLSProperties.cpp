// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
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
#include <ibtk/HierarchyMathOps.h>

#include "SetLSProperties.h"

#include <CartesianGridGeometry.h>

#include <ibamr/app_namespaces.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

void
callSetLSCallbackFunction(int ls_idx,
                          Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                          const int integrator_step,
                          const double current_time,
                          const bool initial_time,
                          const bool regrid_time,
                          void* ctx)
{
    // Set the density from the level set information
    static SetLSProperties* ptr_SetLSProperties = static_cast<SetLSProperties*>(ctx);
    ptr_SetLSProperties->setLSPatchData(
        ls_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;

} // callSetFluidDensityCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetLSProperties::SetLSProperties(const std::string& object_name, Pointer<LSInitStrategy> ls_ops)
    : d_object_name(object_name), d_ls_ops(ls_ops)
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
SetLSProperties::setLSPatchData(int ls_idx,
                                SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                const int integrator_step,
                                const double current_time,
                                const bool initial_time,
                                const bool regrid_time)
{
    // If at the regrid time, force reinitialization
    pout << "Resetting level set data" << std::endl;
    d_ls_ops->setReinitializeLSData(regrid_time);
    d_ls_ops->initializeLSData(ls_idx, hier_math_ops, integrator_step, current_time, initial_time);

    return;
} // setLSPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
