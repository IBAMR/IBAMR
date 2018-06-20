// Filename: LSInitStrategy.cpp
// Created on 27 Sep 2017 by Amneet Bhalla and Nishant Nangia
//
// Copyright (c) 2002-2017, Amneet Bhalla and Nishant Nangia
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/LSInitStrategy.h"
#include "ibamr/namespaces.h"
#include "tbox/RestartManager.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LSInitStrategy::LSInitStrategy(const std::string& object_name, bool register_for_restart)
    : d_object_name(object_name), d_registered_for_restart(register_for_restart)
{
    // Some default values.
    d_ls_order = FIRST_ORDER_LS;
    d_max_its = 100;
    d_abs_tol = 1e-5;
    d_enable_logging = false;
    d_bc_coef = NULL;
    d_reinitialize_ls = false;
    d_reinit_interval = 0;

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    return;
} // LSInitStrategy

LSInitStrategy::~LSInitStrategy()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }
    d_registered_for_restart = false;

} // ~LSInitStrategy

void
LSInitStrategy::registerPhysicalBoundaryCondition(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* robin_bc_coef)
{
    d_bc_coef = robin_bc_coef;
    return;
} // registerPhysicalBoundaryCondition

void
LSInitStrategy::registerInterfaceNeighborhoodLocatingFcn(LocateInterfaceNeighborhoodFcnPtr callback_fcn, void* ctx)
{
    d_locate_interface_fcns.push_back(callback_fcn);
    d_locate_interface_fcns_ctx.push_back(ctx);

    return;
} // registerInterfaceNeighborhoodLocatingFcn

void
LSInitStrategy::setReinitializeLSData(bool reinit_ls_data)
{
    d_reinitialize_ls = reinit_ls_data;
    return;
} // setReinitializeLSData

void LSInitStrategy::putToDatabase(Pointer<Database> /*db*/)
{
    // intentionally blank
    return;
} // putToDatabase

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
