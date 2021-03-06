// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/LSInitStrategy.h"

#include "tbox/Database.h"
#include "tbox/RestartManager.h"

#include <utility>

#include "ibamr/namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LSInitStrategy::LSInitStrategy(std::string object_name, bool register_for_restart)
    : d_object_name(std::move(object_name)), d_registered_for_restart(register_for_restart)
{
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
