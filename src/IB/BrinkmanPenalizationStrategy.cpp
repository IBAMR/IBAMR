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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/BrinkmanPenalizationStrategy.h"

#include "tbox/Database.h"
#include "tbox/RestartManager.h"

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC //////////////////////////////////////
BrinkmanPenalizationStrategy::BrinkmanPenalizationStrategy(std::string object_name, bool register_for_restart)
    : d_object_name(std::move(object_name)), d_registered_for_restart(register_for_restart)
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }
    return;
} // BrinkmanPenalizationStrategy

BrinkmanPenalizationStrategy::~BrinkmanPenalizationStrategy()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~BrinkmanPenalizationStrategy

void
BrinkmanPenalizationStrategy::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

void
BrinkmanPenalizationStrategy::preprocessComputeBrinkmanPenalization(double /*current_time*/,
                                                                    double /*new_time*/,
                                                                    int /*num_cycles*/)
{
    return;
} // preprocessComputeBrinkmanPenalization

void
BrinkmanPenalizationStrategy::postprocessComputeBrinkmanPenalization(double /*current_time*/,
                                                                     double /*new_time*/,
                                                                     int /*num_cycles*/)
{
    d_current_time = std::numeric_limits<double>::quiet_NaN(), d_new_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessComputeBrinkmanPenalization

void BrinkmanPenalizationStrategy::putToDatabase(Pointer<Database> /*db*/)
{
    return;
} // putToDatabase

void
BrinkmanPenalizationStrategy::setBrinkmanCoefficient(double chi)
{
    d_chi = chi;
    return;
} // setBrinkmanCoefficient

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
