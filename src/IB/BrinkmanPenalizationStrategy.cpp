// Filename BrinkmanPenalizationStrategy.cpp
// Created on Dec 5, 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/BrinkmanPenalizationStrategy.h"
#include "ibamr/namespaces.h"
#include "tbox/RestartManager.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC //////////////////////////////////////
BrinkmanPenalizationStrategy::BrinkmanPenalizationStrategy(std::string object_name, bool register_for_restart)
    : d_object_name(std::move(object_name))
{
    d_registered_for_restart = false;
    if (register_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
        d_registered_for_restart = true;
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
