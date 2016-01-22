// Filename: AdvDiffConvectiveOperatorManager.cpp
// Created on 17 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <stddef.h>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "CellVariable.h"
#include "ibamr/AdvDiffCenteredConvectiveOperator.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffPPMConvectiveOperator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

const std::string AdvDiffConvectiveOperatorManager::DEFAULT = "DEFAULT";
const std::string AdvDiffConvectiveOperatorManager::CENTERED = "CENTERED";
const std::string AdvDiffConvectiveOperatorManager::PPM = "PPM";

AdvDiffConvectiveOperatorManager* AdvDiffConvectiveOperatorManager::s_operator_manager_instance = NULL;
bool AdvDiffConvectiveOperatorManager::s_registered_callback = false;
unsigned char AdvDiffConvectiveOperatorManager::s_shutdown_priority = 200;

AdvDiffConvectiveOperatorManager*
AdvDiffConvectiveOperatorManager::getManager()
{
    if (!s_operator_manager_instance)
    {
        s_operator_manager_instance = new AdvDiffConvectiveOperatorManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_operator_manager_instance;
} // getManager

void
AdvDiffConvectiveOperatorManager::freeManager()
{
    delete s_operator_manager_instance;
    s_operator_manager_instance = NULL;
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<ConvectiveOperator>
AdvDiffConvectiveOperatorManager::allocateOperator(const std::string& operator_type,
                                                   const std::string& operator_object_name,
                                                   Pointer<CellVariable<NDIM, double> > Q_var,
                                                   Pointer<Database> input_db,
                                                   ConvectiveDifferencingType difference_form,
                                                   const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs) const
{
    std::map<std::string, OperatorMaker>::const_iterator it = d_operator_maker_map.find(operator_type);
    if (it == d_operator_maker_map.end())
    {
        TBOX_ERROR("AdvDiffConvectiveOperatorManager::allocateOperator():\n"
                   << "  unrecognized operator type: "
                   << operator_type
                   << "\n");
    }
    return (it->second)(operator_object_name, Q_var, input_db, difference_form, bc_coefs);
} // allocateOperator

void
AdvDiffConvectiveOperatorManager::registerOperatorFactoryFunction(const std::string& operator_type,
                                                                  OperatorMaker operator_maker)
{
    if (d_operator_maker_map.find(operator_type) != d_operator_maker_map.end())
    {
        pout << "AdvDiffConvectiveOperatorManager::registerOperatorFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for operator_type = " << operator_type << "\n";
    }
    d_operator_maker_map[operator_type] = operator_maker;
    return;
} // registerOperatorFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

AdvDiffConvectiveOperatorManager::AdvDiffConvectiveOperatorManager() : d_operator_maker_map()
{
    registerOperatorFactoryFunction(DEFAULT, AdvDiffPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(CENTERED, AdvDiffCenteredConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(PPM, AdvDiffPPMConvectiveOperator::allocate_operator);
    return;
} // AdvDiffConvectiveOperatorManager

AdvDiffConvectiveOperatorManager::~AdvDiffConvectiveOperatorManager()
{
    // intentionally blank
    return;
} // ~AdvDiffConvectiveOperatorManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
