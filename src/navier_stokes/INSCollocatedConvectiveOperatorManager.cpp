// Filename: INSCollocatedConvectiveOperatorManager.cpp
// Created on 19 Aug 2012 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include <cstddef>
#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSCollocatedCenteredConvectiveOperator.h"
#include "ibamr/INSCollocatedConvectiveOperatorManager.h"
#include "ibamr/INSCollocatedPPMConvectiveOperator.h"
#include "ibamr/INSCollocatedWavePropConvectiveOperator.h"
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

const std::string INSCollocatedConvectiveOperatorManager::DEFAULT = "DEFAULT";
const std::string INSCollocatedConvectiveOperatorManager::CENTERED = "CENTERED";
const std::string INSCollocatedConvectiveOperatorManager::PPM = "PPM";
const std::string INSCollocatedConvectiveOperatorManager::WAVE_PROP = "WAVE_PROP";

INSCollocatedConvectiveOperatorManager* INSCollocatedConvectiveOperatorManager::s_operator_manager_instance = nullptr;
bool INSCollocatedConvectiveOperatorManager::s_registered_callback = false;
unsigned char INSCollocatedConvectiveOperatorManager::s_shutdown_priority = 200;

INSCollocatedConvectiveOperatorManager*
INSCollocatedConvectiveOperatorManager::getManager()
{
    if (!s_operator_manager_instance)
    {
        s_operator_manager_instance = new INSCollocatedConvectiveOperatorManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_operator_manager_instance;
} // getManager

void
INSCollocatedConvectiveOperatorManager::freeManager()
{
    delete s_operator_manager_instance;
    s_operator_manager_instance = nullptr;
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<ConvectiveOperator>
INSCollocatedConvectiveOperatorManager::allocateOperator(const std::string& operator_type,
                                                         const std::string& operator_object_name,
                                                         Pointer<Database> input_db,
                                                         const ConvectiveDifferencingType difference_form,
                                                         const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs) const
{
    auto it = d_operator_maker_map.find(operator_type);
    if (it == d_operator_maker_map.end())
    {
        TBOX_ERROR("INSCollocatedConvectiveOperatorManager::allocateOperator():\n"
                   << "  unrecognized operator type: "
                   << operator_type
                   << "\n");
    }
    return (it->second)(operator_object_name, input_db, difference_form, bc_coefs);
} // allocateOperator

void
INSCollocatedConvectiveOperatorManager::registerOperatorFactoryFunction(const std::string& operator_type,
                                                                        OperatorMaker operator_maker)
{
    if (d_operator_maker_map.find(operator_type) != d_operator_maker_map.end())
    {
        pout << "INSCollocatedConvectiveOperatorManager::registerOperatorFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for operator_type = " << operator_type << "\n";
    }
    d_operator_maker_map[operator_type] = operator_maker;
    return;
} // registerOperatorFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

INSCollocatedConvectiveOperatorManager::INSCollocatedConvectiveOperatorManager() : d_operator_maker_map()
{
    registerOperatorFactoryFunction(DEFAULT, INSCollocatedPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(CENTERED, INSCollocatedCenteredConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(PPM, INSCollocatedPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(WAVE_PROP, INSCollocatedWavePropConvectiveOperator::allocate_operator);
    return;
} // INSCollocatedConvectiveOperatorManager

INSCollocatedConvectiveOperatorManager::~INSCollocatedConvectiveOperatorManager()
{
    // intentionally blank
    return;
} // ~INSCollocatedConvectiveOperatorManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
