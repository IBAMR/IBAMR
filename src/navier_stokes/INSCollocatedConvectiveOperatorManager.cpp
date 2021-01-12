// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSCollocatedCenteredConvectiveOperator.h"
#include "ibamr/INSCollocatedConvectiveOperatorManager.h"
#include "ibamr/INSCollocatedPPMConvectiveOperator.h"
#include "ibamr/INSCollocatedWavePropConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/ShutdownRegistry.h"
#include "tbox/Utilities.h"

#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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
                   << "  unrecognized operator type: " << operator_type << "\n");
    }
    return (it->second)(operator_object_name, input_db, difference_form, bc_coefs);
} // allocateOperator

void
INSCollocatedConvectiveOperatorManager::registerOperatorFactoryFunction(const std::string& operator_type,
                                                                        OperatorMaker operator_maker)
{
    if (d_operator_maker_map.find(operator_type) != d_operator_maker_map.end())
    {
        pout << "INSCollocatedConvectiveOperatorManager::"
                "registerOperatorFactoryFunction():\n"
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

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
