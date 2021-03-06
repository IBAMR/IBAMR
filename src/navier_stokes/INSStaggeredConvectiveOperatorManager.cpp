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
#include "ibamr/INSStaggeredCUIConvectiveOperator.h"
#include "ibamr/INSStaggeredCenteredConvectiveOperator.h"
#include "ibamr/INSStaggeredConvectiveOperatorManager.h"
#include "ibamr/INSStaggeredPPMConvectiveOperator.h"
#include "ibamr/INSStaggeredStabilizedPPMConvectiveOperator.h"
#include "ibamr/INSStaggeredUpwindConvectiveOperator.h"
#include "ibamr/INSStaggeredWavePropConvectiveOperator.h"
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

const std::string INSStaggeredConvectiveOperatorManager::DEFAULT = "DEFAULT";
const std::string INSStaggeredConvectiveOperatorManager::CENTERED = "CENTERED";
const std::string INSStaggeredConvectiveOperatorManager::PPM = "PPM";
const std::string INSStaggeredConvectiveOperatorManager::UPWIND = "UPWIND";
const std::string INSStaggeredConvectiveOperatorManager::STABILIZED_PPM = "STABILIZED_PPM";
const std::string INSStaggeredConvectiveOperatorManager::WAVE_PROP = "WAVE_PROP";
const std::string INSStaggeredConvectiveOperatorManager::CUI = "CUI";

INSStaggeredConvectiveOperatorManager* INSStaggeredConvectiveOperatorManager::s_operator_manager_instance = nullptr;
bool INSStaggeredConvectiveOperatorManager::s_registered_callback = false;
unsigned char INSStaggeredConvectiveOperatorManager::s_shutdown_priority = 200;

INSStaggeredConvectiveOperatorManager*
INSStaggeredConvectiveOperatorManager::getManager()
{
    if (!s_operator_manager_instance)
    {
        s_operator_manager_instance = new INSStaggeredConvectiveOperatorManager();
    }
    if (!s_registered_callback)
    {
        ShutdownRegistry::registerShutdownRoutine(freeManager, s_shutdown_priority);
        s_registered_callback = true;
    }
    return s_operator_manager_instance;
} // getManager

void
INSStaggeredConvectiveOperatorManager::freeManager()
{
    delete s_operator_manager_instance;
    s_operator_manager_instance = nullptr;
    return;
} // freeManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

Pointer<ConvectiveOperator>
INSStaggeredConvectiveOperatorManager::allocateOperator(const std::string& operator_type,
                                                        const std::string& operator_object_name,
                                                        Pointer<Database> input_db,
                                                        const ConvectiveDifferencingType difference_form,
                                                        const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs) const
{
    auto it = d_operator_maker_map.find(operator_type);
    if (it == d_operator_maker_map.end())
    {
        TBOX_ERROR("INSStaggeredConvectiveOperatorManager::allocateOperator():\n"
                   << "  unrecognized operator type: " << operator_type << "\n");
    }
    return (it->second)(operator_object_name, input_db, difference_form, bc_coefs);
} // allocateOperator

void
INSStaggeredConvectiveOperatorManager::registerOperatorFactoryFunction(const std::string& operator_type,
                                                                       OperatorMaker operator_maker)
{
    if (d_operator_maker_map.find(operator_type) != d_operator_maker_map.end())
    {
        pout << "INSStaggeredConvectiveOperatorManager::"
                "registerOperatorFactoryFunction():\n"
             << "  NOTICE: overriding initialization function for operator_type = " << operator_type << "\n";
    }
    d_operator_maker_map[operator_type] = operator_maker;
    return;
} // registerOperatorFactoryFunction

/////////////////////////////// PROTECTED ////////////////////////////////////

INSStaggeredConvectiveOperatorManager::INSStaggeredConvectiveOperatorManager() : d_operator_maker_map()
{
    registerOperatorFactoryFunction(DEFAULT, INSStaggeredPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(CENTERED, INSStaggeredCenteredConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(PPM, INSStaggeredPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(UPWIND, INSStaggeredUpwindConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(STABILIZED_PPM, INSStaggeredStabilizedPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(WAVE_PROP, INSStaggeredWavePropConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(CUI, INSStaggeredCUIConvectiveOperator::allocate_operator);
    return;
} // INSStaggeredConvectiveOperatorManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
