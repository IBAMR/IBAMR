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

#include "ibamr/AdvDiffCUIConvectiveOperator.h"
#include "ibamr/AdvDiffCenteredConvectiveOperator.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffPPMConvectiveOperator.h"
#include "ibamr/AdvDiffWavePropConvectiveOperator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

#include "CellVariable.h"
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

const std::string AdvDiffConvectiveOperatorManager::DEFAULT = "DEFAULT";
const std::string AdvDiffConvectiveOperatorManager::CENTERED = "CENTERED";
const std::string AdvDiffConvectiveOperatorManager::CUI = "CUI";
const std::string AdvDiffConvectiveOperatorManager::PPM = "PPM";
const std::string AdvDiffConvectiveOperatorManager::WAVE_PROP = "WAVE_PROP";

AdvDiffConvectiveOperatorManager* AdvDiffConvectiveOperatorManager::s_operator_manager_instance = nullptr;
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
    s_operator_manager_instance = nullptr;
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
    auto it = d_operator_maker_map.find(operator_type);
    if (it == d_operator_maker_map.end())
    {
        TBOX_ERROR("AdvDiffConvectiveOperatorManager::allocateOperator():\n"
                   << "  unrecognized operator type: " << operator_type << "\n");
    }
    return (it->second)(operator_object_name, Q_var, input_db, difference_form, bc_coefs);
} // allocateOperator

void
AdvDiffConvectiveOperatorManager::registerOperatorFactoryFunction(const std::string& operator_type,
                                                                  OperatorMaker operator_maker)
{
    if (d_operator_maker_map.find(operator_type) != d_operator_maker_map.end())
    {
        pout << "AdvDiffConvectiveOperatorManager::registerOperatorFactoryFunction("
                "):\n"
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
    registerOperatorFactoryFunction(CUI, AdvDiffCUIConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(PPM, AdvDiffPPMConvectiveOperator::allocate_operator);
    registerOperatorFactoryFunction(WAVE_PROP, AdvDiffWavePropConvectiveOperator::allocate_operator);
    return;
} // AdvDiffConvectiveOperatorManager

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
