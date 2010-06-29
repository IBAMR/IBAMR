// Filename: IBLagrangianSourceStrategy.C
// Last modified: <27.Jun.2010 15:29:53 griffith@griffith-macbook-pro.local>
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)

#include "IBLagrangianSourceStrategy.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBLagrangianSourceStrategy::IBLagrangianSourceStrategy()
{
    // intentionally blank
    return;
}// IBLagrangianSourceStrategy

IBLagrangianSourceStrategy::~IBLagrangianSourceStrategy()
{
    // intentionally blank
    return;
}// ~IBLagrangianSourceStrategy

void
IBLagrangianSourceStrategy::setTimeInterval(
    const double current_time,
    const double new_time)
{
    // intentionally blank
    return;
}// setTimeInterval

void
IBLagrangianSourceStrategy::initializeLevelData(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // intentionally blank
    return;
}// initializeLevelData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBLagrangianSourceStrategy>;

//////////////////////////////////////////////////////////////////////////////
