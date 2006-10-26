// Filename: IBLagrangianSourceStrategy.C
// Created on 18 Jun 2005 by Boyce Griffith (boyce@bigboy.verizon.net)
// Last modified: <25.Oct.2006 18:30:00 boyce@bigboy.nyconnect.com>

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
IBLagrangianSourceStrategy::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    const LDataManager* const lag_manager)
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
template class SAMRAI::tbox::Pointer<IBAMR::IBLagrangianSourceStrategy>;

//////////////////////////////////////////////////////////////////////////////
