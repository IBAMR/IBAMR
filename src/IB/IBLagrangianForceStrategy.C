// Filename: IBLagrangianForceStrategy.C
// Created on 03 May 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)
// Last modified: <03.Oct.2006 10:35:33 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "IBLagrangianForceStrategy.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBLagrangianForceStrategy::IBLagrangianForceStrategy()
{
    // intentionally blank
    return;
}// IBLagrangianForceStrategy

IBLagrangianForceStrategy::~IBLagrangianForceStrategy()
{
    // intentionally blank
    return;
}// ~IBLagrangianForceStrategy

void
IBLagrangianForceStrategy::initializeLevelData(
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
template class SAMRAI::tbox::Pointer<IBAMR::IBLagrangianForceStrategy>;

//////////////////////////////////////////////////////////////////////////////
