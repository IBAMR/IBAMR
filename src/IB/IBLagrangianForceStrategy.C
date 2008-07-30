// Filename: IBLagrangianForceStrategy.C
// Last modified: <29.Jul.2008 15:37:46 griffith@box230.cims.nyu.edu>
// Created on 03 May 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

#include "IBLagrangianForceStrategy.h"

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
IBLagrangianForceStrategy::setTimeInterval(
    const double current_time,
    const double new_time)
{
    // intentionally blank
    return;
}// setTimeInterval

void
IBLagrangianForceStrategy::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    // intentionally blank
    return;
}// initializeLevelData

void
IBLagrangianForceStrategy::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobianNonzeroStructure():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// initializeLevelData

void
IBLagrangianForceStrategy::computeLagrangianForceJacobian(
    Mat& J_mat,
    MatAssemblyType assembly_type,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobian():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForceJacobian

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBLagrangianForceStrategy>;

//////////////////////////////////////////////////////////////////////////////
