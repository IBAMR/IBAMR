// Filename: IBLagrangianForceStrategySet.C
// Last modified: <29.Jul.2008 15:38:02 griffith@box230.cims.nyu.edu>
// Created on 04 April 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBLagrangianForceStrategySet.h"

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

IBLagrangianForceStrategySet::~IBLagrangianForceStrategySet()
{
    // intentionally blank
    return;
}// ~IBLagrangianForceStrategySet

void
IBLagrangianForceStrategySet::setTimeInterval(
    const double current_time,
    const double new_time)
{
    for (std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->setTimeInterval(current_time, new_time);
    }
    return;
}// setTimeInterval

void
IBLagrangianForceStrategySet::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    for (std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->initializeLevelData(
            hierarchy, level_number, init_data_time, initial_time, lag_manager);
    }
    return;
}// initializeLevelData

void
IBLagrangianForceStrategySet::computeLagrangianForce(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    for (std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->computeLagrangianForce(
            F_data, X_data, hierarchy, level_number, data_time, lag_manager);
    }
    return;
}// computeLagrangianForce

void
IBLagrangianForceStrategySet::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    for (std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->computeLagrangianForceJacobianNonzeroStructure(
            d_nnz, o_nnz, hierarchy, level_number, data_time, lag_manager);
    }
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBLagrangianForceStrategySet::computeLagrangianForceJacobian(
    Mat& J_mat,
    MatAssemblyType assembly_type,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    for (std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->computeLagrangianForceJacobian(
            J_mat, MAT_FLUSH_ASSEMBLY, X_data, hierarchy, level_number, data_time, lag_manager);
    }
    if (assembly_type != MAT_FLUSH_ASSEMBLY)
    {
        int ierr;
        ierr = MatAssemblyBegin(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
    }
    return;
}// computeLagrangianForceJacobian

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBLagrangianForceStrategySet>;

//////////////////////////////////////////////////////////////////////////////
