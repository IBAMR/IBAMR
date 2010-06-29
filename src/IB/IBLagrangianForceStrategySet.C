// Filename: IBLagrangianForceStrategySet.C
// Last modified: <27.Jun.2010 15:29:49 griffith@griffith-macbook-pro.local>
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

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

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
    for (std::vector<Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->setTimeInterval(current_time, new_time);
    }
    return;
}// setTimeInterval

void
IBLagrangianForceStrategySet::initializeLevelData(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    for (std::vector<Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->initializeLevelData(
            hierarchy, level_number, init_data_time, initial_time, lag_manager);
    }
    return;
}// initializeLevelData

void
IBLagrangianForceStrategySet::computeLagrangianForce(
    Pointer<LNodeLevelData> F_data,
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    for (std::vector<Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->computeLagrangianForce(
            F_data, X_data, U_data, hierarchy, level_number, data_time, lag_manager);
    }
    return;
}// computeLagrangianForce

void
IBLagrangianForceStrategySet::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    for (std::vector<Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
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
    const double X_coef,
    Pointer<LNodeLevelData> X_data,
    const double U_coef,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    for (std::vector<Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        (*it)->computeLagrangianForceJacobian(
            J_mat, MAT_FLUSH_ASSEMBLY, X_coef, X_data, U_coef, U_data, hierarchy, level_number, data_time, lag_manager);
    }
    if (assembly_type != MAT_FLUSH_ASSEMBLY)
    {
        int ierr;
        ierr = MatAssemblyBegin(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
    }
    return;
}// computeLagrangianForceJacobian

double
IBLagrangianForceStrategySet::computeLagrangianEnergy(
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    double ret_val = 0.0;
    for (std::vector<Pointer<IBLagrangianForceStrategy> >::iterator it = d_strategy_set.begin();
         it != d_strategy_set.end(); ++it)
    {
        ret_val += (*it)->computeLagrangianEnergy(X_data, U_data, hierarchy, level_number, data_time, lag_manager);
    }
    return ret_val;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBLagrangianForceStrategySet>;

//////////////////////////////////////////////////////////////////////////////
