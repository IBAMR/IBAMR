// Filename: IBLagrangianForceStrategy.C
// Last modified: <27.Jun.2010 15:29:45 griffith@griffith-macbook-pro.local>
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

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

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
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    // intentionally blank
    return;
}// initializeLevelData

void
IBLagrangianForceStrategy::computeLagrangianForceAndTorque(
    Pointer<LNodeLevelData> F_data,
    Pointer<LNodeLevelData> N_data,
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> D_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceAndTorque():\n"
               << "  not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForceAndTorque

void
IBLagrangianForceStrategy::computeLagrangianForce(
    Pointer<LNodeLevelData> F_data,
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForce():\n"
               << "  not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForce

void
IBLagrangianForceStrategy::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobianNonzeroStructure():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBLagrangianForceStrategy::computeLagrangianForceJacobian(
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
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobian():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForceJacobian

double
IBLagrangianForceStrategy::computeLagrangianEnergy(
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianEnergy():\n"
               << "  potential energy functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
    return 0.0;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBLagrangianForceStrategy>;

//////////////////////////////////////////////////////////////////////////////
