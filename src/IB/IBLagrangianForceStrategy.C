// Filename: IBLagrangianForceStrategy.C
// Last modified: <22.Jun.2010 17:45:29 griffith@boyce-griffiths-mac-pro.local>
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
IBLagrangianForceStrategy::computeLagrangianForceAndTorque(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> N_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> D_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceAndTorque():\n"
               << "  not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForceAndTorque

void
IBLagrangianForceStrategy::computeLagrangianForce(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForce():\n"
               << "  not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForce

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
}// computeLagrangianForceJacobianNonzeroStructure

void
IBLagrangianForceStrategy::computeLagrangianForceJacobian(
    Mat& J_mat,
    MatAssemblyType assembly_type,
    const double X_coef,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const double U_coef,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    TBOX_ERROR("IBLagrangianForceStrategy::computeLagrangianForceJacobian():\n"
               << "  Jacobian functionality not implemented for this IBLagrangianForceStrategy." << std::endl);
    return;
}// computeLagrangianForceJacobian

double
IBLagrangianForceStrategy::computeLagrangianEnergy(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
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
template class SAMRAI::tbox::Pointer<IBAMR::IBLagrangianForceStrategy>;

//////////////////////////////////////////////////////////////////////////////
