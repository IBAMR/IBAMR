// Filename: IBMovingTargetPointForceGen.C
// Last modified: <01.Mar.2010 15:54:35 griffith@boyce-griffiths-mac-pro.local>
// Created on 14 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

#include "IBMovingTargetPointForceGen.h"

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
#include <ibamr/IBMovingTargetPointForceSpec.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LNodeIndexData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/MathUtilities.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force_jacobian;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force_jacobian_nonzero_structure;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMovingTargetPointForceGen::IBMovingTargetPointForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_spec_fcn_map()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBMovingTargetPointForceGen::computeLagrangianForce()");
        t_compute_lagrangian_force_jacobian = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBMovingTargetPointForceGen::computeLagrangianForceJacobian()");
        t_compute_lagrangian_force_jacobian_nonzero_structure = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBMovingTargetPointForceGen::computeLagrangianForceJacobianNonzeroStructure()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBMovingTargetPointForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// IBMovingTargetPointForceGen

IBMovingTargetPointForceGen::~IBMovingTargetPointForceGen()
{
    // intentionally blank
    return;
}// ~IBMovingTargetPointForceGen

void
IBMovingTargetPointForceGen::registerPositionAndVelocityFunction(
    const int spec_fcn_index,
    void (*spec_fcn)(double X_target[NDIM], double U_target[NDIM], const double X[NDIM], const double U[NDIM], const double& data_time, const int& lag_mastr_idx))
{
    d_spec_fcn_map[spec_fcn_index] = spec_fcn;
    return;
}// registerPositionAndVelocityFunction

void
IBMovingTargetPointForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force->start();

    int ierr;

    // Extract the local arrays.
    Vec F_vec = F_data->getGlobalVec();
    double* F_arr;
    ierr = VecGetArray(F_vec, &F_arr);  IBTK_CHKERRQ(ierr);

    Vec X_vec = X_data->getGlobalVec();
    double* X_arr;
    ierr = VecGetArray(X_vec, &X_arr);  IBTK_CHKERRQ(ierr);

    Vec U_vec = U_data->getGlobalVec();
    double* U_arr;
    ierr = VecGetArray(U_vec, &U_arr);  IBTK_CHKERRQ(ierr);

    // Get the grid geometry object and determine the extents of the physical
    // domain.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    if (!grid_geom->getDomainIsSingleBox())
    {
        TBOX_ERROR("IBMovingTargetPointForceGen::computeLagrangianForce():\n"
                   << "  physical domain must be a single box.\n");
    }

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Compute the penalty force associated with the Lagrangian target points.
    double X_target[NDIM], U_target[NDIM];
    static double max_displacement = 0.0;
    double max_config_displacement = 0.0;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const IBTK::LNodeIndex& node_idx = *it;
            SAMRAI::tbox::Pointer<IBMovingTargetPointForceSpec> force_spec = node_idx.getStashData<IBMovingTargetPointForceSpec>();
            if (!force_spec.isNull())
            {
                const int& mastr_idx = node_idx.getLagrangianIndex();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const double& kappa_target = force_spec->getStiffness();
                const double& eta_target = force_spec->getDamping();
                const int& spec_fcn_idx = force_spec->getPositionAndVelocityFunctionIndex();
                if (!SAMRAI::tbox::MathUtilities<double>::equalEps(kappa_target,0.0))
                {
                    const int& petsc_idx = node_idx.getLocalPETScIndex();
                    const double* const X = &X_arr[NDIM*petsc_idx];
                    const double* const U = &U_arr[NDIM*petsc_idx];
                    d_spec_fcn_map[spec_fcn_idx](X_target, U_target, X, U, data_time, mastr_idx);

                    double* const F = &F_arr[NDIM*petsc_idx];
                    double displacement = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F[d] += kappa_target*(X_target[d] - X[d]) + eta_target*(U_target[d] - U[d]);
                        displacement += pow(X_target[d] - X[d],2.0);
                    }
                    displacement = sqrt(displacement);
                    if (displacement > max_config_displacement)
                    {
                        max_config_displacement = displacement;
                    }
                }
            }
        }
    }

    ierr = VecRestoreArray(F_vec, &F_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(X_vec, &X_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(U_vec, &U_arr);  IBTK_CHKERRQ(ierr);

    max_config_displacement = SAMRAI::tbox::SAMRAI_MPI::maxReduction(max_config_displacement);
    if (max_config_displacement > max_displacement)
    {
        max_displacement = max_config_displacement;
    }

    if (!SAMRAI::tbox::MathUtilities<double>::equalEps(max_config_displacement,0.0))
    {
        SAMRAI::tbox::plog << "IBMovingTargetPointForceGen::computeLagrangianForce():" << std::endl;
        SAMRAI::tbox::plog << "  maximum target point displacement [present configuration] = " << max_config_displacement << std::endl;
        SAMRAI::tbox::plog << "  maximum target point displacement [entire simulation] = " << max_displacement << std::endl;
    }

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

void
IBMovingTargetPointForceGen::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force_jacobian_nonzero_structure->start();

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the PETSc indices of the target point nodes and the
    // corresponding target force spring constant.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const IBTK::LNodeIndex& node_idx = *it;
            SAMRAI::tbox::Pointer<IBMovingTargetPointForceSpec> force_spec = node_idx.getStashData<IBMovingTargetPointForceSpec>();
            if (!force_spec.isNull())
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                const int& mastr_idx = node_idx.getLagrangianIndex();
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const int& local_petsc_idx = node_idx.getLocalPETScIndex();
                ++d_nnz[local_petsc_idx];
            }
        }
    }

    t_compute_lagrangian_force_jacobian_nonzero_structure->stop();
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBMovingTargetPointForceGen::computeLagrangianForceJacobian(
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
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force_jacobian->start();

    int ierr;

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the PETSc indices of the target point nodes and the
    // corresponding target force spring constants.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_number);
    std::vector<int> global_petsc_idxs;
    std::vector<double> spring_stiffnesses, damping_coefficients;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const IBTK::LNodeIndex& node_idx = *it;
            SAMRAI::tbox::Pointer<IBMovingTargetPointForceSpec> force_spec = node_idx.getStashData<IBMovingTargetPointForceSpec>();
            if (!force_spec.isNull())
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                const int& mastr_idx = node_idx.getLagrangianIndex();
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const int& local_petsc_idx = node_idx.getLocalPETScIndex()+global_node_offset;
                const int global_petsc_idx = local_petsc_idx+global_node_offset;
                global_petsc_idxs.push_back(global_petsc_idx);

                const double& spring_stiffness = force_spec->getStiffness();
                spring_stiffnesses.push_back(spring_stiffness);

                const double& damping_coefficient = force_spec->getDamping();
                damping_coefficients.push_back(damping_coefficient);
            }
        }
    }

    // Compute the elements of the Jacobian matrix.
    std::vector<double> dF_dX(NDIM*NDIM,0.0);
    const int num_local_idxs = global_petsc_idxs.size();
    for (int k = 0; k < num_local_idxs; ++k)
    {
        const int& global_petsc_idx = global_petsc_idxs[k];
        const double& spring_stiffness = spring_stiffnesses[k];
        const double& damping_coefficient = damping_coefficients[k];
        for (int alpha = 0; alpha < NDIM; ++alpha)
        {
            dF_dX[alpha+alpha*NDIM] = -X_coef*spring_stiffness-U_coef*damping_coefficient;
        }
        ierr = MatSetValuesBlocked(J_mat,1,&global_petsc_idx,1,&global_petsc_idx,&dF_dX[0],ADD_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_jacobian->stop();
    return;
}// computeLagrangianForceJacobian

double
IBMovingTargetPointForceGen::computeLagrangianEnergy(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    TBOX_WARNING("IBMovingTargetPointForceGen::computeLagrangianEnergy():\n"
               << "  unimplemented; returning 0.0." << std::endl);
    return 0.0;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBMovingTargetPointForceGen::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    if (!db.isNull())
    {
        // intentionally blank
    }
    return;
}// getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBMovingTargetPointForceGen>;

//////////////////////////////////////////////////////////////////////////////
