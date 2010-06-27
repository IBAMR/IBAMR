// Filename: IBTargetPointForceGen.C
// Last modified: <27.Jun.2010 16:03:06 griffith@griffith-macbook-pro.local>
// Created on 21 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBTargetPointForceGen.h"

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
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/namespaces.h>

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

// BLITZ++ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_compute_lagrangian_force;
static Pointer<Timer> t_compute_lagrangian_force_jacobian;
static Pointer<Timer> t_compute_lagrangian_force_jacobian_nonzero_structure;
static Pointer<Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBTargetPointForceGen::IBTargetPointForceGen(
    Pointer<Database> input_db)
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = TimerManager::getManager()->
            getTimer("IBAMR::IBTargetPointForceGen::computeLagrangianForce()");
        t_compute_lagrangian_force_jacobian = TimerManager::getManager()->
            getTimer("IBAMR::IBTargetPointForceGen::computeLagrangianForceJacobian()");
        t_compute_lagrangian_force_jacobian_nonzero_structure = TimerManager::getManager()->
            getTimer("IBAMR::IBTargetPointForceGen::computeLagrangianForceJacobianNonzeroStructure()");
        t_initialize_level_data = TimerManager::getManager()->
            getTimer("IBAMR::IBTargetPointForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// IBTargetPointForceGen

IBTargetPointForceGen::~IBTargetPointForceGen()
{
    // intentionally blank
    return;
}// ~IBTargetPointForceGen

void
IBTargetPointForceGen::computeLagrangianForce(
    Pointer<LNodeLevelData> F_data,
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
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
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    if (!grid_geom->getDomainIsSingleBox())
    {
        TBOX_ERROR("IBTargetPointForceGen::computeLagrangianForce():\n"
                   << "  physical domain must be a single box.\n");
    }

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Compute the penalty force associated with the Lagrangian target points.
    static double max_displacement = 0.0;
    double max_config_displacement = 0.0;
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const LNodeIndex& node_idx = *it;
            Pointer<IBTargetPointForceSpec> force_spec = node_idx.getStashData<IBTargetPointForceSpec>();
            if (!force_spec.isNull())
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                const int& mastr_idx = node_idx.getLagrangianIndex();
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const double& kappa_target = force_spec->getStiffness();
                const double&   eta_target = force_spec->getDamping();
                if (!MathUtilities<double>::equalEps(kappa_target,0.0) ||
                    !MathUtilities<double>::equalEps(  eta_target,0.0))
                {
                    const int& petsc_idx = node_idx.getLocalPETScIndex();
                    const double* const X = &X_arr[NDIM*petsc_idx];
                    const double* const U = &U_arr[NDIM*petsc_idx];
                    const std::vector<double>& X_target = force_spec->getTargetPointPosition();

                    double* const F = &F_arr[NDIM*petsc_idx];
                    double displacement = 0.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F[d] += kappa_target*(X_target[d] - X[d]) - eta_target*(U[d]);
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

    max_config_displacement = SAMRAI_MPI::maxReduction(max_config_displacement);
    if (max_config_displacement > max_displacement)
    {
        max_displacement = max_config_displacement;
    }

    if (!MathUtilities<double>::equalEps(max_config_displacement,0.0))
    {
        plog << "IBTargetPointForceGen::computeLagrangianForce():" << std::endl;
        plog << "  maximum target point displacement [present configuration] = " << max_config_displacement << std::endl;
        plog << "  maximum target point displacement [entire simulation] = " << max_displacement << std::endl;
    }

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

void
IBTargetPointForceGen::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force_jacobian_nonzero_structure->start();

    // Get the patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the PETSc indices of the target point nodes and the
    // corresponding target force spring constant.
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const LNodeIndex& node_idx = *it;
            Pointer<IBTargetPointForceSpec> force_spec = node_idx.getStashData<IBTargetPointForceSpec>();
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
IBTargetPointForceGen::computeLagrangianForceJacobian(
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
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const LNodeIndex& node_idx = *it;
            Pointer<IBTargetPointForceSpec> force_spec = node_idx.getStashData<IBTargetPointForceSpec>();
            if (!force_spec.isNull())
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                const int& mastr_idx = node_idx.getLagrangianIndex();
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const int& local_petsc_idx = node_idx.getLocalPETScIndex();
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
    blitz::Array<double,2> dF_dX(NDIM,NDIM);  dF_dX = 0.0;
    const int num_local_idxs = global_petsc_idxs.size();
    for (int k = 0; k < num_local_idxs; ++k)
    {
        const int& global_petsc_idx = global_petsc_idxs[k];
        const double& spring_stiffness = spring_stiffnesses[k];
        const double& damping_coefficient = damping_coefficients[k];
        for (int alpha = 0; alpha < NDIM; ++alpha)
        {
            dF_dX(alpha,alpha) = -X_coef*spring_stiffness-U_coef*damping_coefficient;
        }
        ierr = MatSetValuesBlocked(J_mat,1,&global_petsc_idx,1,&global_petsc_idx,dF_dX.data(),ADD_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, assembly_type);  IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_jacobian->stop();
    return;
}// computeLagrangianForceJacobian

double
IBTargetPointForceGen::computeLagrangianEnergy(
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    TBOX_WARNING("IBTargetPointForceGen::computeLagrangianEnergy():\n"
                 << "  unimplemented; returning 0.0." << std::endl);
    return 0.0;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBTargetPointForceGen::getFromInput(
    Pointer<Database> db)
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
template class Pointer<IBAMR::IBTargetPointForceGen>;

//////////////////////////////////////////////////////////////////////////////
