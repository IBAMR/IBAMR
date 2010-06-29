// Filename: IBStandardForceGen.C
// Last modified: <27.Jun.2010 15:30:56 griffith@griffith-macbook-pro.local>
// Created on 03 May 2005 by Boyce Griffith (boyce@mstu1.cims.nyu.edu)

#include "IBStandardForceGen.h"

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

// IBTK INCLUDES
#include <ibtk/LNodeIndexData.h>
#include <ibtk/LNodeLevelData.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceGen::IBStandardForceGen(
    Pointer<IBSpringForceGen> spring_force_gen,
    Pointer<IBBeamForceGen> beam_force_gen,
    Pointer<IBTargetPointForceGen> target_point_force_gen)
    : d_force_strategy_set(NULL),
      d_X_orig_vec(),
      d_shift_vec()
{
    std::vector<Pointer<IBLagrangianForceStrategy> > strategy_set;

    if (spring_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen():\n"
                     << "  spring forces disabled." << std::endl);
    }
    else
    {
        strategy_set.push_back(spring_force_gen);
    }

    if (beam_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen():\n"
                     << "  beam forces disabled." << std::endl);
    }
    else
    {
        strategy_set.push_back(beam_force_gen);
    }

    if (target_point_force_gen.isNull())
    {
        TBOX_WARNING("IBStandardForceGen():\n"
                     << "  target point forces disabled." << std::endl);
    }
    else
    {
        strategy_set.push_back(target_point_force_gen);
    }

    d_force_strategy_set = new IBLagrangianForceStrategySet(strategy_set.begin(), strategy_set.end());
    return;
}// IBStandardForceGen

IBStandardForceGen::~IBStandardForceGen()
{
    int ierr;
    for (std::vector<Vec>::iterator it = d_X_orig_vec.begin(); it != d_X_orig_vec.end(); ++it)
    {
        if (*it != PETSC_NULL)
        {
            ierr = VecDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    d_X_orig_vec.clear();
    for (std::vector<Vec>::iterator it = d_shift_vec.begin(); it != d_shift_vec.end(); ++it)
    {
        if (*it != PETSC_NULL)
        {
            ierr = VecDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    d_shift_vec.clear();
    return;
}// ~IBStandardForceGen

void
IBStandardForceGen::initializeLevelData(
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    d_X_orig_vec.resize(std::max(d_X_orig_vec.size(),size_t(level_number+1)),PETSC_NULL);
    d_shift_vec .resize(std::max( d_shift_vec.size(),size_t(level_number+1)),PETSC_NULL);
    if (lag_manager->levelContainsLagrangianData(level_number))
    {
        int ierr;

        // Clean up existing vectors.
        if (d_X_orig_vec[level_number] != PETSC_NULL)
        {
            ierr = VecDestroy(d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
        }
        if (d_shift_vec[level_number] != PETSC_NULL)
        {
            ierr = VecDestroy(d_shift_vec[level_number]);  IBTK_CHKERRQ(ierr);
        }

        // Create a duplicate of the position vector.
        Pointer<LNodeLevelData> X_data = lag_manager->getLNodeLevelData("X", level_number);
        ierr = VecDuplicate(X_data->getGlobalVec(), &d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(X_data->getGlobalVec(), & d_shift_vec[level_number]);  IBTK_CHKERRQ(ierr);

        ierr = VecSet(d_X_orig_vec[level_number], 0.0);  IBTK_CHKERRQ(ierr);
        ierr = VecSet( d_shift_vec[level_number], 0.0);  IBTK_CHKERRQ(ierr);

        // Compute the periodic shifts (if any).
        double* shift_arr;
        ierr = VecGetArray(d_shift_vec[level_number], &shift_arr);  IBTK_CHKERRQ(ierr);
        const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();
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
                if (node_idx.getPeriodicOffset() != IntVector<NDIM>(0))
                {
                    const int petsc_local_idx = node_idx.getLocalPETScIndex();
                    const std::vector<double>& periodic_displacement = node_idx.getPeriodicDisplacement();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        shift_arr[NDIM*petsc_local_idx+d] = -periodic_displacement[d];
                    }
                }
            }
        }
        ierr = VecRestoreArray(d_shift_vec[level_number], &shift_arr);  IBTK_CHKERRQ(ierr);
    }
    d_force_strategy_set->initializeLevelData(hierarchy, level_number, init_data_time, initial_time, lag_manager);
    return;
}// initializeLevelData

void
IBStandardForceGen::computeLagrangianForce(
    Pointer<LNodeLevelData> F_data,
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    int ierr;
    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    ierr = VecWAXPY(X_data->getGlobalVec(), 1.0, d_shift_vec[level_number], d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);

    d_force_strategy_set->computeLagrangianForce(F_data, X_data, U_data, hierarchy, level_number, data_time, lag_manager);

    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    return;
}// computeLagrangianForce

void
IBStandardForceGen::computeLagrangianForceJacobianNonzeroStructure(
    std::vector<int>& d_nnz,
    std::vector<int>& o_nnz,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;
    d_force_strategy_set->computeLagrangianForceJacobianNonzeroStructure(d_nnz, o_nnz, hierarchy, level_number, data_time, lag_manager);
    return;
}// computeLagrangianForceJacobianNonzeroStructure

void
IBStandardForceGen::computeLagrangianForceJacobian(
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

    int ierr;
    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    ierr = VecWAXPY(X_data->getGlobalVec(), 1.0, d_shift_vec[level_number], d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);

    d_force_strategy_set->computeLagrangianForceJacobian(J_mat, assembly_type, X_coef, X_data, U_coef, U_data, hierarchy, level_number, data_time, lag_manager);

    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    return;
}// computeLagrangianForceJacobian

double
IBStandardForceGen::computeLagrangianEnergy(
    Pointer<LNodeLevelData> X_data,
    Pointer<LNodeLevelData> U_data,
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return 0.0;

    int ierr;
    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    ierr = VecWAXPY(X_data->getGlobalVec(), 1.0, d_shift_vec[level_number], d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);

    double ret_val = d_force_strategy_set->computeLagrangianEnergy(X_data, U_data, hierarchy, level_number, data_time, lag_manager);

    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    return ret_val;
}// computeLagrangianEnergy

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBStandardForceGen>;

//////////////////////////////////////////////////////////////////////////////
