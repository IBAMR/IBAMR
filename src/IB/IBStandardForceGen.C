// Filename: IBStandardForceGen.C
// Last modified: <02.Nov.2009 11:12:07 griffith@griffith-macbook-pro.local>
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

// IBTK INCLUDES
#include <ibtk/LNodeIndexData2.h>
#include <ibtk/LNodeLevelData.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceGen::IBStandardForceGen(
    SAMRAI::tbox::Pointer<IBSpringForceGen> spring_force_gen,
    SAMRAI::tbox::Pointer<IBBeamForceGen> beam_force_gen,
    SAMRAI::tbox::Pointer<IBTargetPointForceGen> target_point_force_gen)
    : d_force_strategy_set(NULL),
      d_X_orig_vec(),
      d_shift_vec()
{
    std::vector<SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> > strategy_set;

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
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
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
        SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data = lag_manager->getLNodeLevelData("X", level_number);
        ierr = VecDuplicate(X_data->getGlobalVec(), &d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
        ierr = VecDuplicate(X_data->getGlobalVec(), & d_shift_vec[level_number]);  IBTK_CHKERRQ(ierr);

        ierr = VecSet(d_X_orig_vec[level_number], 0.0);  IBTK_CHKERRQ(ierr);
        ierr = VecSet( d_shift_vec[level_number], 0.0);  IBTK_CHKERRQ(ierr);

        // Compute the periodic shifts (if any).
        double* shift_arr;
        ierr = VecGetArray(d_shift_vec[level_number], &shift_arr);  IBTK_CHKERRQ(ierr);
        const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
            const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
            const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data = patch->getPatchData(lag_node_index_idx);
            for (IBTK::LNodeIndexData2::Iterator it(patch_box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
                const IBTK::LNodeIndexSet& node_set = (*idx_data)(i);
                for (IBTK::LNodeIndexSet::const_iterator n = node_set.begin(); n != node_set.end(); ++n)
                {
                    const IBTK::LNodeIndexSet::value_type& node_idx = *n;
                    if (node_idx->getPeriodicOffset() != SAMRAI::hier::IntVector<NDIM>(0))
                    {
                        const int petsc_local_idx = node_idx->getLocalPETScIndex();
                        const std::vector<double>& periodic_displacement = node_idx->getPeriodicDisplacement();
                        for (int d = 0; d < NDIM; ++d)
                        {
                            shift_arr[NDIM*petsc_local_idx+d] = -periodic_displacement[d];
                        }
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
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
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
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
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
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const double U_coef,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> U_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    int ierr;
    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    ierr = VecWAXPY(X_data->getGlobalVec(), 1.0, d_shift_vec[level_number], d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);

    d_force_strategy_set->computeLagrangianForceJacobian(J_mat, assembly_type, X_coef, X_data, U_coef, U_data, hierarchy, level_number, data_time, lag_manager);

    ierr = VecSwap(X_data->getGlobalVec(), d_X_orig_vec[level_number]);  IBTK_CHKERRQ(ierr);
    return;
}// computeLagrangianForceJacobian

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardForceGen>;

//////////////////////////////////////////////////////////////////////////////
