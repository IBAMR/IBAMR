// Filename: SpringForceGen.C
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <18.Jan.2007 17:29:17 boyce@bigboy.nyconnect.com>

#include "SpringForceGen.h"

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
#include <ibamr/LNodeIndexData.h>
#include <ibamr/SpringForceSpec.h>

// STOOLS INCLUDES
#include <stools/PETSC_SAMRAI_ERROR.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <numeric>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data_0;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data_1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

SpringForceGen::SpringForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_L_mats(),
      d_local_src_ids(),
      d_stiffnesses(),
      d_rest_lengths(),
      d_level_initialized()
{
    // Initialize object with data read from the input database.
    if (!input_db.isNull())
    {
        // intentionally blank
    }

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::SpringForceGen::computeLagrangianForce()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::SpringForceGen::initializeLevelData()");
        t_initialize_level_data_0 = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::SpringForceGen::initializeLevelData()[MatSetValues]");
        t_initialize_level_data_1 = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::SpringForceGen::initializeLevelData()[MatAssembly]");
        timers_need_init = false;
    }
    return;
}// SpringForceGen

SpringForceGen::~SpringForceGen()
{
    for (std::vector<Mat>::iterator it = d_L_mats.begin(); it != d_L_mats.end(); ++it)
    {
        if (*it) MatDestroy(*it);
    }
    return;
}// ~SpringForceGen

void
SpringForceGen::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    const LDataManager* const lag_manager)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
#endif
    (void) init_data_time;
    (void) initial_time;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually
    // maintained for separate levels of the patch hierarchy.
    const int level_num = level->getLevelNumber();
    const int new_size = SAMRAI::tbox::Utilities::imax(
        level_num+1, d_level_initialized.size());

    d_L_mats.resize(new_size);
    d_local_src_ids.resize(new_size);
    d_stiffnesses.resize(new_size);
    d_rest_lengths.resize(new_size);
    d_level_initialized.resize(new_size, false);

    std::vector<int>& local_src_ids = d_local_src_ids[level_num];
    std::vector<double>& stiffnesses = d_stiffnesses[level_num];
    std::vector<double>& rest_lengths = d_rest_lengths[level_num];

    local_src_ids.clear();
    stiffnesses.clear();
    rest_lengths.clear();

    int ierr;
    if (d_L_mats[level_num])
    {
        ierr = MatDestroy(d_L_mats[level_num]);
        PETSC_SAMRAI_ERROR(ierr);
    }

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->
        getLNodeIndexPatchDescriptorIndex();

    // Determine the "source" and "destination" Lagrangian indices.
    std::vector<int> src_ids, dst_ids, num_dst_ids;
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<LNodeIndexData> idx_data =
            patch->getPatchData(lag_node_index_idx);

        for (LNodeIndexData::Iterator it(*idx_data); it; it++)
        {
            if (patch_box.contains(it.getIndex()))
            {
                const LNodeIndexSet& node_set = *it;
                for (LNodeIndexSet::const_iterator n = node_set.begin();
                     n != node_set.end(); ++n)
                {
                    const LNodeIndexSet::value_type& node_idx = *n;
                    const int& lag_idx = node_idx->getLagrangianIndex();
                    SAMRAI::tbox::Pointer<SpringForceSpec> force_spec =
                        node_idx->getStashData()[0];
                    if (!force_spec.isNull())
                    {
                        const std::vector<int>& dst_idxs =
                            force_spec->getDestinationNodeIndices();
                        const std::vector<double>& stiff =
                            force_spec->getStiffnesses();
                        const std::vector<double>& rest =
                            force_spec->getRestingLengths();
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(dst_idxs.size() == stiff.size());
                        assert(dst_idxs.size() == rest.size());
#endif
                        if (!dst_idxs.empty())
                        {
                            src_ids.push_back(lag_idx);
                            dst_ids.insert(dst_ids.end(),
                                           dst_idxs.begin(), dst_idxs.end());
                            num_dst_ids.push_back(static_cast<int>(dst_idxs.size()));

                            local_src_ids.insert(local_src_ids.end(),
                                                 dst_idxs.size(),
                                                 node_idx->getLocalPETScIndex());
                            stiffnesses.insert(stiffnesses.end(),
                                               stiff.begin(), stiff.end());
                            rest_lengths.insert(rest_lengths.end(),
                                                rest.begin(), rest.end());
                        }
                    }
                }
            }
        }
    }

    // Map the Lagrangian source/destination indices to the PETSc
    // indices corresponding to the present data distribution.
    lag_manager->mapLagrangianToPETSc(src_ids, level_num);
    lag_manager->mapLagrangianToPETSc(dst_ids, level_num);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_num);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_num);

    // Determine the non-zero structure for the matrix.
    const int local_sz = static_cast<int>(dst_ids.size());
    std::vector<int> d_nz(local_sz,1), o_nz(local_sz,0);

    for (int k = 0; k < local_sz; ++k)
    {
        const int& dst_id = dst_ids[k];
        if (dst_id >= global_node_offset &&
            dst_id <  global_node_offset+num_local_nodes)
        {
            ++d_nz[k]; // a "local"    destination index
        }
        else
        {
            ++o_nz[k]; // a "nonlocal" destination index
        }
    }

    // Create a new MPI block AIJ matrix and set the values of the
    // non-zero entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, &d_nz[0],
                            PETSC_DEFAULT, &o_nz[0],
                            &d_L_mats[level_num]);  PETSC_SAMRAI_ERROR(ierr);

    std::vector<double> src_vals(NDIM*NDIM,0.0);
    std::vector<double> dst_vals(NDIM*NDIM,0.0);
    for (int d = 0; d < NDIM; ++d)
    {
        src_vals[d+d*NDIM] = -1.0;
        dst_vals[d+d*NDIM] = +1.0;
    }

    int i_offset;
    ierr = MatGetOwnershipRange(d_L_mats[level_num], &i_offset, PETSC_NULL);
    PETSC_SAMRAI_ERROR(ierr);
    i_offset /= NDIM;

    for (int src = 0, dst = -1; src < static_cast<int>(src_ids.size()); ++src)
    {
        int j_src = src_ids[src];
        for (int n_dst = 0; n_dst < num_dst_ids[src]; ++n_dst)
        {
            int i = i_offset + ++dst;
            int j_dst = dst_ids[dst];

            ierr = MatSetValuesBlocked(d_L_mats[level_num],
                                       1,&i,1,&j_src,&(src_vals[0]),
                                       INSERT_VALUES);
            PETSC_SAMRAI_ERROR(ierr);

            ierr = MatSetValuesBlocked(d_L_mats[level_num],
                                       1,&i,1,&j_dst,&(dst_vals[0]),
                                       INSERT_VALUES);
            PETSC_SAMRAI_ERROR(ierr);
        }
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(d_L_mats[level_num], MAT_FINAL_ASSEMBLY);
    PETSC_SAMRAI_ERROR(ierr);
    ierr = MatAssemblyEnd(d_L_mats[level_num], MAT_FINAL_ASSEMBLY);
    PETSC_SAMRAI_ERROR(ierr);

    // Indicate that the level data has been properly initialized.
    d_level_initialized[level_num] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
SpringForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    const LDataManager* const lag_manager)
{
    t_compute_lagrangian_force->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(level_number < static_cast<int>(d_level_initialized.size()));
    assert(d_level_initialized[level_number]);
#endif
    (void) hierarchy;
    (void) lag_manager;

    int ierr;

    // Create an appropriately sized temporary vector to store the
    // node displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_L_mats[level_number], &i_start, &i_stop);
    PETSC_SAMRAI_ERROR(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          PETSC_SAMRAI_ERROR(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_L_mats[level_number], X_data->getGlobalVec(), D_vec);
    PETSC_SAMRAI_ERROR(ierr);

    // Compute the forces acting on each node.
    Vec F_vec = F_data->getGlobalVec();
    double* F_arr, * D_arr;
    ierr = VecGetArray(F_vec, &F_arr);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecGetArray(D_vec, &D_arr);  PETSC_SAMRAI_ERROR(ierr);

    const std::vector<int>& local_src_ids = d_local_src_ids[level_number];
    const std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    const std::vector<double>& rest_lengths = d_rest_lengths[level_number];
    for (int k = 0; k < static_cast<int>(local_src_ids.size()); ++k)
    {
        const int& f_idx = local_src_ids[k];
        double norm_D_sq = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            norm_D_sq += pow(D_arr[d+k*NDIM],2.0);
        }
        const double norm_D = sqrt(norm_D_sq);

        if (!SAMRAI::tbox::Utilities::deq(norm_D,0.0))
        {
            const double& stf = stiffnesses[k];
            const double& rst = rest_lengths[k];
            const double stf_scal = stf*(norm_D-rst)/norm_D;
            for (int d = 0; d < NDIM; ++d)
            {
                F_arr[d+f_idx*NDIM] += stf_scal*D_arr[d+k*NDIM];
            }
        }
    }

    ierr = VecRestoreArray(F_vec, &F_arr);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecRestoreArray(D_vec, &D_arr);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecDestroy(D_vec);               PETSC_SAMRAI_ERROR(ierr);

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::SpringForceGen>;

//////////////////////////////////////////////////////////////////////////////
