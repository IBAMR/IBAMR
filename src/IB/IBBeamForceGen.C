// Filename: IBBeamForceGen.C
// Last modified: <12.Mar.2008 22:33:47 griffith@box221.cims.nyu.edu>
// Created on 22 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBBeamForceGen.h"

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
#include <ibamr/IBBeamForceSpec.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LNodeIndexData2.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBeamForceGen::IBBeamForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_D_next_mats(),
      d_D_prev_mats(),
      d_petsc_mastr_node_idxs(),
      d_petsc_next_node_idxs(),
      d_petsc_prev_node_idxs(),
      d_bend_rigidities(),
      d_is_initialized()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBBeamForceGen::computeLagrangianForce()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBBeamForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// IBBeamForceGen

IBBeamForceGen::~IBBeamForceGen()
{
    int ierr;
    for (std::vector<Mat>::iterator it = d_D_next_mats.begin(); it != d_D_next_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    for (std::vector<Mat>::iterator it = d_D_prev_mats.begin(); it != d_D_prev_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// ~IBBeamForceGen

void
IBBeamForceGen::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
#endif
    (void) init_data_time;
    (void) initial_time;

    int ierr;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually maintained for
    // separate levels of the patch hierarchy.
    const int level_num = level->getLevelNumber();
    const int new_size = std::max(level_num+1, int(d_is_initialized.size()));

    d_D_next_mats.resize(new_size);
    d_D_prev_mats.resize(new_size);
    d_petsc_mastr_node_idxs.resize(new_size);
    d_petsc_next_node_idxs.resize(new_size);
    d_petsc_prev_node_idxs.resize(new_size);
    d_bend_rigidities.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_next_mat = d_D_next_mats[level_num];
    Mat& D_prev_mat = d_D_prev_mats[level_num];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_num];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_num];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_num];
    std::vector<double>& bend_rigidities = d_bend_rigidities[level_num];

    if (D_next_mat)
    {
        ierr = MatDestroy(D_next_mat);  IBTK_CHKERRQ(ierr);
    }
    if (D_prev_mat)
    {
        ierr = MatDestroy(D_prev_mat);  IBTK_CHKERRQ(ierr);
    }
    petsc_mastr_node_idxs.clear();
    petsc_next_node_idxs.clear();
    petsc_prev_node_idxs.clear();
    bend_rigidities.clear();

    // The patch data descriptor index for the LNodeIndexData2.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the "next" and "prev" node indices for all beams associated
    // with the present MPI process.
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData2> idx_data =
            patch->getPatchData(lag_node_index_idx);

        for (IBTK::LNodeIndexData2::Iterator it(patch_box); it; it++)
        {
            const SAMRAI::pdat::CellIndex<NDIM>& i = *it;
            const IBTK::LNodeIndexSet& node_set = (*idx_data)(i);
            for (IBTK::LNodeIndexSet::const_iterator n = node_set.begin();
                 n != node_set.end(); ++n)
            {
                const IBTK::LNodeIndexSet::value_type& node_idx = *n;
                const int& mastr_idx = node_idx->getLagrangianIndex();
                const std::vector<SAMRAI::tbox::Pointer<IBTK::Stashable> >& stash_data =
                    node_idx->getStashData();
                for (unsigned l = 0; l < stash_data.size(); ++l)
                {
                    SAMRAI::tbox::Pointer<IBBeamForceSpec> force_spec = stash_data[l];
                    if (!force_spec.isNull())
                    {
                        const unsigned num_beams = force_spec->getNumberOfBeams();
#ifdef DEBUG_CHECK_ASSERTIONS
                        TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                        const std::vector<std::pair<int,int> >& nghbrs = force_spec->getNeighborNodeIndices();
                        const std::vector<double>& bend = force_spec->getBendingRigidities();
#ifdef DEBUG_CHECK_ASSERTIONS
                        TBOX_ASSERT(num_beams == nghbrs.size());
                        TBOX_ASSERT(num_beams == bend.size());
#endif
                        for (unsigned k = 0; k < num_beams; ++k)
                        {
                            petsc_mastr_node_idxs.push_back(mastr_idx);
                            petsc_next_node_idxs.push_back(nghbrs[k].second);
                            petsc_prev_node_idxs.push_back(nghbrs[k].first );
                            bend_rigidities.push_back(bend[k]);
                        }
                    }
                }
            }
        }
    }

    // Map the Lagrangian node indices to the PETSc indices corresponding to the
    // present data distribution.
    lag_manager->mapLagrangianToPETSc(petsc_mastr_node_idxs, level_num);
    lag_manager->mapLagrangianToPETSc(petsc_next_node_idxs, level_num);
    lag_manager->mapLagrangianToPETSc(petsc_prev_node_idxs, level_num);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_num);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_num);

    // Determine the non-zero structure for the matrices.
    const int local_sz = static_cast<int>(petsc_mastr_node_idxs.size());

    std::vector<int> next_d_nz(local_sz,1), next_o_nz(local_sz,0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& next_idx = petsc_next_node_idxs[k];
        if (next_idx >= global_node_offset &&
            next_idx <  global_node_offset+num_local_nodes)
        {
            ++next_d_nz[k]; // a "local"    next index
        }
        else
        {
            ++next_o_nz[k]; // a "nonlocal" next index
        }
    }

    std::vector<int> prev_d_nz(local_sz,1), prev_o_nz(local_sz,0);
    for (int k = 0; k < local_sz; ++k)
    {
        const int& prev_idx = petsc_prev_node_idxs[k];
        if (prev_idx >= global_node_offset &&
            prev_idx <  global_node_offset+num_local_nodes)
        {
            ++prev_d_nz[k]; // a "local"    prev index
        }
        else
        {
            ++prev_o_nz[k]; // a "nonlocal" prev index
        }
    }

    // Create new MPI block AIJ matrices and set the values of the non-zero
    // entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, &next_d_nz[0],
                            PETSC_DEFAULT, &next_o_nz[0],
                            &D_next_mat);  IBTK_CHKERRQ(ierr);

    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, &prev_d_nz[0],
                            PETSC_DEFAULT, &prev_o_nz[0],
                            &D_prev_mat);  IBTK_CHKERRQ(ierr);

    std::vector<double> mastr_vals(NDIM*NDIM,0.0);
    std::vector<double> slave_vals(NDIM*NDIM,0.0);
    for (int d = 0; d < NDIM; ++d)
    {
        mastr_vals[d+d*NDIM] = -1.0;
        slave_vals[d+d*NDIM] = +1.0;
    }

    int i_offset;

    ierr = MatGetOwnershipRange(D_next_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_next_node_idxs[k];

        ierr = MatSetValuesBlocked(D_next_mat,1,&i,1,&j_mastr,&(mastr_vals[0]),
                                   INSERT_VALUES);  IBTK_CHKERRQ(ierr);

        ierr = MatSetValuesBlocked(D_next_mat,1,&i,1,&j_slave,&(slave_vals[0]),
                                   INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    ierr = MatGetOwnershipRange(D_prev_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_prev_node_idxs[k];

        ierr = MatSetValuesBlocked(D_prev_mat,1,&i,1,&j_mastr,&(mastr_vals[0]),
                                   INSERT_VALUES);  IBTK_CHKERRQ(ierr);

        ierr = MatSetValuesBlocked(D_prev_mat,1,&i,1,&j_slave,&(slave_vals[0]),
                                   INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(D_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(D_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_num] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBBeamForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    t_compute_lagrangian_force->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < static_cast<int>(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Create an appropriately sized temporary vector to store the node
    // displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_next_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    Vec D_prev_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_prev_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_prev_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_next_mats[level_number], X_data->getGlobalVec(), D_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_D_prev_mats[level_number], X_data->getGlobalVec(), D_prev_vec);  IBTK_CHKERRQ(ierr);

    // Compute the beam forces acting on the nodes of the Lagrangian mesh.
    double* D_next_arr;
    ierr = VecGetArray(D_next_vec, &D_next_arr);  IBTK_CHKERRQ(ierr);

    double* D_prev_arr;
    ierr = VecGetArray(D_prev_vec, &D_prev_arr);  IBTK_CHKERRQ(ierr);

    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_number];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_number];
    std::vector<double>& bend_rigidities = d_bend_rigidities[level_number];

    const int local_sz = static_cast<int>(petsc_mastr_node_idxs.size());
    std::vector<double> F_mastr_node_arr(NDIM*local_sz,0.0);
    std::vector<double> F_nghbr_node_arr(NDIM*local_sz,0.0);

    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the forces applied by the beam to the "master" and "slave"
        // nodes.
        double* const F_mastr_node = &F_mastr_node_arr[k*NDIM];
        double* const F_nghbr_node = &F_nghbr_node_arr[k*NDIM];
        const double* const D_next = &D_next_arr[k*NDIM];
        const double* const D_prev = &D_prev_arr[k*NDIM];
        const double& bnd = bend_rigidities[k];

        for (int d = 0; d < NDIM; ++d)
        {
            F_mastr_node[d] = +2.0*bnd*(D_next[d]+D_prev[d]);
            F_nghbr_node[d] = -1.0*bnd*(D_next[d]+D_prev[d]);
        }
    }

    ierr = VecRestoreArray(D_next_vec, &D_next_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_next_vec);                    IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_prev_vec, &D_prev_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_prev_vec);                    IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getGlobalVec();
    ierr = VecSetValuesBlocked(F_vec, petsc_mastr_node_idxs.size(), &petsc_mastr_node_idxs[0], &F_mastr_node_arr[0], ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValuesBlocked(F_vec, petsc_next_node_idxs.size(), &petsc_next_node_idxs[0], &F_nghbr_node_arr[0], ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValuesBlocked(F_vec, petsc_prev_node_idxs.size(), &petsc_prev_node_idxs[0], &F_nghbr_node_arr[0], ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_vec);    IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBBeamForceGen::getFromInput(
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
template class SAMRAI::tbox::Pointer<IBAMR::IBBeamForceGen>;

//////////////////////////////////////////////////////////////////////////////
