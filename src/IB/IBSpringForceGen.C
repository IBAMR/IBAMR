// Filename: IBSpringForceGen.C
// Last modified: <03.Apr.2007 21:14:42 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBSpringForceGen.h"

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
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/LNodeIndexData.h>

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceGen::IBSpringForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_D_mats(),
      d_lag_mastr_node_idxs(),
      d_lag_slave_node_idxs(),
      d_petsc_mastr_node_idxs(),
      d_petsc_slave_node_idxs(),
      d_force_fcn_idxs(),
      d_stiffnesses(),
      d_rest_lengths(),
      d_is_initialized(),
      d_force_fcn_map()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup the default force generation functions.
    registerSpringForceFunction(0, &IBAMR::default_linear_spring_force);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBSpringForceGen::computeLagrangianForce()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBSpringForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// IBSpringForceGen

IBSpringForceGen::~IBSpringForceGen()
{
    int ierr;
    for (std::vector<Mat>::iterator it = d_D_mats.begin(); it != d_D_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  PETSC_SAMRAI_ERROR(ierr);
        }
    }
    return;
}// ~IBSpringForceGen

void
IBSpringForceGen::registerSpringForceFunction(
    const int force_fcn_index,
    void (*force_fcn)(double F[NDIM], const double D[NDIM], const double& stf, const double& rst, const int& lag_mastr_idx, const int& lag_slave_idx))
{
    d_force_fcn_map[force_fcn_index] = force_fcn;
    return;
}// registerSpringForceFunction

void
IBSpringForceGen::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    LDataManager* const lag_manager)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
#endif
    (void) init_data_time;
    (void) initial_time;

    int ierr;

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Resize the vectors corresponding to data individually
    // maintained for separate levels of the patch hierarchy.
    const int level_num = level->getLevelNumber();
    const int new_size = SAMRAI::tbox::Utilities::imax(
        level_num+1, d_is_initialized.size());

    d_D_mats.resize(new_size);
    d_lag_mastr_node_idxs.resize(new_size);
    d_lag_slave_node_idxs.resize(new_size);
    d_petsc_mastr_node_idxs.resize(new_size);
    d_petsc_slave_node_idxs.resize(new_size);
    d_force_fcn_idxs.resize(new_size);
    d_stiffnesses.resize(new_size);
    d_rest_lengths.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_mat = d_D_mats[level_num];
    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_num];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_num];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_num];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_num];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_num];
    std::vector<double>& stiffnesses = d_stiffnesses[level_num];
    std::vector<double>& rest_lengths = d_rest_lengths[level_num];

    if (D_mat)
    {
        ierr = MatDestroy(D_mat);  PETSC_SAMRAI_ERROR(ierr);
    }
    lag_mastr_node_idxs.clear();
    lag_slave_node_idxs.clear();
    petsc_mastr_node_idxs.clear();
    petsc_slave_node_idxs.clear();
    force_fcn_idxs.clear();
    stiffnesses.clear();
    rest_lengths.clear();

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->
        getLNodeIndexPatchDescriptorIndex();

    // Determine the "master" and "slave" node indices for all springs
    // associated with the present MPI process.
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
                    const int& mastr_idx = node_idx->getLagrangianIndex();
                    const std::vector<SAMRAI::tbox::Pointer<Stashable> >& stash_data =
                        node_idx->getStashData();
                    for (unsigned l = 0; l < stash_data.size(); ++l)
                    {
                        SAMRAI::tbox::Pointer<IBSpringForceSpec> force_spec = stash_data[l];
                        if (!force_spec.isNull())
                        {
                            const unsigned num_springs = force_spec->getNumberOfSprings();
#ifdef DEBUG_CHECK_ASSERTIONS
                            assert(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                            const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
                            const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
                            const std::vector<double>& stf = force_spec->getStiffnesses();
                            const std::vector<double>& rst = force_spec->getRestingLengths();
#ifdef DEBUG_CHECK_ASSERTIONS
                            assert(num_springs == slv.size());
                            assert(num_springs == fcn.size());
                            assert(num_springs == stf.size());
                            assert(num_springs == rst.size());
#endif
                            if (num_springs > 0)
                            {
                                lag_mastr_node_idxs.insert(lag_mastr_node_idxs.end(), num_springs, mastr_idx);
                                lag_slave_node_idxs.insert(lag_slave_node_idxs.end(), slv.begin(), slv.end());
                                force_fcn_idxs.insert(force_fcn_idxs.end(), fcn.begin(), fcn.end());
                                stiffnesses   .insert(stiffnesses   .end(), stf.begin(), stf.end());
                                rest_lengths  .insert(rest_lengths  .end(), rst.begin(), rst.end());
                            }
                        }
                    }
                }
            }
        }
    }

    // Map the Lagrangian master/slave node indices to the PETSc
    // indices corresponding to the present data distribution.
    petsc_mastr_node_idxs = lag_mastr_node_idxs;
    petsc_slave_node_idxs = lag_slave_node_idxs;
    lag_manager->mapLagrangianToPETSc(petsc_mastr_node_idxs, level_num);
    lag_manager->mapLagrangianToPETSc(petsc_slave_node_idxs, level_num);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_num);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_num);

    // Determine the non-zero structure for the matrix.
    const int local_sz = static_cast<int>(petsc_mastr_node_idxs.size());
    std::vector<int> d_nz(local_sz,1), o_nz(local_sz,0);

    for (int k = 0; k < local_sz; ++k)
    {
        const int& slave_idx = petsc_slave_node_idxs[k];
        if (slave_idx >= global_node_offset &&
            slave_idx <  global_node_offset+num_local_nodes)
        {
            ++d_nz[k]; // a "local"    slave index
        }
        else
        {
            ++o_nz[k]; // a "nonlocal" slave index
        }
    }

    // Create a new MPI block AIJ matrix and set the values of the
    // non-zero entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, &d_nz[0],
                            PETSC_DEFAULT, &o_nz[0],
                            &D_mat);  PETSC_SAMRAI_ERROR(ierr);

    std::vector<double> mastr_vals(NDIM*NDIM,0.0);
    std::vector<double> slave_vals(NDIM*NDIM,0.0);
    for (int d = 0; d < NDIM; ++d)
    {
        mastr_vals[d+d*NDIM] = -1.0;
        slave_vals[d+d*NDIM] = +1.0;
    }

    int i_offset;
    ierr = MatGetOwnershipRange(D_mat, &i_offset, PETSC_NULL);
    PETSC_SAMRAI_ERROR(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_slave_node_idxs[k];

        ierr = MatSetValuesBlocked(D_mat,1,&i,1,&j_mastr,&(mastr_vals[0]),
                                   INSERT_VALUES);  PETSC_SAMRAI_ERROR(ierr);

        ierr = MatSetValuesBlocked(D_mat,1,&i,1,&j_slave,&(slave_vals[0]),
                                   INSERT_VALUES);  PETSC_SAMRAI_ERROR(ierr);
    }

    // Assemble the matrix.
    ierr = MatAssemblyBegin(D_mat, MAT_FINAL_ASSEMBLY);  PETSC_SAMRAI_ERROR(ierr);
    ierr = MatAssemblyEnd(D_mat, MAT_FINAL_ASSEMBLY);  PETSC_SAMRAI_ERROR(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_num] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBSpringForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    LDataManager* const lag_manager)
{
    t_compute_lagrangian_force->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(level_number < static_cast<int>(d_is_initialized.size()));
    assert(d_is_initialized[level_number]);
#endif

    int ierr;

    // Create an appropriately sized temporary vector to store the
    // node displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_mats[level_number], &i_start, &i_stop);
    PETSC_SAMRAI_ERROR(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          PETSC_SAMRAI_ERROR(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_mats[level_number], X_data->getGlobalVec(), D_vec);
    PETSC_SAMRAI_ERROR(ierr);

    // Compute the spring forces acting on the nodes of the Lagrangian
    // mesh.
    double* D_arr;
    ierr = VecGetArray(D_vec, &D_arr);  PETSC_SAMRAI_ERROR(ierr);

    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    const int local_sz = static_cast<int>(petsc_mastr_node_idxs.size());
    std::vector<double> F_mastr_node_arr(NDIM*local_sz);
    std::vector<double> F_slave_node_arr(NDIM*local_sz);

    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the force applied by the spring to the "master"
        // node.
        double* const F_mastr_node = &F_mastr_node_arr[k*NDIM];
        const double* const D = &D_arr[k*NDIM];
        const int& lag_mastr_idx = lag_mastr_node_idxs[k];
        const int& lag_slave_idx = lag_slave_node_idxs[k];
        const int& force_fcn_id = force_fcn_idxs[k];
        const double& stf = stiffnesses[k];
        const double& rst = rest_lengths[k];
        d_force_fcn_map[force_fcn_id](F_mastr_node,D,stf,rst,lag_mastr_idx,lag_slave_idx);

        // Compute the force applied by the spring to the
        // corresponding "slave" node.
        double* const F_slave_node = &F_slave_node_arr[k*NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            F_slave_node[d] = -F_mastr_node[d];
        }
    }

    ierr = VecRestoreArray(D_vec, &D_arr);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecDestroy(D_vec);               PETSC_SAMRAI_ERROR(ierr);

    Vec F_vec = F_data->getGlobalVec();
    ierr = VecSetValuesBlocked(F_vec, petsc_mastr_node_idxs.size(), &petsc_mastr_node_idxs[0], &F_mastr_node_arr[0], ADD_VALUES);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecSetValuesBlocked(F_vec, petsc_slave_node_idxs.size(), &petsc_slave_node_idxs[0], &F_slave_node_arr[0], ADD_VALUES);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecAssemblyBegin(F_vec);  PETSC_SAMRAI_ERROR(ierr);
    ierr = VecAssemblyEnd(F_vec);    PETSC_SAMRAI_ERROR(ierr);

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBSpringForceGen::getFromInput(
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
template class SAMRAI::tbox::Pointer<IBAMR::IBSpringForceGen>;

//////////////////////////////////////////////////////////////////////////////
