// Filename: IBSpringForceGen.C
// Last modified: <11.May.2008 16:11:57 griffith@box230.cims.nyu.edu>
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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force_jacobian;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceGen::IBSpringForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_D_mats(),
      d_J_mats(),
      d_lag_mastr_node_idxs(),
      d_lag_slave_node_idxs(),
      d_petsc_mastr_node_idxs(),
      d_petsc_slave_node_idxs(),
      d_force_fcn_idxs(),
      d_stiffnesses(),
      d_rest_lengths(),
      d_is_initialized(),
      d_force_fcn_map(),
      d_force_jacobian_fcn_map()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup the default force generation functions.
    registerSpringForceFunction(0, &IBAMR::default_linear_spring_force);
    registerSpringForceJacobianFunction(0, &IBAMR::default_linear_spring_force_jacobian);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBSpringForceGen::computeLagrangianForce()");
        t_compute_lagrangian_force_jacobian = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBSpringForceGen::computeLagrangianForceJacobian()");
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
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    for (std::vector<Mat>::iterator it = d_J_mats.begin(); it != d_J_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
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
IBSpringForceGen::registerSpringForceJacobianFunction(
    const int force_fcn_index,
    void (*force_jacobian_fcn)(double dF_dX[NDIM*NDIM], const double D[NDIM], const double& stf, const double& rst, const int& lag_mastr_idx, const int& lag_slave_idx))
{
    d_force_jacobian_fcn_map[force_fcn_index] = force_jacobian_fcn;
    return;
}// registerSpringForceJacobianFunction

void
IBSpringForceGen::initializeLevelData(
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
    const int new_size = std::max(level_number+1, int(d_is_initialized.size()));

    d_D_mats.resize(new_size);
    d_J_mats.resize(new_size);
    d_lag_mastr_node_idxs.resize(new_size);
    d_lag_slave_node_idxs.resize(new_size);
    d_petsc_mastr_node_idxs.resize(new_size);
    d_petsc_slave_node_idxs.resize(new_size);
    d_force_fcn_idxs.resize(new_size);
    d_stiffnesses.resize(new_size);
    d_rest_lengths.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_mat = d_D_mats[level_number];
    Mat& J_mat = d_J_mats[level_number];
    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    if (D_mat)
    {
        ierr = MatDestroy(D_mat);  IBTK_CHKERRQ(ierr);
    }
    if (J_mat)
    {
        ierr = MatDestroy(J_mat);  IBTK_CHKERRQ(ierr);
    }
    lag_mastr_node_idxs.clear();
    lag_slave_node_idxs.clear();
    petsc_mastr_node_idxs.clear();
    petsc_slave_node_idxs.clear();
    force_fcn_idxs.clear();
    stiffnesses.clear();
    rest_lengths.clear();

    // The patch data descriptor index for the LNodeIndexData2.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the "master" and "slave" node indices for all springs
    // associated with the present MPI process.
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
                    SAMRAI::tbox::Pointer<IBSpringForceSpec> force_spec = stash_data[l];
                    if (!force_spec.isNull())
                    {
                        const unsigned num_springs = force_spec->getNumberOfSprings();
#ifdef DEBUG_CHECK_ASSERTIONS
                        TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                        const std::vector<int>& slv = force_spec->getSlaveNodeIndices();
                        const std::vector<int>& fcn = force_spec->getForceFunctionIndices();
                        const std::vector<double>& stf = force_spec->getStiffnesses();
                        const std::vector<double>& rst = force_spec->getRestingLengths();
#ifdef DEBUG_CHECK_ASSERTIONS
                        TBOX_ASSERT(num_springs == slv.size());
                        TBOX_ASSERT(num_springs == fcn.size());
                        TBOX_ASSERT(num_springs == stf.size());
                        TBOX_ASSERT(num_springs == rst.size());
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

    // Map the Lagrangian master/slave node indices to the PETSc indices
    // corresponding to the present data distribution.
    petsc_mastr_node_idxs = lag_mastr_node_idxs;
    petsc_slave_node_idxs = lag_slave_node_idxs;
    lag_manager->mapLagrangianToPETSc(petsc_mastr_node_idxs, level_number);
    lag_manager->mapLagrangianToPETSc(petsc_slave_node_idxs, level_number);

    // Determine the global node offset and the number of local nodes.
    const int global_node_offset = lag_manager->getGlobalNodeOffset(level_number);
    const int num_local_nodes = lag_manager->getNumberOfLocalNodes(level_number);

    // Determine the non-zero structure for the matrix used to compute nodal
    // displacements.
    const int local_sz = petsc_mastr_node_idxs.size();
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

    // Create a new MPI block AIJ matrix used to compute nodal displacements and
    // set the values of the non-zero entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, &d_nz[0],
                            PETSC_DEFAULT, &o_nz[0],
                            &D_mat);  IBTK_CHKERRQ(ierr);

    std::vector<double> mastr_vals(NDIM*NDIM,0.0);
    std::vector<double> slave_vals(NDIM*NDIM,0.0);
    for (int d = 0; d < NDIM; ++d)
    {
        mastr_vals[d+d*NDIM] = -1.0;
        slave_vals[d+d*NDIM] = +1.0;
    }

    int i_offset;
    ierr = MatGetOwnershipRange(D_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        int i = i_offset + k;
        int j_mastr = petsc_mastr_node_idxs[k];
        int j_slave = petsc_slave_node_idxs[k];

        ierr = MatSetValuesBlocked(D_mat,1,&i,1,&j_mastr,&(mastr_vals[0]),
                                   INSERT_VALUES);  IBTK_CHKERRQ(ierr);

        ierr = MatSetValuesBlocked(D_mat,1,&i,1,&j_slave,&(slave_vals[0]),
                                   INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Determine the non-zero structure for the matrix used to store the
    // Jacobian of the force.
    //
    // NOTE: Each edge is *only* associated with a single node in the mesh.  We
    // must take this into account when determining the non-zero structure of
    // the matrix.
    MPI_Comm comm = SAMRAI::tbox::SAMRAI_MPI::getCommunicator();
    Vec d_nz_vec, o_nz_vec;
    ierr = VecCreateMPI(comm, num_local_nodes, PETSC_DETERMINE, &d_nz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(comm, num_local_nodes, PETSC_DETERMINE, &o_nz_vec);  IBTK_CHKERRQ(ierr);

    ierr = VecSet(d_nz_vec, 1.0);  IBTK_CHKERRQ(ierr);
    ierr = VecSet(o_nz_vec, 0.0);  IBTK_CHKERRQ(ierr);

    TBOX_ASSERT(local_sz == int(petsc_mastr_node_idxs.size()));
    for (int k = 0; k < local_sz; ++k)
    {
        const int& mastr_idx = petsc_mastr_node_idxs[k];
        const int& slave_idx = petsc_slave_node_idxs[k];

        static const int N = 2;
        const int idxs[N] = { mastr_idx , slave_idx };
        const double vals[N] = { 1.0 , 1.0 };
        if (slave_idx >= global_node_offset &&
            slave_idx <  global_node_offset+num_local_nodes)
        {
            ierr = VecSetValues(d_nz_vec, N, idxs, vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
        else
        {
            ierr = VecSetValues(o_nz_vec, N, idxs, vals, ADD_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    ierr = VecAssemblyBegin(d_nz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(o_nz_vec);  IBTK_CHKERRQ(ierr);

    ierr = VecAssemblyEnd(d_nz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(o_nz_vec);  IBTK_CHKERRQ(ierr);

    double* d_nz_vec_arr;
    double* o_nz_vec_arr;

    ierr = VecGetArray(d_nz_vec, &d_nz_vec_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(o_nz_vec, &o_nz_vec_arr);  IBTK_CHKERRQ(ierr);

    d_nz.clear();
    o_nz.clear();
    d_nz.resize(num_local_nodes,1);
    o_nz.resize(num_local_nodes,0);
    for (int k = 0; k < num_local_nodes; ++k)
    {
        d_nz[k] = int(d_nz_vec_arr[k]);
        o_nz[k] = int(o_nz_vec_arr[k]);
    }

    ierr = VecRestoreArray(d_nz_vec, &d_nz_vec_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArray(o_nz_vec, &o_nz_vec_arr);  IBTK_CHKERRQ(ierr);

    ierr = VecDestroy(d_nz_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(o_nz_vec);  IBTK_CHKERRQ(ierr);

    // Create a new MPI block AIJ matrix used to store the Jacobian of the force
    // and set dummy values for the matrix entries.
    ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                            NDIM, NDIM*num_local_nodes, NDIM*num_local_nodes,
                            PETSC_DETERMINE, PETSC_DETERMINE,
                            PETSC_DEFAULT, &d_nz[0],
                            PETSC_DEFAULT, &o_nz[0],
                            &J_mat);  IBTK_CHKERRQ(ierr);

    static const int N = 2;
    std::vector<double> dummy_vals(N*NDIM*NDIM,1.0);
    ierr = MatGetOwnershipRange(J_mat, &i_offset, PETSC_NULL);
    IBTK_CHKERRQ(ierr);
    i_offset /= NDIM;

    for (int k = 0; k < local_sz; ++k)
    {
        const int& i = petsc_mastr_node_idxs[k];
        const int& j = petsc_slave_node_idxs[k];

        // Set entries on row i.
        const int idxs_i[N] = { j , i };
        ierr = MatSetValuesBlocked(J_mat,1,&i,N,idxs_i,&(dummy_vals[0]),INSERT_VALUES);  IBTK_CHKERRQ(ierr);

        // Set entries on row j.
        const int idxs_j[N] = { i , j };
        ierr = MatSetValuesBlocked(J_mat,1,&j,N,idxs_j,&(dummy_vals[0]),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(D_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(J_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_number] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBSpringForceGen::computeLagrangianForce(
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> F_data,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    t_compute_lagrangian_force->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    // Create an appropriately sized temporary vector to store the node
    // displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_mats[level_number], X_data->getGlobalVec(), D_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the spring forces acting on the nodes of the Lagrangian mesh.
    double* D_arr;
    ierr = VecGetArray(D_vec, &D_arr);  IBTK_CHKERRQ(ierr);

    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    const int local_sz = petsc_mastr_node_idxs.size();
    std::vector<double> F_mastr_node_arr(NDIM*local_sz,0.0);
    std::vector<double> F_slave_node_arr(NDIM*local_sz,0.0);

    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the force applied by the spring to the "master" node.
        double* const F_mastr_node = &F_mastr_node_arr[k*NDIM];
        const double* const D = &D_arr[k*NDIM];
        const double& stf = stiffnesses[k];
        const double& rst = rest_lengths[k];
        const int& lag_mastr_idx = lag_mastr_node_idxs[k];
        const int& lag_slave_idx = lag_slave_node_idxs[k];
        const int& force_fcn_id = force_fcn_idxs[k];
        d_force_fcn_map[force_fcn_id](F_mastr_node,D,stf,rst,lag_mastr_idx,lag_slave_idx);

        // Compute the force applied by the spring to the corresponding "slave"
        // node.
        double* const F_slave_node = &F_slave_node_arr[k*NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            F_slave_node[d] = -F_mastr_node[d];
        }
    }

    ierr = VecRestoreArray(D_vec, &D_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_vec);               IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getGlobalVec();
    ierr = VecSetValuesBlocked(F_vec, petsc_mastr_node_idxs.size(), &petsc_mastr_node_idxs[0], &F_mastr_node_arr[0], ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecSetValuesBlocked(F_vec, petsc_slave_node_idxs.size(), &petsc_slave_node_idxs[0], &F_slave_node_arr[0], ADD_VALUES);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(F_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(F_vec);    IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force->stop();
    return;
}// computeLagrangianForce

void
IBSpringForceGen::computeLagrangianForceJacobian(
    Mat& J_mat,
    SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> X_data,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double data_time,
    IBTK::LDataManager* const lag_manager)
{
    t_compute_lagrangian_force_jacobian->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    int ierr;

    TBOX_WARNING("this should ADD to the Jacobian, not destroy it!\n");
    if (J_mat)
    {
        ierr = MatDestroy(J_mat);  IBTK_CHKERRQ(ierr);
    }
    ierr = MatDuplicate(d_J_mats[level_number], MAT_DO_NOT_COPY_VALUES, &J_mat);  IBTK_CHKERRQ(ierr);

    // Create an appropriately sized temporary vector to store the node
    // displacements.
    int i_start, i_stop;
    ierr = MatGetOwnershipRange(d_D_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_mats[level_number], X_data->getGlobalVec(), D_vec);
    IBTK_CHKERRQ(ierr);

    // Compute the force Jacobians and insert them into the Jacobian matrix.
    double* D_arr;
    ierr = VecGetArray(D_vec, &D_arr);  IBTK_CHKERRQ(ierr);

    std::vector<int>& lag_mastr_node_idxs = d_lag_mastr_node_idxs[level_number];
    std::vector<int>& lag_slave_node_idxs = d_lag_slave_node_idxs[level_number];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];
    std::vector<int>& petsc_slave_node_idxs = d_petsc_slave_node_idxs[level_number];
    std::vector<int>& force_fcn_idxs = d_force_fcn_idxs[level_number];
    std::vector<double>& stiffnesses = d_stiffnesses[level_number];
    std::vector<double>& rest_lengths = d_rest_lengths[level_number];

    const int local_sz = petsc_mastr_node_idxs.size();
    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the Jacobian of the force applied by the spring to the
        // "master" node with respect to the position of the "slave" node.
        std::vector<double> dF_dX(NDIM*NDIM,0.0);
        const double* const D = &D_arr[k*NDIM];
        const double& stf = stiffnesses[k];
        const double& rst = rest_lengths[k];
        const int& lag_mastr_idx = lag_mastr_node_idxs[k];
        const int& lag_slave_idx = lag_slave_node_idxs[k];
        const int& force_fcn_id = force_fcn_idxs[k];
        d_force_jacobian_fcn_map[force_fcn_id](&dF_dX[0],D,stf,rst,lag_mastr_idx,lag_slave_idx);

        // Get the PETSc indices corresponding to the "master" and "slave"
        // nodes.
        const int& petsc_mastr_idx = petsc_mastr_node_idxs[k];
        const int& petsc_slave_idx = petsc_slave_node_idxs[k];

        // Accumulate the off-diagonal parts of the matrix.
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_slave_idx,&dF_dX[0],ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_slave_idx,1,&petsc_mastr_idx,&dF_dX[0],ADD_VALUES);  IBTK_CHKERRQ(ierr);

        // Negate dF_dX to obtain the Jacobian of the force applied by the
        // spring to the "master" node with respect to the position of the
        // "master" node.
        for (int beta = 0; beta < NDIM; ++beta)
        {
            for (int alpha = 0; alpha < NDIM; ++alpha)
            {
                dF_dX[alpha+beta*NDIM] = -dF_dX[alpha+beta*NDIM];
            }
        }

        // Accumulate the diagonal parts of the matrix.
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_mastr_idx,1,&petsc_mastr_idx,&dF_dX[0],ADD_VALUES);  IBTK_CHKERRQ(ierr);
        ierr = MatSetValuesBlocked(J_mat,1,&petsc_slave_idx,1,&petsc_slave_idx,&dF_dX[0],ADD_VALUES);  IBTK_CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(D_vec, &D_arr);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_vec);               IBTK_CHKERRQ(ierr);

    // Assemble the matrix.
    ierr = MatAssemblyBegin(J_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(J_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_jacobian->stop();
    return;
}// computeLagrangianForceJacobian

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
