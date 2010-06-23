// Filename: IBKirchhoffRodForceGen.C
// Last modified: <22.Jun.2010 22:05:19 griffith@boyce-griffiths-mac-pro.local>
// Created on 22 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-mac-pro.local)

#include "IBKirchhoffRodForceGen.h"

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
#include <ibtk/LNodeIndexData.h>
#include <ibtk/PETScVecOps.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

// BLITZ++ INCLUDES
#include <blitz/array.h>

// C++ STDLIB INCLUDES
#include <numeric>

// FORTRAN ROUTINES
#define SQRTM_FC FC_FUNC(sqrtm,SQRTM)
extern "C"
{
    void
    SQRTM_FC(
        double* X, const double* A);
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
void
interpolate_directors(
    const blitz::Array<blitz::Array<double,1>*,1>& D_half,
    const blitz::Array<blitz::Array<double,1>*,1>& D,
    const blitz::Array<blitz::Array<double,1>*,1>& D_next)
{
    blitz::Array<double,2> A(3,3,blitz::ColumnMajorArray<2>());
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            A(i,j) = (*D_next(0))(i)*(*D(0))(j) + (*D_next(1))(i)*(*D(1))(j) + (*D_next(2))(i)*(*D(2))(j);
        }
    }

    blitz::Array<double,2> sqrt_A(3,3,blitz::ColumnMajorArray<2>());
    SQRTM_FC(sqrt_A.data(), A.data());

    for (int alpha = 0; alpha < 3; ++alpha)
    {
        blitz::firstIndex i;
        blitz::secondIndex j;
        *D_half(alpha) = blitz::sum(sqrt_A(i,j)*(*D(alpha))(j),j);
    }
    return;
}// interpolate_directors

inline double
dot(
    const blitz::Array<double,1>& x,
    const blitz::Array<double,1>& y)
{
    blitz::firstIndex i;
    return blitz::sum(x(i) * y(i));
}// dot

inline blitz::Array<double,1>
cross(
    const blitz::Array<double,1>& x,
    const blitz::Array<double,1>& y)
{
    blitz::Array<double,1> x_cross_y(3);
    x_cross_y = x(1)*y(2) - y(1)*x(2) , y(0)*x(2) - x(0)*y(2) , x(0)*y(1) - y(0)*x(1);
    return x_cross_y;
}// cross

void
compute_force_and_torque(
    blitz::Array<double,1>& F_half,
    blitz::Array<double,1>& N_half,
    blitz::Array<double,1>& X,
    blitz::Array<double,1>& X_next,
    blitz::Array<double,1>& D1,
    blitz::Array<double,1>& D1_next,
    blitz::Array<double,1>& D2,
    blitz::Array<double,1>& D2_next,
    blitz::Array<double,1>& D3,
    blitz::Array<double,1>& D3_next)
{
    blitz::Array<blitz::Array<double,1>*,1> D(3);
    D = &D1 , &D2 , &D3;

    blitz::Array<blitz::Array<double,1>*,1> D_next(3);
    D_next = &D1_next , &D2_next , &D3_next;

    blitz::Array<blitz::Array<double,1>*,1> D_half(3);
    blitz::Array<double,1> D1_half(3);
    D_half(0) = &D1_half;
    blitz::Array<double,1> D2_half(3);
    D_half(1) = &D2_half;
    blitz::Array<double,1> D3_half(3);
    D_half(2) = &D3_half;

    interpolate_directors(D_half, D, D_next);

    const double a1 = 0.3;
    const double a2 = 0.3;
    const double a3 = 0.2;
    const double b1 = 54.0;
    const double b2 = 54.0;
    const double b3 = 54.0;
    const double ds = 0.0785398163397448;

    const blitz::Array<double,1> dX_ds((X_next-X)/ds);
    const double F1 = b1* dot(D1_half, dX_ds);
    const double F2 = b2* dot(D2_half, dX_ds);
    const double F3 = b3*(dot(D3_half, dX_ds) - 1.0);
    F_half = F1*D1_half + F2*D2_half + F3*D3_half;

    const blitz::Array<double,1> dD1_ds((D1_next-D1)/ds);
    const blitz::Array<double,1> dD2_ds((D2_next-D2)/ds);
    const blitz::Array<double,1> dD3_ds((D3_next-D3)/ds);
    const double N1 = a1*dot(dD2_ds,D3_half);
    const double N2 = a2*dot(dD3_ds,D1_half);
    const double N3 = a3*dot(dD1_ds,D2_half);
    N_half = N1*D1_half + N2*D2_half + N3*D3_half;
    return;
}// compute_force_and_torque

// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_compute_lagrangian_force_and_torque;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBKirchhoffRodForceGen::IBKirchhoffRodForceGen(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_D_next_mats(),
      d_D_prev_mats(),
      d_X_next_mats(),
      d_X_prev_mats(),
      d_petsc_mastr_node_idxs(),
      d_petsc_next_node_idxs(),
      d_petsc_prev_node_idxs(),
      d_is_initialized()
{
    // Initialize object with data read from the input database.
    getFromInput(input_db);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_compute_lagrangian_force_and_torque = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBKirchhoffRodForceGen::computeLagrangianForceAndTorque()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::IBKirchhoffRodForceGen::initializeLevelData()");
        timers_need_init = false;
    }
    return;
}// IBKirchhoffRodForceGen

IBKirchhoffRodForceGen::~IBKirchhoffRodForceGen()
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
    for (std::vector<Mat>::iterator it = d_X_next_mats.begin(); it != d_X_next_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    for (std::vector<Mat>::iterator it = d_X_prev_mats.begin(); it != d_X_prev_mats.end(); ++it)
    {
        if (*it)
        {
            ierr = MatDestroy(*it);  IBTK_CHKERRQ(ierr);
        }
    }
    return;
}// ~IBKirchhoffRodForceGen

void
IBKirchhoffRodForceGen::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool initial_time,
    IBTK::LDataManager* const lag_manager)
{
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

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
    d_X_next_mats.resize(new_size);
    d_X_prev_mats.resize(new_size);
    d_petsc_mastr_node_idxs.resize(new_size);
    d_petsc_next_node_idxs.resize(new_size);
    d_petsc_prev_node_idxs.resize(new_size);
    d_is_initialized.resize(new_size, false);

    Mat& D_next_mat = d_D_next_mats[level_num];
    Mat& D_prev_mat = d_D_prev_mats[level_num];
    Mat& X_next_mat = d_X_next_mats[level_num];
    Mat& X_prev_mat = d_X_prev_mats[level_num];
    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_num];
    std::vector<int>& petsc_next_node_idxs = d_petsc_next_node_idxs[level_num];
    std::vector<int>& petsc_prev_node_idxs = d_petsc_prev_node_idxs[level_num];

    if (D_next_mat)
    {
        ierr = MatDestroy(D_next_mat);  IBTK_CHKERRQ(ierr);
    }
    if (D_prev_mat)
    {
        ierr = MatDestroy(D_prev_mat);  IBTK_CHKERRQ(ierr);
    }
    if (X_next_mat)
    {
        ierr = MatDestroy(X_next_mat);  IBTK_CHKERRQ(ierr);
    }
    if (X_prev_mat)
    {
        ierr = MatDestroy(X_prev_mat);  IBTK_CHKERRQ(ierr);
    }
    petsc_mastr_node_idxs.clear();
    petsc_next_node_idxs.clear();
    petsc_prev_node_idxs.clear();

    // The patch data descriptor index for the LNodeIndexData.
    const int lag_node_index_idx = lag_manager->getLNodeIndexPatchDescriptorIndex();

    // Determine the "next" and "prev" node indices for all beams associated
    // with the present MPI process.
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        const SAMRAI::hier::Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::tbox::Pointer<IBTK::LNodeIndexData> idx_data = patch->getPatchData(lag_node_index_idx);
        for (IBTK::LNodeIndexData::LNodeIndexIterator it = idx_data->lnode_index_begin(patch_box);
             it != idx_data->lnode_index_end(); ++it)
        {
            const IBTK::LNodeIndex& node_idx = *it;
            const SAMRAI::tbox::Pointer<IBBeamForceSpec> force_spec = node_idx.getStashData<IBBeamForceSpec>();
            if (!force_spec.isNull())
            {
                const int& mastr_idx = node_idx.getLagrangianIndex();
                const unsigned num_beams = force_spec->getNumberOfBeams();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(mastr_idx == force_spec->getMasterNodeIndex());
#endif
                const std::vector<std::pair<int,int> >& nghbrs = force_spec->getNeighborNodeIndices();
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(num_beams == nghbrs.size());
#endif
                for (unsigned k = 0; k < num_beams; ++k)
                {
                    petsc_mastr_node_idxs.push_back(mastr_idx);
                    petsc_next_node_idxs.push_back(nghbrs[k].first );
                    petsc_prev_node_idxs.push_back(nghbrs[k].second);
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
    const int local_sz = petsc_mastr_node_idxs.size();

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
    {
        ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                                3*3, 3*3*local_sz, 3*3*num_local_nodes,
                                PETSC_DETERMINE, PETSC_DETERMINE,
                                PETSC_DEFAULT, local_sz > 0 ? &next_d_nz[0] : PETSC_NULL,
                                PETSC_DEFAULT, local_sz > 0 ? &next_o_nz[0] : PETSC_NULL,
                                &D_next_mat);  IBTK_CHKERRQ(ierr);

        ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                                3*3, 3*3*local_sz, 3*3*num_local_nodes,
                                PETSC_DETERMINE, PETSC_DETERMINE,
                                PETSC_DEFAULT, local_sz > 0 ? &prev_d_nz[0] : PETSC_NULL,
                                PETSC_DEFAULT, local_sz > 0 ? &prev_o_nz[0] : PETSC_NULL,
                                &D_prev_mat);  IBTK_CHKERRQ(ierr);

        blitz::Array<double,2> mastr_vals(3*3,3*3);  mastr_vals = 0.0;
        blitz::Array<double,2> slave_vals(3*3,3*3);  slave_vals = 0.0;
        for (int d = 0; d < 3*3; ++d)
        {
            mastr_vals(d,d) =  0.0;
            slave_vals(d,d) = +1.0;
        }

        int i_offset;

        ierr = MatGetOwnershipRange(D_next_mat, &i_offset, PETSC_NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= 3*3;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_mastr = petsc_mastr_node_idxs[k];
            int j_slave = petsc_next_node_idxs[k];
            ierr = MatSetValuesBlocked(D_next_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(D_next_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }

        ierr = MatGetOwnershipRange(D_prev_mat, &i_offset, PETSC_NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= 3*3;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_mastr = petsc_mastr_node_idxs[k];
            int j_slave = petsc_prev_node_idxs[k];
            ierr = MatSetValuesBlocked(D_prev_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(D_prev_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }

    }

    {
        ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                                NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                                PETSC_DETERMINE, PETSC_DETERMINE,
                                PETSC_DEFAULT, local_sz > 0 ? &next_d_nz[0] : PETSC_NULL,
                                PETSC_DEFAULT, local_sz > 0 ? &next_o_nz[0] : PETSC_NULL,
                                &X_next_mat);  IBTK_CHKERRQ(ierr);

        ierr = MatCreateMPIBAIJ(PETSC_COMM_WORLD,
                                NDIM, NDIM*local_sz, NDIM*num_local_nodes,
                                PETSC_DETERMINE, PETSC_DETERMINE,
                                PETSC_DEFAULT, local_sz > 0 ? &prev_d_nz[0] : PETSC_NULL,
                                PETSC_DEFAULT, local_sz > 0 ? &prev_o_nz[0] : PETSC_NULL,
                                &X_prev_mat);  IBTK_CHKERRQ(ierr);

        blitz::Array<double,2> mastr_vals(NDIM,NDIM);  mastr_vals = 0.0;
        blitz::Array<double,2> slave_vals(NDIM,NDIM);  slave_vals = 0.0;
        for (int d = 0; d < NDIM; ++d)
        {
            mastr_vals(d,d) =  0.0;
            slave_vals(d,d) = +1.0;
        }

        int i_offset;

        ierr = MatGetOwnershipRange(X_next_mat, &i_offset, PETSC_NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= NDIM;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_mastr = petsc_mastr_node_idxs[k];
            int j_slave = petsc_next_node_idxs[k];
            ierr = MatSetValuesBlocked(X_next_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(X_next_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }

        ierr = MatGetOwnershipRange(X_prev_mat, &i_offset, PETSC_NULL);
        IBTK_CHKERRQ(ierr);
        i_offset /= NDIM;

        for (int k = 0; k < local_sz; ++k)
        {
            int i = i_offset + k;
            int j_mastr = petsc_mastr_node_idxs[k];
            int j_slave = petsc_prev_node_idxs[k];
            ierr = MatSetValuesBlocked(X_prev_mat,1,&i,1,&j_mastr,mastr_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
            ierr = MatSetValuesBlocked(X_prev_mat,1,&i,1,&j_slave,slave_vals.data(),INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(D_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(D_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(X_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyBegin(X_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(D_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(X_next_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(X_prev_mat, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);

    // Indicate that the level data has been initialized.
    d_is_initialized[level_num] = true;

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
IBKirchhoffRodForceGen::computeLagrangianForceAndTorque(
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
    if (!lag_manager->levelContainsLagrangianData(level_number)) return;

    t_compute_lagrangian_force_and_torque->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(level_number < int(d_is_initialized.size()));
    TBOX_ASSERT(d_is_initialized[level_number]);
#endif

    const int global_offset = lag_manager->getGlobalNodeOffset(level_number);

    int ierr;

    // Create appropriately sized temporary vectors.
    int i_start, i_stop;

    ierr = MatGetOwnershipRange(d_D_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec D_vec = D_data->getGlobalVec();

    Vec D_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_next_vec, 3*3);                                           IBTK_CHKERRQ(ierr);

    Vec D_prev_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &D_prev_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(D_prev_vec, 3*3);                                           IBTK_CHKERRQ(ierr);

    ierr = MatGetOwnershipRange(d_X_next_mats[level_number], &i_start, &i_stop);
    IBTK_CHKERRQ(ierr);

    Vec X_vec = X_data->getGlobalVec();

    Vec X_next_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &X_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(X_next_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    Vec X_prev_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD, i_stop-i_start, PETSC_DECIDE, &X_prev_vec);  IBTK_CHKERRQ(ierr);
    ierr = VecSetBlockSize(X_prev_vec, NDIM);                                          IBTK_CHKERRQ(ierr);

    // Compute the node displacements.
    ierr = MatMult(d_D_next_mats[level_number], D_data->getGlobalVec(), D_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_D_prev_mats[level_number], D_data->getGlobalVec(), D_prev_vec);  IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_X_next_mats[level_number], X_data->getGlobalVec(), X_next_vec);  IBTK_CHKERRQ(ierr);
    ierr = MatMult(d_X_prev_mats[level_number], X_data->getGlobalVec(), X_prev_vec);  IBTK_CHKERRQ(ierr);

    // Compute the beam forces acting on the nodes of the Lagrangian mesh.
    double* D_vals;
    ierr = VecGetArray(D_vec, &D_vals);  IBTK_CHKERRQ(ierr);

    double* D_next_vals;
    ierr = VecGetArray(D_next_vec, &D_next_vals);  IBTK_CHKERRQ(ierr);

    double* D_prev_vals;
    ierr = VecGetArray(D_prev_vec, &D_prev_vals);  IBTK_CHKERRQ(ierr);

    double* X_vals;
    ierr = VecGetArray(X_vec, &X_vals);  IBTK_CHKERRQ(ierr);

    double* X_next_vals;
    ierr = VecGetArray(X_next_vec, &X_next_vals);  IBTK_CHKERRQ(ierr);

    double* X_prev_vals;
    ierr = VecGetArray(X_prev_vec, &X_prev_vals);  IBTK_CHKERRQ(ierr);

    std::vector<int>& petsc_mastr_node_idxs = d_petsc_mastr_node_idxs[level_number];

    const int local_sz = petsc_mastr_node_idxs.size();
    std::vector<double> F_node_vals(NDIM*local_sz,0.0);
    std::vector<double> N_node_vals(NDIM*local_sz,0.0);

    for (int k = 0; k < local_sz; ++k)
    {
        // Compute the forces applied by the beam to the "master" and "slave"
        // nodes.
        blitz::Array<double,1> F(&F_node_vals[k*NDIM],blitz::shape(NDIM),blitz::neverDeleteData);
        blitz::Array<double,1> N(&N_node_vals[k*NDIM],blitz::shape(NDIM),blitz::neverDeleteData);

        blitz::Array<double,1> D1(&D_vals[(petsc_mastr_node_idxs[k]-global_offset)*3*3],blitz::shape(3),blitz::neverDeleteData);
        blitz::Array<double,1> D1_next(&D_next_vals[k*3*3],blitz::shape(3),blitz::neverDeleteData);
        blitz::Array<double,1> D1_prev(&D_prev_vals[k*3*3],blitz::shape(3),blitz::neverDeleteData);

        const int D2_offset = 3;
        blitz::Array<double,1> D2(&D_vals[(petsc_mastr_node_idxs[k]-global_offset)*3*3+D2_offset],blitz::shape(3),blitz::neverDeleteData);
        blitz::Array<double,1> D2_next(&D_next_vals[k*3*3+D2_offset],blitz::shape(3),blitz::neverDeleteData);
        blitz::Array<double,1> D2_prev(&D_prev_vals[k*3*3+D2_offset],blitz::shape(3),blitz::neverDeleteData);

        const int D3_offset = 6;
        blitz::Array<double,1> D3(&D_vals[(petsc_mastr_node_idxs[k]-global_offset)*3*3+D3_offset],blitz::shape(3),blitz::neverDeleteData);
        blitz::Array<double,1> D3_next(&D_next_vals[k*3*3+D3_offset],blitz::shape(3),blitz::neverDeleteData);
        blitz::Array<double,1> D3_prev(&D_prev_vals[k*3*3+D3_offset],blitz::shape(3),blitz::neverDeleteData);

        blitz::Array<double,1> X(&X_vals[(petsc_mastr_node_idxs[k]-global_offset)*NDIM],blitz::shape(NDIM),blitz::neverDeleteData);
        blitz::Array<double,1> X_next(&X_next_vals[k*NDIM],blitz::shape(NDIM),blitz::neverDeleteData);
        blitz::Array<double,1> X_prev(&X_prev_vals[k*NDIM],blitz::shape(NDIM),blitz::neverDeleteData);

        blitz::Array<double,1> F_plus_half(3), N_plus_half(3);
        compute_force_and_torque(F_plus_half, N_plus_half, X, X_next, D1, D1_next, D2, D2_next, D3, D3_next);

        blitz::Array<double,1> F_minus_half(3), N_minus_half(3);
        compute_force_and_torque(F_minus_half, N_minus_half, X_prev, X, D1_prev, D1, D2_prev, D2, D3_prev, D3);

        F = F_plus_half-F_minus_half;
        N = N_plus_half-N_minus_half + 0.5*(cross(blitz::Array<double,1>(X_next-X),F_plus_half) + cross(blitz::Array<double,1>(X-X_prev),F_minus_half));
    }

    ierr = VecRestoreArray(D_vec, &D_vals);            IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_next_vec, &D_next_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_next_vec);                     IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(D_prev_vec, &D_prev_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(D_prev_vec);                     IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_vec, &X_vals);            IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_next_vec, &X_next_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(X_next_vec);                     IBTK_CHKERRQ(ierr);

    ierr = VecRestoreArray(X_prev_vec, &X_prev_vals);  IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(X_prev_vec);                     IBTK_CHKERRQ(ierr);

    Vec F_vec = F_data->getGlobalVec();
    ierr = IBTK::PETScVecOps::VecSetValuesBlocked(F_vec,
                                                  petsc_mastr_node_idxs.size(),
                                                  !petsc_mastr_node_idxs.empty() ? &petsc_mastr_node_idxs[0] : PETSC_NULL,
                                                  !petsc_mastr_node_idxs.empty() ? &          F_node_vals[0] : PETSC_NULL,
                                                  ADD_VALUES);  IBTK_CHKERRQ(ierr);

    Vec N_vec = N_data->getGlobalVec();
    ierr = IBTK::PETScVecOps::VecSetValuesBlocked(N_vec,
                                                  petsc_mastr_node_idxs.size(),
                                                  !petsc_mastr_node_idxs.empty() ? &petsc_mastr_node_idxs[0] : PETSC_NULL,
                                                  !petsc_mastr_node_idxs.empty() ? &          N_node_vals[0] : PETSC_NULL,
                                                  ADD_VALUES);  IBTK_CHKERRQ(ierr);

    ierr = IBTK::PETScVecOps::VecAssemblyBegin(F_vec);  IBTK_CHKERRQ(ierr);
    ierr = IBTK::PETScVecOps::VecAssemblyBegin(N_vec);  IBTK_CHKERRQ(ierr);
    ierr = IBTK::PETScVecOps::VecAssemblyEnd(F_vec);    IBTK_CHKERRQ(ierr);
    ierr = IBTK::PETScVecOps::VecAssemblyEnd(N_vec);    IBTK_CHKERRQ(ierr);

    t_compute_lagrangian_force_and_torque->stop();
    return;
}// computeLagrangianForceAndTorque

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBKirchhoffRodForceGen::getFromInput(
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
template class SAMRAI::tbox::Pointer<IBAMR::IBKirchhoffRodForceGen>;

//////////////////////////////////////////////////////////////////////////////
