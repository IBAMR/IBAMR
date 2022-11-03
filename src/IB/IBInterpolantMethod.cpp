// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBInterpolantMethod.h"

#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/LInitStrategy.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "CartesianGridGeometry.h"
#include "CellVariable.h"
#include "CoarsenSchedule.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include "petscvec.h"
#include <petscsys.h>

#include "ibamr/namespaces.h" // IWYU pragma: keep

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

IBTK_DISABLE_EXTRA_WARNINGS
#include <Eigen/Core>
#include <Eigen/Geometry>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace SAMRAI
{
namespace xfer
{
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBInterpolantMethod restart file data.
static const int IB_INTERPOLANT_METHOD_VERSION = 1;

inline void
set_rotation_matrix(const EigenAlignedVector<Eigen::Vector3d>& rot_vel,
                    const EigenAlignedVector<Eigen::Quaterniond>& q_old,
                    EigenAlignedVector<Eigen::Quaterniond>& q_new,
                    EigenAlignedVector<Eigen::Matrix3d>& rot_mat,
                    const double dt)
{
    unsigned n_structs = (unsigned)rot_mat.size();
    for (unsigned struct_no = 0; struct_no < n_structs; ++struct_no)
    {
        const double norm = rot_vel[struct_no].norm();
        if (!IBTK::abs_equal_eps(norm, 0.0))
        {
            Eigen::Vector3d rot_axis = rot_vel[struct_no] / norm;
            Eigen::Quaterniond q(Eigen::AngleAxisd(norm * dt, rot_axis));
            q_new[struct_no] = (q.normalized() * q_old[struct_no]).normalized();
        }
        else
        {
            q_new[struct_no] = q_old[struct_no];
        }

        rot_mat[struct_no] = q_new[struct_no].toRotationMatrix();
    }
    return;
} // set_rotation_matrix

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInterpolantMethod::IBInterpolantMethod(std::string object_name,
                                         Pointer<Database> input_db,
                                         int no_structures,
                                         bool register_for_restart)
    : d_num_rigid_parts(no_structures),
      d_object_name(std::move(object_name)),
      d_registered_for_restart(register_for_restart)
{
    // Register object with the restart manager.
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Set some default values.
    d_struct_lag_idx_range.resize(d_num_rigid_parts);
    d_quaternion_current.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());
    d_quaternion_new.resize(d_num_rigid_parts, Eigen::Quaterniond::Identity());
    d_center_of_mass_initial.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_center_of_mass_current.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());
    d_center_of_mass_new.resize(d_num_rigid_parts, Eigen::Vector3d::Zero());

    // Set some default values.
    d_ghosts = std::max(LEInteractor::getMinimumGhostWidth(d_interp_kernel_fcn),
                        LEInteractor::getMinimumGhostWidth(d_spread_kernel_fcn));

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Check the choices for the kernel function.
    if (d_interp_kernel_fcn != d_spread_kernel_fcn)
    {
        pout << "WARNING: different kernel functions are being used for  "
                "interpolation and "
                "spreading.\n";
    }

    // Get the Lagrangian Data Manager.
    d_l_data_manager = LDataManager::getManager(d_object_name + "::LDataManager",
                                                d_interp_kernel_fcn,
                                                d_spread_kernel_fcn,
                                                d_error_if_points_leave_domain,
                                                d_ghosts,
                                                d_registered_for_restart);
    d_ghosts = d_l_data_manager->getGhostCellWidth();

    return;
} // IBInterpolantMethod

IBInterpolantMethod::~IBInterpolantMethod()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~IBInterpolantMethod

void
IBInterpolantMethod::registerEulerianVariables()
{
    const IntVector<NDIM> ib_ghosts = getMinimumGhostCellWidth();
    for (auto& q_pair : d_q_interp_idx)
    {
        const std::string& var_name = q_pair.first;
        q_pair.second = -1;
        registerVariable(q_pair.second, d_q_var[var_name], ib_ghosts, NULL /*d_ib_solver->getScratchContext()*/);
    }

    return;
} // registerEulerianVariables

void
IBInterpolantMethod::registerEulerianCommunicationAlgorithms()
{
    Pointer<RefineAlgorithm<NDIM> > ghost_fill_alg = new RefineAlgorithm<NDIM>();
    Pointer<RefineOperator<NDIM> > refine_op = NULL;

    for (const auto& q_pair : d_q_interp_idx)
    {
        const int q_interp_idx = q_pair.second;
        ghost_fill_alg->registerRefine(q_interp_idx, q_interp_idx, q_interp_idx, refine_op);
    }
    registerGhostfillRefineAlgorithm(d_object_name + "::ghost_fill_alg", ghost_fill_alg);

    return;
} // registerEulerianCommunicationAlgorithms

void
IBInterpolantMethod::registerVariableAndHierarchyIntegrator(const std::string& var_name,
                                                            const int var_depth,
                                                            Pointer<Variable<NDIM> > var,
                                                            Pointer<HierarchyIntegrator> hier_integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_q_var.find(var_name) == d_q_var.end());
    TBOX_ASSERT(d_q_depth.find(var_name) == d_q_depth.end());
    TBOX_ASSERT(d_q_hier_integrator.find(var_name) == d_q_hier_integrator.end());
#endif
    d_q_var[var_name] = var;
    d_q_depth[var_name] = var_depth;
    d_q_hier_integrator[var_name] = hier_integrator;
    d_q_interp_idx[var_name] = -1;
    d_Q_current_data[var_name] = {};
    d_Q_new_data[var_name] = {};

    return;
} // registerVariableAndHierarchyIntegrator

void
IBInterpolantMethod::registerLInitStrategy(Pointer<LInitStrategy> l_initializer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(l_initializer);
#endif
    d_l_initializer = l_initializer;
    d_l_data_manager->registerLInitStrategy(d_l_initializer);
    return;
} // registerLInitStrategy

void
IBInterpolantMethod::freeLInitStrategy()
{
    d_l_initializer.setNull();
    d_l_data_manager->freeLInitStrategy();
    return;
} // freeLInitStrategy

LDataManager*
IBInterpolantMethod::getLDataManager() const
{
    return d_l_data_manager;
} // getLDataManager

int
IBInterpolantMethod::getStructuresLevelNumber() const
{
    return d_hierarchy->getFinestLevelNumber();

} // getStructuresLevelNumber

int
IBInterpolantMethod::getStructureHandle(const int lag_idx) const
{
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        const std::pair<int, int>& lag_idx_range = d_struct_lag_idx_range[struct_no];
        if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second) return struct_no;
    }

    return -1;
} // getStructureHandle

void
IBInterpolantMethod::registerLSiloDataWriter(Pointer<LSiloDataWriter> silo_writer)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(silo_writer);
#endif
    d_silo_writer = silo_writer;
    d_l_data_manager->registerLSiloDataWriter(d_silo_writer);
    return;
} // registerLSiloDataWriter

const IntVector<NDIM>&
IBInterpolantMethod::getMinimumGhostCellWidth() const
{
    return d_ghosts;
} // getMinimumGhostCellWidth

void
IBInterpolantMethod::setupTagBuffer(Array<int>& tag_buffer, Pointer<GriddingAlgorithm<NDIM> > gridding_alg) const
{
    const int finest_hier_ln = gridding_alg->getMaxLevels() - 1;
    const int tsize = tag_buffer.size();
    tag_buffer.resizeArray(finest_hier_ln);
    for (int i = tsize; i < finest_hier_ln; ++i) tag_buffer[i] = 0;
    const int gcw = d_ghosts.max();
    for (int tag_ln = 0; tag_ln < finest_hier_ln; ++tag_ln)
    {
        const int data_ln = tag_ln + 1;
        const int can_be_refined = data_ln < finest_hier_ln;
        if (!d_l_initializer->getLevelHasLagrangianData(data_ln, can_be_refined)) continue;
        tag_buffer[tag_ln] = std::max(tag_buffer[tag_ln], gcw);
    }
    for (int ln = finest_hier_ln - 2; ln >= 0; --ln)
    {
        tag_buffer[ln] =
            std::max(tag_buffer[ln], tag_buffer[ln + 1] / gridding_alg->getRatioToCoarserLevel(ln + 1).max() + 1);
    }
    return;
} // setupTagBuffer

void
IBInterpolantMethod::preprocessIntegrateData(double current_time, double new_time, int /*num_cycles*/)
{
    d_current_time = current_time;
    d_new_time = new_time;
    d_half_time = current_time + 0.5 * (new_time - current_time);

    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Look-up or allocate Lagangian data.
    d_X_current_data.resize(finest_ln + 1);
    d_X_new_data.resize(finest_ln + 1);
    for (auto& Q_pair : d_Q_current_data)
    {
        const std::string& name = Q_pair.first;
        Q_pair.second.resize(finest_ln + 1);
        d_Q_new_data[name].resize(finest_ln + 1);
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_X_current_data[ln] = d_l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
        d_X_new_data[ln] = d_l_data_manager->createLData("X_new", ln, NDIM);

        for (auto& Q_pair : d_Q_current_data)
        {
            const std::string& name = Q_pair.first;
            Q_pair.second[ln] = d_l_data_manager->getLData(name, ln);
            d_Q_new_data[name][ln] = d_l_data_manager->createLData(name + "_new", ln, d_q_depth[name]);
        }

        // Initialize X^{n+1} to equal X^{n}, and initialize Q^{n+1} to equal Q^{n}.
        ierr = VecCopy(d_X_current_data[ln]->getVec(), d_X_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        for (auto& Q_pair : d_Q_current_data)
        {
            const std::string& name = Q_pair.first;
            ierr = VecCopy(Q_pair.second[ln]->getVec(), d_Q_new_data[name][ln]->getVec());
            IBTK_CHKERRQ(ierr);
        }
    }

    return;
} // preprocessIntegrateData

void
IBInterpolantMethod::postprocessIntegrateData(double /*current_time*/, double /*new_time*/, int /*num_cycles*/)
{
    int ierr;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Reset time-dependent Lagrangian data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecSwap(d_X_current_data[ln]->getVec(), d_X_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
        for (auto& Q_pair : d_Q_current_data)
        {
            const std::string& name = Q_pair.first;
            ierr = VecSwap(Q_pair.second[ln]->getVec(), d_Q_new_data[name][ln]->getVec());
            IBTK_CHKERRQ(ierr);
        }
    }

    // Deallocate Lagrangian scratch data.
    d_X_current_data.clear();
    d_X_new_data.clear();
    for (auto& Q_pair : d_Q_current_data)
    {
        const std::string& name = Q_pair.first;
        Q_pair.second.clear();
        d_Q_new_data[name].clear();
    }

    // New state becomes current state for the next timestep.
    d_center_of_mass_current = d_center_of_mass_new;
    d_quaternion_current = d_quaternion_new;

    // Reset the current time step interval.
    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();
    d_half_time = std::numeric_limits<double>::quiet_NaN();
    return;
} // postprocessIntegrateData

void
IBInterpolantMethod::interpolateVelocity(const int /*u_data_idx*/,
                                         const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
                                         const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
                                         const double /*data_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::interpolateVelocity(). This method is not implemented." << std::endl);
    return;
} // interpolateVelocity

void
IBInterpolantMethod::interpolateQ(const double data_time)
{
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> >* X_data;
    getPositionData(&X_data, data_time);
    std::vector<Pointer<LData> >* Q_data;

    for (auto& Q_pair : d_Q_current_data)
    {
        const std::string& name = Q_pair.first;
        int q_data_idx = d_q_interp_idx[name];
        copyEulerianDataFromIntegrator(name, q_data_idx, data_time);
        getQData(name, &Q_data, data_time);
        d_l_data_manager->interp(q_data_idx,
                                 *Q_data,
                                 *X_data,
                                 std::vector<Pointer<CoarsenSchedule<NDIM> > >(finest_ln + 1, NULL),
                                 getGhostfillRefineSchedules(d_object_name + "::ghost_fill_alg"),
                                 data_time);
    }

    return;
} // interpolateQ

void
IBInterpolantMethod::interpolateQ()
{
    interpolateQ(d_current_time);
    interpolateQ(d_new_time);

    return;
} // interpolateQ

void
IBInterpolantMethod::forwardEulerStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::forwardEulerStep(). This method is not implemented." << std::endl);
    return;
} // forwardEulerStep

void
IBInterpolantMethod::backwardEulerStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::backwardEulerStep(). This method is not implemented." << std::endl);
    return;
} // backwardEulerStep

void
IBInterpolantMethod::midpointStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::midpointStep(). This method is not implemented." << std::endl);
    return;
} // midpointStep

void
IBInterpolantMethod::trapezoidalStep(const double /*current_time*/, const double /*new_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::trapezoidalStep(). This method is not implemented." << std::endl);
    return;
} // trapezoidalStep

void
IBInterpolantMethod::updateMeshPosition(double current_time,
                                        double new_time,
                                        const EigenAlignedVector<Eigen::Vector3d>& U,
                                        const EigenAlignedVector<Eigen::Vector3d>& W)
{
    const double dt = new_time - current_time;

    // Fill the rotation matrix of structures with rotation angle (W^n+1)*dt.
    EigenAlignedVector<Eigen::Matrix3d> rotation_mat(d_num_rigid_parts, Eigen::Matrix3d::Identity(3, 3));
    set_rotation_matrix(W, d_quaternion_current, d_quaternion_new, rotation_mat, dt);

    // Rotate the body with new rotational velocity about origin
    // and translate the body to newer position.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        boost::multi_array_ref<double, 2>& X_new_array = *d_X_new_data[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X0_array =
            *(d_l_data_manager->getLData("X0", ln)->getLocalFormVecArray());
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const unsigned structs_on_this_ln = (unsigned)structIDs.size();
#if !defined(NDEBUG)
        TBOX_ASSERT(structs_on_this_ln == d_num_rigid_parts);
#endif

        for (const auto& node_idx : local_nodes)
        {
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            double* const X_new = &X_new_array[local_idx][0];
            const double* const X0 = &X0_array[local_idx][0];
            Eigen::Vector3d dr = Eigen::Vector3d::Zero();

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = getStructureHandle(lag_idx);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                dr[d] = X0[d] - d_center_of_mass_initial[struct_handle][d];
            }

            // Rotate dr vector using the rotation matrix.
            const Eigen::Vector3d R_dr = rotation_mat[struct_handle] * dr;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X_new[d] = d_center_of_mass_current[struct_handle][d] + R_dr[d] + dt * U[struct_handle][d];
            }
        }
        d_X_new_data[ln]->restoreArrays();
        d_l_data_manager->getLData("X0", ln)->restoreArrays();
    }

    // Compute the new center of mass.
    for (unsigned struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_center_of_mass_new[struct_no][d] = d_center_of_mass_current[struct_no][d] + dt * U[struct_no][d];
        }
    }

    return;
} // updateMeshPosition

void
IBInterpolantMethod::computeLagrangianForce(double /*data_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::computeLagrangianForce(). This method is not implemented." << std::endl);
    return;
} // computeLagrangianForce

void
IBInterpolantMethod::spreadForce(const int /*f_data_idx*/,
                                 RobinPhysBdryPatchStrategy* /*f_phys_bdry_op*/,
                                 const std::vector<Pointer<RefineSchedule<NDIM> > >& /*f_prolongation_scheds*/,
                                 const double /*data_time*/)
{
    TBOX_ERROR("IBInterpolantMethod::spreadForce(). This method is not implemented." << std::endl);
    return;
} // spreadForce

void
IBInterpolantMethod::spreadQ(double data_time)
{
    int finest_ln = d_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = d_hierarchy->getPatchLevel(finest_ln);
    const IntVector<NDIM>& ratio = finest_level->getRatio();
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* dx0 = grid_geom->getDx();
    std::array<double, NDIM> dx;
    PetscScalar vol = 1.0;
    for (unsigned i = 0; i < NDIM; ++i)
    {
        dx[i] = dx0[i] / ratio(i);
        vol *= dx[i];
    }

    std::vector<Pointer<LData> >* X_data;
    getPositionData(&X_data, data_time);
    std::vector<Pointer<LData> >* Q_data;

    for (auto& Q_pair : d_Q_current_data)
    {
        const std::string& name = Q_pair.first;
        int q_data_idx = d_q_interp_idx[name];
        zeroOutEulerianData(name, q_data_idx);

        getQData(name, &Q_data, data_time);
        Vec l_data_vec = (*Q_data)[finest_ln]->getVec();
        VecScale(l_data_vec, vol);
        d_l_data_manager->spread(q_data_idx, *Q_data, *X_data, (RobinPhysBdryPatchStrategy*)NULL);
        VecScale(l_data_vec, 1.0 / vol);
    }

    return;
} // spreadQ

void
IBInterpolantMethod::copyEulerianDataToIntegrator(double data_time)
{
    for (auto& Q_pair : d_Q_current_data)
    {
        const std::string& name = Q_pair.first;
        int q_data_idx = d_q_interp_idx[name];
        copyEulerianDataToIntegrator(name, q_data_idx, data_time);
    }
    return;
} // copyEulerianDataToIntegrator

void
IBInterpolantMethod::initializePatchHierarchy(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg,
    int /*u_data_idx*/,
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& /*u_synch_scheds*/,
    const std::vector<Pointer<RefineSchedule<NDIM> > >& /*u_ghost_fill_scheds*/,
    int /*integrator_step*/,
    double /*init_data_time*/,
    bool initial_time)
{
    // Cache pointers to the patch hierarchy and gridding algorithm.
    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;

    // Set structure index info.
    const int struct_ln = getStructuresLevelNumber();
    std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(struct_ln);
    std::sort(structIDs.begin(), structIDs.end());
    const auto structs_on_this_ln = static_cast<unsigned>(structIDs.size());

    for (unsigned struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
    {
        d_struct_lag_idx_range[struct_no] =
            d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], struct_ln);
    }

    // Initialize initial center of mass of structures.
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<Pointer<LData> > X0_data_vec(finest_ln + 1, Pointer<LData>(NULL));
    X0_data_vec[finest_ln] = d_l_data_manager->getLData("X0", finest_ln);
    computeCenterOfMass(d_center_of_mass_initial, X0_data_vec);

    if (initial_time)
    {
        d_center_of_mass_current = d_center_of_mass_initial;
    }

    return;
} // initializePatchHierarchy

void
IBInterpolantMethod::addWorkloadEstimate(Pointer<PatchHierarchy<NDIM> > hierarchy, const int workload_data_idx)
{
    d_l_data_manager->addWorkloadEstimate(hierarchy, workload_data_idx);
    return;
} // addWorkloadEstimate

void IBInterpolantMethod::beginDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                  Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->beginDataRedistribution();
    return;
} // beginDataRedistribution

void IBInterpolantMethod::endDataRedistribution(Pointer<PatchHierarchy<NDIM> > /*hierarchy*/,
                                                Pointer<GriddingAlgorithm<NDIM> > /*gridding_alg*/)
{
    d_l_data_manager->endDataRedistribution();
    return;
} // endDataRedistribution

void
IBInterpolantMethod::initializeLevelData(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                         int level_number,
                                         double init_data_time,
                                         bool can_be_refined,
                                         bool initial_time,
                                         Pointer<BasePatchLevel<NDIM> > old_level,
                                         bool allocate_data)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->setPatchLevels(0, finest_hier_level);
    d_l_data_manager->initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);
    if (initial_time && d_l_data_manager->levelContainsLagrangianData(level_number))
    {
        for (const auto& Q_pair : d_Q_current_data)
        {
            const std::string& Q_name = Q_pair.first;
            const int Q_depth = d_q_depth[Q_name];
            Pointer<LData> Q_data = d_l_data_manager->createLData(Q_name, level_number, Q_depth, /*manage_data*/ true);
        }
    }

    return;
} // initializeLevelData

void
IBInterpolantMethod::resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
                                                 int coarsest_level,
                                                 int finest_level)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_l_data_manager->setPatchHierarchy(hierarchy);
    d_l_data_manager->setPatchLevels(0, finest_hier_level);
    d_l_data_manager->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    return;
} // resetHierarchyConfiguration

void
IBInterpolantMethod::applyGradientDetector(Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
                                           int level_number,
                                           double error_data_time,
                                           int tag_index,
                                           bool initial_time,
                                           bool uses_richardson_extrapolation_too)
{
    Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif

    // Tag cells that contain Lagrangian nodes.
    d_l_data_manager->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
} // applyGradientDetector

void
IBInterpolantMethod::putToDatabase(Pointer<Database> db)
{
    db->putInteger("IB_INTERPOLANT_METHOD_VERSION", IB_INTERPOLANT_METHOD_VERSION);
    db->putString("d_interp_kernel_fcn", d_interp_kernel_fcn);
    db->putString("d_spread_kernel_fcn", d_spread_kernel_fcn);
    db->putIntegerArray("d_ghosts", d_ghosts, NDIM);

    const int finest_hier_level = d_hierarchy->getFinestLevelNumber();
    db->putInteger("finest_hier_level", finest_hier_level);

    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        std::ostringstream C, Q;
        C << "C_" << struct_no;
        Q << "Q_" << struct_no;

        double Q_coeffs[4] = { d_quaternion_current[struct_no].w(),
                               d_quaternion_current[struct_no].x(),
                               d_quaternion_current[struct_no].y(),
                               d_quaternion_current[struct_no].z() };

        db->putDoubleArray(C.str(), &d_center_of_mass_current[struct_no][0], 3);
        db->putDoubleArray(Q.str(), &Q_coeffs[0], 4);
    }

    return;
} // putToDatabase

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBInterpolantMethod::getPositionData(std::vector<Pointer<LData> >** X_data, double data_time)
{
    if (IBTK::rel_equal_eps(data_time, d_current_time))
    {
        *X_data = &d_X_current_data;
    }
    else if (IBTK::rel_equal_eps(data_time, d_new_time))
    {
        *X_data = &d_X_new_data;
    }
    else
    {
        TBOX_ERROR(
            "IBInterpolantMethod::getPositionData() Structure position inquired at times other than current and new. "
            "\n");
    }
    return;
} // getPositionData

void
IBInterpolantMethod::copyEulerianDataFromIntegrator(const std::string& var_name, int q_data_idx, double data_time)
{
    int q_hier_idx = invalid_index;

    Pointer<HierarchyIntegrator> hier_integrator = d_q_hier_integrator[var_name];
    Pointer<Variable<NDIM> > var = d_q_var[var_name];

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    if (IBTK::rel_equal_eps(data_time, d_current_time))
    {
        q_hier_idx = var_db->mapVariableAndContextToIndex(var, hier_integrator->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(data_time, d_new_time))
    {
        q_hier_idx = var_db->mapVariableAndContextToIndex(var, hier_integrator->getNewContext());
    }
    else
    {
        TBOX_ERROR(
            "IBInterpolantMethod::getEulerianDataForInterpolation() Interpolation not supported at times other than "
            "current and new.\n");
    }

    // Copy integrator data into q_data_idx
    Pointer<CellVariable<NDIM, double> > cc_var = var;
    Pointer<SideVariable<NDIM, double> > sc_var = var;
    if (cc_var)
    {
        HierarchyCellDataOpsReal<NDIM, double> hier_data_ops(d_hierarchy);
        hier_data_ops.copyData(q_data_idx, q_hier_idx);
    }

    if (sc_var)
    {
        HierarchySideDataOpsReal<NDIM, double> hier_data_ops(d_hierarchy);
        hier_data_ops.copyData(q_data_idx, q_hier_idx);
    }

    return;
} // copyEulerianDataFromIntegrator

void
IBInterpolantMethod::zeroOutEulerianData(const std::string& var_name, int q_data_idx)
{
    Pointer<Variable<NDIM> > var = d_q_var[var_name];

    Pointer<CellVariable<NDIM, double> > cc_var = var;
    Pointer<SideVariable<NDIM, double> > sc_var = var;
    if (cc_var)
    {
        HierarchyCellDataOpsReal<NDIM, double> hier_data_ops(d_hierarchy);
        hier_data_ops.setToScalar(q_data_idx, 0.0);
    }

    if (sc_var)
    {
        HierarchySideDataOpsReal<NDIM, double> hier_data_ops(d_hierarchy);
        hier_data_ops.setToScalar(q_data_idx, 0.0);
    }

    return;
} // zeroOutEulerianData

void
IBInterpolantMethod::copyEulerianDataToIntegrator(const std::string& var_name, int q_data_idx, double data_time)
{
    int q_hier_idx = invalid_index;

    Pointer<HierarchyIntegrator> hier_integrator = d_q_hier_integrator[var_name];
    Pointer<Variable<NDIM> > var = d_q_var[var_name];

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    if (IBTK::rel_equal_eps(data_time, d_current_time))
    {
        q_hier_idx = var_db->mapVariableAndContextToIndex(var, hier_integrator->getCurrentContext());
    }
    else if (IBTK::rel_equal_eps(data_time, d_new_time))
    {
        q_hier_idx = var_db->mapVariableAndContextToIndex(var, hier_integrator->getNewContext());
    }
    else
    {
        TBOX_ERROR(
            "IBInterpolantMethod::getEulerianDataForInterpolation() Interpolation not supported at times other than "
            "current and new.\n");
    }

    // Copy integrator data into q_data_idx
    Pointer<CellVariable<NDIM, double> > cc_var = var;
    Pointer<SideVariable<NDIM, double> > sc_var = var;
    if (cc_var)
    {
        HierarchyCellDataOpsReal<NDIM, double> hier_data_ops(d_hierarchy);
        hier_data_ops.copyData(q_hier_idx, q_data_idx);
    }

    if (sc_var)
    {
        HierarchySideDataOpsReal<NDIM, double> hier_data_ops(d_hierarchy);
        hier_data_ops.copyData(q_hier_idx, q_data_idx);
    }

    return;
} // copyEulerianDataToIntegrator

void
IBInterpolantMethod::getQData(const std::string& var_name, std::vector<Pointer<LData> >** Q_data, double data_time)
{
    if (IBTK::rel_equal_eps(data_time, d_current_time))
    {
        *Q_data = &d_Q_current_data[var_name];
    }
    else if (IBTK::rel_equal_eps(data_time, d_new_time))
    {
        *Q_data = &d_Q_new_data[var_name];
    }
    else
    {
        TBOX_ERROR("IBInterpolantMethod::getQData() Q data inquired at times other than current and new.\n");
    }
    return;
} // getQData

void
IBInterpolantMethod::computeCenterOfMass(EigenAlignedVector<Eigen::Vector3d>& center_of_mass,
                                         std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Zero out the COM vector.
    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        center_of_mass[struct_no].setZero();
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& X_data_array = *X_data[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const unsigned structs_on_this_ln = static_cast<unsigned>(structIDs.size());
#if !defined(NDEBUG)
        TBOX_ASSERT(structs_on_this_ln == d_num_rigid_parts);
#endif

        for (const auto& node_idx : local_nodes)
        {
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_data_array[local_idx][0];

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = getStructureHandle(lag_idx);

            for (unsigned int d = 0; d < NDIM; ++d) center_of_mass[struct_handle][d] += X[d];
        }

        for (unsigned struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            IBTK_MPI::sumReduction(&center_of_mass[struct_no][0], NDIM);
            const int total_nodes = getNumberOfNodes(struct_no);
            center_of_mass[struct_no] /= total_nodes;
        }

        X_data[ln]->restoreArrays();
    }
    return;
} // computeCenterOfMass

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBInterpolantMethod::getFromInput(Pointer<Database> db, bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->isString("interp_kernel_fcn") && db->isString("spread_kernel_fcn"))
        {
            d_interp_kernel_fcn = db->getString("interp_kernel_fcn");
            d_spread_kernel_fcn = db->getString("spread_kernel_fcn");
        }
        if (db->isString("interp_delta_fcn") && db->isString("spread_delta_fcn"))
        {
            d_interp_kernel_fcn = db->getString("interp_delta_fcn");
            d_spread_kernel_fcn = db->getString("spread_delta_fcn");
        }
        else if (db->keyExists("delta_fcn"))
        {
            d_interp_kernel_fcn = db->getString("delta_fcn");
            d_spread_kernel_fcn = db->getString("delta_fcn");
        }
        else if (db->keyExists("kernel_fcn"))
        {
            d_interp_kernel_fcn = db->getString("kernel_fcn");
            d_spread_kernel_fcn = db->getString("kernel_fcn");
        }
        else if (db->keyExists("IB_delta_fcn"))
        {
            d_interp_kernel_fcn = db->getString("IB_delta_fcn");
            d_spread_kernel_fcn = db->getString("IB_delta_fcn");
        }
        else if (db->keyExists("IB_kernel_fcn"))
        {
            d_interp_kernel_fcn = db->getString("IB_kernel_fcn");
            d_spread_kernel_fcn = db->getString("IB_kernel_fcn");
        }

        if (db->isInteger("min_ghost_cell_width"))
        {
            d_ghosts = db->getInteger("min_ghost_cell_width");
        }
        else if (db->isDouble("min_ghost_cell_width"))
        {
            d_ghosts = static_cast<int>(std::ceil(db->getDouble("min_ghost_cell_width")));
        }
    }
    TBOX_ASSERT(LEInteractor::isKnownKernel(d_interp_kernel_fcn));
    TBOX_ASSERT(LEInteractor::isKnownKernel(d_spread_kernel_fcn));
    if (db->keyExists("error_if_points_leave_domain"))
        d_error_if_points_leave_domain = db->getBool("error_if_points_leave_domain");
    if (db->keyExists("do_log"))
        d_do_log = db->getBool("do_log");
    else if (db->keyExists("enable_logging"))
        d_do_log = db->getBool("enable_logging");
    return;
} // getFromInput

void
IBInterpolantMethod::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_INTERPOLANT_METHOD_VERSION");
    if (ver != IB_INTERPOLANT_METHOD_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    if (db->keyExists("d_interp_kernel_fcn"))
        d_interp_kernel_fcn = db->getString("d_interp_kernel_fcn");
    else if (db->keyExists("d_interp_delta_fcn"))
        d_interp_kernel_fcn = db->getString("d_interp_delta_fcn");
    if (db->keyExists("d_spread_kernel_fcn"))
        d_spread_kernel_fcn = db->getString("d_spread_kernel_fcn");
    else if (db->keyExists("d_spread_delta_fcn"))
        d_spread_kernel_fcn = db->getString("d_spread_delta_fcn");
    TBOX_ASSERT(LEInteractor::isKnownKernel(d_interp_kernel_fcn));
    TBOX_ASSERT(LEInteractor::isKnownKernel(d_spread_kernel_fcn));
    db->getIntegerArray("d_ghosts", d_ghosts, NDIM);

    for (unsigned int struct_no = 0; struct_no < d_num_rigid_parts; ++struct_no)
    {
        std::ostringstream C, Q;
        C << "C_" << struct_no;
        Q << "Q_" << struct_no;

        double Q_coeffs[4];
        db->getDoubleArray(C.str(), &d_center_of_mass_current[struct_no][0], 3);
        db->getDoubleArray(Q.str(), &Q_coeffs[0], 4);

        d_quaternion_current[struct_no].w() = Q_coeffs[0];
        d_quaternion_current[struct_no].x() = Q_coeffs[1];
        d_quaternion_current[struct_no].y() = Q_coeffs[2];
        d_quaternion_current[struct_no].z() = Q_coeffs[3];
        d_quaternion_current[struct_no].normalized();
    }

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
