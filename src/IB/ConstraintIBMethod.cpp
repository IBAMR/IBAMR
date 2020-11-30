// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include "ibamr/ConstraintIBKinematics.h"
#include "ibamr/ConstraintIBMethod.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/app_namespaces.h" // IWYU pragma: keep

#include "ibtk/CCLaplaceOperator.h"
#include "ibtk/CCPoissonPointRelaxationFACOperator.h"
#include "ibtk/FACPreconditioner.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"
#include "ibtk/LIndexSetData.h"
#include "ibtk/LMesh.h"
#include "ibtk/LNode.h"
#include "ibtk/LNodeSetData.h"
#include "ibtk/LSet.h"
#include "ibtk/LSetData.h"
#include "ibtk/LSetDataIterator.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_utilities.h"

#include "ArrayData.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellIterator.h"
#include "ComponentSelector.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchyDataOpsReal.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "VariableDatabase.h"
#include "VariableFillPattern.h"
#include "tbox/Array.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "petscvec.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <algorithm>
#include <limits>
#include <utility>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_postprocessSolveFluidEquation;
static Pointer<Timer> t_calculateCOMandMOIOfStructures;
static Pointer<Timer> t_interpolateFluidSolveVelocity;
static Pointer<Timer> t_calculateKinematicsVelocity;
static Pointer<Timer> t_calculateRigidMomentum;
static Pointer<Timer> t_correctVelocityOnLagrangianMesh;
static Pointer<Timer> t_spreadCorrectedLagrangianVelocity;
static Pointer<Timer> t_applyProjection;
static Pointer<Timer> t_eulerStep;
static Pointer<Timer> t_midpointStep;

// Number of ghost cells used for each variable quantity.
static const int CELLG = 1;

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string CELL_DATA_COARSEN_TYPE = "CUBIC_COARSEN";
static const std::string SIDE_DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

class find_struct_handle
{
private:
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    std::pair<int, int> struct_to_find_range;

public:
    find_struct_handle(std::pair<int, int> struct_range) : struct_to_find_range(std::move(struct_range))
    {
    }

    inline bool operator()(Pointer<ConstraintIBKinematics> ib_kinematics_ptr)
    {
        const StructureParameters& struct_param = ib_kinematics_ptr->getStructureParameters();
        const std::vector<std::pair<int, int> >& range = struct_param.getLagIdxRange();

        bool is_in_range = false;
        for (const auto& idx : range)
        {
            if (struct_to_find_range.first == idx.first && struct_to_find_range.second == idx.second)
            {
                is_in_range = true;
                break;
            }
        }

        return is_in_range;
    }
};

template <typename itr, typename T>
inline int
find_struct_handle_position(itr begin, itr end, const T& value)
{
    int position = 0;
    while (begin != end)
    {
        if (*begin++ == value)
            return position;
        else
            ++position;
    }

    TBOX_ASSERT(false);
    return -1;
}

#if (NDIM == 3)
// Routine to solve 3X3 equation to get rigid body rotational velocity.
inline void
solveSystemOfEqns(std::vector<double>& ang_mom, const Eigen::Matrix3d& inertiaTensor)
{
    const double a1 = inertiaTensor(0, 0), a2 = inertiaTensor(0, 1), a3 = inertiaTensor(0, 2), b1 = inertiaTensor(1, 0),
                 b2 = inertiaTensor(1, 1), b3 = inertiaTensor(1, 2), c1 = inertiaTensor(2, 0), c2 = inertiaTensor(2, 1),
                 c3 = inertiaTensor(2, 2), d1 = ang_mom[0], d2 = ang_mom[1], d3 = ang_mom[2];
    const double Dnr = (a3 * b2 * c1 - a2 * b3 * c1 - a3 * b1 * c2 + a1 * b3 * c2 + a2 * b1 * c3 - a1 * b2 * c3);
    ang_mom[0] = (b3 * c2 * d1 - b2 * c3 * d1 - a3 * c2 * d2 + a2 * c3 * d2 + a3 * b2 * d3 - a2 * b3 * d3) / Dnr;
    ang_mom[1] = -(b3 * c1 * d1 - b1 * c3 * d1 - a3 * c1 * d2 + a1 * c3 * d2 + a3 * b1 * d3 - a1 * b3 * d3) / Dnr;
    ang_mom[2] = (b2 * c1 * d1 - b1 * c2 * d1 - a2 * c1 * d2 + a1 * c2 * d2 + a2 * b1 * d3 - a1 * b2 * d3) / Dnr;
    return;
}
#endif
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

ConstraintIBMethod::ConstraintIBMethod(std::string object_name,
                                       Pointer<Database> input_db,
                                       const int no_structures,
                                       bool register_for_restart)
    : IBMethod(std::move(object_name), input_db, register_for_restart),
      d_no_structures(no_structures),
      d_ib_kinematics(d_no_structures, Pointer<ConstraintIBKinematics>(nullptr)),
      d_vol_element(d_no_structures, 0.0),
      d_vol_element_is_set(d_no_structures, false),
      d_structure_vol(d_no_structures, 0.0),
      d_structure_mom(d_no_structures, std::vector<double>(3, 0.0)),
      d_structure_rotational_mom(d_no_structures, std::vector<double>(3, 0.0)),
      d_rigid_trans_vel_current(d_no_structures, std::vector<double>(3, 0.0)),
      d_rigid_trans_vel_new(d_no_structures, std::vector<double>(3, 0.0)),
      d_rigid_rot_vel_current(d_no_structures, std::vector<double>(3, 0.0)),
      d_rigid_rot_vel_new(d_no_structures, std::vector<double>(3, 0.0)),
      d_incremented_angle_from_reference_axis(d_no_structures, std::vector<double>(3, 0.0)),
      d_vel_com_def_current(d_no_structures, std::vector<double>(3, 0.0)),
      d_vel_com_def_new(d_no_structures, std::vector<double>(3, 0.0)),
      d_omega_com_def_current(d_no_structures, std::vector<double>(3, 0.0)),
      d_omega_com_def_new(d_no_structures, std::vector<double>(3, 0.0)),
      d_center_of_mass_unshifted_current(d_no_structures, std::vector<double>(3, 0.0)),
      d_center_of_mass_unshifted_new(d_no_structures, std::vector<double>(3, 0.0)),
      d_center_of_mass_current(d_no_structures, std::vector<double>(3, 0.0)),
      d_center_of_mass_new(d_no_structures, std::vector<double>(3, 0.0)),
      d_moment_of_inertia_current(d_no_structures, Eigen::Matrix3d::Zero()),
      d_moment_of_inertia_new(d_no_structures, Eigen::Matrix3d::Zero()),
      d_tagged_pt_lag_idx(d_no_structures, 0),
      d_tagged_pt_position(d_no_structures, std::vector<double>(3, 0.0)),
      d_rho_solid(d_no_structures, std::numeric_limits<double>::quiet_NaN())
{
    // NOTE: Parent class constructor registers class with the restart manager,
    // sets object name.

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Setup the cell centered Poisson Solver needed for projection.
    if (d_needs_div_free_projection)
    {
        const std::string velcorrection_projection_prefix = "cIB_";
        // Setup the various solver components.
        for (int d = 0; d < NDIM; ++d)
        {
            d_velcorrection_projection_bc_coef.setBoundarySlope(2 * d, 0.0);
            d_velcorrection_projection_bc_coef.setBoundarySlope(2 * d + 1, 0.0);
        }

        d_velcorrection_projection_spec.reset(
            new PoissonSpecifications(d_object_name + "::ConstraintIBMethodProjection::Spec"));
        d_velcorrection_projection_op =
            new CCLaplaceOperator(d_object_name + "ConstraintIBMethodProjection::PoissonOperator", true);
        d_velcorrection_projection_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
        d_velcorrection_projection_op->setPhysicalBcCoef(&d_velcorrection_projection_bc_coef);

        d_velcorrection_projection_solver =
            new PETScKrylovPoissonSolver(d_object_name + "ConstraintIBMethodProjection::PoissonKrylovSolver",
                                         Pointer<Database>(nullptr),
                                         velcorrection_projection_prefix);
        d_velcorrection_projection_solver->setInitialGuessNonzero(false);
        d_velcorrection_projection_solver->setOperator(d_velcorrection_projection_op);

        if (d_velcorrection_projection_fac_pc_db.isNull())
        {
            TBOX_WARNING(d_object_name << "::ConstraintIBMethod():\n"
                                       << " ConstraintIBMethodProjection:: Poisson "
                                          "FAC PC solver database is null."
                                       << std::endl);
        }

        d_velcorrection_projection_fac_op = new CCPoissonPointRelaxationFACOperator(
            d_object_name + ":: ConstraintIBMethodProjection::PoissonFACOperator",
            d_velcorrection_projection_fac_pc_db,
            "");
        d_velcorrection_projection_fac_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
        d_velcorrection_projection_fac_pc =
            new IBTK::FACPreconditioner(d_object_name + "::ConstraintIBMethodProjection::PoissonPreconditioner",
                                        d_velcorrection_projection_fac_op,
                                        d_velcorrection_projection_fac_pc_db,
                                        "");
        d_velcorrection_projection_solver->setPreconditioner(d_velcorrection_projection_fac_pc);

        // Set some default options.
        d_velcorrection_projection_solver->setKSPType("gmres");
        d_velcorrection_projection_solver->setAbsoluteTolerance(1.0e-12);
        d_velcorrection_projection_solver->setRelativeTolerance(1.0e-08);
        d_velcorrection_projection_solver->setMaxIterations(25);

        // NOTE: We always use homogeneous Neumann boundary conditions for the
        // velocity correction projection Poisson solver.
        d_velcorrection_projection_solver->setNullspace(true);
    }
    else
    {
        d_velcorrection_projection_spec = nullptr;
        d_velcorrection_projection_op = nullptr;
        d_velcorrection_projection_fac_op = nullptr;
        d_velcorrection_projection_fac_pc = nullptr;
        d_velcorrection_projection_solver = nullptr;
    }

    // Do printing operation for processor 0 only.
    if (!IBTK_MPI::getRank() && d_print_output)
    {
        d_trans_vel_stream.resize(d_no_structures);
        d_rot_vel_stream.resize(d_no_structures);
        d_drag_force_stream.resize(d_no_structures);
        d_moment_of_inertia_stream.resize(d_no_structures);
        d_torque_stream.resize(d_no_structures);
        d_position_COM_stream.resize(d_no_structures);
        d_power_spent_stream.resize(d_no_structures);

        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            const std::string struct_no_str = std::to_string(struct_no);

            if (from_restart)
                d_trans_vel_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Trans_vel_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_trans_vel_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Trans_vel_struct_no_" + struct_no_str, std::fstream::out));

            if (from_restart)
                d_rot_vel_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Rot_vel_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_rot_vel_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Rot_vel_struct_no_" + struct_no_str, std::fstream::out));

            if (from_restart)
                d_drag_force_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Drag_force_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_drag_force_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Drag_force_struct_no_" + struct_no_str, std::fstream::out));

            if (from_restart)
                d_moment_of_inertia_stream[struct_no].reset(
                    new std::ofstream(d_base_output_filename + "_MOI_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_moment_of_inertia_stream[struct_no].reset(
                    new std::ofstream(d_base_output_filename + "_MOI_struct_no_" + struct_no_str, std::fstream::out));

            if (from_restart)
                d_torque_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Torque_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_torque_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Torque_struct_no_" + struct_no_str, std::fstream::out));

            if (from_restart)
                d_position_COM_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_COM_coordinates_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_position_COM_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_COM_coordinates_struct_no_" + struct_no_str, std::fstream::out));

            if (from_restart)
                d_power_spent_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Power_spent_struct_no_" + struct_no_str, std::fstream::app));
            else
                d_power_spent_stream[struct_no].reset(new std::ofstream(
                    d_base_output_filename + "_Power_spent_struct_no_" + struct_no_str, std::fstream::out));
        }

        // Output Eulerian momentum.
        std::string eul_mom_str = d_base_output_filename + "_EULERIAN_MOMENTUM.TXT";
        if (from_restart)
            d_eulerian_mom_stream.open(eul_mom_str.c_str(), std::fstream::app | std::fstream::out);
        else
            d_eulerian_mom_stream.open(eul_mom_str.c_str(), std::fstream::out);
    }

    // Setup the Timers.
    IBTK_DO_ONCE(
        t_postprocessSolveFluidEquation =
            TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::postprocessSolveFluidEquation()", true);
        t_calculateCOMandMOIOfStructures =
            TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::calculateCOMandMOIOfStructures()", true);
        t_interpolateFluidSolveVelocity =
            TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::interpolateFluidSolveVelocity()", true);
        t_calculateKinematicsVelocity =
            TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::calculateKinematicsVelocity()", true);
        t_calculateRigidMomentum =
            TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::calculateRigidMomentum()", true);
        t_correctVelocityOnLagrangianMesh =
            TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::correctVelocityOnLagrangianMesh()", true);
        t_spreadCorrectedLagrangianVelocity = TimerManager::getManager()->getTimer(
            "IBAMR::ConstraintIBMethod::spreadCorrectedLagrangianVelocity()", true);
        t_applyProjection = TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::applyProjection()", true);
        t_eulerStep = TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::eulerStep()", true);
        t_midpointStep = TimerManager::getManager()->getTimer("IBAMR::ConstraintIBMethod::midpointStep()", true););

    return;
} // ConstraintIBMethod

ConstraintIBMethod::~ConstraintIBMethod()
{
    // Deallocate the scratch fluid solve variable.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_u_fluidSolve_cib_idx)) level->deallocatePatchData(d_u_fluidSolve_cib_idx);
        if (!d_rho_is_const && level->checkAllocated(d_rho_scratch_idx)) level->deallocatePatchData(d_rho_scratch_idx);
    }

    return;
} //~ConstraintIBMethod

void
ConstraintIBMethod::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    IBMethod::preprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Call any registered pre fluid solve callback functions.
    for (unsigned i = 0; i < d_prefluidsolve_callback_fns.size(); ++i)
        d_prefluidsolve_callback_fns[i](current_time, new_time, cycle_num, d_prefluidsolve_callback_fns_ctx[i]);
    return;
}

void
ConstraintIBMethod::postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    IBMethod::postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    IBTK_TIMER_START(t_postprocessSolveFluidEquation);

    setCounter();

    IBTK_TIMER_START(t_calculateCOMandMOIOfStructures);
    calculateCOMandMOIOfStructures();
    IBTK_TIMER_STOP(t_calculateCOMandMOIOfStructures);

    IBTK_TIMER_START(t_interpolateFluidSolveVelocity);
    copyFluidVariable(d_u_fluidSolve_idx, d_u_fluidSolve_cib_idx);
    interpolateFluidSolveVelocity();
    IBTK_TIMER_STOP(t_interpolateFluidSolveVelocity);

    IBTK_TIMER_START(t_calculateKinematicsVelocity);
    calculateKinematicsVelocity();
    IBTK_TIMER_STOP(t_calculateKinematicsVelocity);

    IBTK_TIMER_START(t_calculateRigidMomentum);
    calculateRigidTranslationalMomentum();
    calculateRigidRotationalMomentum();
    IBTK_TIMER_STOP(t_calculateRigidMomentum);

    IBTK_TIMER_START(t_correctVelocityOnLagrangianMesh);
    correctVelocityOnLagrangianMesh();
    IBTK_TIMER_STOP(t_correctVelocityOnLagrangianMesh);

    IBTK_TIMER_START(t_spreadCorrectedLagrangianVelocity);
    spreadCorrectedLagrangianVelocity();
    IBTK_TIMER_STOP(t_spreadCorrectedLagrangianVelocity);

    if (d_needs_div_free_projection)
    {
        IBTK_TIMER_START(t_applyProjection);
        applyProjection();
        IBTK_TIMER_STOP(t_applyProjection);
    }

    if (d_output_drag) calculateDrag();
    if (d_output_torque) calculateTorque();
    if (d_output_eul_mom) calculateEulerianMomentum();
    if (d_output_power) calculatePower();
    if (d_calculate_structure_linear_mom) calculateStructureMomentum();
    if (d_calculate_structure_rotational_mom) calculateStructureRotationalMomentum();

    IBTK_TIMER_STOP(t_postprocessSolveFluidEquation);

    // call any other registered post fluid solve callback functions.
    for (unsigned i = 0; i < d_postfluidsolve_callback_fns.size(); ++i)
        d_postfluidsolve_callback_fns[i](current_time, new_time, cycle_num, d_postfluidsolve_callback_fns_ctx[i]);

    return;
}

void
ConstraintIBMethod::calculateEulerianMomentum()
{
    // Compute Eulerian momentum.
    std::vector<double> momentum(3, 0.0);
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    SAMRAIVectorReal<NDIM, double> wgt_sc("sc_wgt_original", d_hierarchy, coarsest_ln, finest_ln);
    wgt_sc.addComponent(getHierarchyMathOps()->getSideWeightVariable(),
                        getHierarchyMathOps()->getSideWeightPatchDescriptorIndex());

    for (int active = 0; active < NDIM; ++active)
    {
        Pointer<SAMRAIVectorReal<NDIM, double> > wgt_sc_active = wgt_sc.cloneVector("");
        wgt_sc_active->allocateVectorData();
        wgt_sc_active->copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&wgt_sc, false));

        // Zero out components other than active dimension.
        const int wgt_sc_active_idx = wgt_sc_active->getComponentDescriptorIndex(0);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > wgt_sc_active_data = patch->getPatchData(wgt_sc_active_idx);
                for (int d = 0; d < NDIM; ++d)
                {
                    if (d != active)
                    {
                        ArrayData<NDIM, double>& arraydata = wgt_sc_active_data->getArrayData(d);
                        arraydata.fill(0.0);
                    }
                }
            }
        }
        if (d_rho_is_const)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_rho_fluid >= 0.0);
#endif
            d_hier_sc_data_ops->scale(wgt_sc_active_idx, d_rho_fluid, wgt_sc_active_idx);
        }
        else
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_rho_ins_idx > 0);
#endif
            d_hier_sc_data_ops->multiply(wgt_sc_active_idx, d_rho_ins_idx, wgt_sc_active_idx);
        }

        momentum[active] = d_hier_sc_data_ops->dot(d_u_fluidSolve_idx, wgt_sc_active_idx);

        wgt_sc_active->freeVectorComponents();
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_eul_mom && (d_timestep_counter % d_output_interval) == 0)
    {
        d_eulerian_mom_stream << d_FuRMoRP_new_time << '\t' << momentum[0] << '\t' << momentum[1] << '\t' << momentum[2]
                              << '\t' << std::endl;
    }

    return;
} // calculateEulerianMomentum

void
ConstraintIBMethod::registerEulerianVariables()
{
    IBMethod::registerEulerianVariables();

    // Register a scratch fluid velocity variable with appropriate IB-width.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_u_fluidSolve_var = d_ib_solver->getVelocityVariable();
    Pointer<VariableContext> u_new_ctx = d_ib_solver->getNewContext();
    d_scratch_context = var_db->getContext(d_object_name + "::SCRATCH");
    d_u_fluidSolve_idx = var_db->mapVariableAndContextToIndex(d_u_fluidSolve_var, u_new_ctx);
    d_u_fluidSolve_cib_idx =
        var_db->registerVariableAndContext(d_u_fluidSolve_var, d_scratch_context, getMinimumGhostCellWidth());

    // Initialize  variables & variable contexts associated with projection step.
    if (d_needs_div_free_projection)
    {
        d_u_var = new SideVariable<NDIM, double>(d_object_name + "::u");
        d_Div_u_var = new CellVariable<NDIM, double>(d_object_name + "::Div_u");
        d_phi_var = new CellVariable<NDIM, double>(d_object_name + "::phi");
        const IntVector<NDIM> cell_ghosts = CELLG;
        d_phi_idx = var_db->registerVariableAndContext(d_phi_var, d_scratch_context, cell_ghosts);
        d_Div_u_scratch_idx = var_db->registerVariableAndContext(d_Div_u_var, d_scratch_context, cell_ghosts);
    }

    auto p_vc_ins_hier_integrator =
        dynamic_cast<INSVCStaggeredHierarchyIntegrator*>(IBStrategy::getINSHierarchyIntegrator());
    // If using constant rho INS solver,
    // then assert rho_fluid == rho_solid
    if (!p_vc_ins_hier_integrator)
    {
        d_rho_is_const = true;
        INSHierarchyIntegrator* p_ins_hier_integrator = IBStrategy::getINSHierarchyIntegrator();
        d_rho_fluid = p_ins_hier_integrator->getStokesSpecifications()->getRho();
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            if (d_rho_solid[struct_no] != d_rho_fluid)
            {
                TBOX_ERROR(d_object_name << "::registerEulerianVariables():\n"
                                         << "  for constant density cases, rho_solid[struct_no]\n"
                                         << "  must equal rho_fluid");
            }
        }
    }
    else
    {
        d_rho_is_const = p_vc_ins_hier_integrator->rhoIsConstant();
        if (d_rho_is_const)
        {
            d_rho_fluid = p_vc_ins_hier_integrator->getStokesSpecifications()->getRho();
            for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
            {
                if (d_rho_solid[struct_no] != d_rho_fluid)
                {
                    TBOX_ERROR(d_object_name << "::registerEulerianVariables():\n"
                                             << "  for constant density cases, rho_solid[struct_no]\n"
                                             << "  must equal rho_fluid");
                }
            }
        }
        else
        {
            // Get the density maintained by the integrator
            d_rho_ins_idx = p_vc_ins_hier_integrator->getLinearOperatorRhoPatchDataIndex();
#if !defined(NDEBUG)
            TBOX_ASSERT(d_rho_ins_idx >= 0);
#endif
            d_rho_var = new SideVariable<NDIM, double>(d_object_name + "::rho");
            d_rho_scratch_idx =
                var_db->registerVariableAndContext(d_rho_var, d_scratch_context, getMinimumGhostCellWidth());
        }
    }

    return;
} // registerEulerianVariables

void
ConstraintIBMethod::initializeHierarchyOperatorsandData()
{
    // Obtain the Hierarchy data operations objects
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM, double> > cc_var = new CellVariable<NDIM, double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, d_hierarchy, true);
    Pointer<SideVariable<NDIM, double> > sc_var = new SideVariable<NDIM, double>("sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(sc_var, d_hierarchy, true);
    d_wgt_cc_idx = getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();
    d_wgt_sc_idx = getHierarchyMathOps()->getSideWeightPatchDescriptorIndex();
    d_volume = getHierarchyMathOps()->getVolumeOfPhysicalDomain();

    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (!from_restart) calculateVolumeElement();
    setInitialLagrangianVelocity();

    return;
} // initializeHierarchyOperatorsandData

void
ConstraintIBMethod::registerConstraintIBKinematics(const std::vector<Pointer<ConstraintIBKinematics> >& ib_kinematics)
{
    if (ib_kinematics.size() != static_cast<unsigned int>(d_no_structures))
    {
        TBOX_ERROR(
            "ConstraintIBMethod::registerConstraintIBKinematics(). No of "
            "structures "
            << ib_kinematics.size() << " in vector passed to this method is not equal to no. of structures "
            << d_no_structures << " registered with this class" << std::endl);
    }
    else
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            if (ib_kinematics[struct_no].isNull())
                TBOX_ERROR("NULL ConstraintIBKinematics encountered in vector at " << struct_no << std::endl);
            else
                d_ib_kinematics[struct_no] = ib_kinematics[struct_no];
        }
    }

    // Get tagged point index info from objects.
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();
        d_tagged_pt_lag_idx[struct_no] = struct_param.getTaggedPtIdx();
    }
    return;

} // registerConstraintIBKinematics

void
ConstraintIBMethod::putToDatabase(Pointer<Database> db)
{
    IBMethod::putToDatabase(db);

    // Put the following quantities to restart database.
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const std::string struct_no_str = std::to_string(struct_no);

        db->putDoubleArray("POSN_COM_STRUCT_" + struct_no_str, &d_center_of_mass_current[struct_no][0], 3);

        db->putDoubleArray("VEL_COM_DEF_STRUCT_" + struct_no_str, &d_vel_com_def_current[struct_no][0], 3);
        db->putDoubleArray("OMEGA_COM_DEF_STRUCT_" + struct_no_str, &d_omega_com_def_current[struct_no][0], 3);

        db->putDoubleArray("VOL_ELEMENT_STRUCT_" + struct_no_str, &d_vol_element[0], d_no_structures);

        db->putDoubleArray("VOL_STRUCT_" + struct_no_str, &d_structure_vol[0], d_no_structures);

        db->putDoubleArray("VEL_COM_RIG_STRUCT_" + struct_no_str, &d_rigid_trans_vel_current[struct_no][0], 3);
        db->putDoubleArray("OMEGA_COM_RIG_STRUCT_" + struct_no_str, &d_rigid_rot_vel_current[struct_no][0], 3);

        db->putDoubleArray(
            "DELTA_THETA_STRUCT_" + struct_no_str, &d_incremented_angle_from_reference_axis[struct_no][0], 3);
    }
    db->putDouble("TIMESTEP_COUNTER", d_timestep_counter);
    db->putDouble("FuRMoRP_CURRENT_TIME", d_FuRMoRP_new_time);

    return;

} // putToDatabase

void
ConstraintIBMethod::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::preprocessIntegrateData(current_time, new_time, num_cycles);

    // Allocate memory for Lagrangian data.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    d_l_data_U_interp.resize(finest_ln + 1);
    d_l_data_U_correction.resize(finest_ln + 1);
    d_l_data_U_new.resize(finest_ln + 1);
    d_l_data_U_half.resize(finest_ln + 1);
    d_l_data_X_half_Euler.resize(finest_ln + 1);
    d_l_data_X_new_MidPoint.resize(finest_ln + 1);
    d_l_data_U_current.resize(finest_ln + 1);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        d_l_data_U_interp[ln] = d_l_data_manager->createLData(d_object_name + "interp_lag_vel", ln, NDIM, false);
        d_l_data_U_correction[ln] = d_l_data_manager->createLData(d_object_name + "correct_lag_vel", ln, NDIM, false);
        d_l_data_U_new[ln] = d_l_data_manager->createLData(d_object_name + "new_lag_vel", ln, NDIM, false);
        d_l_data_U_half[ln] = d_l_data_manager->createLData(d_object_name + "half_lag_vel", ln, NDIM, false);
        d_l_data_X_half_Euler[ln] = d_l_data_manager->createLData(d_object_name + "X_Euler", ln, NDIM, false);
        d_l_data_X_new_MidPoint[ln] = d_l_data_manager->createLData(d_object_name + "X_MidPoint", ln, NDIM, false);
        d_l_data_U_current[ln] = d_l_data_manager->createLData(d_object_name + "current_lag_vel", ln, NDIM, false);
    }

    // Compue the current Lagrangian velocity according to constraint for the
    // predictor Euler step.
    calculateCurrentLagrangianVelocity();
    return;

} // preprocessIntegrateData

void
ConstraintIBMethod::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    IBMethod::postprocessIntegrateData(current_time, new_time, num_cycles);

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        d_vel_com_def_current[struct_no] = d_vel_com_def_new[struct_no];
        d_omega_com_def_current[struct_no] = d_omega_com_def_new[struct_no];
        d_rigid_trans_vel_current[struct_no] = d_rigid_trans_vel_new[struct_no];
        d_rigid_rot_vel_current[struct_no] = d_rigid_rot_vel_new[struct_no];
    }

    // Deallocate memory for Lagrangian data.
    d_l_data_U_interp.clear();
    d_l_data_U_correction.clear();
    d_l_data_U_new.clear();
    d_l_data_U_half.clear();
    d_l_data_X_half_Euler.clear();
    d_l_data_X_new_MidPoint.clear();
    d_l_data_U_current.clear();

    return;
} // postprocessIntegrateData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
ConstraintIBMethod::getFromInput(Pointer<Database> input_db, const bool from_restart)
{
    // Read in control parameters from input database.
    d_needs_div_free_projection = input_db->getBoolWithDefault("needs_divfree_projection", d_needs_div_free_projection);
    input_db->getDoubleArray("rho_solid", &d_rho_solid[0], d_no_structures);
    d_rho_fluid = input_db->getDoubleWithDefault("rho_fluid", d_rho_fluid);
    d_calculate_structure_linear_mom =
        input_db->getBoolWithDefault("calculate_structure_linear_mom", d_calculate_structure_linear_mom);
    d_calculate_structure_rotational_mom =
        input_db->getBoolWithDefault("calculate_structure_rotational_mom", d_calculate_structure_rotational_mom);

    // Printing stuff to files.
    Pointer<Database> output_db = input_db->getDatabase("PrintOutput");
    d_print_output = output_db->getBoolWithDefault("print_output", d_print_output);
    d_output_interval = output_db->getIntegerWithDefault("output_interval", d_output_interval);
    d_output_drag = output_db->getBoolWithDefault("output_drag", d_output_drag);
    d_output_torque = output_db->getBoolWithDefault("output_torque", d_output_torque);
    d_output_power = output_db->getBoolWithDefault("output_power", d_output_power);
    d_output_trans_vel = output_db->getBoolWithDefault("output_rig_transvel", d_output_trans_vel);
    d_output_rot_vel = output_db->getBoolWithDefault("output_rig_rotvel", d_output_rot_vel);
    d_output_COM_coordinates = output_db->getBoolWithDefault("output_com_coords", d_output_COM_coordinates);
    d_output_MOI = output_db->getBoolWithDefault("output_moment_inertia", d_output_MOI);
    d_output_eul_mom = output_db->getBoolWithDefault("output_eulerian_mom", d_output_eul_mom);
    d_dir_name = output_db->getStringWithDefault("output_dirname", d_dir_name) + "/";
    d_base_output_filename = d_dir_name + output_db->getStringWithDefault("base_filename", d_base_output_filename);

    // Sanity check.
    if (d_output_eul_mom && !d_needs_div_free_projection)
        TBOX_WARNING(
            "WARNING ConstraintIBMethod::getFromInput() Eulerian momentum "
            "is calculated but divergence free projection "
            "is not active"
            << std::endl);

    if (!from_restart)
        tbox::Utilities::recursiveMkdir(d_dir_name);
    else
    {
        const bool restart_output_dump_dir_exists = output_db->keyExists("restart_output_dump_dir");
        if (restart_output_dump_dir_exists)
        {
            std::string restart_output_dump_dir = output_db->getString("restart_output_dump_dir");
            tbox::Utilities::recursiveMkdir(restart_output_dump_dir);
            d_base_output_filename = restart_output_dump_dir + "/" +
                                     output_db->getStringWithDefault("base_filename", d_base_output_filename);
        }
    }

    return;
} // getFromInput

void
ConstraintIBMethod::getFromRestart()
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

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const std::string struct_no_str = std::to_string(struct_no);

        db->getDoubleArray("POSN_COM_STRUCT_" + struct_no_str, &d_center_of_mass_current[struct_no][0], 3);

        db->getDoubleArray("VEL_COM_DEF_STRUCT_" + struct_no_str, &d_vel_com_def_current[struct_no][0], 3);
        db->getDoubleArray("OMEGA_COM_DEF_STRUCT_" + struct_no_str, &d_omega_com_def_current[struct_no][0], 3);

        db->getDoubleArray("VOL_ELEMENT_STRUCT_" + struct_no_str, &d_vol_element[0], d_no_structures);

        db->getDoubleArray("VOL_STRUCT_" + struct_no_str, &d_structure_vol[0], d_no_structures);

        db->getDoubleArray("VEL_COM_RIG_STRUCT_" + struct_no_str, &d_rigid_trans_vel_current[struct_no][0], 3);
        db->getDoubleArray("OMEGA_COM_RIG_STRUCT_" + struct_no_str, &d_rigid_rot_vel_current[struct_no][0], 3);

        db->getDoubleArray(
            "DELTA_THETA_STRUCT_" + struct_no_str, &d_incremented_angle_from_reference_axis[struct_no][0], 3);
    }
    d_timestep_counter = db->getDouble("TIMESTEP_COUNTER");
    d_FuRMoRP_current_time = db->getDouble("FuRMoRP_CURRENT_TIME");

    return;
} // getFromRestart

void
ConstraintIBMethod::setInitialLagrangianVelocity()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;

    const bool from_restart = RestartManager::getManager()->isFromRestart();
    if (!from_restart) calculateCOMandMOIOfStructures();

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        d_ib_kinematics[struct_no]->setKinematicsVelocity(d_FuRMoRP_current_time,
                                                          d_incremented_angle_from_reference_axis[struct_no],
                                                          d_center_of_mass_current[struct_no],
                                                          d_tagged_pt_position[struct_no]);
        d_ib_kinematics[struct_no]->setShape(d_FuRMoRP_current_time,
                                             d_incremented_angle_from_reference_axis[struct_no]);

        if (!from_restart)
        {
            const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();
            if (struct_param.getStructureIsSelfTranslating()) calculateMomentumOfKinematicsVelocity(struct_no);

            d_vel_com_def_current[struct_no] = d_vel_com_def_new[struct_no];
            d_omega_com_def_current[struct_no] = d_omega_com_def_new[struct_no];
        }
    }
    return;
} // setInitialLagrangianVelocity

void
ConstraintIBMethod::calculateCOMandMOIOfStructures()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();
    const double* const domain_x_upper = grid_geom->getXUpper();
    double domain_length[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        domain_length[d] = domain_x_upper[d] - domain_x_lower[d];
    }
    const IntVector<NDIM>& periodic_shift = grid_geom->getPeriodicShift();

    // Zero out the COM vector.
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        for (int d = 0; d < 3; ++d)
        {
            d_center_of_mass_unshifted_current[struct_no][d] = 0.0;
            d_center_of_mass_unshifted_new[struct_no][d] = 0.0;
            d_center_of_mass_current[struct_no][d] = 0.0;
            d_center_of_mass_new[struct_no][d] = 0.0;
        }
    }
    std::vector<std::vector<double> > tagged_position(d_no_structures, std::vector<double>(3, 0.0));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get LData corresponding to the present and new position of the
        // structures.
        Pointer<LData> ptr_x_lag_data_current(nullptr), ptr_x_lag_data_new(nullptr);
        ptr_x_lag_data_current = d_l_data_manager->getLData("X", ln);
        if (tbox::MathUtilities<double>::equalEps(d_FuRMoRP_current_time, 0.0))
        {
            ptr_x_lag_data_new = d_l_data_manager->getLData("X", ln);
        }
        else
        {
            ptr_x_lag_data_new = d_l_data_X_half_Euler[ln];
        }

        const boost::multi_array_ref<double, 2>& X_data_current = *ptr_x_lag_data_current->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data_new = *ptr_x_lag_data_new->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            double X_com_current[NDIM] = { 0.0 }, X_com_new[NDIM] = { 0.0 };
            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    const double* const X_current = &X_data_current[local_idx][0];
                    const double* const X_new = &X_data_new[local_idx][0];
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        X_com_current[d] += X_current[d] + displacement[d];
                        X_com_new[d] += X_new[d] + displacement[d];
                    }
                    if (lag_idx == d_tagged_pt_lag_idx[location_struct_handle])
                    {
                        for (unsigned int d = 0; d < NDIM; ++d) tagged_position[location_struct_handle][d] = X_new[d];
                    }
                }
            }
            for (int d = 0; d < NDIM; ++d)
            {
                d_center_of_mass_unshifted_current[location_struct_handle][d] += X_com_current[d];
                d_center_of_mass_unshifted_new[location_struct_handle][d] += X_com_new[d];
            }
        }
        ptr_x_lag_data_current->restoreArrays();
        ptr_x_lag_data_new->restoreArrays();
    }

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();
        const int total_nodes = struct_param.getTotalNodes();

        IBTK_MPI::sumReduction(d_center_of_mass_unshifted_current[struct_no].data(),
                               d_center_of_mass_unshifted_current[struct_no].size());
        IBTK_MPI::sumReduction(d_center_of_mass_unshifted_new[struct_no].data(),
                               d_center_of_mass_unshifted_new[struct_no].size());

        for (int i = 0; i < 3; ++i)
        {
            d_center_of_mass_unshifted_current[struct_no][i] /= total_nodes;
            d_center_of_mass_unshifted_new[struct_no][i] /= total_nodes;

            d_center_of_mass_current[struct_no][i] = d_center_of_mass_unshifted_current[struct_no][i];
            d_center_of_mass_new[struct_no][i] = d_center_of_mass_unshifted_new[struct_no][i];
        }
    }

    // now apply displacement
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        for (unsigned int i = 0; i < NDIM; ++i)
        {
            if (periodic_shift[i])
            {
                while (d_center_of_mass_current[struct_no][i] < domain_x_lower[i])
                {
                    d_center_of_mass_current[struct_no][i] += domain_length[i];
                }
                while (d_center_of_mass_current[struct_no][i] >= domain_x_upper[i])
                {
                    d_center_of_mass_current[struct_no][i] -= domain_length[i];
                }

                while (d_center_of_mass_new[struct_no][i] < domain_x_lower[i])
                {
                    d_center_of_mass_new[struct_no][i] += domain_length[i];
                }
                while (d_center_of_mass_new[struct_no][i] >= domain_x_upper[i])
                {
                    d_center_of_mass_new[struct_no][i] -= domain_length[i];
                }
            }
        }
    }

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        IBTK_MPI::sumReduction(&tagged_position[struct_no][0], 3);
        d_tagged_pt_position[struct_no] = tagged_position[struct_no];
    }

    // Zero out the moment of inertia tensor.
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        d_moment_of_inertia_current[struct_no].setZero();
        d_moment_of_inertia_new[struct_no].setZero();
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get LData corresponding to the present position of the structures.
        Pointer<LData> ptr_x_lag_data_current, ptr_x_lag_data_new;
        ptr_x_lag_data_current = d_l_data_manager->getLData("X", ln);
        if (tbox::MathUtilities<double>::equalEps(d_FuRMoRP_current_time, 0.0))
        {
            ptr_x_lag_data_new = d_l_data_manager->getLData("X", ln);
        }
        else
        {
            ptr_x_lag_data_new = d_l_data_X_half_Euler[ln];
        }

        const boost::multi_array_ref<double, 2>& X_data_current = *ptr_x_lag_data_current->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data_new = *ptr_x_lag_data_new->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            if (!struct_param.getStructureIsSelfRotating()) continue;

            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);
            const std::vector<double>& X_com_current = d_center_of_mass_unshifted_current[location_struct_handle];
            const std::vector<double>& X_com_new = d_center_of_mass_unshifted_new[location_struct_handle];

            Eigen::Matrix3d Inertia_current(3, 3), Inertia_new(3, 3);
            Inertia_current.setZero();
            Inertia_new.setZero();

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    const double* const X_current = &X_data_current[local_idx][0];
                    const double* const X_new = &X_data_new[local_idx][0];
#if (NDIM == 2)
                    Inertia_current(0, 0) += std::pow(displacement[1] + X_current[1] - X_com_current[1], 2);
                    Inertia_current(0, 1) += -(displacement[0] + X_current[0] - X_com_current[0]) *
                                             (displacement[1] + X_current[1] - X_com_current[1]);
                    Inertia_current(1, 1) += std::pow(displacement[0] + X_current[0] - X_com_current[0], 2);
                    Inertia_current(2, 2) += std::pow(displacement[0] + X_current[0] - X_com_current[0], 2) +
                                             std::pow(displacement[1] + X_current[1] - X_com_current[1], 2);

                    Inertia_new(0, 0) += std::pow(displacement[1] + X_new[1] - X_com_new[1], 2);
                    Inertia_new(0, 1) +=
                        -(displacement[0] + X_new[0] - X_com_new[0]) * (displacement[1] + X_new[1] - X_com_new[1]);
                    Inertia_new(1, 1) += std::pow(displacement[0] + X_new[0] - X_com_new[0], 2);
                    Inertia_new(2, 2) += std::pow(displacement[0] + X_new[0] - X_com_new[0], 2) +
                                         std::pow(displacement[1] + X_new[1] - X_com_new[1], 2);
#endif

#if (NDIM == 3)
                    Inertia_current(0, 0) += std::pow(displacement[1] + X_current[1] - X_com_current[1], 2) +
                                             std::pow(displacement[2] + X_current[2] - X_com_current[2], 2);
                    Inertia_current(0, 1) += -(displacement[0] + X_current[0] - X_com_current[0]) *
                                             (displacement[1] + X_current[1] - X_com_current[1]);
                    Inertia_current(0, 2) += -(displacement[0] + X_current[0] - X_com_current[0]) *
                                             (displacement[2] + X_current[2] - X_com_current[2]);
                    Inertia_current(1, 1) += std::pow(displacement[0] + X_current[0] - X_com_current[0], 2) +
                                             std::pow(displacement[2] + X_current[2] - X_com_current[2], 2);
                    Inertia_current(1, 2) += -(displacement[1] + X_current[1] - X_com_current[1]) *
                                             (displacement[2] + X_current[2] - X_com_current[2]);
                    Inertia_current(2, 2) += std::pow(displacement[0] + X_current[0] - X_com_current[0], 2) +
                                             std::pow(displacement[1] + X_current[1] - X_com_current[1], 2);

                    Inertia_new(0, 0) += std::pow(displacement[1] + X_new[1] - X_com_new[1], 2) +
                                         std::pow(displacement[2] + X_new[2] - X_com_new[2], 2);
                    Inertia_new(0, 1) +=
                        -(displacement[0] + X_new[0] - X_com_new[0]) * (displacement[1] + X_new[1] - X_com_new[1]);
                    Inertia_new(0, 2) +=
                        -(displacement[0] + X_new[0] - X_com_new[0]) * (displacement[2] + X_new[2] - X_com_new[2]);
                    Inertia_new(1, 1) += std::pow(displacement[0] + X_new[0] - X_com_new[0], 2) +
                                         std::pow(displacement[2] + X_new[2] - X_com_new[2], 2);
                    Inertia_new(1, 2) +=
                        -(displacement[1] + X_new[1] - X_com_new[1]) * (displacement[2] + X_new[2] - X_com_new[2]);
                    Inertia_new(2, 2) += std::pow(displacement[0] + X_new[0] - X_com_new[0], 2) +
                                         std::pow(displacement[1] + X_new[1] - X_com_new[1], 2);
#endif
                }
            }
            d_moment_of_inertia_current[location_struct_handle] += Inertia_current;
            d_moment_of_inertia_new[location_struct_handle] += Inertia_new;
        } // all structs
        ptr_x_lag_data_current->restoreArrays();
        ptr_x_lag_data_new->restoreArrays();
    } // all levels

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();
        if (struct_param.getStructureIsSelfRotating())
        {
            IBTK_MPI::sumReduction(&d_moment_of_inertia_current[struct_no](0, 0), 9);
            IBTK_MPI::sumReduction(&d_moment_of_inertia_new[struct_no](0, 0), 9);
        }
    }

    // Fill-in symmetric part of inertia tensor.
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        d_moment_of_inertia_current[struct_no](1, 0) = d_moment_of_inertia_current[struct_no](0, 1);
        d_moment_of_inertia_current[struct_no](2, 0) = d_moment_of_inertia_current[struct_no](0, 2);
        d_moment_of_inertia_current[struct_no](2, 1) = d_moment_of_inertia_current[struct_no](1, 2);

        d_moment_of_inertia_new[struct_no](1, 0) = d_moment_of_inertia_new[struct_no](0, 1);
        d_moment_of_inertia_new[struct_no](2, 0) = d_moment_of_inertia_new[struct_no](0, 2);
        d_moment_of_inertia_new[struct_no](2, 1) = d_moment_of_inertia_new[struct_no](1, 2);
    }

    // write the COM and MOI to the output file
    if (!IBTK_MPI::getRank() && d_print_output && d_output_COM_coordinates &&
        (d_timestep_counter % d_output_interval) == 0 && !MathUtilities<double>::equalEps(d_FuRMoRP_current_time, 0.0))
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_position_COM_stream[struct_no]
                << d_FuRMoRP_current_time << '\t' << d_center_of_mass_current[struct_no][0] << '\t'
                << d_center_of_mass_current[struct_no][1] << '\t' << d_center_of_mass_current[struct_no][2] << '\t'
                << d_center_of_mass_unshifted_current[struct_no][0] << '\t'
                << d_center_of_mass_unshifted_current[struct_no][1] << '\t'
                << d_center_of_mass_unshifted_current[struct_no][2] << std::endl;
        }
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_MOI && (d_timestep_counter % d_output_interval) == 0 &&
        !MathUtilities<double>::equalEps(d_FuRMoRP_current_time, 0.0))
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_moment_of_inertia_stream[struct_no]
                << d_FuRMoRP_current_time << '\t' << d_moment_of_inertia_current[struct_no](0, 0) << '\t'
                << d_moment_of_inertia_current[struct_no](0, 1) << '\t' << d_moment_of_inertia_current[struct_no](0, 2)
                << '\t' << d_moment_of_inertia_current[struct_no](1, 0) << '\t'
                << d_moment_of_inertia_current[struct_no](1, 1) << '\t' << d_moment_of_inertia_current[struct_no](1, 2)
                << '\t' << d_moment_of_inertia_current[struct_no](2, 0) << '\t'
                << d_moment_of_inertia_current[struct_no](2, 1) << '\t' << d_moment_of_inertia_current[struct_no](2, 2)
                << std::endl;
        }
    }

    return;
} // calculateCOMandMOIOfStructures

void
ConstraintIBMethod::calculateKinematicsVelocity()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;
    // Theta_new = Theta_old + Omega_old*dt
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();

        for (int d = 0; d < 3; ++d)
            d_incremented_angle_from_reference_axis[struct_no][d] +=
                (d_rigid_rot_vel_current[struct_no][d] - d_omega_com_def_current[struct_no][d]) * dt;

        d_ib_kinematics[struct_no]->setKinematicsVelocity(d_FuRMoRP_new_time,
                                                          d_incremented_angle_from_reference_axis[struct_no],
                                                          d_center_of_mass_new[struct_no],
                                                          d_tagged_pt_position[struct_no]);

        d_ib_kinematics[struct_no]->setShape(d_FuRMoRP_new_time, d_incremented_angle_from_reference_axis[struct_no]);

        if (struct_param.getStructureIsSelfTranslating()) calculateMomentumOfKinematicsVelocity(struct_no);
    }

    return;
} // calculateKinematicsVelocity

void
ConstraintIBMethod::calculateMomentumOfKinematicsVelocity(const int position_handle)
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    Pointer<ConstraintIBKinematics> ptr_ib_kinematics = d_ib_kinematics[position_handle];
    const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
    tbox::Array<int> calculate_trans_mom = struct_param.getCalculateTranslationalMomentum();
    tbox::Array<int> calculate_rot_mom = struct_param.getCalculateRotationalMomentum();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    const std::vector<std::pair<int, int> >& range = struct_param.getLagIdxRange();
    const int total_nodes = struct_param.getTotalNodes();

    // Zero out linear momentum of kinematics velocity of the structure.
    for (int d = 0; d < 3; ++d) d_vel_com_def_new[position_handle][d] = 0.0;

    // Calculate linear momentum
    for (int ln = coarsest_ln, itr = 0; ln <= finest_ln && static_cast<unsigned int>(itr) < range.size(); ++ln, ++itr)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_l_data_manager->levelContainsLagrangianData(ln));
#endif

        std::pair<int, int> lag_idx_range = range[itr];
        const int offset = lag_idx_range.first;
        double U_com_def[NDIM] = { 0.0 };

        // Get LMesh corresponding to the present position of the structures
        // on this level.
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        const std::vector<std::vector<double> >& def_vel = ptr_ib_kinematics->getKinematicsVelocity(ln);

        for (const auto& node_idx : local_nodes)
        {
            const int lag_idx = node_idx->getLagrangianIndex();
            if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_com_def[d] += def_vel[d][lag_idx - offset];
                }
            }
        }
        for (int d = 0; d < NDIM; ++d)
        {
            d_vel_com_def_new[position_handle][d] += U_com_def[d];
        }
    }
    IBTK_MPI::sumReduction(d_vel_com_def_new[position_handle].data(), d_vel_com_def_new[position_handle].size());

    for (int d = 0; d < 3; ++d)
    {
        if (calculate_trans_mom[d])
            d_vel_com_def_new[position_handle][d] /= total_nodes;
        else
            d_vel_com_def_new[position_handle][d] = 0.0;
    }

    // Calculate angular momentum.
    if (struct_param.getStructureIsSelfRotating())
    {
        // Zero out angular momentum of kinematics velocity of the structure.
        for (int d = 0; d < 3; ++d) d_omega_com_def_new[position_handle][d] = 0.0;

        for (int ln = coarsest_ln, itr = 0; ln <= finest_ln && static_cast<unsigned int>(itr) < range.size();
             ++ln, ++itr)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(d_l_data_manager->levelContainsLagrangianData(ln));
#endif

            std::pair<int, int> lag_idx_range = range[itr];
            const int offset = lag_idx_range.first;
            double R_cross_U_def[3] = { 0.0 };

            // Get LData corresponding to the present position of the structures.
            Pointer<LData> ptr_x_lag_data;
            if (MathUtilities<double>::equalEps(d_FuRMoRP_current_time, 0.0))
            {
                ptr_x_lag_data = d_l_data_manager->getLData("X", ln);
            }
            else
            {
                ptr_x_lag_data = d_l_data_X_half_Euler[ln];
            }

            const boost::multi_array_ref<double, 2>& X_data = *ptr_x_lag_data->getLocalFormVecArray();
            const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
            const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
            const std::vector<std::vector<double> >& def_vel = ptr_ib_kinematics->getKinematicsVelocity(ln);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const X = &X_data[local_idx][0];
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
#if (NDIM == 2)
                    double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[position_handle][0];
                    double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[position_handle][1];
                    R_cross_U_def[2] += (x * (def_vel[1][lag_idx - offset]) - y * (def_vel[0][lag_idx - offset]));

#endif

#if (NDIM == 3)
                    double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[position_handle][0];
                    double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[position_handle][1];
                    double z = displacement[2] + X[2] - d_center_of_mass_unshifted_new[position_handle][2];

                    R_cross_U_def[0] += (y * (def_vel[2][lag_idx - offset]) - z * (def_vel[1][lag_idx - offset]));

                    R_cross_U_def[1] += (-x * (def_vel[2][lag_idx - offset]) + z * (def_vel[0][lag_idx - offset]));

                    R_cross_U_def[2] += (x * (def_vel[1][lag_idx - offset]) - y * (def_vel[0][lag_idx - offset]));
#endif
                }
            }
            for (int d = 0; d < 3; ++d)
            {
                d_omega_com_def_new[position_handle][d] += R_cross_U_def[d];
            }
            ptr_x_lag_data->restoreArrays();
        } // all levels
        IBTK_MPI::sumReduction(&d_omega_com_def_new[position_handle][0], 3);

// Find angular velocity of deformational velocity.
#if (NDIM == 2)
        d_omega_com_def_new[position_handle][2] /= d_moment_of_inertia_new[position_handle](2, 2);
#endif

#if (NDIM == 3)
        solveSystemOfEqns(d_omega_com_def_new[position_handle], d_moment_of_inertia_new[position_handle]);
        for (int d = 0; d < 3; ++d)
            if (!calculate_rot_mom[d]) d_omega_com_def_new[position_handle][d] = 0.0;
#endif
    } // if struct is rotating

    return;
} // calculateMomentumOfKinematicsVelocity

void
ConstraintIBMethod::calculateVolumeElement()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;

    // Initialize variables and variable contexts associated with Eulerian
    // tracking of the Lagrangian points.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> cell_ghosts = 0;
    Pointer<CellVariable<NDIM, int> > vol_cc_var = new CellVariable<NDIM, int>(d_object_name + "::vol_cc_var");
    const int vol_cc_scratch_idx = var_db->registerVariableAndContext(vol_cc_var, d_scratch_context, cell_ghosts);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(vol_cc_scratch_idx, 0.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CellData<NDIM, int> > vol_cc_scratch_idx_data = patch->getPatchData(vol_cc_scratch_idx);
            vol_cc_scratch_idx_data->fill(0, 0);
        }
    }

    const int lag_node_index_idx = d_l_data_manager->getLNodePatchDescriptorIndex();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get LData corresponding to the present position of the structures.
        const boost::multi_array_ref<double, 2>& X_data = *d_l_data_manager->getLData("X", ln)->getLocalFormVecArray();
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();
        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            // If the volume element has already been set, then skip the volume
            // computation
            if (d_vol_element_is_set[struct_no])
            {
                tbox::plog << "Skipping volume element computation for structure no. " << struct_no << std::endl;
                tbox::pout << "Skipping volume element computation for structure no. " << struct_no << std::endl;
                continue;
            }
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<CellData<NDIM, int> > vol_cc_scratch_idx_data = patch->getPatchData(vol_cc_scratch_idx);
                Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                Pointer<CartesianGridGeometry<NDIM> > ggeom = level->getGridGeometry();
                const IntVector<NDIM>& ratio = level->getRatio();
                const double* const dx = pgeom->getDx();
                const Box<NDIM>& patch_box = patch->getBox();
#if (NDIM == 2)
                const double patch_cell_vol = dx[0] * dx[1];
#endif

#if (NDIM == 3)
                const double patch_cell_vol = dx[0] * dx[1] * dx[2];
#endif
                const Pointer<LNodeSetData> lag_node_index_data = patch->getPatchData(lag_node_index_idx);
                for (LNodeSetData::DataIterator it = lag_node_index_data->data_begin(patch_box);
                     it != lag_node_index_data->data_end();
                     ++it)
                {
                    LNode* const node_idx = *it;
                    const int lag_idx = node_idx->getLagrangianIndex();
                    if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                    {
                        const int local_idx = node_idx->getLocalPETScIndex();
                        const double* const X = &X_data[local_idx][0];
                        const CellIndex<NDIM> Lag2Eul_cellindex = IndexUtilities::getCellIndex(X, ggeom, ratio);

                        (*vol_cc_scratch_idx_data)(Lag2Eul_cellindex)++;
                    }
                } // on a patch

                for (CellData<NDIM, int>::Iterator it(patch_box); it; it++)
                {
                    if ((*vol_cc_scratch_idx_data)(*it)) d_structure_vol[location_struct_handle] += patch_cell_vol;

                } // on the same patch
                vol_cc_scratch_idx_data->fill(0, patch_box, 0);

            } // all patches
        }     // all structs
        d_l_data_manager->getLData("X", ln)->restoreArrays();
    } // all levels
    IBTK_MPI::sumReduction(&d_structure_vol[0], d_no_structures);

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        Pointer<ConstraintIBKinematics> ptr_ib_kinematics = d_ib_kinematics[struct_no];
        const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();

        // If the volume element has already been set, then no need to compute it
        if (d_vol_element_is_set[struct_no])
        {
            d_structure_vol[struct_no] = struct_param.getTotalNodes() * d_vol_element[struct_no];
        }
        else
        {
            d_vol_element[struct_no] = d_structure_vol[struct_no] / struct_param.getTotalNodes();
            d_vol_element_is_set[struct_no] = true;
        }

        tbox::plog << " ++++++++++++++++ "
                   << " STRUCTURE NO. " << struct_no << "  ++++++++++++++++++++++++++ \n\n\n"
                   << " VOLUME OF THE MATERIAL ELEMENT           = " << d_vol_element[struct_no] << "\n"
                   << " VOLUME OF THE STRUCTURE                  = " << d_structure_vol[struct_no] << "\n"
                   << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                      "+++++++\n"
                   << std::endl;

        tbox::pout << " ++++++++++++++++ "
                   << " STRUCTURE NO. " << struct_no << "  ++++++++++++++++++++++++++ \n\n\n"
                   << " VOLUME OF THE MATERIAL ELEMENT           = " << d_vol_element[struct_no] << "\n"
                   << " VOLUME OF THE STRUCTURE                  = " << d_structure_vol[struct_no] << "\n"
                   << " +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                      "+++++++\n"
                   << std::endl;
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(vol_cc_scratch_idx)) level->deallocatePatchData(vol_cc_scratch_idx);
    }
    var_db->removePatchDataIndex(vol_cc_scratch_idx);

    return;

} // calculateVolumeElement

void
ConstraintIBMethod::calculateRigidTranslationalMomentum()
{
    // Zero out new rigid momentum.
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        for (int d = 0; d < 3; ++d) d_rigid_trans_vel_new[struct_no][d] = 0.0;
    }

    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Calculate rigid translational velocity.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get LData corresponding to the present position of the structures.
        const boost::multi_array_ref<double, 2>& U_interp_data = *d_l_data_U_interp[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            if (!struct_param.getStructureIsSelfTranslating()) continue;

            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            double U_rigid[NDIM] = { 0.0 };
            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U = &U_interp_data[local_idx][0];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        U_rigid[d] += U[d];
                    }
                }
            }
            for (int d = 0; d < NDIM; ++d) d_rigid_trans_vel_new[location_struct_handle][d] += U_rigid[d];
        } // all structs
        d_l_data_U_interp[ln]->restoreArrays();
    } // all levels

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();
        if (struct_param.getStructureIsSelfTranslating())
        {
            IBTK_MPI::sumReduction(d_rigid_trans_vel_new[struct_no].data(), d_rigid_trans_vel_new[struct_no].size());
            tbox::Array<int> calculate_trans_mom = struct_param.getCalculateTranslationalMomentum();
            for (int d = 0; d < NDIM; ++d)
            {
                if (calculate_trans_mom[d])
                    d_rigid_trans_vel_new[struct_no][d] /= struct_param.getTotalNodes();
                else
                    d_rigid_trans_vel_new[struct_no][d] = 0.0;
            }
        }
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_trans_vel && (d_timestep_counter % d_output_interval) == 0)
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_trans_vel_stream[struct_no]
                << d_FuRMoRP_new_time << '\t' << d_rigid_trans_vel_new[struct_no][0] << '\t'
                << d_rigid_trans_vel_new[struct_no][1] << '\t' << d_rigid_trans_vel_new[struct_no][2] << '\t'
                << d_vel_com_def_new[struct_no][0] << '\t' << d_vel_com_def_new[struct_no][1] << '\t'
                << d_vel_com_def_new[struct_no][2] << std::endl;
        }
    }

    return;

} // calculateRigidTranslationalMomentum

void
ConstraintIBMethod::calculateRigidRotationalMomentum()
{
    // Zero out new rigid momentum.
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        for (int d = 0; d < 3; ++d) d_rigid_rot_vel_new[struct_no][d] = 0.0;
    }

    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Calculate rigid rotational velocity.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get ponter to LData.
        const boost::multi_array_ref<double, 2>& U_interp_data = *d_l_data_U_interp[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data = *d_l_data_X_half_Euler[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            if (!struct_param.getStructureIsSelfRotating()) continue;

            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            double Omega_rigid[3] = { 0.0 };
            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    const double* const U = &U_interp_data[local_idx][0];
                    const double* const X = &X_data[local_idx][0];
#if (NDIM == 2)
                    const double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[location_struct_handle][0];
                    const double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[location_struct_handle][1];
                    Omega_rigid[2] += x * U[1] - y * U[0];
#endif

#if (NDIM == 3)
                    const double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[location_struct_handle][0];
                    const double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[location_struct_handle][1];
                    const double z = displacement[2] + X[2] - d_center_of_mass_unshifted_new[location_struct_handle][2];
                    Omega_rigid[0] += y * U[2] - z * U[1];
                    Omega_rigid[1] += -x * U[2] + z * U[0];
                    Omega_rigid[2] += x * U[1] - y * U[0];
#endif
                }
            }
            for (int d = 0; d < 3; ++d) d_rigid_rot_vel_new[location_struct_handle][d] += Omega_rigid[d];
        } // all structs
        d_l_data_U_interp[ln]->restoreArrays();
        d_l_data_X_half_Euler[ln]->restoreArrays();
    } // all levels

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        const StructureParameters& struct_param = d_ib_kinematics[struct_no]->getStructureParameters();
        if (struct_param.getStructureIsSelfRotating())
        {
            IBTK_MPI::sumReduction(&d_rigid_rot_vel_new[struct_no][0], 3);
#if (NDIM == 2)
            d_rigid_rot_vel_new[struct_no][2] /= d_moment_of_inertia_new[struct_no](2, 2);
#endif

#if (NDIM == 3)
            solveSystemOfEqns(d_rigid_rot_vel_new[struct_no], d_moment_of_inertia_new[struct_no]);
            tbox::Array<int> calculate_rot_mom = struct_param.getCalculateRotationalMomentum();
            for (int d = 0; d < NDIM; ++d)
            {
                if (!calculate_rot_mom[d]) d_rigid_rot_vel_new[struct_no][d] = 0.0;
            }
#endif
        }
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_rot_vel && (d_timestep_counter % d_output_interval) == 0)
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_rot_vel_stream[struct_no] << d_FuRMoRP_new_time << '\t' << d_rigid_rot_vel_new[struct_no][0] << '\t'
                                         << d_rigid_rot_vel_new[struct_no][1] << '\t'
                                         << d_rigid_rot_vel_new[struct_no][2] << '\t'
                                         << d_omega_com_def_new[struct_no][0] << '\t'
                                         << d_omega_com_def_new[struct_no][1] << '\t'
                                         << d_omega_com_def_new[struct_no][2] << std::endl;
        }
    }

    return;

} // calculateRigidRotationalMomentum

void
ConstraintIBMethod::calculateCurrentLagrangianVelocity()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<double> WxR(3, 0.0), R(3, 0.0);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get pointer to LData.
        boost::multi_array_ref<double, 2>& U_current_data = *d_l_data_U_current[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data = *d_l_data_manager->getLData("X", ln)->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            const int offset = lag_idx_range.first;
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            const std::vector<std::vector<double> >& current_vel = ptr_ib_kinematics->getKinematicsVelocity(ln);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();

                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    double* const U_current = &U_current_data[local_idx][0];
                    const double* const X = &X_data[local_idx][0];

                    // Imposed velocity
                    for (int d = 0; d < NDIM; ++d)
                    {
                        U_current[d] = current_vel[d][lag_idx - offset];
                    }

                    // Translational velocity
                    if (struct_param.getStructureIsSelfTranslating())
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            U_current[d] += d_rigid_trans_vel_current[location_struct_handle][d] -
                                            d_vel_com_def_current[location_struct_handle][d];
                        }
                    }

                    // Rotational velocity
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    if (struct_param.getStructureIsSelfRotating())
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            R[d] =
                                displacement[d] + X[d] - d_center_of_mass_unshifted_current[location_struct_handle][d];
                        }

                        WxR[0] = R[2] * (d_rigid_rot_vel_current[location_struct_handle][1] -
                                         d_omega_com_def_current[location_struct_handle][1]) -
                                 R[1] * (d_rigid_rot_vel_current[location_struct_handle][2] -
                                         d_omega_com_def_current[location_struct_handle][2]);

                        WxR[1] = -R[2] * (d_rigid_rot_vel_current[location_struct_handle][0] -
                                          d_omega_com_def_current[location_struct_handle][0]) +
                                 R[0] * (d_rigid_rot_vel_current[location_struct_handle][2] -
                                         d_omega_com_def_current[location_struct_handle][2]);

                        WxR[2] = R[1] * (d_rigid_rot_vel_current[location_struct_handle][0] -
                                         d_omega_com_def_current[location_struct_handle][0]) -
                                 R[0] * (d_rigid_rot_vel_current[location_struct_handle][1] -
                                         d_omega_com_def_current[location_struct_handle][1]);
                        for (int d = 0; d < NDIM; ++d)
                        {
                            U_current[d] += WxR[d];
                        }
                    }

                } // choose a struct
            }     // all nodes on a level
        }         // all structs
        d_l_data_U_current[ln]->restoreArrays();
        d_l_data_manager->getLData("X", ln)->restoreArrays();
    } // all levels

    return;

} // calculateCurrentLagrangianVelocity

void
ConstraintIBMethod::correctVelocityOnLagrangianMesh()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    std::vector<double> WxR(3, 0.0), R(3, 0.0);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get pointer to LData.
        const boost::multi_array_ref<double, 2>& U_interp_data = *d_l_data_U_interp[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& U_corr_data = *d_l_data_U_correction[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& U_new_data = *d_l_data_U_new[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data = *d_l_data_X_half_Euler[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            const int offset = lag_idx_range.first;
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            const std::vector<std::vector<double> >& new_vel = ptr_ib_kinematics->getKinematicsVelocity(ln);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U = &U_interp_data[local_idx][0];
                    double* const U_corr = &U_corr_data[local_idx][0];
                    double* const U_new = &U_new_data[local_idx][0];
                    const double* const X = &X_data[local_idx][0];

                    // Imposed velocity
                    for (int d = 0; d < NDIM; ++d)
                    {
                        U_new[d] = new_vel[d][lag_idx - offset];
                    }

                    // Translational velocity
                    if (struct_param.getStructureIsSelfTranslating())
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            U_new[d] += d_rigid_trans_vel_new[location_struct_handle][d] -
                                        d_vel_com_def_new[location_struct_handle][d];
                        }
                    }

                    // Rotational velocity
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    if (struct_param.getStructureIsSelfRotating())
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            R[d] = displacement[d] + X[d] - d_center_of_mass_unshifted_new[location_struct_handle][d];
                        }
                        WxR[0] = R[2] * (d_rigid_rot_vel_new[location_struct_handle][1] -
                                         d_omega_com_def_new[location_struct_handle][1]) -
                                 R[1] * (d_rigid_rot_vel_new[location_struct_handle][2] -
                                         d_omega_com_def_new[location_struct_handle][2]);

                        WxR[1] = -R[2] * (d_rigid_rot_vel_new[location_struct_handle][0] -
                                          d_omega_com_def_new[location_struct_handle][0]) +
                                 R[0] * (d_rigid_rot_vel_new[location_struct_handle][2] -
                                         d_omega_com_def_new[location_struct_handle][2]);

                        WxR[2] = R[1] * (d_rigid_rot_vel_new[location_struct_handle][0] -
                                         d_omega_com_def_new[location_struct_handle][0]) -
                                 R[0] * (d_rigid_rot_vel_new[location_struct_handle][1] -
                                         d_omega_com_def_new[location_struct_handle][1]);

                        for (int d = 0; d < NDIM; ++d)
                        {
                            U_new[d] += WxR[d];
                        }
                    }

                    for (int d = 0; d < NDIM; ++d)
                    {
                        U_corr[d] = (U_new[d] - U[d]) * d_vol_element[location_struct_handle];
                    }

                } // choose a struct
            }     // all nodes on a level
        }         // all structs
        d_l_data_U_interp[ln]->restoreArrays();
        d_l_data_U_correction[ln]->restoreArrays();
        d_l_data_U_new[ln]->restoreArrays();
        d_l_data_X_half_Euler[ln]->restoreArrays();
    } // all levels

    return;

} // correctVelocityOnLagrangianMesh

void
ConstraintIBMethod::applyProjection()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate temporary data.
    ComponentSelector scratch_idxs;
    scratch_idxs.setFlag(d_phi_idx);
    scratch_idxs.setFlag(d_Div_u_scratch_idx);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(scratch_idxs, d_FuRMoRP_new_time);
    }

    // Compute div U before applying the projection operator.
    const bool U_current_cf_bdry_synch = true;
    getHierarchyMathOps()->div(d_Div_u_scratch_idx,
                               d_Div_u_var, // dst
                               +1.0,        // alpha
                               d_u_fluidSolve_idx,
                               Pointer<SideVariable<NDIM, double> >(d_u_fluidSolve_var), // src
                               d_no_fill_op,                                             // src_bdry_fill
                               d_FuRMoRP_new_time,                                       // src_bdry_fill_time
                               U_current_cf_bdry_synch);                                 // src_cf_bdry_synch

    if (d_do_log)
    {
        const double Div_u_norm_1 = d_hier_cc_data_ops->L1Norm(d_Div_u_scratch_idx, d_wgt_cc_idx);
        const double Div_u_norm_2 = d_hier_cc_data_ops->L2Norm(d_Div_u_scratch_idx, d_wgt_cc_idx);
        const double Div_u_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_u_scratch_idx, d_wgt_cc_idx);
        tbox::plog << d_object_name << "::applyProjection():\n"
                   << "  performing velocity correction projection\n"
                   << "  before projection:\n"
                   << "    ||Div U||_1  = " << Div_u_norm_1 << "\n"
                   << "    ||Div U||_2  = " << Div_u_norm_2 << "\n"
                   << "    ||Div U||_oo = " << Div_u_norm_oo << "\n";
    }

    // Setup the solver vectors.
    d_hier_cc_data_ops->setToScalar(d_phi_idx, 0.0, false);
    d_hier_cc_data_ops->scale(d_Div_u_scratch_idx, -1.0, d_Div_u_scratch_idx);
    const double Div_u_mean = (1.0 / d_volume) * d_hier_cc_data_ops->integral(d_Div_u_scratch_idx, d_wgt_cc_idx);
    d_hier_cc_data_ops->addScalar(d_Div_u_scratch_idx, d_Div_u_scratch_idx, -Div_u_mean);

    SAMRAIVectorReal<NDIM, double> sol_vec(d_object_name + "::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(d_phi_var, d_phi_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAIVectorReal<NDIM, double> rhs_vec(d_object_name + "::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_Div_u_var, d_Div_u_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the Poisson solver.
    // Note that here, the delta t is absorbed in the phi term
    d_velcorrection_projection_spec->setCZero();
    if (!d_rho_is_const)
    {
        // Copy rho from INS integrator and take reciprocal and scale
        copyDensityVariable(d_rho_ins_idx, d_rho_scratch_idx);
        d_hier_sc_data_ops->reciprocal(d_rho_scratch_idx, d_rho_scratch_idx);
        d_hier_sc_data_ops->scale(d_rho_scratch_idx, -1.0, d_rho_scratch_idx);

        // Synchronize the coefficient patch data
        using SynchronizationTransactionComponent = SideDataSynchronization::SynchronizationTransactionComponent;
        SynchronizationTransactionComponent coef_synch_transaction =
            SynchronizationTransactionComponent(d_rho_scratch_idx, "CONSERVATIVE_COARSEN");
        Pointer<SideDataSynchronization> side_synch_op = new SideDataSynchronization();
        side_synch_op->initializeOperatorState(coef_synch_transaction, d_hierarchy);
        side_synch_op->synchronizeData(d_FuRMoRP_new_time);
        d_velcorrection_projection_spec->setDPatchDataId(d_rho_scratch_idx);
    }
    else
    {
        d_velcorrection_projection_spec->setDConstant(-1.0 / d_rho_fluid);
    }

    d_velcorrection_projection_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
    d_velcorrection_projection_op->setPhysicalBcCoef(&d_velcorrection_projection_bc_coef);
    d_velcorrection_projection_op->setHomogeneousBc(true);
    d_velcorrection_projection_op->setHierarchyMathOps(getHierarchyMathOps());

    d_velcorrection_projection_fac_op->setPoissonSpecifications(*d_velcorrection_projection_spec);
    d_velcorrection_projection_fac_op->setPhysicalBcCoef(&d_velcorrection_projection_bc_coef);

    d_velcorrection_projection_solver->setInitialGuessNonzero(false);
    d_velcorrection_projection_solver->setOperator(d_velcorrection_projection_op);

    // NOTE: We always use homogeneous Neumann boundary conditions for the
    // velocity correction projection Poisson solver.
    d_velcorrection_projection_solver->setNullspace(true);

    // Solve the projection Poisson problem.
    d_velcorrection_projection_solver->initializeSolverState(sol_vec, rhs_vec);
    d_velcorrection_projection_solver->solveSystem(sol_vec, rhs_vec);
    d_velcorrection_projection_solver->deallocateSolverState();

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    InterpolationTransactionComponent Phi_bc_component(
        d_phi_idx, "LINEAR_REFINE", true, "CUBIC_COARSEN", "LINEAR", false, &d_velcorrection_projection_bc_coef);
    Pointer<HierarchyGhostCellInterpolation> Phi_bdry_bc_fill_op = new HierarchyGhostCellInterpolation();
    Phi_bdry_bc_fill_op->initializeOperatorState(Phi_bc_component, d_hierarchy);

    // Fill the physical boundary conditions for Phi.
    Phi_bdry_bc_fill_op->setHomogeneousBc(true);
    Phi_bdry_bc_fill_op->fillData(d_FuRMoRP_new_time);

    // Set U := U - 1/rho * grad Phi.
    const bool U_scratch_cf_bdry_synch = true;
    if (!d_rho_is_const)
    {
        getHierarchyMathOps()->grad(d_u_fluidSolve_idx,
                                    Pointer<SideVariable<NDIM, double> >(d_u_var),
                                    U_scratch_cf_bdry_synch,
                                    d_rho_scratch_idx,
                                    Pointer<SideVariable<NDIM, double> >(d_rho_var),
                                    d_phi_idx,
                                    d_phi_var,
                                    d_no_fill_op,
                                    d_FuRMoRP_new_time,
                                    1.0,
                                    d_u_fluidSolve_idx,
                                    Pointer<SideVariable<NDIM, double> >(d_u_var));
    }
    else
    {
        getHierarchyMathOps()->grad(d_u_fluidSolve_idx,
                                    Pointer<SideVariable<NDIM, double> >(d_u_var),
                                    U_scratch_cf_bdry_synch,
                                    -1.0 / d_rho_fluid,
                                    d_phi_idx,
                                    d_phi_var,
                                    d_no_fill_op,
                                    d_FuRMoRP_new_time,
                                    1.0,
                                    d_u_fluidSolve_idx,
                                    Pointer<SideVariable<NDIM, double> >(d_u_var));
    }

    // Update pressure p = p + phi/dt
    const Pointer<Variable<NDIM> > p_var = d_ib_solver->getPressureVariable();
    const Pointer<VariableContext> p_ctx = d_ib_solver->getNewContext();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;
    d_hier_cc_data_ops->axpy(p_idx, 1.0 / dt, d_phi_idx, p_idx);

    // Compute div U after applying the projection operator
    if (d_do_log)
    {
        // Compute div U before applying the projection operator.
        const bool U_current_cf_bdry_synch = true;
        getHierarchyMathOps()->div(d_Div_u_scratch_idx,
                                   d_Div_u_var, // dst
                                   +1.0,        // alpha
                                   d_u_fluidSolve_idx,
                                   Pointer<SideVariable<NDIM, double> >(d_u_fluidSolve_var), // src
                                   d_no_fill_op,                                             // src_bdry_fill
                                   d_FuRMoRP_new_time,                                       // src_bdry_fill_time
                                   U_current_cf_bdry_synch);                                 // src_cf_bdry_synch

        const double Div_u_norm_1 = d_hier_cc_data_ops->L1Norm(d_Div_u_scratch_idx, d_wgt_cc_idx);
        const double Div_u_norm_2 = d_hier_cc_data_ops->L2Norm(d_Div_u_scratch_idx, d_wgt_cc_idx);
        const double Div_u_norm_oo = d_hier_cc_data_ops->maxNorm(d_Div_u_scratch_idx, d_wgt_cc_idx);
        tbox::plog << "  after projection:\n"
                   << "    ||Div U||_1  = " << Div_u_norm_1 << "\n"
                   << "    ||Div U||_2  = " << Div_u_norm_2 << "\n"
                   << "    ||Div U||_oo = " << Div_u_norm_oo << "\n";
    }

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(scratch_idxs);
    }

    return;

} // applyProjection

void
ConstraintIBMethod::updateStructurePositionEulerStep()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        boost::multi_array_ref<double, 2>& X_half_Euler_data = *d_l_data_X_half_Euler[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_current_data =
            *d_l_data_manager->getLData("X", ln)->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_current_data = *d_l_data_U_current[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            const int offset = lag_idx_range.first;
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            const std::string position_update_method = struct_param.getPositionUpdateMethod();
            const std::vector<std::vector<double> >& current_shape = ptr_ib_kinematics->getShape(ln);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U_current = &U_current_data[local_idx][0];
                    const double* const X_current = &X_current_data[local_idx][0];
                    double* const X_half = &X_half_Euler_data[local_idx][0];
                    if (position_update_method == "CONSTRAINT_VELOCITY")
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_half[d] = X_current[d] + 0.5 * dt * U_current[d];
                        }
                    }
                    else if (position_update_method == "CONSTRAINT_POSITION")
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_half[d] = d_center_of_mass_current[location_struct_handle][d] +
                                        current_shape[d][lag_idx - offset] +
                                        0.5 * dt * (d_rigid_trans_vel_current[location_struct_handle][d]);
                        }
                    }
                    else if (position_update_method == "CONSTRAINT_EXPT_POSITION")
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_half[d] = current_shape[d][lag_idx - offset];
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "ConstraintIBMethod::updateStructurePositionEulerStep():"
                            ": Unknown position update method encountered\n"
                            << "Supported methods are : CONSTRAINT_VELOCITY, "
                               "CONSTRAINT_POSITION AND "
                               "CONSTRAINT_EXPT_POSITION "
                            << std::endl);
                    }
                }
            }
        } // all structs
        d_l_data_X_half_Euler[ln]->restoreArrays();
        d_l_data_U_current[ln]->restoreArrays();
        d_l_data_manager->getLData("X", ln)->restoreArrays();
    }
    return;
} // updateStructurePositionEulerStep

void
ConstraintIBMethod::forwardEulerStep(double current_time, double new_time)
{
    IBMethod::forwardEulerStep(current_time, new_time);

    setFuRMoRPTime(current_time, new_time);

    IBTK_TIMER_START(t_eulerStep);
    updateStructurePositionEulerStep();
    IBTK_TIMER_STOP(t_eulerStep);

    return;
} // eulerStep

void
ConstraintIBMethod::updateStructurePositionMidPointStep()
{
    using StructureParameters = ConstraintIBKinematics::StructureParameters;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;

    calculateMidPointVelocity();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        boost::multi_array_ref<double, 2>& X_new_MidPoint_data = *d_l_data_X_new_MidPoint[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_current_data =
            *d_l_data_manager->getLData("X", ln)->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_half_data = *d_l_data_U_half[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            const int offset = lag_idx_range.first;
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);
            const StructureParameters& struct_param = ptr_ib_kinematics->getStructureParameters();
            const std::string position_update_method = struct_param.getPositionUpdateMethod();
            const std::vector<std::vector<double> >& new_shape = ptr_ib_kinematics->getShape(ln);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U_half = &U_half_data[local_idx][0];
                    const double* const X_current = &X_current_data[local_idx][0];
                    double* const X_new = &X_new_MidPoint_data[local_idx][0];

                    if (position_update_method == "CONSTRAINT_VELOCITY")
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_new[d] = X_current[d] + dt * U_half[d];
                        }
                    }
                    else if (position_update_method == "CONSTRAINT_POSITION")
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_new[d] = d_center_of_mass_current[location_struct_handle][d] +
                                       new_shape[d][lag_idx - offset] +
                                       dt * 0.5 *
                                           (d_rigid_trans_vel_current[location_struct_handle][d] +
                                            d_rigid_trans_vel_new[location_struct_handle][d]);
                        }
                    }
                    else if (position_update_method == "CONSTRAINT_EXPT_POSITION")
                    {
                        for (int d = 0; d < NDIM; ++d)
                        {
                            X_new[d] = new_shape[d][lag_idx - offset];
                        }
                    }
                    else
                    {
                        TBOX_ERROR(
                            "ConstraintIBMethod::"
                            "updateStructurePositionMidPointStep():: Unknown "
                            "position update method encountered\n"
                            << "Supported methods are : CONSTRAINT_VELOCITY, "
                               "CONSTRAINT_POSITION AND "
                               "CONSTRAINT_EXPT_POSITION "
                            << std::endl);
                    }
                }
            }
        } // all structs
        d_l_data_X_new_MidPoint[ln]->restoreArrays();
        d_l_data_manager->getLData("X", ln)->restoreArrays();
        d_l_data_U_half[ln]->restoreArrays();
    }
    return;

} // updateStructurePositionMidPointStep

void
ConstraintIBMethod::midpointStep(double current_time, double new_time)
{
    IBMethod::midpointStep(current_time, new_time);

    IBTK_TIMER_START(t_midpointStep);

    updateStructurePositionMidPointStep();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    int ierr;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecCopy(d_l_data_X_new_MidPoint[ln]->getVec(), d_X_new_data[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }

    IBTK_TIMER_STOP(t_midpointStep);
    return;

} // midpointStep

void
ConstraintIBMethod::copyFluidVariable(int copy_from_idx, int copy_to_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(copy_to_idx)) level->allocatePatchData(copy_to_idx);
    }

    SAMRAIVectorReal<NDIM, double> u_from(d_object_name + "from", d_hierarchy, coarsest_ln, finest_ln);
    SAMRAIVectorReal<NDIM, double> u_to(d_object_name + "to", d_hierarchy, coarsest_ln, finest_ln);

    u_from.addComponent(d_u_fluidSolve_var, copy_from_idx, d_wgt_sc_idx);
    u_to.addComponent(d_u_fluidSolve_var, copy_to_idx, d_wgt_sc_idx);

    u_to.copyVector(Pointer<SAMRAIVectorReal<NDIM, double> >(&u_from, false));

    using InterpolationTransactionComponent = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    InterpolationTransactionComponent component(copy_to_idx,
                                                DATA_REFINE_TYPE,
                                                USE_CF_INTERPOLATION,
                                                SIDE_DATA_COARSEN_TYPE,
                                                BDRY_EXTRAP_TYPE,
                                                CONSISTENT_TYPE_2_BDRY,
                                                std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr),
                                                nullptr);
    transaction_comps.push_back(component);

    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(transaction_comps, d_hierarchy, coarsest_ln, finest_ln);
    const bool homogeneous_bc = true;
    hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    hier_bdry_fill->fillData(0.0);

    return;
} // copyFluidVariable

void
ConstraintIBMethod::copyDensityVariable(int copy_from_idx, int copy_to_idx)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(copy_to_idx)) level->allocatePatchData(copy_to_idx);
    }

    d_hier_sc_data_ops->copyData(copy_to_idx,
                                 copy_from_idx,
                                 /*interior_only*/ true);

    return;
} // copyDensityVariable

void
ConstraintIBMethod::interpolateFluidSolveVelocity()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > F_data(finest_ln + 1, SAMRAI::tbox::Pointer<IBTK::LData>(nullptr));
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > X_data(finest_ln + 1, SAMRAI::tbox::Pointer<IBTK::LData>(nullptr));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        F_data[ln] = d_l_data_U_interp[ln];
        X_data[ln] = d_l_data_X_half_Euler[ln];
    }

    d_l_data_manager->interp(d_u_fluidSolve_cib_idx, F_data, X_data);

    return;
} // interpolateFluidSolveVelocity

void
ConstraintIBMethod::spreadCorrectedLagrangianVelocity()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > F_data(finest_ln + 1, SAMRAI::tbox::Pointer<IBTK::LData>(nullptr));
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > X_data(finest_ln + 1, SAMRAI::tbox::Pointer<IBTK::LData>(nullptr));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        F_data[ln] = d_l_data_U_correction[ln];
        X_data[ln] = d_l_data_X_half_Euler[ln];
    }

    // Since we do not want to mess up the boundary values of u_ins, we
    // zero-out the scratch variable, spread to it and then add the correction
    // to u_ins. This assumes that the structure is away from the physical
    // domain.
    SAMRAIVectorReal<NDIM, double> u_cib(d_object_name + "cib", d_hierarchy, coarsest_ln, finest_ln);
    SAMRAIVectorReal<NDIM, double> u_ins(d_object_name + "ins", d_hierarchy, coarsest_ln, finest_ln);

    u_cib.addComponent(d_u_fluidSolve_var, d_u_fluidSolve_cib_idx, d_wgt_sc_idx);
    u_ins.addComponent(d_u_fluidSolve_var, d_u_fluidSolve_idx, d_wgt_sc_idx);

    u_cib.setToScalar(0.0);
    d_l_data_manager->spread(d_u_fluidSolve_cib_idx, F_data, X_data, d_u_phys_bdry_op);

    u_ins.add(Pointer<SAMRAIVectorReal<NDIM, double> >(&u_ins, false),
              Pointer<SAMRAIVectorReal<NDIM, double> >(&u_cib, false));

    return;
} // spreadCorrectedLagrangianVelocity

void
ConstraintIBMethod::calculateMidPointVelocity()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    int ierr;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        ierr = VecAXPBYPCZ(d_l_data_U_half[ln]->getVec(),
                           0.5,
                           0.5,
                           0.0,
                           d_l_data_U_current[ln]->getVec(),
                           d_l_data_U_new[ln]->getVec());
        IBTK_CHKERRQ(ierr);
    }
    return;

} // calculateMidPointVelocity

void
ConstraintIBMethod::calculateDrag()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;

    std::vector<std::vector<double> > inertia_force(d_no_structures, std::vector<double>(3, 0.0));
    std::vector<std::vector<double> > constraint_force(d_no_structures, std::vector<double>(3, 0.0));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& U_new_data = *d_l_data_U_new[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_current_data = *d_l_data_U_current[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_correction_data = *d_l_data_U_correction[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U_new = &U_new_data[local_idx][0];
                    const double* const U_current = &U_current_data[local_idx][0];
                    const double* const U_correction = &U_correction_data[local_idx][0];

                    for (int d = 0; d < NDIM; ++d)
                    {
                        inertia_force[location_struct_handle][d] += U_new[d] - U_current[d];
                        constraint_force[location_struct_handle][d] += U_correction[d];
                    }
                }
            }
        } // all structs
        d_l_data_U_new[ln]->restoreArrays();
        d_l_data_U_current[ln]->restoreArrays();
        d_l_data_U_correction[ln]->restoreArrays();
    }

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        IBTK_MPI::sumReduction(&inertia_force[struct_no][0], 3);
        IBTK_MPI::sumReduction(&constraint_force[struct_no][0], 3);
        for (int d = 0; d < NDIM; ++d)
        {
            inertia_force[struct_no][d] *= (d_rho_solid[struct_no] / dt) * d_vol_element[struct_no];
            constraint_force[struct_no][d] *= (d_rho_solid[struct_no] / dt);
        }
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_drag && (d_timestep_counter % d_output_interval) == 0)
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_drag_force_stream[struct_no]
                << d_FuRMoRP_new_time << '\t' << inertia_force[struct_no][0] << '\t' << inertia_force[struct_no][1]
                << '\t' << inertia_force[struct_no][2] << '\t' << constraint_force[struct_no][0] << '\t'
                << constraint_force[struct_no][1] << '\t' << constraint_force[struct_no][2] << std::endl;
        }
    }

    return;
} // calculateDrag

void
ConstraintIBMethod::calculateTorque()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;

    std::vector<std::vector<double> > inertia_torque(d_no_structures, std::vector<double>(3, 0.0));
    std::vector<std::vector<double> > constraint_torque(d_no_structures, std::vector<double>(3, 0.0));
    double R_cross_U_inertia[3] = { 0.0 }, R_cross_U_constraint[3] = { 0.0 };

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& U_new_data = *d_l_data_U_new[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_current_data = *d_l_data_U_current[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_correction_data = *d_l_data_U_correction[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data = *d_X_new_data[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    const double* const U_new = &U_new_data[local_idx][0];
                    const double* const U_current = &U_current_data[local_idx][0];
                    const double* const U_correction = &U_correction_data[local_idx][0];
                    const double* const X = &X_data[local_idx][0];
#if (NDIM == 2)
                    double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[location_struct_handle][0];
                    double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[location_struct_handle][1];

                    R_cross_U_inertia[2] = (x * (U_new[1] - U_current[1]) - y * (U_new[0] - U_current[0]));
                    R_cross_U_constraint[2] = (x * (U_correction[1]) - y * (U_correction[0]));
#endif

#if (NDIM == 3)
                    double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[location_struct_handle][0];
                    double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[location_struct_handle][1];
                    double z = displacement[2] + X[2] - d_center_of_mass_unshifted_new[location_struct_handle][2];

                    R_cross_U_inertia[0] = (y * (U_new[2] - U_current[2]) - z * (U_new[1] - U_current[1]));

                    R_cross_U_inertia[1] = (-x * (U_new[2] - U_current[2]) + z * (U_new[0] - U_current[0]));

                    R_cross_U_inertia[2] = (x * (U_new[1] - U_current[1]) - y * (U_new[0] - U_current[0]));

                    R_cross_U_constraint[0] = (y * (U_correction[2]) - z * (U_correction[1]));

                    R_cross_U_constraint[1] = (-x * (U_correction[2]) + z * (U_correction[0]));

                    R_cross_U_constraint[2] = (x * (U_correction[1]) - y * (U_correction[0]));
#endif

                    for (int d = 0; d < 3; ++d)
                    {
                        inertia_torque[location_struct_handle][d] += R_cross_U_inertia[d];
                        constraint_torque[location_struct_handle][d] += R_cross_U_constraint[d];
                    }
                }
            }
        } // all structs
        d_l_data_U_new[ln]->restoreArrays();
        d_l_data_U_current[ln]->restoreArrays();
        d_l_data_U_correction[ln]->restoreArrays();
        d_X_new_data[ln]->restoreArrays();
    }
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        IBTK_MPI::sumReduction(&inertia_torque[struct_no][0], 3);
        IBTK_MPI::sumReduction(&constraint_torque[struct_no][0], 3);
        for (int d = 0; d < 3; ++d)
        {
            inertia_torque[struct_no][d] *= (d_rho_solid[struct_no] / dt) * d_vol_element[struct_no];
            constraint_torque[struct_no][d] *= (d_rho_solid[struct_no] / dt);
        }
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_torque && (d_timestep_counter % d_output_interval) == 0)
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_torque_stream[struct_no] << d_FuRMoRP_new_time << '\t' << inertia_torque[struct_no][0] << '\t'
                                        << inertia_torque[struct_no][1] << '\t' << inertia_torque[struct_no][2] << '\t'
                                        << constraint_torque[struct_no][0] << '\t' << constraint_torque[struct_no][1]
                                        << '\t' << constraint_torque[struct_no][2] << std::endl;
        }
    }

    return;
} // calculateTorque

void
ConstraintIBMethod::calculatePower()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = d_FuRMoRP_new_time - d_FuRMoRP_current_time;

    std::vector<std::vector<double> > inertia_power(d_no_structures, std::vector<double>(3, 0.0));
    std::vector<std::vector<double> > constraint_power(d_no_structures, std::vector<double>(3, 0.0));

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& U_new_data = *d_l_data_U_new[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_current_data = *d_l_data_U_current[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& U_correction_data = *d_l_data_U_correction[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U_new = &U_new_data[local_idx][0];
                    const double* const U_current = &U_current_data[local_idx][0];
                    const double* const U_correction = &U_correction_data[local_idx][0];

                    for (int d = 0; d < NDIM; ++d)
                    {
                        inertia_power[location_struct_handle][d] += (U_new[d] - U_current[d]) * U_new[d];
                        constraint_power[location_struct_handle][d] += U_correction[d] * U_new[d];
                    }
                }
            }
        } // all structs
        d_l_data_U_new[ln]->restoreArrays();
        d_l_data_U_current[ln]->restoreArrays();
        d_l_data_U_correction[ln]->restoreArrays();
    }

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        IBTK_MPI::sumReduction(&inertia_power[struct_no][0], 3);
        IBTK_MPI::sumReduction(&constraint_power[struct_no][0], 3);
        for (int d = 0; d < NDIM; ++d)
        {
            inertia_power[struct_no][d] *= (d_rho_solid[struct_no] / dt) * d_vol_element[struct_no];
            constraint_power[struct_no][d] *= (d_rho_solid[struct_no] / dt);
        }
    }

    if (!IBTK_MPI::getRank() && d_print_output && d_output_drag && (d_timestep_counter % d_output_interval) == 0)
    {
        for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
        {
            *d_power_spent_stream[struct_no]
                << d_FuRMoRP_new_time << '\t' << inertia_power[struct_no][0] << '\t' << inertia_power[struct_no][1]
                << '\t' << inertia_power[struct_no][2] << '\t' << constraint_power[struct_no][0] << '\t'
                << constraint_power[struct_no][1] << '\t' << constraint_power[struct_no][2] << std::endl;
        }
    }

    return;
} // calculatePower

void
ConstraintIBMethod::calculateStructureMomentum()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& U_new_data = *d_l_data_U_new[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const double* const U_new = &U_new_data[local_idx][0];

                    for (int d = 0; d < NDIM; ++d)
                    {
                        d_structure_mom[location_struct_handle][d] += U_new[d];
                    }
                }
            }
        } // all structs
        d_l_data_U_new[ln]->restoreArrays();
    }

    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        IBTK_MPI::sumReduction(&d_structure_mom[struct_no][0], 3);
        for (int d = 0; d < NDIM; ++d)
        {
            d_structure_mom[struct_no][d] *= d_rho_solid[struct_no] * d_vol_element[struct_no];
        }
    }

    return;
} // calculateStructureMomentum

void
ConstraintIBMethod::calculateStructureRotationalMomentum()
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    double R_cross_U[3] = { 0.0 };

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;

        const boost::multi_array_ref<double, 2>& U_new_data = *d_l_data_U_new[ln]->getLocalFormVecArray();
        const boost::multi_array_ref<double, 2>& X_data = *d_X_new_data[ln]->getLocalFormVecArray();

        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const size_t structs_on_this_ln = structIDs.size();

        for (size_t struct_no = 0; struct_no < structs_on_this_ln; ++struct_no)
        {
            std::pair<int, int> lag_idx_range =
                d_l_data_manager->getLagrangianStructureIndexRange(structIDs[struct_no], ln);
            Pointer<ConstraintIBKinematics> ptr_ib_kinematics =
                *std::find_if(d_ib_kinematics.begin(), d_ib_kinematics.end(), find_struct_handle(lag_idx_range));
            const int location_struct_handle =
                find_struct_handle_position(d_ib_kinematics.begin(), d_ib_kinematics.end(), ptr_ib_kinematics);

            for (const auto& node_idx : local_nodes)
            {
                const int lag_idx = node_idx->getLagrangianIndex();
                if (lag_idx_range.first <= lag_idx && lag_idx < lag_idx_range.second)
                {
                    const int local_idx = node_idx->getLocalPETScIndex();
                    const Vector& displacement = node_idx->getPeriodicDisplacement();
                    const double* const U_new = &U_new_data[local_idx][0];
                    const double* const X = &X_data[local_idx][0];
#if (NDIM == 2)
                    double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[location_struct_handle][0];
                    double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[location_struct_handle][1];
                    R_cross_U[2] = (x * (U_new[1]) - y * (U_new[0]));
#endif

#if (NDIM == 3)
                    double x = displacement[0] + X[0] - d_center_of_mass_unshifted_new[location_struct_handle][0];
                    double y = displacement[1] + X[1] - d_center_of_mass_unshifted_new[location_struct_handle][1];
                    double z = displacement[2] + X[2] - d_center_of_mass_unshifted_new[location_struct_handle][2];

                    R_cross_U[0] = (y * (U_new[2]) - z * (U_new[1]));

                    R_cross_U[1] = (-x * (U_new[2]) + z * (U_new[0]));

                    R_cross_U[2] = (x * (U_new[1]) - y * (U_new[0]));
#endif

                    for (int d = 0; d < 3; ++d)
                    {
                        d_structure_rotational_mom[location_struct_handle][d] += R_cross_U[d];
                    }
                }
            }
        } // all structs
        d_l_data_U_new[ln]->restoreArrays();
        d_X_new_data[ln]->restoreArrays();
    }
    for (int struct_no = 0; struct_no < d_no_structures; ++struct_no)
    {
        IBTK_MPI::sumReduction(&d_structure_rotational_mom[struct_no][0], 3);
        for (int d = 0; d < 3; ++d)
        {
            d_structure_rotational_mom[struct_no][d] *= d_rho_solid[struct_no] * d_vol_element[struct_no];
        }
    }

    return;
} // calculateStructureRotationalMomentum

} // namespace IBAMR
