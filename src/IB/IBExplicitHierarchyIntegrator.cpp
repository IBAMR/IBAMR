// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#include "ibamr/IBExplicitHierarchyIntegrator.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/LEInteractor.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_enums.h"

#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsManager.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <ostream>
#include <string>
#include <utility>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBExplicitHierarchyIntegrator restart file data.
static const int IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION = 2;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBExplicitHierarchyIntegrator::IBExplicitHierarchyIntegrator(std::string object_name,
                                                             Pointer<Database> input_db,
                                                             Pointer<IBStrategy> ib_method_ops,
                                                             Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                                             bool register_for_restart)
    : IBHierarchyIntegrator(std::move(object_name), input_db, ib_method_ops, ins_hier_integrator, register_for_restart)
{
    // Set default configuration options.
    d_use_structure_predictor = true;

    // Set options from input.
    if (input_db)
    {
        if (input_db->keyExists("use_structure_predictor"))
            d_use_structure_predictor = input_db->getBool("use_structure_predictor");
        if (input_db->keyExists("IB_delta_fcn")) d_marker_kernel = input_db->getString("IB_delta_fcn");
        if (input_db->keyExists("viz_dump_dirname"))
        {
            d_viz_dump_dirname = input_db->getString("viz_dump_dirname");
            Utilities::recursiveMkdir(d_viz_dump_dirname);
        }
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBExplicitHierarchyIntegrator

void
IBExplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                            const double new_time,
                                                            const int num_cycles)
{
    // preprocess our dependencies...
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(current_time);
        if (d_enable_logging)
            plog << d_object_name
                 << "::preprocessIntegrateHierarchy(): spreading Lagrangian force "
                    "to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), current_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_hier_velocity_data_ops->copyData(d_f_current_idx, d_f_idx);
        break;
    case MIDPOINT_RULE:
        // intentionally blank
        break;
    default:
        TBOX_ERROR(
            d_object_name
            << "::preprocessIntegrateHierarchy():\n"
            << "  unsupported time stepping type: " << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    if (d_use_structure_predictor)
    {
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
        d_ib_method_ops->forwardEulerStep(current_time, new_time);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IBExplicitHierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double half_time = current_time + 0.5 * (new_time - current_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());
    const int p_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case BACKWARD_EULER:
        // intentionally blank
        break;
    case MIDPOINT_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(half_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case TRAPEZOIDAL_RULE:
        if (d_use_structure_predictor || cycle_num > 0)
        {
            // NOTE: We do not re-compute the force unless it could have changed.
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
            d_ib_method_ops->computeLagrangianForce(new_time);
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
            d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
            d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
            d_u_phys_bdry_op->setHomogeneousBc(true);
            d_ib_method_ops->spreadForce(
                d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
            d_u_phys_bdry_op->setHomogeneousBc(false);
            d_hier_velocity_data_ops->linearSum(d_f_idx, 0.5, d_f_current_idx, 0.5, d_f_idx);
        }
        break;
    default:
        TBOX_ERROR(
            d_object_name
            << "::integrateHierarchy():\n"
            << "  unsupported time stepping type: " << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute the Lagrangian source/sink strengths and spread them to the
    // Eulerian grid.
    if (d_ib_method_ops->hasFluidSources())
    {
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): computing Lagrangian fluid source strength\n";
        d_ib_method_ops->computeLagrangianFluidSource(half_time);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): spreading Lagrangian fluid source "
                    "strength to the Eulerian grid\n";
        d_hier_pressure_data_ops->setToScalar(d_q_idx, 0.0);
        // NOTE: This does not correctly treat the case in which the structure
        // is close to the physical boundary.
        d_ib_method_ops->spreadFluidSource(
            d_q_idx, nullptr, getProlongRefineSchedules(d_object_name + "::q"), half_time);
    }

    // Solve the incompressible Navier-Stokes equations.
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): solving the incompressible Navier-Stokes equations\n";
    if (d_current_num_cycles > 1)
    {
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_current_num_cycles == 1);
#endif
        const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
        for (int ins_cycle_num = 0; ins_cycle_num < ins_num_cycles; ++ins_cycle_num)
        {
            d_ins_hier_integrator->integrateHierarchy(current_time, new_time, ins_cycle_num);
        }
    }
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case BACKWARD_EULER:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    case MIDPOINT_RULE:
    {
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             half_time);
    }
    break;
    case TRAPEZOIDAL_RULE:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    default:
        TBOX_ERROR(
            d_object_name
            << "::integrateHierarchy():\n"
            << "  unsupported time stepping type: " << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
            << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    if (d_current_num_cycles > 1 && d_current_cycle_num == 0)
    {
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): performing Lagrangian forward-Euler step\n";
        d_ib_method_ops->forwardEulerStep(current_time, new_time);
    }
    else
    {
        switch (d_time_stepping_type)
        {
        case FORWARD_EULER:
            break;
        case BACKWARD_EULER:
            d_ib_method_ops->backwardEulerStep(current_time, new_time);
            break;
        case MIDPOINT_RULE:
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): performing Lagrangian midpoint-rule step\n";
            d_ib_method_ops->midpointStep(current_time, new_time);
            break;
        case TRAPEZOIDAL_RULE:
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): performing Lagrangian trapezoidal-rule step\n";
            d_ib_method_ops->trapezoidalStep(current_time, new_time);
            break;
        default:
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                     << "  unsupported time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                     << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, "
                                        "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    // Compute the pressure at the updated locations of any distributed internal
    // fluid sources or sinks.
    if (d_ib_method_ops->hasFluidSources())
    {
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian fluid "
                    "pressure to the Lagrangian mesh\n";
        d_hier_pressure_data_ops->copyData(d_p_idx, p_new_idx);
        d_p_phys_bdry_op->setPatchDataIndex(d_p_idx);
        d_p_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolatePressure(d_p_idx,
                                             getCoarsenSchedules(d_object_name + "::p::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::p"),
                                             half_time);
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
IBExplicitHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                             const double new_time,
                                                             const bool skip_synchronize_new_state_data,
                                                             const int num_cycles)
{
    auto ops = HierarchyDataOpsManager<NDIM>::getManager()->getOperationsDouble(d_u_var, d_hierarchy, true);

    auto velocity_ghost_update = [&](const std::vector<int>& indices)
    {
        using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghostfills;
        for (const int& idx : indices)
        {
            ghostfills.emplace_back(idx,
                                    "CONSERVATIVE_LINEAR_REFINE",
                                    /*use_cf_bdry_interpolation*/ true,
                                    "CONSERVATIVE_COARSEN",
                                    "LINEAR",
                                    false,
                                    d_ins_hier_integrator->getVelocityBoundaryConditions());
        }
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghostfills, d_hierarchy);
        ghost_fill_op.fillData(current_time);
    };

    // Update the marker points, should they exist:
    if (d_markers && !d_marker_velocities_set)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                       d_ins_hier_integrator->getCurrentContext());
        // Clear any ghost data outside the domain:
#ifndef NDEBUG
        ops->setToScalar(d_u_idx, std::numeric_limits<double>::quiet_NaN(), false);
#endif
        ops->copyData(d_u_idx, u_current_idx);
        velocity_ghost_update({ d_u_idx });
        d_markers->setVelocities(d_u_idx, d_marker_kernel);
        d_marker_velocities_set = true;
    }

    // The last thing we need to do (before we really postprocess) is update the structure velocity:
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());
#ifndef NDEBUG
    ops->setToScalar(d_u_idx, std::numeric_limits<double>::quiet_NaN(), false);
#endif
    ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging)
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                "velocity to the Lagrangian mesh\n";
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ib_method_ops->interpolateVelocity(d_u_idx,
                                         getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                         getGhostfillRefineSchedules(d_object_name + "::u"),
                                         new_time);

    if (d_markers)
    {
#ifndef NDEBUG
        ops->setToScalar(d_u_idx, std::numeric_limits<double>::quiet_NaN(), false);
#endif
        ops->copyData(d_u_idx, u_new_idx);
        d_hier_velocity_data_ops->linearSum(d_u_half_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        velocity_ghost_update({ d_u_idx, d_u_half_idx });

        const double dt = new_time - current_time;
        switch (d_time_stepping_type)
        {
        case FORWARD_EULER:
            d_markers->forwardEulerStep(dt, d_u_idx, d_marker_kernel);
            break;
        case BACKWARD_EULER:
            d_markers->backwardEulerStep(dt, d_u_idx, d_marker_kernel);
            break;
        case MIDPOINT_RULE:
            d_markers->midpointStep(dt, d_u_half_idx, d_u_idx, d_marker_kernel);
            break;
        case TRAPEZOIDAL_RULE:
            d_markers->trapezoidalStep(dt, d_u_idx, d_marker_kernel);
            break;
        default:
            TBOX_ERROR(d_object_name << "::postprocessIntegrateHierarchy():\n"
                                     << "  unsupported time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                     << "  supported time stepping types are: FORWARD_EULER, BACKWARD_EULER, "
                                        "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
IBExplicitHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                             Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Setup the fluid solver for explicit coupling.
    d_ins_hier_integrator->registerBodyForceFunction(new IBEulerianForceFunction(this));

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    return;
} // initializeHierarchyIntegrator

void
IBExplicitHierarchyIntegrator::setMarkers(const EigenAlignedVector<IBTK::Point>& markers)
{
    if (d_marker_kernel.size() == 0)
    {
        TBOX_ERROR(d_object_name << "::setMarkers():\n To use marker points the IB kernel must be specified in "
                                    "the input database via IB_kernel_fcn.");
    }
    if (!d_hierarchy)
    {
        TBOX_ERROR(
            d_object_name << "::setMarkers():\n Markers may only be set after calling initializePatchHierarchy().");
    }

    // Ensure that whichever patch data indices we need to exist are present.
    if (d_u_half_idx == IBTK::invalid_index)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        d_u_half_idx = var_db->registerClonedPatchDataIndex(getVelocityVariable(), d_u_idx);
        d_ib_data.setFlag(d_u_half_idx);
    }

    EigenAlignedVector<IBTK::Vector> velocities(markers.size());
    // Eigen 'bug': no default initialization
    IBTK::Vector v;
    v.fill(0.0);
    std::fill(velocities.begin(), velocities.end(), v);
    // explicitly destroy the present markers (should they exist) since the
    // new one will have the same name. This is necessary because otherwise we
    // would have two objects with the same name in the restart database.
    d_markers.setNull();
    d_markers = new MarkerPatchHierarchy(d_object_name + "::markers", d_hierarchy, markers, velocities);
    d_marker_velocities_set = false;
} // setMarkers

std::pair<EigenAlignedVector<IBTK::Point>, EigenAlignedVector<IBTK::Vector> >
IBExplicitHierarchyIntegrator::collectAllMarkers() const
{
    if (d_markers)
        return d_markers->collectAllMarkers();
    else
    {
        // Some older compilers cannot convert just '{}' into a std::pair so help it out
        EigenAlignedVector<IBTK::Point> positions;
        EigenAlignedVector<IBTK::Vector> velocities;
        return std::make_pair(positions, velocities);
    }
} // collectAllMarkers

void
IBExplicitHierarchyIntegrator::writeMarkerPlotData(const int time_step,
                                                   const double simulation_time,
                                                   const bool save_velocites) const
{
    (void)simulation_time;
    if (d_markers)
    {
        if (d_viz_dump_dirname.size() == 0)
        {
            TBOX_ERROR(
                d_object_name << "::writeMarkerPlotData():\n This function requires that viz_dump_dirname was set "
                                 "in the input database.");
        }
        d_markers->writeH5Part(d_viz_dump_dirname + "/markerpoints-" + std::to_string(time_step) + ".h5part",
                               time_step,
                               simulation_time,
                               save_velocites);
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBExplicitHierarchyIntegrator::regridHierarchyBeginSpecialized()
{
    IBHierarchyIntegrator::regridHierarchyBeginSpecialized();

    if (d_markers)
    {
        d_regrid_temporary_data = new IBExplicitHierarchyIntegrator::RegridData();
        auto pair = d_markers->collectAllMarkers();
        d_regrid_temporary_data->d_marker_positions = std::move(pair.first);
        d_regrid_temporary_data->d_marker_velocities = std::move(pair.second);
        d_markers->reinit({}, {});
    }
} // regridHierarchyBeginSpecialized

void
IBExplicitHierarchyIntegrator::regridHierarchyEndSpecialized()
{
    if (d_regrid_temporary_data)
    {
        d_markers->reinit(d_regrid_temporary_data->d_marker_positions, d_regrid_temporary_data->d_marker_velocities);
        d_regrid_temporary_data = nullptr;
    }

    IBHierarchyIntegrator::regridHierarchyEndSpecialized();
} // regridHierarchyEndSpecialized

void
IBExplicitHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION", IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBExplicitHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_EXPLICIT_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
