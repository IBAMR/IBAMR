// Filename: IBExplicitHierarchyIntegrator.cpp
// Created on 12 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <deque>
#include <ostream>
#include <string>

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/IBExplicitHierarchyIntegrator.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_enums.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBExplicitHierarchyIntegrator::IBExplicitHierarchyIntegrator(const std::string& object_name,
                                                             Pointer<Database> input_db,
                                                             Pointer<IBStrategy> ib_method_ops,
                                                             Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                                                             bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_method_ops, ins_hier_integrator, register_for_restart)
{
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBExplicitHierarchyIntegrator

IBExplicitHierarchyIntegrator::~IBExplicitHierarchyIntegrator()
{
    // intentionally blank
    return;
} // ~IBExplicitHierarchyIntegrator

void
IBExplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                            const double new_time,
                                                            const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate Eulerian scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        if (d_f_current_idx != -1) level->allocatePatchData(d_f_current_idx, current_time);
        if (d_ib_method_ops->hasFluidSources())
        {
            level->allocatePatchData(d_p_idx, current_time);
            level->allocatePatchData(d_q_idx, current_time);
        }
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Initialize IB data.
    d_ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    if (ins_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  attempting to perform "
                                 << d_current_num_cycles
                                 << " cycles of fixed point iteration.\n"
                                 << "  number of cycles required by Navier-Stokes solver = "
                                 << ins_num_cycles
                                 << ".\n"
                                 << "  current implementation requires either that both solvers "
                                    "use the same number of cycles,\n"
                                 << "  or that the IB solver use only a single cycle.\n");
    }
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, ins_num_cycles);

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(current_time);
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): spreading Lagrangian force "
                                     "to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), current_time);
        d_hier_velocity_data_ops->copyData(d_f_current_idx, d_f_idx);
        break;
    case MIDPOINT_RULE:
        // intentionally blank
        break;
    default:
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                                 << "\n"
                                 << "  supported time stepping types are: FORWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
        // intentionally blank
        break;
    case MIDPOINT_RULE:
    case TRAPEZOIDAL_RULE:
        if (num_cycles == 1)
        {
            if (d_enable_logging)
                plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian "
                                         "forward Euler step\n";
            d_ib_method_ops->eulerStep(current_time, new_time);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                                 << "\n"
                                 << "  supported time stepping types are: FORWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
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
        if (cycle_num > 0)
        {
            IBAMR_DO_ONCE(
                {
                    pout << "IBExplicitHierarchyIntegrator::integrateHierarchy():\n"
                         << "  WARNING: time_stepping_type = " << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                         << " but num_cycles = " << d_current_num_cycles << " > 1.\n";
                });
        }
        break;
    case MIDPOINT_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(half_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
        break;
    case TRAPEZOIDAL_RULE:
        if (d_current_num_cycles == 1 || cycle_num > 0)
        {
            // NOTE: If (current_num_cycles > 1 && cycle_num == 0), the
            // force computed here would be the same as that computed above
            // in preprocessIntegrateHierarchy(), so we don't bother to
            // recompute it.
            if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
            d_ib_method_ops->computeLagrangianForce(new_time);
            if (d_enable_logging)
                plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force "
                                         "to the Eulerian grid\n";
            d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
            d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
            d_ib_method_ops->spreadForce(
                d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
            d_hier_velocity_data_ops->linearSum(d_f_idx, 0.5, d_f_current_idx, 0.5, d_f_idx);
        }
        break;
    default:
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                                 << "\n"
                                 << "  supported time stepping types are: FORWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute the Lagrangian source/sink strengths and spread them to the
    // Eulerian grid.
    if (d_ib_method_ops->hasFluidSources())
    {
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): computing Lagrangian fluid source strength\n";
        d_ib_method_ops->computeLagrangianFluidSource(half_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian fluid source "
                                     "strength to the Eulerian grid\n";
        d_hier_pressure_data_ops->setToScalar(d_q_idx, 0.0);
        d_ib_method_ops->spreadFluidSource(d_q_idx, getProlongRefineSchedules(d_object_name + "::q"), half_time);
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
        // intentionally blank
        break;
    case MIDPOINT_RULE:
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian velocity to "
                                     "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             half_time);
        break;
    case TRAPEZOIDAL_RULE:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian velocity to "
                                     "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    default:
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                                 << "\n"
                                 << "  supported time stepping types are: FORWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    if (d_current_num_cycles > 1 && d_current_cycle_num == 0)
    {
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): performing Lagrangian forward-Euler step\n";
        d_ib_method_ops->eulerStep(current_time, new_time);
    }
    else
    {
        switch (d_time_stepping_type)
        {
        case FORWARD_EULER:
            // intentionally blank
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
                                     << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                                     << "\n"
                                     << "  supported time stepping types are: FORWARD_EULER, "
                                        "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }

    // Compute the pressure at the updated locations of any distributed internal
    // fluid sources or sinks.
    if (d_ib_method_ops->hasFluidSources())
    {
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian fluid "
                                     "pressure to the Lagrangian mesh\n";
        d_hier_pressure_data_ops->copyData(d_p_idx, p_new_idx);
        d_p_phys_bdry_op->setPatchDataIndex(d_p_idx);
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
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                                 "velocity to the Lagrangian mesh\n";
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_ib_method_ops->interpolateVelocity(d_u_idx,
                                         getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                         getGhostfillRefineSchedules(d_object_name + "::u"),
                                         new_time);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
    double cfl_max = 0.0;
    PatchCellDataOpsReal<NDIM, double> patch_cc_ops;
    PatchSideDataOpsReal<NDIM, double> patch_sc_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx, dx + NDIM));
            Pointer<CellData<NDIM, double> > u_cc_new_data = patch->getPatchData(u_new_idx);
            Pointer<SideData<NDIM, double> > u_sc_new_data = patch->getPatchData(u_new_idx);
            double u_max = 0.0;
            if (u_cc_new_data) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (u_sc_new_data) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max * dt / dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): estimated upper bound on IB "
                                 "point displacement since last regrid = "
             << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, ins_num_cycles);

    // Deallocate IB data.
    d_ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        if (d_f_current_idx != -1) level->deallocatePatchData(d_f_current_idx);
        if (d_ib_method_ops->hasFluidSources())
        {
            level->deallocatePatchData(d_p_idx);
            level->deallocatePatchData(d_q_idx);
        }
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

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

/////////////////////////////// PROTECTED ////////////////////////////////////

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
                                 << " not found in restart file."
                                 << std::endl);
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
