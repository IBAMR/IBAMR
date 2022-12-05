// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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
#include "ibamr/MassIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"

#include "BasePatchHierarchy.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CoarseFineBoundary.h"
#include "FaceData.h"
#include "HierarchyDataOpsManager.h"
#include "HierarchySideDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include <array>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
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
// NOTE: The number of ghost cells required by the convection scheme depends
// on the chosen convective limiter, which will be set via input file
static const int GPPMG = 4;
static const int GCUIG = 3;

// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_integrate;
static Timer* t_initialize_integrator;
static Timer* t_deallocate_integrator;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////
MassIntegrator::MassIntegrator(std::string object_name, Pointer<Database> input_db)
    : d_object_name(std::move(object_name)), d_u_bc_coefs(NDIM), d_rho_bc_coefs(NDIM)
{
    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type"))
        {
            d_density_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        }
        if (input_db->keyExists("density_bdry_extrap_type"))
        {
            d_density_bdry_extrap_type = input_db->getString("density_bdry_extrap_type");
        }
        if (input_db->keyExists("convective_limiter"))
        {
            d_density_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("convective_limiter"));
        }
        if (input_db->keyExists("density_convective_limiter"))
        {
            d_density_convective_limiter =
                IBAMR::string_to_enum<LimiterType>(input_db->getString("density_convective_limiter"));
        }

        if (input_db->keyExists("density_time_stepping_type"))
        {
            d_density_time_stepping_type =
                IBAMR::string_to_enum<TimeSteppingType>(input_db->getString("density_time_stepping_type"));
        }
        if (input_db->keyExists("enable_logging"))
        {
            d_enable_logging = input_db->getBool("enable_logging");
        }
    }

    switch (d_density_convective_limiter)
    {
    case PPM:
        d_density_limiter_gcw = GPPMG;
        break;
    case CUI:
        d_density_limiter_gcw = GCUIG;
        break;
    default:
        TBOX_ERROR(
            "MassIntegrator::"
            "MassIntegrator():\n"
            << "  unsupported density convective limiter: "
            << IBAMR::enum_to_string<LimiterType>(d_density_convective_limiter) << " \n"
            << "  valid choices are: PPM, CUI\n");
    }

    switch (d_density_time_stepping_type)
    {
    case FORWARD_EULER:
    case SSPRK2:
        d_num_steps = 2;
        break;
    case SSPRK3:
        d_num_steps = 3;
        break;
    default:
        TBOX_ERROR(
            "MassIntegrator::"
            "MassIntegrator():\n"
            << "  unsupported density time stepping type: "
            << IBAMR::enum_to_string<TimeSteppingType>(d_density_time_stepping_type) << " \n"
            << "  valid choices are: FORWARD_EULER, SSPRK2, SSPRK3\n");
    }

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::"
                                                                                     "applyConvectiveOperator()");
                  t_integrate = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::integrate("
                                                                     ")");
                  t_initialize_integrator = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::"
                                                                                 "initializeTimeIntegrator()");
                  t_deallocate_integrator = TimerManager::getManager()->getTimer("IBAMR::MassIntegrator::"
                                                                                 "deallocateTimeIntegrator()"););
    return;
} // MassIntegrator

void
MassIntegrator::setDensityPatchDataIndex(int rho_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rho_idx >= 0);
#endif
    d_rho_current_idx = rho_idx;
} // setDensityPatchDataIndex

void
MassIntegrator::setConvectiveDerivativePatchDataIndex(int N_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(N_idx >= 0);
#endif
    d_N_idx = N_idx;
} // setConvectiveDerivativePatchDataIndex

void
MassIntegrator::setDensityBoundaryConditions(const std::vector<RobinBcCoefStrategy<NDIM>*>& rho_bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(rho_bc_coefs.size() == NDIM);
#endif
    d_rho_bc_coefs = rho_bc_coefs;
    return;
} // setSideCenteredDensityBoundaryConditions

int
MassIntegrator::getUpdatedDensityPatchDataIndex()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_rho_new_idx >= 0);
#endif
    return d_rho_new_idx;
} // getUpdatedDensityPatchDataIndex

void
MassIntegrator::setFluidVelocityPatchDataIndices(int V_old_idx, int V_current_idx, int V_new_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(V_current_idx >= 0);
#endif

    // Set the old velocity if it has been set, otherwise set to current.
    if (V_old_idx >= 0)
    {
        d_V_old_idx = V_old_idx;
    }
    else
    {
        d_V_old_idx = V_current_idx;
    }

    // Set the current velocity
    d_V_current_idx = V_current_idx;

    // Set the new velocity if it has been set, otherwise set to current.
    if (V_new_idx >= 0)
    {
        d_V_new_idx = V_new_idx;
    }
    else
    {
        d_V_new_idx = V_current_idx;
    }
    return;
} // setFluidVelocityPatchDataIndices

void
MassIntegrator::setCycleNumber(int cycle_num)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(cycle_num >= 0);
#endif
    d_cycle_num = cycle_num;
    return;
} // setCycleNumber

void
MassIntegrator::setSolutionTime(double solution_time)
{
    d_solution_time = solution_time;
} // setSolutionTime

void
MassIntegrator::setTimeInterval(double current_time, double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
} // setTimeInterval

std::pair<double, double>
MassIntegrator::getTimeInterval() const
{
    return std::make_pair(d_current_time, d_new_time);
} // getTimeInterval

double
MassIntegrator::getDt() const
{
    return d_new_time - d_current_time;
} // getDt

void
MassIntegrator::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    d_hier_math_ops = hier_math_ops;
    d_hier_math_ops_external = d_hier_math_ops;
    return;
} // setHierarchyMathOps

Pointer<HierarchyMathOps>
MassIntegrator::getHierarchyMathOps() const
{
    return d_hier_math_ops;
} // getHierarchyMathOps

void
MassIntegrator::setPreviousTimeStepSize(double dt_prev)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(dt_prev > 0.0);
#endif
    d_dt_prev = dt_prev;
    return;
} // setPreviousTimeStepSize

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
