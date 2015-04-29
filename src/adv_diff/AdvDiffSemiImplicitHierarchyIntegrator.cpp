// Filename: AdvDiffSemiImplicitHierarchyIntegrator.cpp
// Created on 22 May 2012 by Boyce Griffith
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

#include <stddef.h>
#include <algorithm>
#include <deque>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyDataOpsManager.h"
#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchFaceDataOpsReal.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibamr/AdvDiffConvectiveOperatorManager.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/PoissonSolver.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
class Box;
} // namespace hier
} // namespace SAMRAI

namespace SAMRAI
{
namespace solv
{
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffSemiImplicitHierarchyIntegrator::AdvDiffSemiImplicitHierarchyIntegrator(
    const std::string& object_name,
    const boost::shared_ptr<Database>& input_db,
    bool register_for_restart)
    : AdvDiffHierarchyIntegrator(object_name, input_db, register_for_restart)
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);

    // Default values.
    d_default_convective_time_stepping_type = MIDPOINT_RULE;
    d_default_init_convective_time_stepping_type = MIDPOINT_RULE;
    d_default_convective_op_type = AdvDiffConvectiveOperatorManager::DEFAULT;
    d_default_convective_op_input_db =
        boost::make_shared<MemoryDatabase>(d_object_name + "::default_convective_op_input_db");

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Check to make sure the time stepping types are supported.
    switch (d_default_diffusion_time_stepping_type)
    {
    case BACKWARD_EULER:
    case FORWARD_EULER:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::AdvDiffSemiImplicitHierarchyIntegrator():\n"
                                 << "  unsupported default diffusion time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_default_diffusion_time_stepping_type) << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }

    switch (d_default_convective_time_stepping_type)
    {
    case ADAMS_BASHFORTH:
    case FORWARD_EULER:
    case MIDPOINT_RULE:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::AdvDiffSemiImplicitHierarchyIntegrator():\n"
                                 << "  unsupported default convective time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_default_convective_time_stepping_type) << " \n"
                                 << "  valid choices are: ADAMS_BASHFORTH, FORWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    if (is_multistep_time_stepping_type(d_default_convective_time_stepping_type))
    {
        switch (d_default_init_convective_time_stepping_type)
        {
        case FORWARD_EULER:
        case MIDPOINT_RULE:
        case TRAPEZOIDAL_RULE:
            break;
        default:
            TBOX_ERROR(d_object_name << "::AdvDiffSemiImplicitHierarchyIntegrator():\n"
                                     << "  unsupported default initial convective time stepping type: "
                                     << enum_to_string<TimeSteppingType>(d_default_init_convective_time_stepping_type)
                                     << " \n"
                                     << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }
    return;
}

AdvDiffSemiImplicitHierarchyIntegrator::~AdvDiffSemiImplicitHierarchyIntegrator()
{
    // intentionally blank
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setDefaultConvectiveTimeSteppingType(
    TimeSteppingType default_convective_time_stepping_type)
{
    d_default_convective_time_stepping_type = default_convective_time_stepping_type;
    return;
}

TimeSteppingType AdvDiffSemiImplicitHierarchyIntegrator::getDefaultConvectiveTimeSteppingType() const
{
    return d_default_convective_time_stepping_type;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setDefaultInitialConvectiveTimeSteppingType(
    TimeSteppingType default_init_convective_time_stepping_type)
{
    d_default_init_convective_time_stepping_type = default_init_convective_time_stepping_type;
    return;
}

TimeSteppingType AdvDiffSemiImplicitHierarchyIntegrator::getDefaultInitialConvectiveTimeSteppingType() const
{
    return d_default_init_convective_time_stepping_type;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setDefaultConvectiveOperatorType(const std::string& op_type)
{
    d_default_convective_op_type = op_type;
    return;
}

const std::string& AdvDiffSemiImplicitHierarchyIntegrator::getDefaultConvectiveOperatorType() const
{
    return d_default_convective_op_type;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setDefaultConvectiveOperatorInputDatabase(
    const boost::shared_ptr<Database>& input_db)
{
    d_default_convective_op_input_db = input_db;
    return;
}

boost::shared_ptr<Database> AdvDiffSemiImplicitHierarchyIntegrator::getDefaultConvectiveOperatorInputDatabase() const
{
    return d_default_convective_op_input_db;
}

void AdvDiffSemiImplicitHierarchyIntegrator::registerTransportedQuantity(
    const boost::shared_ptr<CellVariable<double>>& Q_var)
{
    AdvDiffHierarchyIntegrator::registerTransportedQuantity(Q_var);

    // Set default values.
    d_Q_convective_time_stepping_type[Q_var] = d_default_convective_time_stepping_type;
    d_Q_init_convective_time_stepping_type[Q_var] = d_default_init_convective_time_stepping_type;
    d_Q_convective_op_type[Q_var] = d_default_convective_op_type;
    d_Q_convective_op_input_db[Q_var] = d_default_convective_op_input_db;
    d_Q_convective_op[Q_var] = NULL;
    d_Q_convective_op_needs_init[Q_var] = false;
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveTimeSteppingType(
    const boost::shared_ptr<CellVariable<double>>& Q_var,
    TimeSteppingType convective_time_stepping_type)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    d_Q_convective_time_stepping_type[Q_var] = convective_time_stepping_type;
    return;
}

TimeSteppingType AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveTimeSteppingType(
    const boost::shared_ptr<CellVariable<double>>& Q_var) const
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    return d_Q_convective_time_stepping_type.find(Q_var)->second;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setInitialConvectiveTimeSteppingType(
    const boost::shared_ptr<CellVariable<double>>& Q_var,
    TimeSteppingType init_convective_time_stepping_type)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    d_Q_init_convective_time_stepping_type[Q_var] = init_convective_time_stepping_type;
    return;
}

TimeSteppingType AdvDiffSemiImplicitHierarchyIntegrator::getInitialConvectiveTimeSteppingType(
    const boost::shared_ptr<CellVariable<double>>& Q_var) const
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    return d_Q_init_convective_time_stepping_type.find(Q_var)->second;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperatorType(
    const boost::shared_ptr<CellVariable<double>>& Q_var,
    const std::string& op_type)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    d_Q_convective_op_type[Q_var] = op_type;
    return;
}

const std::string& AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveOperatorType(
    const boost::shared_ptr<CellVariable<double>>& Q_var) const
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    return d_Q_convective_op_type.find(Q_var)->second;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperatorInputDatabase(
    const boost::shared_ptr<CellVariable<double>>& Q_var,
    const boost::shared_ptr<Database>& input_db)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    d_Q_convective_op_input_db[Q_var] = input_db;
    return;
}

boost::shared_ptr<Database> AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveOperatorInputDatabase(
    const boost::shared_ptr<CellVariable<double>>& Q_var) const
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    return d_Q_convective_op_input_db.find(Q_var)->second;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperator(
    const boost::shared_ptr<CellVariable<double>>& Q_var,
    const boost::shared_ptr<ConvectiveOperator>& convective_op)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    TBOX_ASSERT(!d_integrator_is_initialized);
    d_Q_convective_op[Q_var] = convective_op;
    d_Q_convective_op_needs_init[Q_var] = true;
    return;
}

boost::shared_ptr<ConvectiveOperator>
AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveOperator(const boost::shared_ptr<CellVariable<double>>& Q_var)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    if (!d_Q_convective_op[Q_var])
    {
        auto convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
        d_Q_convective_op[Q_var] = convective_op_manager->allocateOperator(
            d_Q_convective_op_type[Q_var], d_object_name + "::" + Q_var->getName() + "::ConvectiveOperator", Q_var,
            d_Q_convective_op_input_db[Q_var], d_Q_difference_form[Q_var], d_Q_bc_coef[Q_var]);
        d_Q_convective_op_needs_init[Q_var] = true;
    }
    return d_Q_convective_op[Q_var];
}

void AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperatorsNeedInit()
{
    for (auto it = d_Q_var.begin(); it != d_Q_var.end(); ++it)
    {
        setConvectiveOperatorNeedsInit(*it);
    }
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperatorNeedsInit(
    const boost::shared_ptr<CellVariable<double>>& Q_var)
{
    TBOX_ASSERT(std::find(d_Q_var.begin(), d_Q_var.end(), Q_var) != d_Q_var.end());
    d_Q_convective_op_needs_init[Q_var] = true;
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(
    const boost::shared_ptr<PatchHierarchy>& hierarchy,
    const boost::shared_ptr<GriddingAlgorithm>& gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    auto grid_geom = BOOST_CAST<CartesianGridGeometry>(d_hierarchy->getGridGeometry());

    // Obtain the Hierarchy data operations objects.
    auto hier_ops_manager = HierarchyDataOpsManager::getManager();
    auto fc_var = boost::make_shared<FaceVariable<double>>(DIM, "fc_var");
    d_hier_fc_data_ops =
        BOOST_CAST<HierarchyFaceDataOpsReal<double>>(hier_ops_manager->getOperationsDouble(fc_var, hierarchy, true));

    // Register variables using the default variable registration routine.
    AdvDiffHierarchyIntegrator::registerVariables();

    // Setup the convective operators.
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        auto Q_var = *cit;
        getConvectiveOperator(Q_var);
    }

    // Register additional variables required for present time stepping algorithm.
    const IntVector cell_ghosts(DIM, CELLG);
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        auto Q_var = *cit;
        const int Q_depth = Q_var->getDepth();

        auto N_var = boost::make_shared<CellVariable<double>>(DIM, Q_var->getName() + "::N", Q_depth);
        d_N_var.insert(N_var);
        d_Q_N_map[Q_var] = N_var;
        int N_scratch_idx;
        registerVariable(N_scratch_idx, N_var, cell_ghosts, getScratchContext());

        auto N_old_var = boost::make_shared<CellVariable<double>>(DIM, Q_var->getName() + "::N_old", Q_depth);
        d_N_old_var.insert(N_old_var);
        d_Q_N_old_map[Q_var] = N_old_var;
        int N_old_current_idx, N_old_new_idx, N_old_scratch_idx;
        registerVariable(N_old_current_idx, N_old_new_idx, N_old_scratch_idx, N_old_var, cell_ghosts,
                         "CONSERVATIVE_COARSEN", "CONSERVATIVE_LINEAR_REFINE");
    }

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffHierarchyIntegrator.
    AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}

int AdvDiffSemiImplicitHierarchyIntegrator::getNumberOfCycles() const
{
    int num_cycles = d_num_cycles;
    if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
    {
        for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
        {
            auto Q_var = *cit;
            if (!d_Q_u_map.find(Q_var)->second) continue;
            if (is_multistep_time_stepping_type(d_Q_convective_time_stepping_type.find(Q_var)->second) &&
                d_Q_init_convective_time_stepping_type.find(Q_var)->second != FORWARD_EULER)
            {
                num_cycles = std::max(2, num_cycles);
            }
        }
    }
    return num_cycles;
}

void AdvDiffSemiImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                          const double new_time,
                                                                          const int num_cycles)
{
    AdvDiffHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    auto var_db = VariableDatabase::getDatabase();

    // Indicate that all solvers need to be reinitialized if the current
    // timestep size is different from the previous one.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0]);
    if (dt_change)
    {
        std::fill(d_helmholtz_solvers_need_init.begin(), d_helmholtz_solvers_need_init.end(), true);
        std::fill(d_helmholtz_rhs_ops_need_init.begin(), d_helmholtz_rhs_ops_need_init.end(), true);
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Update the advection velocity.
    for (auto cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        auto u_var = *cit;
        const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
        const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
        const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
        if (d_u_fcn[u_var])
        {
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_current_idx, u_var, d_hierarchy, current_time);
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_new_idx, u_var, d_hierarchy, new_time);
        }
        else
        {
            d_hier_fc_data_ops->copyData(u_new_idx, u_current_idx);
        }
        d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    }

    // Update the diffusion coefficient
    for (auto cit = d_diffusion_coef_var.begin(); cit != d_diffusion_coef_var.end(); ++cit)
    {
        auto D_var = *cit;
        auto D_fcn = d_diffusion_coef_fcn[D_var];
        if (D_fcn)
        {
            const int D_current_idx = var_db->mapVariableAndContextToIndex(D_var, getCurrentContext());
            D_fcn->setDataOnPatchHierarchy(D_current_idx, D_var, d_hierarchy, current_time);
        }
    }

    // Setup the operators and solvers and compute the right-hand-side terms.
    unsigned int l = 0;
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit, ++l)
    {
        auto Q_var = *cit;
        auto Q_rhs_var = d_Q_Q_rhs_map[Q_var];
        auto D_var = d_Q_diffusion_coef_variable[Q_var];
        auto D_rhs_var = d_diffusion_coef_rhs_map[D_var];
        TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];
        const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& Q_bc_coef = d_Q_bc_coef[Q_var];

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());
        const int D_current_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getCurrentContext()) : -1);
        const int D_scratch_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getScratchContext()) : -1);
        const int D_rhs_scratch_idx =
            (D_rhs_var ? var_db->mapVariableAndContextToIndex(D_rhs_var, getScratchContext()) : -1);

        // Setup the problem coefficients for the linear solve for Q(n+1).
        double K = 0.0;
        switch (diffusion_time_stepping_type)
        {
        case BACKWARD_EULER:
            K = 1.0;
            break;
        case FORWARD_EULER:
            K = 0.0;
            break;
        case TRAPEZOIDAL_RULE:
            K = 0.5;
            break;
        default:
            TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                     << "  unsupported diffusion time stepping type: "
                                     << enum_to_string<TimeSteppingType>(diffusion_time_stepping_type) << " \n"
                                     << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
        }
        PoissonSpecifications solver_spec(d_object_name + "::solver_spec::" + Q_var->getName());
        PoissonSpecifications rhs_op_spec(d_object_name + "::rhs_op_spec::" + Q_var->getName());
        solver_spec.setCConstant(1.0 / dt + K * lambda);
        rhs_op_spec.setCConstant(1.0 / dt - (1.0 - K) * lambda);
        if (isDiffusionCoefficientVariable(Q_var))
        {
            // set -K*kappa in solver_spec
            d_hier_sc_data_ops->scale(D_scratch_idx, -K, D_current_idx);
            solver_spec.setDPatchDataId(D_scratch_idx);
            // set (1.0-K)*kappa in rhs_op_spec
            d_hier_sc_data_ops->scale(D_rhs_scratch_idx, (1.0 - K), D_current_idx);
            rhs_op_spec.setDPatchDataId(D_rhs_scratch_idx);
        }
        else
        {
            const double kappa = d_Q_diffusion_coef[Q_var];
            solver_spec.setDConstant(-K * kappa);
            rhs_op_spec.setDConstant(+(1.0 - K) * kappa);
        }

        // Initialize the RHS operator and compute the RHS vector.
        auto helmholtz_rhs_op = d_helmholtz_rhs_ops[l];
        helmholtz_rhs_op->setPoissonSpecifications(rhs_op_spec);
        helmholtz_rhs_op->setPhysicalBcCoefs(Q_bc_coef);
        helmholtz_rhs_op->setHomogeneousBc(false);
        helmholtz_rhs_op->setSolutionTime(current_time);
        helmholtz_rhs_op->setTimeInterval(current_time, new_time);
        if (d_helmholtz_rhs_ops_need_init[l])
        {
            if (d_enable_logging)
            {
                plog << d_object_name << ": "
                     << "Initializing Helmholtz RHS operator for variable number " << l << "\n";
            }
            helmholtz_rhs_op->initializeOperatorState(*d_sol_vecs[l], *d_rhs_vecs[l]);
            d_helmholtz_rhs_ops_need_init[l] = false;
        }
        d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx, false);
        helmholtz_rhs_op->apply(*d_sol_vecs[l], *d_rhs_vecs[l]);

        // Initialize the linear solver.
        auto helmholtz_solver = d_helmholtz_solvers[l];
        helmholtz_solver->setPoissonSpecifications(solver_spec);
        helmholtz_solver->setPhysicalBcCoefs(Q_bc_coef);
        helmholtz_solver->setHomogeneousBc(false);
        helmholtz_solver->setSolutionTime(new_time);
        helmholtz_solver->setTimeInterval(current_time, new_time);
        if (d_helmholtz_solvers_need_init[l])
        {
            if (d_enable_logging)
            {
                plog << d_object_name << ": "
                     << "Initializing Helmholtz solvers for variable number " << l << "\n";
            }
            helmholtz_solver->initializeSolverState(*d_sol_vecs[l], *d_rhs_vecs[l]);
            d_helmholtz_solvers_need_init[l] = false;
        }

        // Account for the convective difference term.
        auto u_var = d_Q_u_map[Q_var];
        if (u_var)
        {
            auto N_var = d_Q_N_map[Q_var];
            auto N_old_var = d_Q_N_old_map[Q_var];
            TimeSteppingType convective_time_stepping_type = d_Q_convective_time_stepping_type[Q_var];
            if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(convective_time_stepping_type))
            {
                convective_time_stepping_type = d_Q_init_convective_time_stepping_type[Q_var];
            }
            if ((num_cycles == 1) &&
                (convective_time_stepping_type == MIDPOINT_RULE || convective_time_stepping_type == TRAPEZOIDAL_RULE))
            {
                TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                         << "  time stepping type: "
                                         << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                                         << " requires num_cycles > 1.\n"
                                         << "  at current time step, num_cycles = " << num_cycles << "\n");
            }
            if (d_Q_convective_op_needs_init[Q_var])
            {
                d_Q_convective_op[Q_var]->initializeOperatorState(*d_sol_vecs[l], *d_rhs_vecs[l]);
                d_Q_convective_op_needs_init[Q_var] = false;
            }
            const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
            d_Q_convective_op[Q_var]->setAdvectionVelocity(u_current_idx);
            const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
            const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx);
            d_Q_convective_op[Q_var]->setSolutionTime(current_time);
            d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
            const int N_old_new_idx = var_db->mapVariableAndContextToIndex(N_old_var, getNewContext());
            d_hier_cc_data_ops->copyData(N_old_new_idx, N_scratch_idx);
            if (convective_time_stepping_type == FORWARD_EULER)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }

        // Set the initial guess.
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_current_idx);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                                const double new_time,
                                                                const int cycle_num)
{
    AdvDiffHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    auto var_db = VariableDatabase::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles
                     << " but expected num_cycles = " << expected_num_cycles << ".\n";
            });
    }

    // Perform a single step of fixed point iteration.
    unsigned int l = 0;
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit, ++l)
    {
        auto Q_var = *cit;
        auto F_var = d_Q_F_map[Q_var];
        auto Q_rhs_var = d_Q_Q_rhs_map[Q_var];

        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int F_scratch_idx =
            d_F_fcn[F_var] ? var_db->mapVariableAndContextToIndex(F_var, getScratchContext()) : -1;
        const int F_new_idx = d_F_fcn[F_var] ? var_db->mapVariableAndContextToIndex(F_var, getNewContext()) : -1;
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());

        // Update the advection velocity.
        if (cycle_num > 0)
        {
            for (auto cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
            {
                auto u_var = *cit;
                const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
                const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
                const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
                if (d_u_fcn[u_var])
                {
                    d_u_fcn[u_var]->setDataOnPatchHierarchy(u_new_idx, u_var, d_hierarchy, new_time);
                }
                else
                {
                    d_hier_fc_data_ops->copyData(u_new_idx, u_current_idx);
                }
                d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
            }
        }

        // Account for the convective difference term.
        auto u_var = d_Q_u_map[Q_var];
        auto N_var = d_Q_N_map[Q_var];
        TimeSteppingType convective_time_stepping_type = UNKNOWN_TIME_STEPPING_TYPE;
        if (u_var)
        {
            auto N_old_var = d_Q_N_old_map[Q_var];
            convective_time_stepping_type = d_Q_convective_time_stepping_type[Q_var];
            if (is_multistep_time_stepping_type(convective_time_stepping_type))
            {
                TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
                if (getIntegratorStep() == 0)
                {
                    convective_time_stepping_type = d_Q_init_convective_time_stepping_type[Q_var];
                }
                else if (cycle_num > 0)
                {
                    convective_time_stepping_type = MIDPOINT_RULE;
                    IBAMR_DO_ONCE(
                        {
                            pout << "AdvDiffSemiImplicitHierarchyIntegrator::"
                                    "integrateHierarchy():"
                                    "\n"
                                 << "  WARNING: convective_time_stepping_type = "
                                 << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var])
                                 << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                                 << "           using "
                                 << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var])
                                 << " only for the first cycle in each time step;\n"
                                 << "           using "
                                 << enum_to_string<TimeSteppingType>(convective_time_stepping_type)
                                 << " for subsequent cycles.\n";
                        });
                }
            }
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            if (cycle_num > 0)
            {
                if (convective_time_stepping_type == MIDPOINT_RULE)
                {
                    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
                    d_Q_convective_op[Q_var]->setAdvectionVelocity(u_scratch_idx);
                    const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
                    const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
                    const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
                    d_hier_cc_data_ops->linearSum(Q_scratch_idx, 0.5, Q_current_idx, 0.5, Q_new_idx);
                    d_Q_convective_op[Q_var]->setSolutionTime(half_time);
                    d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
                }
                else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
                {
                    const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
                    d_Q_convective_op[Q_var]->setAdvectionVelocity(u_new_idx);
                    const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
                    const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
                    d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_new_idx);
                    d_Q_convective_op[Q_var]->setSolutionTime(new_time);
                    d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
                }
            }
            if (convective_time_stepping_type == ADAMS_BASHFORTH)
            {
                TBOX_ASSERT(cycle_num == 0);
                const int N_old_current_idx = var_db->mapVariableAndContextToIndex(N_old_var, getCurrentContext());
                const double omega = dt / d_dt_previous[0];
                d_hier_cc_data_ops->linearSum(N_scratch_idx, 1.0 + 0.5 * omega, N_scratch_idx, -0.5 * omega,
                                              N_old_current_idx);
            }
            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }

        // Account for forcing terms.
        if (d_F_fcn[F_var])
        {
            d_F_fcn[F_var]->setDataOnPatchHierarchy(F_scratch_idx, F_var, d_hierarchy, half_time);
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, 1.0, F_scratch_idx, Q_rhs_scratch_idx);
        }

        // Solve for Q(n+1).
        auto helmholtz_solver = d_helmholtz_solvers[l];
        helmholtz_solver->solveSystem(*d_sol_vecs[l], *d_rhs_vecs[l]);
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_scratch_idx);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): diffusion solve number of iterations = "
                 << helmholtz_solver->getNumIterations() << "\n";
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): diffusion solve residual norm        = "
                 << helmholtz_solver->getResidualNorm() << "\n";
        if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
        {
            pout << d_object_name << "::integrateHierarchy():"
                 << "  WARNING: linear solver iterations == max iterations\n";
        }

        // Reset the right-hand side vector.
        if (u_var)
        {
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +1.0, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +0.5, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }
        if (d_F_fcn[F_var])
        {
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0, F_scratch_idx, Q_rhs_scratch_idx);
            d_hier_cc_data_ops->copyData(F_new_idx, F_scratch_idx);
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                           const double new_time,
                                                                           const bool skip_synchronize_new_state_data,
                                                                           const int num_cycles)
{
    AdvDiffHierarchyIntegrator::postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data,
                                                              num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Determine the CFL number.
    double cfl_max = 0.0;
    auto var_db = VariableDatabase::getDatabase();
    for (unsigned int k = 0; k < d_u_var.size(); ++k)
    {
        const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_var[k], getNewContext());
        PatchFaceDataOpsReal<double> patch_fc_ops;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            auto level = d_hierarchy->getPatchLevel(ln);
            for (auto p = level->begin(); p != level->end(); ++p)
            {
                auto patch = *p;
                const Box& patch_box = patch->getBox();
                const auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
                const double* const dx = pgeom->getDx();
                const double dx_min = *(std::min_element(dx, dx + NDIM));
                auto u_fc_new_data = BOOST_CAST<FaceData<double>>(patch->getPatchData(u_new_idx));
                double u_max = 0.0;
                u_max = patch_fc_ops.maxNorm(u_fc_new_data, patch_box);
                cfl_max = std::max(cfl_max, u_max * dt / dx_min);
            }
        }
    }
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.AllReduce(&cfl_max, 1, MPI_MAX);
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(current_time, new_time, skip_synchronize_new_state_data,
                                                     num_cycles);
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void AdvDiffSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const boost::shared_ptr<PatchHierarchy>& hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);
    for (auto cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        auto Q_var = *cit;
        d_Q_convective_op_needs_init[Q_var] = true;
    }
    AdvDiffHierarchyIntegrator::resetHierarchyConfigurationSpecialized(hierarchy, coarsest_level, finest_level);
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::putToRestartSpecialized(const boost::shared_ptr<Database>& db) const
{
    db->putString("d_default_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_default_convective_time_stepping_type));
    db->putString("d_default_init_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_default_init_convective_time_stepping_type));
    db->putString("d_default_convective_op_type", d_default_convective_op_type);
    AdvDiffHierarchyIntegrator::putToRestartSpecialized(db);
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void AdvDiffSemiImplicitHierarchyIntegrator::getFromInput(const boost::shared_ptr<Database>& db, bool is_from_restart)
{
    TBOX_ASSERT(db);

    if (!is_from_restart)
    {
        if (db->keyExists("convective_time_stepping_type"))
            d_default_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("convective_time_stepping_type"));
        else if (db->keyExists("convective_timestepping_type"))
            d_default_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("convective_timestepping_type"));
        else if (db->keyExists("default_convective_time_stepping_type"))
            d_default_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("default_convective_time_stepping_type"));
        else if (db->keyExists("default_convective_timestepping_type"))
            d_default_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("default_convective_timestepping_type"));

        if (db->keyExists("init_convective_time_stepping_type"))
            d_default_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("init_convective_time_stepping_type"));
        else if (db->keyExists("init_convective_timestepping_type"))
            d_default_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("init_convective_timestepping_type"));
        else if (db->keyExists("default_init_convective_time_stepping_type"))
            d_default_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("default_init_convective_time_stepping_type"));
        else if (db->keyExists("default_init_convective_timestepping_type"))
            d_default_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("default_init_convective_timestepping_type"));

        if (db->keyExists("convective_op_type"))
            d_default_convective_op_type = db->getString("convective_op_type");
        else if (db->keyExists("convective_operator_type"))
            d_default_convective_op_type = db->getString("convective_operator_type");
        else if (db->keyExists("default_convective_op_type"))
            d_default_convective_op_type = db->getString("default_convective_op_type");
        else if (db->keyExists("default_convective_operator_type"))
            d_default_convective_op_type = db->getString("default_convective_operator_type");

        if (db->keyExists("convective_op_db"))
            d_default_convective_op_input_db = db->getDatabase("convective_op_db");
        else if (db->keyExists("default_convective_op_db"))
            d_default_convective_op_input_db = db->getDatabase("default_convective_op_db");
    }
    return;
}

void AdvDiffSemiImplicitHierarchyIntegrator::getFromRestart()
{
    auto restart_db = RestartManager::getManager()->getRootDatabase();
    boost::shared_ptr<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    d_default_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_default_convective_time_stepping_type"));
    d_default_init_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_default_init_convective_time_stepping_type"));
    d_default_convective_op_type = db->getString("d_default_convective_op_type");
    return;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
