// Filename: AdvDiffPredictorCorrectorHierarchyIntegrator.cpp
// Created on 17 Mar 2004 by Boyce Griffith
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
#include <string>
#include <vector>

#include "BaseGriddingAlgorithm.h"
#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellDataFactory.h"
#include "CellVariable.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "Geometry.h"
#include "GriddingAlgorithm.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "HyperbolicLevelIntegrator.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchFaceDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PoissonSpecifications.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "VisItDataWriter.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h"
#include "ibamr/AdvDiffPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/PoissonSolver.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/NullDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
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
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffPredictorCorrectorHierarchyIntegrator::AdvDiffPredictorCorrectorHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor,
    bool register_for_restart)
    : AdvDiffHierarchyIntegrator(object_name, input_db, register_for_restart),
      d_hyp_level_integrator(NULL),
      d_hyp_level_integrator_db(NULL),
      d_hyp_patch_ops(NULL),
      d_hyp_patch_ops_db(NULL),
      d_explicit_predictor(explicit_predictor)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);
    TBOX_ASSERT(explicit_predictor);
#endif
    // Get initialization data for the hyperbolic patch strategy objects.
    if (input_db->keyExists("HyperbolicLevelIntegrator"))
    {
        d_hyp_level_integrator_db = input_db->getDatabase("HyperbolicLevelIntegrator");
    }
    else
    {
        d_hyp_level_integrator_db = new NullDatabase();
    }
    if (input_db->keyExists("AdvDiffPredictorCorrectorHyperbolicPatchOps"))
    {
        d_hyp_patch_ops_db = input_db->getDatabase("AdvDiffPredictorCorrectorHyperbolicPatchOps");
    }
    else
    {
        d_hyp_patch_ops_db = new NullDatabase();
    }

    // Check to make sure the time stepping types are supported.
    switch (d_default_diffusion_time_stepping_type)
    {
    case BACKWARD_EULER:
    case FORWARD_EULER:
    case TRAPEZOIDAL_RULE:
        break;
    default:
        TBOX_ERROR(d_object_name << "::AdvDiffPredictorCorrectorHierarchyIntegrator():\n"
                                 << "  unsupported default diffusion time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_default_diffusion_time_stepping_type)
                                 << " \n"
                                 << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }
    return;
} // AdvDiffPredictorCorrectorHierarchyIntegrator

AdvDiffPredictorCorrectorHierarchyIntegrator::~AdvDiffPredictorCorrectorHierarchyIntegrator()
{
    // intentionally blank
    return;
} // ~AdvDiffPredictorCorrectorHierarchyIntegrator

Pointer<HyperbolicLevelIntegrator<NDIM> >
AdvDiffPredictorCorrectorHierarchyIntegrator::getHyperbolicLevelIntegrator() const
{
    return d_hyp_level_integrator;
} // getHyperbolicLevelIntegrator

Pointer<AdvDiffPredictorCorrectorHyperbolicPatchOps>
AdvDiffPredictorCorrectorHierarchyIntegrator::getHyperbolicPatchStrategy() const
{
    return d_hyp_patch_ops;
} // getHyperbolicPatchStrategy

void
AdvDiffPredictorCorrectorHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                           const double new_time,
                                                                           const int num_cycles)
{
    AdvDiffHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
AdvDiffPredictorCorrectorHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Initialize the HyperbolicPatchStrategy and HyperbolicLevelIntegrator
    // objects that provide numerical routines for explicitly integrating the
    // advective terms.
    d_hyp_patch_ops =
        new AdvDiffPredictorCorrectorHyperbolicPatchOps(d_object_name + "::AdvDiffPredictorCorrectorHyperbolicPatchOps",
                                                        d_hyp_patch_ops_db,
                                                        d_explicit_predictor,
                                                        grid_geom,
                                                        d_registered_for_restart);
    d_hyp_level_integrator = new HyperbolicLevelIntegrator<NDIM>(d_object_name + "::HyperbolicLevelIntegrator",
                                                                 d_hyp_level_integrator_db,
                                                                 d_hyp_patch_ops,
                                                                 d_registered_for_restart,
                                                                 /*using_time_refinement*/ false);

    // Setup variable contexts.
    d_current_context = d_hyp_level_integrator->getCurrentContext();
    d_new_context = d_hyp_level_integrator->getNewContext();

    // Register the VisIt data writer with the patch strategy object, which is
    // the object that actually registers variables for plotting.
    if (d_visit_writer)
    {
        d_hyp_patch_ops->registerVisItDataWriter(d_visit_writer);
    }

    // Register variables with the hyperbolic level integrator.
    for (std::vector<Pointer<FaceVariable<NDIM, double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end();
         ++cit)
    {
        Pointer<FaceVariable<NDIM, double> > u_var = *cit;
        d_hyp_patch_ops->registerAdvectionVelocity(u_var);
        d_hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var, d_u_is_div_free[u_var]);
        if (d_u_fcn[u_var]) d_hyp_patch_ops->setAdvectionVelocityFunction(u_var, d_u_fcn[u_var]);
    }

    const IntVector<NDIM> cell_ghosts = CELLG;
    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_F_var.begin(); cit != d_F_var.end();
         ++cit)
    {
        Pointer<CellVariable<NDIM, double> > F_var = *cit;
        d_hyp_level_integrator->registerVariable(F_var,
                                                 cell_ghosts,
                                                 HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                                 d_hierarchy->getGridGeometry(),
                                                 "CONSERVATIVE_COARSEN",
                                                 "CONSERVATIVE_LINEAR_REFINE");
    }

    for (std::vector<Pointer<SideVariable<NDIM, double> > >::const_iterator cit = d_diffusion_coef_var.begin();
         cit != d_diffusion_coef_var.end();
         ++cit)
    {
        Pointer<SideVariable<NDIM, double> > D_var = *cit;
        d_hyp_level_integrator->registerVariable(D_var,
                                                 cell_ghosts,
                                                 HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                                                 d_hierarchy->getGridGeometry(),
                                                 "CONSERVATIVE_COARSEN",
                                                 "CONSERVATIVE_LINEAR_REFINE");
        int D_scratch_idx;
        registerVariable(D_scratch_idx, D_var, cell_ghosts, getScratchContext());
    }

    for (std::vector<Pointer<SideVariable<NDIM, double> > >::const_iterator cit = d_diffusion_coef_rhs_var.begin();
         cit != d_diffusion_coef_rhs_var.end();
         ++cit)
    {
        Pointer<SideVariable<NDIM, double> > D_rhs_var = *cit;
        int D_rhs_scratch_idx;
        registerVariable(D_rhs_scratch_idx, D_rhs_var, cell_ghosts, getScratchContext());
    }

    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_Q_rhs_var.begin();
         cit != d_Q_rhs_var.end();
         ++cit)
    {
        Pointer<CellVariable<NDIM, double> > Q_rhs_var = *cit;
        int Q_rhs_scratch_idx;
        registerVariable(Q_rhs_scratch_idx, Q_rhs_var, cell_ghosts, getScratchContext());
        d_hyp_patch_ops->registerSourceTerm(Q_rhs_var);
    }

    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end();
         ++cit)
    {
        Pointer<CellVariable<NDIM, double> > Q_var = *cit;
        int Q_scratch_idx;
        registerVariable(Q_scratch_idx, Q_var, cell_ghosts, getScratchContext());
        d_hyp_patch_ops->registerTransportedQuantity(Q_var);
        if (d_Q_u_map[Q_var]) d_hyp_patch_ops->setAdvectionVelocity(Q_var, d_Q_u_map[Q_var]);
        d_hyp_patch_ops->setSourceTerm(Q_var, d_Q_Q_rhs_map[Q_var]);
        d_hyp_patch_ops->setConvectiveDifferencingType(Q_var, d_Q_difference_form[Q_var]);
        if (d_Q_init[Q_var]) d_hyp_patch_ops->setInitialConditions(Q_var, d_Q_init[Q_var]);
        if (!d_Q_bc_coef[Q_var].empty()) d_hyp_patch_ops->setPhysicalBcCoefs(Q_var, d_Q_bc_coef[Q_var]);
    }

    // Initialize the HyperbolicLevelIntegrator.
    //
    // WARNING: This must be done AFTER all variables have been registered with
    // the level integrator.
    d_hyp_level_integrator->initializeLevelIntegrator(d_gridding_alg);

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffHierarchyIntegrator.
    AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
} // initializeHierarchyIntegrator

void
AdvDiffPredictorCorrectorHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                                 const double new_time,
                                                                 const int cycle_num)
{
    AdvDiffHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt = new_time - current_time;
    const double half_time = current_time + 0.5 * dt;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "AdvDiffPredictorCorrectorHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles
                     << " but expected num_cycles = " << expected_num_cycles << ".\n";
            });
    }

    // Reset time-dependent data when necessary.
    if (cycle_num > 0)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            d_hyp_level_integrator->resetDataToPreadvanceState(d_hierarchy->getPatchLevel(ln));
        }
    }

    // Compute any time-dependent source terms at time-level n.
    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_F_var.begin(); cit != d_F_var.end();
         ++cit)
    {
        Pointer<CellVariable<NDIM, double> > F_var = *cit;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        if (F_fcn && F_fcn->isTimeDependent())
        {
            const int F_current_idx = var_db->mapVariableAndContextToIndex(F_var, getCurrentContext());
            F_fcn->setDataOnPatchHierarchy(F_current_idx, F_var, d_hierarchy, current_time);
        }
    }

    // Compute any time-dependent variable diffusion coefficients at time-level n.
    for (std::vector<Pointer<SideVariable<NDIM, double> > >::const_iterator cit = d_diffusion_coef_var.begin();
         cit != d_diffusion_coef_var.end();
         ++cit)
    {
        Pointer<SideVariable<NDIM, double> > D_var = *cit;
        Pointer<CartGridFunction> D_fcn = d_diffusion_coef_fcn[D_var];
        if (D_fcn)
        {
            const int D_current_idx = var_db->mapVariableAndContextToIndex(D_var, getCurrentContext());
            D_fcn->setDataOnPatchHierarchy(D_current_idx, D_var, d_hierarchy, current_time);
        }
    }

    // Predict the advective terms and synchronize them across all levels of the
    // patch hierarchy.
    unsigned int l = 0;
    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end();
         ++cit, ++l)
    {
        Pointer<CellVariable<NDIM, double> > Q_var = *cit;
        Pointer<CellVariable<NDIM, double> > F_var = d_Q_F_map[Q_var];
        Pointer<SideVariable<NDIM, double> > D_var = d_Q_diffusion_coef_variable[Q_var];
        Pointer<CellVariable<NDIM, double> > Q_rhs_var = d_Q_Q_rhs_map[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];

        Pointer<CellDataFactory<NDIM, double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int F_current_idx = (F_var ? var_db->mapVariableAndContextToIndex(F_var, getCurrentContext()) : -1);
        const int D_current_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getCurrentContext()) : -1);
        const int Q_rhs_current_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getCurrentContext());

        // Allocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_scratch_idx, current_time);
        }

        // Setup the right hand side for the advective flux prediction.
        d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx, false);
        d_hier_bdry_fill_ops[l]->setHomogeneousBc(false);
        d_hier_bdry_fill_ops[l]->fillData(current_time);

        PoissonSpecifications kappa_spec("kappa_spec");
        kappa_spec.setCConstant(-lambda);
        if (isDiffusionCoefficientVariable(Q_var))
        {
            kappa_spec.setDPatchDataId(D_current_idx);
        }
        else
        {
            const double kappa = d_Q_diffusion_coef[Q_var];
            kappa_spec.setDConstant(kappa);
        }

        for (int depth = 0; depth < Q_depth; ++depth)
        {
            d_hier_math_ops->laplace(Q_rhs_current_idx,
                                     Q_rhs_var,  // Q_rhs(n)
                                     kappa_spec, // Poisson spec
                                     Q_scratch_idx,
                                     Q_var,        // Q(n)
                                     d_no_fill_op, // don't need to re-fill Q(n) data
                                     current_time, // Q(n) bdry fill time
                                     1.0,          // gamma
                                     F_current_idx,
                                     F_var, // F(n)
                                     depth,
                                     depth,
                                     depth); // dst_depth, src1_depth, src2_depth
        }

        // Deallocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_scratch_idx);
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        static const bool first_step = true;
        static const bool last_step = false;
        d_hyp_level_integrator->advanceLevel(
            d_hierarchy->getPatchLevel(ln), d_hierarchy, current_time, new_time, first_step, last_step);
    }

    if (finest_ln > 0)
    {
        d_hyp_level_integrator->standardLevelSynchronization(
            d_hierarchy, coarsest_ln, finest_ln, new_time, current_time);
    }

    // Compute any time-dependent source terms at time-level n+1/2.
    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_F_var.begin(); cit != d_F_var.end();
         ++cit, ++l)
    {
        Pointer<CellVariable<NDIM, double> > F_var = *cit;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        if (F_fcn && F_fcn->isTimeDependent())
        {
            const int F_current_idx = var_db->mapVariableAndContextToIndex(F_var, getCurrentContext());
            F_fcn->setDataOnPatchHierarchy(F_current_idx, F_var, d_hierarchy, half_time);
        }
    }

    // Compute any time-dependent variable diffusion coefficients at time-level n+1/2.
    for (std::vector<Pointer<SideVariable<NDIM, double> > >::const_iterator cit = d_diffusion_coef_var.begin();
         cit != d_diffusion_coef_var.end();
         ++cit, ++l)
    {
        Pointer<SideVariable<NDIM, double> > D_var = *cit;
        Pointer<CartGridFunction> D_fcn = d_diffusion_coef_fcn[D_var];
        if (D_fcn)
        {
            const int D_current_idx = var_db->mapVariableAndContextToIndex(D_var, getCurrentContext());
            D_fcn->setDataOnPatchHierarchy(D_current_idx, D_var, d_hierarchy, half_time);
        }
    }

    // Indicate that all linear solvers need to be reinitialized if the current
    // timestep size is different from the previous one.
    if (cycle_num == 0 && (initial_time || !MathUtilities<double>::equalEps(dt, d_dt_previous[0])))
    {
        std::fill(d_helmholtz_solvers_need_init.begin(), d_helmholtz_solvers_need_init.end(), true);
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Solve for Q(n+1).
    l = 0;
    for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end();
         ++cit, ++l)
    {
        Pointer<CellVariable<NDIM, double> > Q_var = *cit;
        Pointer<CellVariable<NDIM, double> > F_var = d_Q_F_map[Q_var];
        Pointer<SideVariable<NDIM, double> > D_var = d_Q_diffusion_coef_variable[Q_var];
        Pointer<SideVariable<NDIM, double> > D_rhs_var = d_diffusion_coef_rhs_map[D_var];
        Pointer<CellVariable<NDIM, double> > Q_rhs_var = d_Q_Q_rhs_map[Q_var];
        TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int F_current_idx = (F_var ? var_db->mapVariableAndContextToIndex(F_var, getCurrentContext()) : -1);
        const int D_current_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getCurrentContext()) : -1);
        const int D_scratch_idx = (D_var ? var_db->mapVariableAndContextToIndex(D_var, getScratchContext()) : -1);
        const int D_rhs_scratch_idx =
            (D_rhs_var ? var_db->mapVariableAndContextToIndex(D_rhs_var, getScratchContext()) : -1);
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());

        // Allocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_scratch_idx, current_time);
            level->allocatePatchData(Q_rhs_scratch_idx, new_time);
            if (isDiffusionCoefficientVariable(Q_var))
            {
                level->allocatePatchData(D_scratch_idx, new_time);
                level->allocatePatchData(D_rhs_scratch_idx, new_time);
            }
        }

        if (F_var)
        {
            d_hier_cc_data_ops->add(Q_new_idx, F_current_idx, Q_new_idx);
        }

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
                                     << enum_to_string<TimeSteppingType>(diffusion_time_stepping_type)
                                     << " \n"
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
        Pointer<LaplaceOperator> helmholtz_rhs_op = d_helmholtz_rhs_ops[l];
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
        d_hier_cc_data_ops->add(Q_rhs_scratch_idx, Q_rhs_scratch_idx, Q_new_idx);

        // Initialize the linear solver.
        Pointer<PoissonSolver> helmholtz_solver = d_helmholtz_solvers[l];
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

        // Solve for Q(n+1).
        helmholtz_solver->solveSystem(*d_sol_vecs[l], *d_rhs_vecs[l]);
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_scratch_idx);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): linear solve number of iterations = "
                 << helmholtz_solver->getNumIterations() << "\n";
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): linear solve residual norm        = "
                 << helmholtz_solver->getResidualNorm() << "\n";
        if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
        {
            pout << d_object_name << "::integrateHierarchy():"
                 << "  WARNING: linear solver iterations == max iterations\n";
        }

        // Deallocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_scratch_idx);
            level->deallocatePatchData(Q_rhs_scratch_idx);
            if (isDiffusionCoefficientVariable(Q_var))
            {
                level->deallocatePatchData(D_scratch_idx);
                level->deallocatePatchData(D_rhs_scratch_idx);
            }
        }
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy

void
AdvDiffPredictorCorrectorHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                            const double new_time,
                                                                            const bool skip_synchronize_new_state_data,
                                                                            const int num_cycles)
{
    AdvDiffHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;

    // Determine the CFL number.
    double cfl_max = 0.0;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    for (unsigned int k = 0; k < d_u_var.size(); ++k)
    {
        const int u_new_idx = var_db->mapVariableAndContextToIndex(d_u_var[k], getNewContext());
        PatchFaceDataOpsReal<NDIM, double> patch_fc_ops;
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
                Pointer<FaceData<NDIM, double> > u_fc_new_data = patch->getPatchData(u_new_idx);
                double u_max = 0.0;
                u_max = patch_fc_ops.maxNorm(u_fc_new_data, patch_box);
                cfl_max = std::max(cfl_max, u_max * dt / dx_min);
            }
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

double
AdvDiffPredictorCorrectorHierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = HierarchyIntegrator::getMaximumTimeStepSizeSpecialized();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt = std::min(dt, d_hyp_level_integrator->getLevelDt(level, d_integrator_time, initial_time));
    }
    return dt;
} // getMaximumTimeStepSizeSpecialized

void
AdvDiffPredictorCorrectorHierarchyIntegrator::resetTimeDependentHierarchyDataSpecialized(const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_dt_previous.push_front(new_time - d_integrator_time);
    static const unsigned int MAX_DT_PREVIOUS_SIZE = 31;
    if (d_dt_previous.size() > MAX_DT_PREVIOUS_SIZE) d_dt_previous.pop_back();
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Reset the time dependent data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->resetTimeDependentData(
            d_hierarchy->getPatchLevel(ln), d_integrator_time, d_gridding_alg->levelCanBeRefined(ln));
    }
    return;
} // resetTimeDependentHierarchyDataSpecialized

void
AdvDiffPredictorCorrectorHierarchyIntegrator::resetIntegratorToPreadvanceStateSpecialized()
{
    // We use the HyperbolicLevelIntegrator to handle most data management.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->resetDataToPreadvanceState(d_hierarchy->getPatchLevel(ln));
    }
    return;
} // resetIntegratorToPreadvanceStateSpecialized

void
AdvDiffPredictorCorrectorHierarchyIntegrator::initializeLevelDataSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
#if !defined(NDEBUG)
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
    if (old_level)
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
    // We use the HyperbolicLevelIntegrator to handle as much data management as
    // possible.
    d_hyp_level_integrator->initializeLevelData(
        hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // Set the initial values of any forcing terms and variable-coefficient
    // diffusion coefficient variables.  All other variables are initialized by
    // the hyperbolic level integrator.
    if (initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        for (std::vector<Pointer<CellVariable<NDIM, double> > >::const_iterator cit = d_F_var.begin();
             cit != d_F_var.end();
             ++cit)
        {
            Pointer<CellVariable<NDIM, double> > F_var = *cit;
            const int F_idx = var_db->mapVariableAndContextToIndex(F_var, getCurrentContext());
            Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
            if (F_fcn)
            {
                F_fcn->setDataOnPatchLevel(F_idx, F_var, level, init_data_time, initial_time);
            }
            else
            {
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double> > F_data = patch->getPatchData(F_idx);
#if !defined(NDEBUG)
                    TBOX_ASSERT(F_data);
#endif
                    F_data->fillAll(0.0);
                }
            }
        }

        // Set the initial value of any variable diffusion coefficient
        for (std::vector<Pointer<SideVariable<NDIM, double> > >::const_iterator cit = d_diffusion_coef_var.begin();
             cit != d_diffusion_coef_var.end();
             ++cit)
        {
            Pointer<SideVariable<NDIM, double> > D_var = *cit;
            const int D_idx = var_db->mapVariableAndContextToIndex(D_var, getCurrentContext());
            Pointer<CartGridFunction> D_fcn = d_diffusion_coef_fcn[D_var];
            if (D_fcn)
            {
                D_fcn->setDataOnPatchLevel(D_idx, D_var, level, init_data_time, initial_time);
            }
            else
            {
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<SideData<NDIM, double> > D_data = patch->getPatchData(D_idx);
#if !defined(NDEBUG)
                    TBOX_ASSERT(D_data);
#endif
                    D_data->fillAll(0.0);
                }
            }
        }
    }
    return;
} // initializeLevelDataSpecialized

void
AdvDiffPredictorCorrectorHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    d_hyp_level_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    AdvDiffHierarchyIntegrator::resetHierarchyConfigurationSpecialized(base_hierarchy, coarsest_level, finest_level);
    return;
} // resetHierarchyConfigurationSpecialized

void
AdvDiffPredictorCorrectorHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    // Tag cells for refinement according to the criteria specified by the
    // criteria specified by the level integrator.
    d_hyp_level_integrator->applyGradientDetector(
        hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
} // applyGradientDetectorSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
