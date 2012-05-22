// Filename: AdvDiffCenteredHierarchyIntegrator.C
// Created on 22 May 2012 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "AdvDiffCenteredHierarchyIntegrator.h"

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
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <tbox/NullDatabase.h>

// C++ STDLIB INCLUDES
#include <limits>

// FORTRAN ROUTINES
#if (NDIM == 2)
#define ADVECT_STABLEDT_FC FC_FUNC_(advect_stabledt2d, ADVECT_STABLEDT2D)
#endif

#if (NDIM == 3)
#define ADVECT_STABLEDT_FC FC_FUNC_(advect_stabledt3d, ADVECT_STABLEDT3D)
#endif

extern "C"
{
    void
    ADVECT_STABLEDT_FC(
        const double*,
#if (NDIM == 2)
        const int& , const int& , const int& , const int& ,
        const int& , const int& ,
        const double* , const double* ,
#endif
#if (NDIM == 3)
        const int& , const int& , const int& , const int& , const int& , const int& ,
        const int& , const int& , const int& ,
        const double* , const double* , const double* ,
#endif
        double&
                       );
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Version of AdvDiffCenteredHierarchyIntegrator restart file data.
static const int ADV_DIFF_CENTERED_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffCenteredHierarchyIntegrator::AdvDiffCenteredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : AdvDiffHierarchyIntegrator(object_name, input_db, register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif
    // Default values.
    d_cfl_max = 0.5;

    // Initialize object with data read from the input database.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);
    return;
}// AdvDiffCenteredHierarchyIntegrator

AdvDiffCenteredHierarchyIntegrator::~AdvDiffCenteredHierarchyIntegrator()
{
    // intentionally blank
    return;
}// ~AdvDiffCenteredHierarchyIntegrator

void
AdvDiffCenteredHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Register variables using the default variable registration routine.
    AdvDiffHierarchyIntegrator::registerVariables();

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffHierarchyIntegrator.
    AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

int
AdvDiffCenteredHierarchyIntegrator::getNumberOfCycles() const
{
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (initial_time)
    {
        return std::max(2,d_num_cycles);
    }
    else
    {
        return d_num_cycles;
    }
}// getNumberOfCycles

void
AdvDiffCenteredHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Indicate that all solvers need to be reinitialized if the current
    // timestep size is different from the previous one.
    if (initial_time || !MathUtilities<double>::equalEps(dt,d_dt_previous[0]))
    {
        std::fill(d_helmholtz_solvers_need_init.begin(),d_helmholtz_solvers_need_init.end(), true);
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Keep track of the number of cycles to be used for the present integration
    // step.
    d_num_cycles_step = num_cycles;

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    unsigned int l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var   = *cit;
        Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];
        const double kappa  = d_Q_diffusion_coef[Q_var];
        const double lambda = d_Q_damping_coef  [Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];

        Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int Psi_scratch_idx = var_db->mapVariableAndContextToIndex(Psi_var, getScratchContext());

        // Setup the operators and solvers and compute the right-hand-side terms.
        d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx, false);
        d_hier_bdry_fill_ops[l]->setHomogeneousBc(false);
        d_hier_bdry_fill_ops[l]->fillData(current_time);
        PoissonSpecifications& helmholtz_spec = d_helmholtz_specs[l];
        Pointer<CCLaplaceOperator> helmholtz_op = d_helmholtz_ops[l];
        switch (d_viscous_timestepping_type)
        {
            case BACKWARD_EULER:
            {
                // The backward Euler discretization is:
                //
                //     (I-dt*kappa*L(t_new)) Q(n+1) = Q(n) + F(t_avg) dt
                //
                // where
                //
                //    t_new = (n+1) dt
                //    t_avg = (t_new+t_old)/2
                helmholtz_spec.setCConstant(1.0+dt*lambda);
                helmholtz_spec.setDConstant(   -dt*kappa );

                PoissonSpecifications rhs_spec("rhs_spec");
                rhs_spec.setCConstant(1.0);
                rhs_spec.setDConstant(0.0);

                for (int depth = 0; depth < Q_depth; ++depth)
                {
                    d_hier_math_ops->laplace(
                        Psi_scratch_idx, Psi_var,  // Psi(n+1/2)
                        rhs_spec,                  // Poisson spec
                        Q_scratch_idx  , Q_var  ,  // Q(n)
                        d_no_fill_op,              // don't need to re-fill Q(n) data
                        current_time,              // Q(n) bdry fill time
                        0.0, -1, NULL,
                        depth, depth, depth);      // dst_depth, src1_depth, src2_depth
                }
                break;
            }
            case CRANK_NICOLSON:
            {
                // The Crank-Nicolson discretization is:
                //
                //     (I-0.5*dt*kappa*L(t_new)) Q(n+1) = (I+0.5*dt*kappa*L(t_old)) Q(n) + F(t_avg) dt
                //
                // where
                //
                //    t_old = n dt
                //    t_new = (n+1) dt
                //    t_avg = (t_new+t_old)/2
                helmholtz_spec.setCConstant(1.0+0.5*dt*lambda);
                helmholtz_spec.setDConstant(   -0.5*dt*kappa );

                PoissonSpecifications rhs_spec("rhs_spec");
                rhs_spec.setCConstant(1.0-0.5*dt*lambda);
                rhs_spec.setDConstant(   +0.5*dt*kappa );

                for (int depth = 0; depth < Q_depth; ++depth)
                {
                    d_hier_math_ops->laplace(
                        Psi_scratch_idx, Psi_var,  // Psi(n+1/2)
                        rhs_spec,                  // Poisson spec
                        Q_scratch_idx  , Q_var  ,  // Q(n)
                        d_no_fill_op,              // don't need to re-fill Q(n) data
                        current_time,              // Q(n) bdry fill time
                        0.0, -1, NULL,
                        depth, depth, depth);      // dst_depth, src1_depth, src2_depth
                }
                break;
            }
            default:
                TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                           << "  unrecognized viscous timestepping type: " << enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type) << "." << std::endl);
        }

        helmholtz_op->setPoissonSpecifications(helmholtz_spec);
        helmholtz_op->setPhysicalBcCoefs(Q_bc_coef);
        helmholtz_op->setHomogeneousBc(false);
        helmholtz_op->setTime(new_time);
        helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        Pointer<CCPoissonPointRelaxationFACOperator> helmholtz_fac_op = d_helmholtz_fac_ops[l];
        Pointer<FACPreconditioner>                   helmholtz_fac_pc = d_helmholtz_fac_pcs[l];
        Pointer<KrylovLinearSolver>                  helmholtz_solver = d_helmholtz_solvers[l];

        if (d_helmholtz_solvers_need_init[l])
        {
            if (d_do_log) plog << d_object_name << "::preprocessIntegrateHierarchy(): initializing Helmholtz solvers for variable number " << l << ", dt = " << dt << "\n";
            if (d_using_FAC)
            {
                helmholtz_fac_op->setPoissonSpecifications(helmholtz_spec);
                helmholtz_fac_op->setTime(new_time);
                helmholtz_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
            }
            helmholtz_solver->initializeSolverState(*d_sol_vecs[l],*d_rhs_vecs[l]);

            // Indicate that the solvers do not presently require
            // re-initialization.
            d_helmholtz_solvers_need_init[l] = false;
        }

        // Set the initial guess.
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_current_idx);
    }

    // Update the advection velocity.
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_idx, u_var, d_hierarchy, current_time);
        }
    }
    return;
}// preprocessIntegrateHierarchy

void
AdvDiffCenteredHierarchyIntegrator::integrateHierarchy(
    const double /*current_time*/,
    const double new_time,
    const int /*cycle_num*/)
{
//  const int coarsest_ln = 0;
//  const int finest_ln = d_hierarchy->getFinestLevelNumber();
//  const double dt  = new_time-current_time;
//  const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Perform a single step of fixed point iteration.
    unsigned int l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var   = *cit;
//      Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];

//      const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
//      const int Psi_scratch_idx = var_db->mapVariableAndContextToIndex(Psi_var, getScratchContext());

        Pointer<CCLaplaceOperator>                   helmholtz_op     = d_helmholtz_ops    [l];
        Pointer<CCPoissonPointRelaxationFACOperator> helmholtz_fac_op = d_helmholtz_fac_ops[l];
        Pointer<KrylovLinearSolver>                  helmholtz_solver = d_helmholtz_solvers[l];

        // Solve for Q(n+1).
        switch (d_viscous_timestepping_type)
        {
            case BACKWARD_EULER:
            case CRANK_NICOLSON:
                helmholtz_op->setTime(new_time);
                if (d_using_FAC) helmholtz_fac_op->setTime(new_time);
                helmholtz_solver->solveSystem(*d_sol_vecs[l],*d_rhs_vecs[l]);
                d_hier_cc_data_ops->copyData(Q_new_idx, Q_scratch_idx);

                if (d_do_log) plog << d_object_name << "::integrateHierarchy(): linear solve number of iterations = " << helmholtz_solver->getNumIterations() << "\n";
                if (d_do_log) plog << d_object_name << "::integrateHierarchy(): linear solve residual norm        = " << helmholtz_solver->getResidualNorm()  << "\n";

                if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
                {
                    pout << d_object_name << "::integrateHierarchy():"
                         <<"  WARNING: linear solver iterations == max iterations\n";
                }
                break;
            default:
                TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                           << "  unrecognized viscous timestepping type: " << enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type) << "." << std::endl);
        }
    }
    return;
}// integrateHierarchy

void
AdvDiffCenteredHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double /*current_time*/,
    const double new_time,
    const bool /*skip_synchronize_new_state_data*/,
    const int /*num_cycles*/)
{
    // Update the advection velocity.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_idx, u_var, d_hierarchy, new_time);
        }
    }
    return;
}// postprocessIntegrateHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

double
AdvDiffCenteredHierarchyIntegrator::getTimeStepSizeSpecialized()
{
    double dt = d_dt_max;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Index<NDIM>& ilower = patch_box.lower();
            const Index<NDIM>& iupper = patch_box.upper();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
            {
                Pointer<FaceVariable<NDIM,double> > u_var = *cit;
                Pointer<FaceData<NDIM,double> > u_data = patch->getPatchData(u_var, getCurrentContext());
                const IntVector<NDIM>& u_ghost_cells = u_data->getGhostCellWidth();
                double stable_dt = std::numeric_limits<double>::max();
#if (NDIM == 2)
                ADVECT_STABLEDT_FC(
                    dx,
                    ilower(0),iupper(0),ilower(1),iupper(1),
                    u_ghost_cells(0),u_ghost_cells(1),
                    u_data->getPointer(0),u_data->getPointer(1),
                    stable_dt);
#endif
#if (NDIM == 3)
                ADVECT_STABLEDT_FC(
                    dx,
                    ilower(0),iupper(0),ilower(1),iupper(1),ilower(2),iupper(2),
                    u_ghost_cells(0),u_ghost_cells(1),u_ghost_cells(2),
                    u_data->getPointer(0),u_data->getPointer(1),u_data->getPointer(2),
                    stable_dt);
#endif
                dt = std::min(dt,stable_dt);
            }
        }
    }
    if (!initial_time && d_dt_growth_factor >= 1.0)
    {
        dt = std::min(dt,d_dt_growth_factor*d_dt_previous[0]);
    }
    return dt;
}// getTimeStepSizeSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvDiffCenteredHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool /*is_from_restart*/)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    if      (db->keyExists("cfl_max")) d_cfl_max = db->getDouble("cfl_max");
    else if (db->keyExists("CFL_max")) d_cfl_max = db->getDouble("CFL_max");
    if      (db->keyExists("cfl"    )) d_cfl_max = db->getDouble("cfl"    );
    else if (db->keyExists("CFL"    )) d_cfl_max = db->getDouble("CFL"    );
    return;
}// getFromInput

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
