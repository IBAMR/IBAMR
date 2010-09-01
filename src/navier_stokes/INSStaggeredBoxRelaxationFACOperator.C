// Filename: INSStaggeredBoxRelaxationFACOperator.C
// Created on 11 Jun 2010 by Boyce Griffith
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

#include "INSStaggeredBoxRelaxationFACOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/NormOps.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/SideSynchCopyFillPattern.h>
#include <ibtk/SideSynchCopyTransactionFactory.h>
#include <ibtk/compiler_hints.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>

// C++ STDLIB INCLUDES
#include <algorithm>
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_restrict_solution;
static Pointer<Timer> t_restrict_residual;
static Pointer<Timer> t_prolong_error_and_correct;
static Pointer<Timer> t_smooth_error;
static Pointer<Timer> t_solve_coarsest_level;
static Pointer<Timer> t_compute_composite_residual_on_level;
static Pointer<Timer> t_compute_residual_norm;
static Pointer<Timer> t_initialize_operator_state;
static Pointer<Timer> t_deallocate_operator_state;

// Number of ghosts cells used for each variable quantity.
static const int GHOSTS = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Size of box to use in the box relaxation.
static const int BOX_GHOSTS = 1;

inline int
compute_side_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box,
    const int axis)
{
    int offset = 0;
    for (int d = 0; d < axis; ++d)
    {
        offset += SideGeometry<NDIM>::toSideBox(box,d).size();
    }
    return offset + SideGeometry<NDIM>::toSideBox(box,axis).offset(i);
}// compute_side_index

inline int
compute_cell_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box)
{
    int offset = 0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        offset += SideGeometry<NDIM>::toSideBox(box,axis).size();
    }
    return box.offset(i) + offset;
}// compute_cell_index

void
buildBoxOperator(
    Mat& A,
    const INSCoefs& problem_coefs,
    const double dt,
    const Box<NDIM>& box,
    const Box<NDIM>& ghost_box,
    const double* const dx)
{
    int ierr;

    const double rho = problem_coefs.getRho();
    const double mu = problem_coefs.getMu();
    const double lambda = problem_coefs.getLambda();

    // Allocate a PETSc matrix for the patch operator.
    Box<NDIM> side_boxes[NDIM];
    BoxList<NDIM> side_ghost_boxes[NDIM];
    for (int axis = 0; axis < NDIM; ++axis)
    {
        side_boxes[axis] = SideGeometry<NDIM>::toSideBox(box, axis);
        side_ghost_boxes[axis] = SideGeometry<NDIM>::toSideBox(ghost_box, axis);
        side_ghost_boxes[axis].removeIntersections(side_boxes[axis]);
    }
    BoxList<NDIM> cell_ghost_boxes(ghost_box);
    cell_ghost_boxes.removeIntersections(box);

    int size = ghost_box.size();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        size += SideGeometry<NDIM>::toSideBox(ghost_box, axis).size();
    }

    static const int U_stencil_sz = 2*NDIM+3;
    static const int P_stencil_sz = 2*NDIM+1;
    std::vector<int> nnz(size, 0);

    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            nnz[compute_side_index(b(), ghost_box, axis)] = U_stencil_sz;
        }
        for (BoxList<NDIM>::Iterator bl(side_ghost_boxes[axis]); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                nnz[compute_side_index(b(), ghost_box, axis)] = 1;
            }
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        nnz[compute_cell_index(b(), ghost_box)] = P_stencil_sz;
    }

    for (BoxList<NDIM>::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            nnz[compute_cell_index(b(), ghost_box)] = 1;
        }
    }

    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, size, size, PETSC_DEFAULT, &nnz[0], &A);  IBTK_CHKERRQ(ierr);

    // Set some general matrix options.
#ifdef DEBUG_CHECK_ASSERTIONS
    ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);  IBTK_CHKERRQ(ierr);
#endif

    // Set the matrix coefficients to correspond to the standard finite
    // difference approximation to the time-dependent incompressible Stokes
    // operator.
    //
    // Note that boundary conditions at both physical boundaries and at
    // coarse-fine interfaces are implicitly treated by setting ghost cell
    // values appropriately.  Thus the matrix coefficients are independent of
    // any boundary conditions.
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            Index<NDIM> i = b();
            const int mat_row = compute_side_index(i, ghost_box, axis);

            std::vector<int> mat_cols(U_stencil_sz,-1);
            std::vector<double> mat_vals(U_stencil_sz,0.0);

            mat_cols[0] = mat_row;
            mat_vals[0] = (rho/dt) + 0.5*lambda;

            for (int d = 0; d < NDIM; ++d)
            {
                Index<NDIM> shift = 0;
                shift(d) = 1;
                const Index<NDIM> i_left = i - shift;
                const Index<NDIM> i_rght = i + shift;
                mat_cols[2*d+1] = compute_side_index(i_left, ghost_box, axis);
                mat_cols[2*d+2] = compute_side_index(i_rght, ghost_box, axis);

                mat_vals[    0] +=      mu/(dx[d]*dx[d]);
                mat_vals[2*d+1]  = -0.5*mu/(dx[d]*dx[d]);
                mat_vals[2*d+2]  = -0.5*mu/(dx[d]*dx[d]);
            }

            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> p_left = i - shift;
            const Index<NDIM> p_rght = i;
            mat_cols[2*NDIM+1] = compute_cell_index(p_left, ghost_box);
            mat_cols[2*NDIM+2] = compute_cell_index(p_rght, ghost_box);

            mat_vals[2*NDIM+1] = -1.0/dx[axis];
            mat_vals[2*NDIM+2] =  1.0/dx[axis];

            static const int m = 1;
            static const int n = U_stencil_sz;
            ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        Index<NDIM> i = b();
        const int mat_row = compute_cell_index(i, ghost_box);

        std::vector<int> mat_cols(P_stencil_sz,-1);
        std::vector<double> mat_vals(P_stencil_sz,0.0);

        mat_cols[0] = mat_row;
        mat_vals[0] = 0.0;

        for (int axis = 0; axis < NDIM; ++axis)
        {
            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> u_left = i;
            const Index<NDIM> u_rght = i + shift;
            mat_cols[2*axis+1] = compute_side_index(u_left, ghost_box, axis);
            mat_cols[2*axis+2] = compute_side_index(u_rght, ghost_box, axis);

            mat_vals[2*axis+1] =  1.0/dx[axis];
            mat_vals[2*axis+2] = -1.0/dx[axis];
        }

        static const int m = 1;
        static const int n = P_stencil_sz;
        ierr = MatSetValues(A, m, &mat_row, n, &mat_cols[0], &mat_vals[0], INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    // Set the entries in the ghost cell region so that ghost cell values are
    // not modified by the smoother.
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (BoxList<NDIM>::Iterator bl(side_ghost_boxes[axis]); bl; bl++)
        {
            for (Box<NDIM>::Iterator b(bl()); b; b++)
            {
                const int i = compute_side_index(b(), ghost_box, axis);
                ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
            }
        }
    }

    for (BoxList<NDIM>::Iterator bl(cell_ghost_boxes); bl; bl++)
    {
        for (Box<NDIM>::Iterator b(bl()); b; b++)
        {
            const int i = compute_cell_index(b(), ghost_box);
            ierr = MatSetValue(A, i, i, 1.0, INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    // Assemble the matrices.
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);  IBTK_CHKERRQ(ierr);
    return;
}// buildBoxOperator

inline void
copyToVec(
    Vec& v,
    const SideData<NDIM,double>& U_data,
    const CellData<NDIM,double>& P_data,
    const Box<NDIM>& box,
    const Box<NDIM>& ghost_box)
{
    int ierr;

    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i = b();
            const SideIndex<NDIM> s_i(i, axis, 0);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecSetValue(v, idx, U_data(s_i), INSERT_VALUES);  IBTK_CHKERRQ(ierr);
        }
    }

    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecSetValue(v, idx, P_data(i), INSERT_VALUES);  IBTK_CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(v);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);  IBTK_CHKERRQ(ierr);
    return;
}// copyToVec

inline void
copyFromVec(
    Vec& v,
    SideData<NDIM,double>& U_data,
    CellData<NDIM,double>& P_data,
    const Box<NDIM>& box,
    const Box<NDIM>& ghost_box)
{
    int ierr;

    static const double omega = 0.5;

    double U;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i = b();
            const SideIndex<NDIM> s_i(i, axis, 0);
            const int idx = compute_side_index(i, ghost_box, axis);
            ierr = VecGetValues(v, 1, &idx, &U);  IBTK_CHKERRQ(ierr);
            U_data(s_i) = (1.0-omega)*U_data(s_i) + omega*U;
        }
    }

    double P;
    for (Box<NDIM>::Iterator b(box); b; b++)
    {
        const Index<NDIM>& i = b();
        const int idx = compute_cell_index(i, ghost_box);
        ierr = VecGetValues(v, 1, &idx, &P);  IBTK_CHKERRQ(ierr);
        P_data(i) = (1.0-omega)*P_data(i) + omega*P;
    }
    return;
}// copyFromVec
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredBoxRelaxationFACOperator::INSStaggeredBoxRelaxationFACOperator(
    const std::string& object_name,
    const INSCoefs& problem_coefs,
    const double dt,
    const Pointer<Database>& input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_solution(NULL),
      d_rhs(NULL),
      d_gcw(GHOSTS),
      d_patch_side_bc_box_overlap(),
      d_patch_cell_bc_box_overlap(),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_problem_coefs(problem_coefs),
      d_dt(dt),
      d_smoother_choice("additive"),
      d_U_prolongation_method("CONSTANT_REFINE"),
      d_P_prolongation_method("LINEAR_REFINE"),
      d_U_restriction_method("CUBIC_COARSEN"),
      d_P_restriction_method("CONSERVATIVE_COARSEN"),
      d_preconditioner(NULL),
      d_fac_max_cycles(1),
      d_fac_uses_presmoothing(false),
      d_fac_initial_guess_nonzero(false),
      d_skip_restrict_sol(true),
      d_skip_restrict_residual(false),
      d_coarse_solver_choice("block_jacobi"),
      d_coarse_solver_tol(1.0e-6),
      d_coarse_solver_max_its(10),
      d_context(NULL),
      d_side_scratch_idx(-1),
      d_cell_scratch_idx(-1),
      d_U_bc_op(NULL),
      d_default_U_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_U_bc_coef", Pointer<Database>(NULL))),
      d_U_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM,d_default_U_bc_coef)),
      d_P_bc_op(NULL),
      d_default_P_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_P_bc_coef", Pointer<Database>(NULL))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_homogeneous_bc(false),
      d_current_time(0.0),
      d_new_time(0.0),
      d_U_op_stencil_fill_pattern(),
      d_P_op_stencil_fill_pattern(),
      d_U_synch_fill_pattern(),
      d_U_prolongation_refine_operator(),
      d_P_prolongation_refine_operator(),
      d_U_cf_bdry_op(),
      d_P_cf_bdry_op(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_U_urestriction_coarsen_operator(),
      d_P_urestriction_coarsen_operator(),
      d_urestriction_coarsen_algorithm(),
      d_urestriction_coarsen_schedules(),
      d_U_rrestriction_coarsen_operator(),
      d_P_rrestriction_coarsen_operator(),
      d_rrestriction_coarsen_algorithm(),
      d_rrestriction_coarsen_schedules(),
      d_U_ghostfill_refine_operator(),
      d_P_ghostfill_refine_operator(),
      d_ghostfill_refine_algorithm(),
      d_ghostfill_refine_schedules(),
      d_U_ghostfill_nocoarse_refine_operator(),
      d_P_ghostfill_nocoarse_refine_operator(),
      d_ghostfill_nocoarse_refine_algorithm(),
      d_ghostfill_nocoarse_refine_schedules(),
      d_U_synch_refine_operator(),
      d_U_synch_refine_algorithm(),
      d_U_synch_refine_schedules()
{
    // Get values from the input database.
    if (!input_db.isNull())
    {
        d_smoother_choice = input_db->getStringWithDefault("smoother_choice", d_smoother_choice);
        d_U_prolongation_method = input_db->getStringWithDefault("U_prolongation_method", d_U_prolongation_method);
        d_P_prolongation_method = input_db->getStringWithDefault("P_prolongation_method", d_P_prolongation_method);
        d_U_restriction_method = input_db->getStringWithDefault("U_restriction_method", d_U_restriction_method);
        d_P_restriction_method = input_db->getStringWithDefault("P_restriction_method", d_P_restriction_method);
        d_fac_max_cycles = input_db->getIntegerWithDefault("fac_max_cycles", d_fac_max_cycles);
        d_fac_uses_presmoothing = input_db->getBoolWithDefault("fac_uses_presmoothing", d_fac_uses_presmoothing);
        d_fac_initial_guess_nonzero = input_db->getBoolWithDefault("fac_initial_guess_nonzero", d_fac_initial_guess_nonzero);
        d_skip_restrict_sol = input_db->getBoolWithDefault("skip_restrict_sol", d_skip_restrict_sol);
        d_skip_restrict_residual = input_db->getBoolWithDefault("skip_restrict_residual", d_skip_restrict_residual);
        d_coarse_solver_choice = input_db->getStringWithDefault("coarse_solver_choice", d_coarse_solver_choice);
        d_coarse_solver_tol = input_db->getDoubleWithDefault("coarse_solver_tolerance", d_coarse_solver_tol);
        d_coarse_solver_max_its = input_db->getIntegerWithDefault("coarse_solver_max_iterations", d_coarse_solver_max_its);
    }
    sanityCheck();

    // Create the hypre solver, if needed.
    setCoarsestLevelSolverChoice(d_coarse_solver_choice);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (int d = 0; d < NDIM; ++d)
    {
        d_default_U_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_U_bc_coef->setBoundaryValue(2*d+1,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d  ,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setHomogeneousBc(d_homogeneous_bc);
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM,d_default_U_bc_coef),d_default_P_bc_coef);

    // Setup scratch variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    const IntVector<NDIM> side_ghosts = d_gcw;
    Pointer<SideVariable<NDIM,double> > side_scratch_var =
        new SideVariable<NDIM,double>(d_object_name+"::side_scratch");
    if (var_db->checkVariableExists(side_scratch_var->getName()))
    {
        side_scratch_var = var_db->getVariable(side_scratch_var->getName());
        d_side_scratch_idx = var_db->mapVariableAndContextToIndex(side_scratch_var, d_context);
        var_db->removePatchDataIndex(d_side_scratch_idx);
    }
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);

    const IntVector<NDIM> cell_ghosts = d_gcw;
    Pointer<CellVariable<NDIM,double> > cell_scratch_var =
        new CellVariable<NDIM,double>(d_object_name+"::cell_scratch");
    if (var_db->checkVariableExists(cell_scratch_var->getName()))
    {
        cell_scratch_var = var_db->getVariable(cell_scratch_var->getName());
        d_cell_scratch_idx = var_db->mapVariableAndContextToIndex(cell_scratch_var, d_context);
        var_db->removePatchDataIndex(d_cell_scratch_idx);
    }
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_restrict_solution                   = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::restrictSolution()");
        t_restrict_residual                   = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::restrictResidual()");
        t_prolong_error_and_correct           = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::prolongErrorAndCorrect()");
        t_smooth_error                        = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::smoothError()");
        t_solve_coarsest_level                = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel()");
        t_compute_composite_residual_on_level = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::computeCompositeResidualOnLevel()");
        t_compute_residual_norm               = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::computeResidualNorm()");
        t_initialize_operator_state           = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::initializeOperatorState()");
        t_deallocate_operator_state           = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::deallocateOperatorState()");
        timers_need_init = false;
    }
    return;
}// INSStaggeredBoxRelaxationFACOperator

INSStaggeredBoxRelaxationFACOperator::~INSStaggeredBoxRelaxationFACOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    delete d_default_U_bc_coef;
    delete d_default_P_bc_coef;
    return;
}// ~INSStaggeredBoxRelaxationFACOperator

void
INSStaggeredBoxRelaxationFACOperator::setProblemCoefficients(
    const INSCoefs& problem_coefs,
    const double dt)
{
    d_problem_coefs = problem_coefs;
    d_dt = dt;
    return;
}// setProblemCoefficients

void
INSStaggeredBoxRelaxationFACOperator::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    d_U_bc_coefs.resize(U_bc_coefs.size());
    for (unsigned l = 0; l < U_bc_coefs.size(); ++l)
    {
        if (U_bc_coefs[l] != NULL)
        {
            d_U_bc_coefs[l] = U_bc_coefs[l];
        }
        else
        {
            d_U_bc_coefs[l] = d_default_U_bc_coef;
        }
    }

    if (P_bc_coef != NULL)
    {
        d_P_bc_coef = P_bc_coef;
    }
    else
    {
        d_P_bc_coef = d_default_P_bc_coef;
    }
    return;
}// setPhysicalBcCoefs

void
INSStaggeredBoxRelaxationFACOperator::setHomogeneousBc(
    const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}// setHomogeneousBc

void
INSStaggeredBoxRelaxationFACOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

void
INSStaggeredBoxRelaxationFACOperator::setPreconditioner(
    const FACPreconditioner<NDIM>* preconditioner)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setPreconditioner()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_preconditioner = preconditioner;
    sanityCheck();
    return;
}// setPreconditioner

void
INSStaggeredBoxRelaxationFACOperator::setResetLevels(
    const int coarsest_ln,
    const int finest_ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT((coarsest_ln == -1 && finest_ln == -1) ||
                (coarsest_ln >=  0 && finest_ln >= coarsest_ln));
#endif
    if (d_is_initialized)
    {
        d_coarsest_reset_ln = coarsest_ln;
        d_finest_reset_ln = finest_ln;
    }
    return;
}// setResetLevels

void
INSStaggeredBoxRelaxationFACOperator::setGhostCellWidth(
    const IntVector<NDIM>& ghost_cell_width)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setGhostCellWidth()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    if (ghost_cell_width.min() == 0)
    {
        TBOX_ERROR(d_object_name << "::setGhostCellWidth()\n"
                   << "  ghost_cell_width.min() must be greater than zero" << std::endl);
    }
    d_gcw = ghost_cell_width;
    sanityCheck();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    Pointer<SideVariable<NDIM,double> > side_scratch_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!side_scratch_var.isNull());
#endif
    var_db->removePatchDataIndex(d_side_scratch_idx);
    const IntVector<NDIM> side_ghosts = d_gcw;
    d_side_scratch_idx = var_db->registerVariableAndContext(side_scratch_var, d_context, side_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<SideDataFactory<NDIM,double> > side_scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_side_scratch_idx);
    TBOX_ASSERT(!side_scratch_pdat_fac.isNull());
    TBOX_ASSERT(side_scratch_pdat_fac->getGhostCellWidth() == d_gcw);
#endif

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    Pointer<CellVariable<NDIM,double> > cell_scratch_var = var;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!cell_scratch_var.isNull());
#endif
    var_db->removePatchDataIndex(d_cell_scratch_idx);
    const IntVector<NDIM> cell_ghosts = d_gcw;
    d_cell_scratch_idx = var_db->registerVariableAndContext(cell_scratch_var, d_context, cell_ghosts);
#ifdef DEBUG_CHECK_ASSERTIONS
    Pointer<CellDataFactory<NDIM,double> > cell_scratch_pdat_fac = var_db->getPatchDescriptor()->getPatchDataFactory(d_cell_scratch_idx);
    TBOX_ASSERT(!cell_scratch_pdat_fac.isNull());
    TBOX_ASSERT(cell_scratch_pdat_fac->getGhostCellWidth() == d_gcw);
#endif
    return;
}// setGhostCellWidth

void
INSStaggeredBoxRelaxationFACOperator::setSmootherChoice(
    const std::string& smoother_choice)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSmootherChoice()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_smoother_choice = smoother_choice;
    sanityCheck();
    return;
}// setSmootherChoice

void
INSStaggeredBoxRelaxationFACOperator::setCoarsestLevelSolverChoice(
    const std::string& coarse_solver_choice)
{
    d_coarse_solver_choice = coarse_solver_choice;
    sanityCheck();
    return;
}// setCoarsestLevelSolverChoice

void
INSStaggeredBoxRelaxationFACOperator::setCoarsestLevelSolverTolerance(
    double coarse_solver_tol)
{
    d_coarse_solver_tol = coarse_solver_tol;
    sanityCheck();
    return;
}// setCoarsestLevelSolverTolerance

void
INSStaggeredBoxRelaxationFACOperator::setCoarsestLevelSolverMaxIterations(
    int coarse_solver_max_its)
{
    d_coarse_solver_max_its = coarse_solver_max_its;
    sanityCheck();
    return;
}// setCoarsestLevelSolverMaxIterations

void
INSStaggeredBoxRelaxationFACOperator::setProlongationMethods(
    const std::string& U_prolongation_method,
    const std::string& P_prolongation_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setProlongationMethods()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_U_prolongation_method = U_prolongation_method;
    d_P_prolongation_method = P_prolongation_method;
    sanityCheck();
    return;
}// setProlongationMethods

void
INSStaggeredBoxRelaxationFACOperator::setRestrictionMethods(
    const std::string& U_restriction_method,
    const std::string& P_restriction_method)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setRestrictionMethods()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_U_restriction_method = U_restriction_method;
    d_P_restriction_method = P_restriction_method;
    sanityCheck();
    return;
}// setRestrictionMethods

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditionerMaxCycles(
    int fac_max_cycles)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditionerMaxCycles()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_fac_max_cycles = fac_max_cycles;
    sanityCheck();
    return;
}// setFACPreconditionerMaxCycles

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditionerUsesPresmoothing(
    bool fac_uses_presmoothing)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditionerUsesPresmoothing()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_fac_uses_presmoothing = fac_uses_presmoothing;
    sanityCheck();
    return;
}// setFACPreconditionerUsesPresmoothing

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditionerInitialGuessNonzero(
    bool fac_initial_guess_nonzero)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditionerInitialGuessNonzero()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_fac_initial_guess_nonzero = fac_initial_guess_nonzero;
    sanityCheck();
    return;
}// setFACPreconditionerInitialGuessNonzero

void
INSStaggeredBoxRelaxationFACOperator::setSkipRestrictSolution(
    bool skip_restrict_sol)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSkipRestrictSolution()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_skip_restrict_sol = skip_restrict_sol;
    sanityCheck();
    return;
}// setSkipRestrictSolution

void
INSStaggeredBoxRelaxationFACOperator::setSkipRestrictResidual(
    bool skip_restrict_residual)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setSkipRestrictResidual()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_skip_restrict_residual = skip_restrict_residual;
    sanityCheck();
    return;
}// setSkipRestrictResidual

///
///  The following routines:
///
///      restrictSolution(),
///      restrictResidual(),
///      prolongErrorAndCorrect(),
///      smoothError(),
///      solveCoarsestLevel(),
///      computeCompositeResidualOnLevel(),
///      computeResidualNorm(),
///      initializeOperatorState(),
///      deallocateOperatorState()
///
///  are concrete implementations of functions declared in the
///  FACOperatorStrategy abstract base class.
///

void
INSStaggeredBoxRelaxationFACOperator::restrictSolution(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    if (d_skip_restrict_sol) return;

    t_restrict_solution->start();

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    xeqScheduleURestriction(dst_idxs, src_idxs, dst_ln);

    static const bool homogeneous_bc = true;
    if (dst_ln == d_coarsest_ln)
    {
        xeqScheduleGhostFillNoCoarse(dst_idxs, dst_ln, homogeneous_bc);
    }
    else
    {
        xeqScheduleGhostFill(dst_idxs, dst_ln, homogeneous_bc);
    }

    t_restrict_solution->stop();
    return;
}// restrictSolution

void
INSStaggeredBoxRelaxationFACOperator::restrictResidual(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    if (d_skip_restrict_residual) return;

    t_restrict_residual->start();

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    xeqScheduleRRestriction(dst_idxs, src_idxs, dst_ln);

    t_restrict_residual->stop();
    return;
}// restrictResidual

void
INSStaggeredBoxRelaxationFACOperator::prolongErrorAndCorrect(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    t_prolong_error_and_correct->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    const bool dst_data_is_prolonged_src_data = !d_fac_uses_presmoothing && !d_fac_initial_guess_nonzero && d_fac_max_cycles == 1;
    TBOX_ASSERT(dst_data_is_prolonged_src_data);
#endif

    // Refine the correction from the coarse level src data directly into the
    // fine level error.
    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    static const bool homogeneous_bc = true;
    xeqScheduleProlongation(dst_idxs, src_idxs, dst_ln, homogeneous_bc);

    t_prolong_error_and_correct->stop();
    return;
}// prolongErrorAndCorrect

void
INSStaggeredBoxRelaxationFACOperator::smoothError(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int level_num,
    int num_sweeps)
{
    if (num_sweeps == 0) return;

    t_smooth_error->start();

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln)
    {
        HierarchySideDataOpsReal<NDIM,double> hierarchy_sc_data_ops(d_hierarchy, level_num, level_num);
        hierarchy_sc_data_ops.copyData(d_side_scratch_idx, U_error_idx, false);

        HierarchyCellDataOpsReal<NDIM,double> hierarchy_cc_data_ops(d_hierarchy, level_num, level_num);
        hierarchy_cc_data_ops.copyData(d_cell_scratch_idx, P_error_idx, false);
    }

    // Smooth the error by the specified number of sweeps.
    for (int isweep = 0; isweep < num_sweeps; ++isweep)
    {
        // Re-fill ghost cell data as needed.
        if (level_num > d_coarsest_ln)
        {
            if (isweep > 0)
            {
                // Copy the coarse-fine interface ghost cell values which are
                // cached in the scratch data into the error data.
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const int patch_num = patch->getPatchNumber();

                    Pointer<SideData<NDIM,double> >   U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM,double> > U_scratch_data = patch->getPatchData(d_side_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
                    TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
                    TBOX_ASSERT(  U_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis).copy(
                            U_scratch_data->getArrayData(axis),
                            d_patch_side_bc_box_overlap[level_num][patch_num][axis],
                            IntVector<NDIM>(0));
                    }

                    Pointer<CellData<NDIM,double> >   P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellData<NDIM,double> > P_scratch_data = patch->getPatchData(d_cell_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
                    TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
                    TBOX_ASSERT(  P_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    P_error_data->getArrayData().copy(
                        P_scratch_data->getArrayData(),
                        d_patch_cell_bc_box_overlap[level_num][patch_num],
                        IntVector<NDIM>(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                const std::pair<int,int> error_idxs = std::make_pair(U_error_idx,P_error_idx);
                static const bool homogeneous_bc = true;
                xeqScheduleGhostFillNoCoarse(error_idxs, level_num, homogeneous_bc);
            }

            // Complete the coarse-fine interface interpolation by computing the
            // normal extension.
            d_U_cf_bdry_op->setPatchDataIndex(U_error_idx);
            d_P_cf_bdry_op->setPatchDataIndex(P_error_idx);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const IntVector<NDIM>& ghost_width_to_fill = d_gcw;
                d_U_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
                d_P_cf_bdry_op->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
            }
        }
        else if (isweep > 0)
        {
            const std::pair<int,int> error_idxs = std::make_pair(U_error_idx,P_error_idx);
            static const bool homogeneous_bc = true;
            xeqScheduleGhostFillNoCoarse(error_idxs, level_num, homogeneous_bc);
        }

        // Smooth the error on the patches.
        Vec& e = d_box_e[level_num];
        Vec& r = d_box_r[level_num];
        KSP& ksp = d_box_ksp[level_num];
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const int patch_num = patch->getPatchNumber();

            // Reset ghost cell values in the copy of the residual data so that
            // patch boundary conditions are properly handled.
            Pointer<SideData<NDIM,double> >    U_error_data = error   .getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > U_residual_data = residual.getComponentPatchData(0, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_residual_data->getGhostBox());
            TBOX_ASSERT(   U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_residual_data->getGhostCellWidth() == d_gcw);
#endif
            for (int axis = 0; axis < NDIM; ++axis)
            {
                U_residual_data->getArrayData(axis).copy(
                    U_error_data->getArrayData(axis),
                    d_patch_side_bc_box_overlap[level_num][patch_num][axis],
                    IntVector<NDIM>(0));
            }
            Pointer<CellData<NDIM,double> >    P_error_data = error   .getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM,double> > P_residual_data = residual.getComponentPatchData(1, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_residual_data->getGhostBox());
            TBOX_ASSERT(   P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_residual_data->getGhostCellWidth() == d_gcw);
#endif
            P_residual_data->getArrayData().copy(
                P_error_data->getArrayData(),
                d_patch_cell_bc_box_overlap[level_num][patch_num],
                IntVector<NDIM>(0));

            // Smooth the error on the patch.
            const Box<NDIM>& patch_box = patch->getBox();
#if 1
            for (int m = 0; m <= 1; ++m)
            {
                for (Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    if ((i(0) + i(1)) % 2 == m)
                    {
                        const Box<NDIM> box(i,i);
                        const Box<NDIM> ghost_box = Box<NDIM>::grow(box,BOX_GHOSTS);

                        copyToVec(e, *U_error_data   , *P_error_data   , ghost_box, ghost_box);
                        copyToVec(r, *U_residual_data, *P_residual_data, ghost_box, ghost_box);
                        int ierr = KSPSolve(ksp, r, e);  IBTK_CHKERRQ(ierr);
                        copyFromVec(e, *U_error_data, *P_error_data, box, ghost_box);
                    }
                }
            }
#else
            int ierr;

            const Box<NDIM> ghost_box = Box<NDIM>::grow(patch_box,d_gcw);
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            Mat A;
            buildBoxOperator(A, d_problem_coefs, d_dt, patch_box, ghost_box, dx);

            Vec e, f;
            ierr = MatGetVecs(A, &e, &f);  IBTK_CHKERRQ(ierr);
            KSP ksp;
            ierr = KSPCreate(PETSC_COMM_SELF, &ksp);  IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, A, A, SAME_PRECONDITIONER);  IBTK_CHKERRQ(ierr);
            ierr = KSPSetOptionsPrefix(ksp, "vanka_");  IBTK_CHKERRQ(ierr);
            ierr = KSPSetFromOptions(ksp);  IBTK_CHKERRQ(ierr);
            copyToVec(e, *U_error_data   , *P_error_data   , ghost_box, ghost_box);
            copyToVec(f, *U_residual_data, *P_residual_data, ghost_box, ghost_box);
            ierr = KSPSolve(ksp, f, e);
            copyFromVec(e, *U_error_data, *P_error_data, patch_box, ghost_box);

            MatDestroy(A);
            VecDestroy(e);
            VecDestroy(f);
            KSPDestroy(ksp);
#endif
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleSideDataSynch(U_error_idx, level_num);

    t_smooth_error->stop();
    return;
}// smoothError

int
INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int coarsest_ln)
{
    t_solve_coarsest_level->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif

    smoothError(error, residual, coarsest_ln, d_coarse_solver_max_its);

    // Re-fill ghost values.
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const std::pair<int,int> error_idxs = std::make_pair(U_error_idx,P_error_idx);

    static const bool homogeneous_bc = true;
    xeqScheduleGhostFillNoCoarse(error_idxs, coarsest_ln, homogeneous_bc);

    t_solve_coarsest_level->stop();
    static const int ret_val = 1;
    return ret_val;
}// solveCoarsestLevel

void
INSStaggeredBoxRelaxationFACOperator::computeCompositeResidualOnLevel(
    SAMRAIVectorReal<NDIM,double>& residual,
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs,
    int level_num,
    bool error_equation_indicator)
{
    t_compute_composite_residual_on_level->start();

    if (!d_fac_uses_presmoothing && !d_fac_initial_guess_nonzero && d_fac_max_cycles == 1)
    {
        // The residual needs to be computed in two different cases:
        //
        // - before each FAC sweep commences, and
        // - after performing any presmoothing sweeps
        //
        // If the FAC preconditioner (a) does not use presmoothing, (b) uses a
        // zero initial guess, and (c) only employs one FAC sweep (only one FAC
        // sweep is needed, for instance, when the preconditioner is being used
        // in conjunction with a Krylov subspace method), then we can simply set
        // the residual equal to the right hand side.
        static const bool interior_only = false;

        HierarchySideDataOpsReal<NDIM,double> hierarchy_sc_data_ops(d_hierarchy, level_num, level_num);
        const int U_dst_idx = residual.getComponentDescriptorIndex(0);
        const int U_src_idx = rhs.getComponentDescriptorIndex(0);
        if (U_dst_idx != U_src_idx)
        {
            hierarchy_sc_data_ops.copyData(U_dst_idx, U_src_idx, interior_only);
        }

        HierarchyCellDataOpsReal<NDIM,double> hierarchy_cc_data_ops(d_hierarchy, level_num, level_num);
        const int P_dst_idx = residual.getComponentDescriptorIndex(1);
        const int P_src_idx = rhs.getComponentDescriptorIndex(1);
        if (P_dst_idx != P_src_idx)
        {
            hierarchy_cc_data_ops.copyData(P_dst_idx, P_src_idx, interior_only);
        }
    }
    else
    {
        TBOX_ERROR("INSStaggeredBoxRelaxationFACOperator::computeResidualOnLevel()\n"
                   << "  current implementation cannot compute residuals,\n"
                   << "  consequently, we require that the FAC algorithm:\n"
                   << "     (1) does not use pre-smoothing,\n"
                   << "     (2) does not use a non-zero initial guess, and\n"
                   << "     (3) only uses one cycle\n"
                   << "  thus the implemented solver is really only suitable for use as a preconditioner." << std::endl);
    }

    t_compute_composite_residual_on_level->stop();
    return;
}// computeCompositeResidualOnLevel

double
INSStaggeredBoxRelaxationFACOperator::computeResidualNorm(
    const SAMRAIVectorReal<NDIM,double>& residual,
    int fine_ln,
    int coarse_ln)
{
    t_compute_residual_norm->start();

    if (coarse_ln != residual.getCoarsestLevelNumber() ||
        fine_ln   != residual.getFinestLevelNumber())
    {
        TBOX_ERROR("INSStaggeredBoxRelaxationFACOperator::computeResidualNorm()\n"
                   << "  residual can only be computed over the range of levels\n"
                   << "  defined in the SAMRAIVectorReal residual vector" << std::endl);
    }

    t_compute_residual_norm->stop();
    return NormOps::L2Norm(&residual);
}// computeResidualNorm

void
INSStaggeredBoxRelaxationFACOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs)
{
    t_initialize_operator_state->start();
    d_in_initialize_operator_state = true;

    // Cache the level range to be reset.
    //
    // NOTE: We cannot use d_coarsest_reset_ln and d_finest_reset_ln since those
    // values are reset by deallocateOperatorState().
    const int coarsest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_coarsest_reset_ln
         : solution.getCoarsestLevelNumber());
    const int finest_reset_ln =
        (d_coarsest_reset_ln != -1 && d_finest_reset_ln != -1
         ? d_finest_reset_ln
         : solution.getFinestLevelNumber());

    // Deallocate the solver state if the solver is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_solution = solution.cloneVector(solution.getName());
    d_solution->allocateVectorData();

    d_rhs = rhs.cloneVector(rhs.getName());
    d_rhs->allocateVectorData();

    // Reset the hierarchy configuration.
    d_hierarchy   = solution.getPatchHierarchy();
    d_coarsest_ln = solution.getCoarsestLevelNumber();
    d_finest_ln   = solution.getFinestLevelNumber();

    // Allocate scratch data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_side_scratch_idx)) level->allocatePatchData(d_side_scratch_idx);
        if (!level->checkAllocated(d_cell_scratch_idx)) level->allocatePatchData(d_cell_scratch_idx);
    }

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    static bool need_to_add_cubic_coarsen = true;
    if (need_to_add_cubic_coarsen)
    {
        geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
        geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
        need_to_add_cubic_coarsen = false;
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > var;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_U_prolongation_method);

    d_U_cf_bdry_op = new CartSideDoubleQuadraticCFInterpolation();
    d_U_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_U_cf_bdry_op->setPatchDataIndex(d_side_scratch_idx);
    d_U_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_prolongation_refine_operator = geometry->lookupRefineOperator(var, d_P_prolongation_method);

    d_P_cf_bdry_op = new CartCellDoubleQuadraticCFInterpolation();
    d_P_cf_bdry_op->setConsistentInterpolationScheme(false);
    d_P_cf_bdry_op->setPatchDataIndex(d_cell_scratch_idx);
    d_P_cf_bdry_op->setPatchHierarchy(d_hierarchy);

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_urestriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_U_restriction_method);
    d_U_rrestriction_coarsen_operator = d_U_urestriction_coarsen_operator;

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_urestriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_P_restriction_method);
    d_P_rrestriction_coarsen_operator = d_P_urestriction_coarsen_operator;

    var_db->mapIndexToVariable(d_side_scratch_idx, var);
    d_U_ghostfill_refine_operator = geometry->lookupRefineOperator(var, d_U_prolongation_method);
    d_U_ghostfill_nocoarse_refine_operator = NULL;

    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_ghostfill_refine_operator = geometry->lookupRefineOperator(var, d_P_prolongation_method);
    d_P_ghostfill_nocoarse_refine_operator = NULL;

    d_U_synch_refine_operator = NULL;

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    d_U_bc_op = new CartSideRobinPhysBdryOp(d_side_scratch_idx, d_U_bc_coefs, false);
    d_P_bc_op = new CartCellRobinPhysBdryOp(d_cell_scratch_idx, d_P_bc_coef , false);

#if (NDIM == 2)
    d_U_op_stencil_fill_pattern = new SideNoCornersFillPattern(GHOSTS);
    d_P_op_stencil_fill_pattern = new CellNoCornersFillPattern(GHOSTS);
#endif
#if (NDIM == 3)
    d_U_op_stencil_fill_pattern = new SideNoCornersFillPattern(GHOSTS,false);
    d_P_op_stencil_fill_pattern = new CellNoCornersFillPattern(GHOSTS,false);
#endif
    d_U_synch_fill_pattern = new SideSynchCopyFillPattern();

    std::vector<RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_U_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_P_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_U_bc_op);
    prolongation_refine_patch_strategies.push_back(d_P_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln+1);
    d_urestriction_coarsen_schedules.resize(d_finest_ln+1);
    d_rrestriction_coarsen_schedules.resize(d_finest_ln+1);
    d_ghostfill_refine_schedules.resize(d_finest_ln+1);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln+1);
    d_U_synch_refine_schedules.resize(d_finest_ln+1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_urestriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_rrestriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    d_ghostfill_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_ghostfill_nocoarse_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_U_synch_refine_algorithm = new RefineAlgorithm<NDIM>();

    d_prolongation_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        solution.getComponentDescriptorIndex(0),
        d_side_scratch_idx,
        d_U_prolongation_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_prolongation_refine_algorithm->registerRefine(
        d_cell_scratch_idx,
        solution.getComponentDescriptorIndex(1),
        d_cell_scratch_idx,
        d_P_prolongation_refine_operator,
        d_P_op_stencil_fill_pattern);

    d_urestriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_urestriction_coarsen_operator);
    d_urestriction_coarsen_algorithm->registerCoarsen(
        d_cell_scratch_idx,
        d_cell_scratch_idx,
        d_P_urestriction_coarsen_operator);

    d_rrestriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_rrestriction_coarsen_operator);
    d_rrestriction_coarsen_algorithm->registerCoarsen(
        d_cell_scratch_idx,
        d_cell_scratch_idx,
        d_P_rrestriction_coarsen_operator);

    d_ghostfill_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_ghostfill_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_ghostfill_refine_algorithm->registerRefine(
        d_cell_scratch_idx,
        d_cell_scratch_idx,
        d_cell_scratch_idx,
        d_P_ghostfill_refine_operator,
        d_P_op_stencil_fill_pattern);

    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_side_scratch_idx,
        d_U_ghostfill_nocoarse_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        d_cell_scratch_idx,
        d_cell_scratch_idx,
        d_cell_scratch_idx,
        d_P_ghostfill_nocoarse_refine_operator,
        d_P_op_stencil_fill_pattern);

    d_U_synch_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        d_U_synch_refine_operator,
        d_U_synch_fill_pattern);

    std::vector<RefinePatchStrategy<NDIM>*> bc_op_ptrs(2);
    bc_op_ptrs[0] = d_U_bc_op.getPointer();
    bc_op_ptrs[1] = d_P_bc_op.getPointer();
    d_U_P_bc_op = new RefinePatchStrategySet(bc_op_ptrs.begin(), bc_op_ptrs.end(), false);

    for (int dst_ln = d_coarsest_ln+1; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_prolongation_refine_schedules[dst_ln] =
            d_prolongation_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln),
                Pointer<PatchLevel<NDIM> >(),
                dst_ln-1, d_hierarchy, d_prolongation_refine_patch_strategy.getPointer());

        d_ghostfill_refine_schedules[dst_ln] =
            d_ghostfill_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln),
                dst_ln-1, d_hierarchy, d_U_P_bc_op);

        d_ghostfill_nocoarse_refine_schedules[dst_ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), d_U_P_bc_op);

        d_U_synch_refine_schedules[dst_ln] =
            d_U_synch_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), NULL, new SideSynchCopyTransactionFactory());
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), d_U_P_bc_op);

    d_U_synch_refine_schedules[d_coarsest_ln] =
        d_U_synch_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), NULL, new SideSynchCopyTransactionFactory());

    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        if (!d_skip_restrict_sol)
        {
            d_urestriction_coarsen_schedules[dst_ln] =
                d_urestriction_coarsen_algorithm->createSchedule(
                    d_hierarchy->getPatchLevel(dst_ln  ),
                    d_hierarchy->getPatchLevel(dst_ln+1));
        }
        if (!d_skip_restrict_residual)
        {
            d_rrestriction_coarsen_schedules[dst_ln] =
                d_rrestriction_coarsen_algorithm->createSchedule(
                    d_hierarchy->getPatchLevel(dst_ln  ),
                    d_hierarchy->getPatchLevel(dst_ln+1));
        }
    }

    // Initialize the box relaxation data on each level of the patch hierarchy.
    d_box_op.resize(d_finest_ln+1);
    d_box_e.resize(d_finest_ln+1);
    d_box_r.resize(d_finest_ln+1);
    d_box_ksp.resize(d_finest_ln+1);

    const Box<NDIM> box(Index<NDIM>(0),Index<NDIM>(0));
    const Box<NDIM> ghost_box = Box<NDIM>::grow(box,BOX_GHOSTS);

    const double* const dx_coarsest = geometry->getDx();
    double dx[NDIM];

    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        const IntVector<NDIM>& ratio = d_hierarchy->getPatchLevel(ln)->getRatio();
        for (int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d]/double(ratio(d));
        }
        buildBoxOperator(d_box_op[ln], d_problem_coefs, d_dt, box, ghost_box, dx);
        int ierr;
        ierr = MatGetVecs(d_box_op[ln], &d_box_e[ln], &d_box_r[ln]);  IBTK_CHKERRQ(ierr);
        ierr = KSPCreate(PETSC_COMM_SELF, &d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
        ierr = KSPSetOperators(d_box_ksp[ln], d_box_op[ln], d_box_op[ln], SAME_PRECONDITIONER);  IBTK_CHKERRQ(ierr);
        ierr = KSPSetType(d_box_ksp[ln], KSPPREONLY);  IBTK_CHKERRQ(ierr);
        PC box_pc;
        ierr = KSPGetPC(d_box_ksp[ln], &box_pc);  IBTK_CHKERRQ(ierr);
        ierr = PCSetType(box_pc, PCLU);  IBTK_CHKERRQ(ierr);
        ierr = PCFactorReorderForNonzeroDiagonal(box_pc, std::numeric_limits<double>::epsilon());  IBTK_CHKERRQ(ierr);
//      ierr = PCFactorSetMatOrderingType(box_pc, MATORDERING_ND);  IBTK_CHKERRQ(ierr);
        ierr = KSPSetUp(d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(d_finest_ln+1);
    d_patch_cell_bc_box_overlap.resize(d_finest_ln+1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const BoxArray<NDIM>& box_array = level->getBoxes();
        const Array<int>& local_ids = level->getProcessorMapping().getLocalIndices();
        for (int k = 0; k < local_ids.size(); ++k)
        {
            const int local_id = local_ids[k];
            d_patch_side_bc_box_overlap[ln][local_id].resize(NDIM);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box_array[local_id],axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][local_id][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][local_id][axis].removeIntersections(side_box);
            }

            const Box<NDIM>& cell_box = box_array[local_id];
            const Box<NDIM> cell_ghost_box = Box<NDIM>::grow(cell_box, 1);
            d_patch_cell_bc_box_overlap[ln][local_id] = BoxList<NDIM>(cell_ghost_box);
            d_patch_cell_bc_box_overlap[ln][local_id].removeIntersections(cell_box);
        }
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;

    d_in_initialize_operator_state = false;
    t_initialize_operator_state->stop();
    return;
}// initializeOperatorState

void
INSStaggeredBoxRelaxationFACOperator::deallocateOperatorState()
{
    if (d_is_initialized && !d_in_initialize_operator_state &&
        (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
    {
        return;
    }

    t_deallocate_operator_state->start();

    if (d_is_initialized)
    {
        const int coarsest_reset_ln =
            (d_in_initialize_operator_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_coarsest_reset_ln : d_coarsest_ln;
        const int finest_reset_ln =
            (d_in_initialize_operator_state &&
             (d_coarsest_reset_ln != -1) && (d_finest_reset_ln != -1))
            ? d_finest_reset_ln : d_finest_ln;

        for (int ln = coarsest_reset_ln; ln <= std::min(d_finest_ln, finest_reset_ln); ++ln)
        {
            if (ln <= d_hierarchy->getFinestLevelNumber())
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(d_side_scratch_idx)) level->deallocatePatchData(d_side_scratch_idx);
            }

            int ierr;
            ierr = MatDestroy(d_box_op [ln]);  IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(d_box_e  [ln]);  IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(d_box_r  [ln]);  IBTK_CHKERRQ(ierr);
            ierr = KSPDestroy(d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
        }

        // Delete the solution and rhs vectors.
        d_solution->resetLevels(d_solution->getCoarsestLevelNumber(), std::min(d_solution->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
        d_solution->freeVectorComponents();
        d_solution.setNull();

        d_rhs->resetLevels(d_rhs->getCoarsestLevelNumber(), std::min(d_rhs->getFinestLevelNumber(),d_hierarchy->getFinestLevelNumber()));
        d_rhs->freeVectorComponents();
        d_rhs.setNull();

        if (!d_in_initialize_operator_state ||
            (d_coarsest_reset_ln == -1) || (d_finest_reset_ln == -1))
        {
            d_patch_side_bc_box_overlap.resize(0);
            d_patch_cell_bc_box_overlap.resize(0);

            d_hierarchy.setNull();
            d_coarsest_ln = -1;
            d_finest_ln   = -1;

            d_U_prolongation_refine_operator    .setNull();
            d_U_cf_bdry_op                      .setNull();
            d_prolongation_refine_patch_strategy.setNull();
            d_prolongation_refine_algorithm     .setNull();
            d_prolongation_refine_schedules     .resize(0);

            d_U_urestriction_coarsen_operator.setNull();
            d_urestriction_coarsen_algorithm .setNull();
            d_urestriction_coarsen_schedules .resize(0);

            d_U_rrestriction_coarsen_operator.setNull();
            d_rrestriction_coarsen_algorithm .setNull();
            d_rrestriction_coarsen_schedules .resize(0);

            d_U_ghostfill_refine_operator.setNull();
            d_ghostfill_refine_algorithm .setNull();
            d_ghostfill_refine_schedules .resize(0);

            d_U_ghostfill_nocoarse_refine_operator.setNull();
            d_ghostfill_nocoarse_refine_algorithm .setNull();
            d_ghostfill_nocoarse_refine_schedules .resize(0);

            d_U_synch_refine_operator .setNull();
            d_U_synch_refine_algorithm.setNull();
            d_U_synch_refine_schedules.resize(0);
        }

        delete d_U_P_bc_op;
    }

    // Clear the "reset level" range.
    d_coarsest_reset_ln = -1;
    d_finest_reset_ln   = -1;

    // Indicate that the operator is not initialized.
    d_is_initialized = false;

    t_deallocate_operator_state->stop();
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleProlongation(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln,
    const bool homogeneous_bc)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(homogeneous_bc);
    d_U_cf_bdry_op->setPatchDataIndex(U_dst_idx);

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setPhysicalBcCoef(d_P_bc_coef);
    d_P_bc_op->setHomogeneousBc(homogeneous_bc);
    d_P_cf_bdry_op->setPatchDataIndex(P_dst_idx);

    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(U_dst_idx, U_src_idx, U_dst_idx, d_U_prolongation_refine_operator, d_U_op_stencil_fill_pattern);
    refiner.registerRefine(P_dst_idx, P_src_idx, P_dst_idx, d_P_prolongation_refine_operator, d_P_op_stencil_fill_pattern);
    refiner.resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    d_prolongation_refine_schedules[dst_ln]->fillData(d_new_time);
    d_prolongation_refine_algorithm->resetSchedule(d_prolongation_refine_schedules[dst_ln]);
    return;
}// xeqScheduleProlongation

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleURestriction(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;

    CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(U_dst_idx, U_src_idx, d_U_urestriction_coarsen_operator);
    coarsener.registerCoarsen(P_dst_idx, P_src_idx, d_P_urestriction_coarsen_operator);
    coarsener.resetSchedule(d_urestriction_coarsen_schedules[dst_ln]);
    d_urestriction_coarsen_schedules[dst_ln]->coarsenData();
    d_urestriction_coarsen_algorithm->resetSchedule(d_urestriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleURestriction

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleRRestriction(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;

    CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(U_dst_idx, U_src_idx, d_U_rrestriction_coarsen_operator);
    coarsener.registerCoarsen(P_dst_idx, P_src_idx, d_P_rrestriction_coarsen_operator);
    coarsener.resetSchedule(d_rrestriction_coarsen_schedules[dst_ln]);
    d_rrestriction_coarsen_schedules[dst_ln]->coarsenData();
    d_rrestriction_coarsen_algorithm->resetSchedule(d_rrestriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleRRestriction

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleGhostFill(
    const std::pair<int,int>& dst_idxs,
    const int dst_ln,
    const bool homogeneous_bc)
{
    const int U_dst_idx = dst_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(homogeneous_bc);

    const int P_dst_idx = dst_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setPhysicalBcCoef(d_P_bc_coef);
    d_P_bc_op->setHomogeneousBc(homogeneous_bc);

    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, d_U_ghostfill_refine_operator, d_U_op_stencil_fill_pattern);
    refiner.registerRefine(P_dst_idx, P_dst_idx, P_dst_idx, d_P_ghostfill_refine_operator, d_P_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_refine_schedules[dst_ln]);
    d_ghostfill_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_refine_algorithm->resetSchedule(d_ghostfill_refine_schedules[dst_ln]);
    return;
}// xeqScheduleGhostFill

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleGhostFillNoCoarse(
    const std::pair<int,int>& dst_idxs,
    const int dst_ln,
    const bool homogeneous_bc)
{
    const int U_dst_idx = dst_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(homogeneous_bc);

    const int P_dst_idx = dst_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setPhysicalBcCoef(d_P_bc_coef);
    d_P_bc_op->setHomogeneousBc(homogeneous_bc);

    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, d_U_ghostfill_nocoarse_refine_operator, d_U_op_stencil_fill_pattern);
    refiner.registerRefine(P_dst_idx, P_dst_idx, P_dst_idx, d_P_ghostfill_nocoarse_refine_operator, d_P_op_stencil_fill_pattern);
    refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    d_ghostfill_nocoarse_refine_schedules[dst_ln]->fillData(d_new_time);
    d_ghostfill_nocoarse_refine_algorithm->resetSchedule(d_ghostfill_nocoarse_refine_schedules[dst_ln]);
    return;
}// xeqScheduleGhostFillNoCoarse

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleSideDataSynch(
    const int U_dst_idx,
    const int dst_ln)
{
    RefineAlgorithm<NDIM> refiner;
    refiner.registerRefine(U_dst_idx, U_dst_idx, U_dst_idx, d_U_synch_refine_operator, d_U_synch_fill_pattern);
    refiner.resetSchedule(d_U_synch_refine_schedules[dst_ln]);
    d_U_synch_refine_schedules[dst_ln]->fillData(d_new_time);
    d_U_synch_refine_algorithm->resetSchedule(d_U_synch_refine_schedules[dst_ln]);
    return;
}// xeqScheduleSideDataSynch

void
INSStaggeredBoxRelaxationFACOperator::sanityCheck()
{
    if (d_gcw.min() <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  ghost_cell_width.min() must be greater than zero" << std::endl);
    }

    if (d_smoother_choice != "additive")
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  unknown smoother type: " << d_smoother_choice << "\n"
                   << "  valid choices are: additive" << std::endl);
    }

    if (d_coarse_solver_choice != "block_jacobi" &&
        d_coarse_solver_choice != "hypre")
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  unknown coarse solver type: " << d_coarse_solver_choice << "\n"
                   << "  valid choices are: block_jacobi, hypre" << std::endl);
    }

    if (d_coarse_solver_tol < 0.0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid coarse solver tolerance: " << d_coarse_solver_tol << std::endl);
    }

    if (d_coarse_solver_max_its <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid coarse solver maximum iterations: " << d_coarse_solver_max_its << std::endl);
    }

    for (unsigned l = 0; l < d_U_bc_coefs.size(); ++l)
    {
        if (d_U_bc_coefs[l] == NULL)
        {
            TBOX_ERROR(d_object_name << ":\n"
                       << "  invalid velocity physical bc object at depth = " << l << std::endl);
        }
    }

    if (d_P_bc_coef == NULL)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid pressure physical bc object" << std::endl);
    }

    if (d_fac_max_cycles <= 0)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid value for fac_max_cycles: " << d_fac_max_cycles << std::endl);
    }
    return;
}// sanityCheck

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredBoxRelaxationFACOperator>;

//////////////////////////////////////////////////////////////////////////////
