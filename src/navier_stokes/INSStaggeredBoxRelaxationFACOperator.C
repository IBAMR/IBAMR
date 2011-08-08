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
#include <ibamr/ibamr_utilities.h>
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
static Timer* t_restrict_residual;
static Timer* t_prolong_error;
static Timer* t_prolong_error_and_correct;
static Timer* t_smooth_error;
static Timer* t_solve_coarsest_level;
static Timer* t_compute_residual;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;

// Number of ghosts cells used for each variable quantity.
static const int GHOSTS = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Number of ghost cells used for box operator.
static const int BOX_GHOSTS = 0;

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values; used only to evaluate composite grid
// residuals.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries; used only to evaluate
// composite grid residuals.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells; used only to evaluate composite grid residuals.
static const bool CONSISTENT_TYPE_2_BDRY = false;

inline int
compute_side_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box,
    const unsigned int axis)
{
    const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box,axis);
    if (!side_box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int d = 0; d < axis; ++d)
    {
        offset += SideGeometry<NDIM>::toSideBox(box,d).size();
    }
    return offset + side_box.offset(i);
}// compute_side_index

inline int
compute_cell_index(
    const Index<NDIM>& i,
    const Box<NDIM>& box)
{
    if (!box.contains(i)) return -1;
    int offset = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
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
    const blitz::TinyVector<double,NDIM>& dx)
{
    int ierr;

    const double rho = problem_coefs.getRho();
    const double mu = problem_coefs.getMu();
    const double lambda = problem_coefs.getLambda();

    // Allocate a PETSc matrix for the box operator.
    blitz::TinyVector<Box<NDIM>,NDIM> side_boxes;
    blitz::TinyVector<BoxList<NDIM>,NDIM> side_ghost_boxes;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        side_boxes[axis] = SideGeometry<NDIM>::toSideBox(box, axis);
        side_ghost_boxes[axis] = SideGeometry<NDIM>::toSideBox(ghost_box, axis);
        side_ghost_boxes[axis].removeIntersections(side_boxes[axis]);
    }
    BoxList<NDIM> cell_ghost_boxes(ghost_box);
    cell_ghost_boxes.removeIntersections(box);

    int size = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        size += SideGeometry<NDIM>::toSideBox(ghost_box, axis).size();
    }
    size += ghost_box.size();

    static const int U_stencil_sz = 2*NDIM+3;
    static const int P_stencil_sz = 2*NDIM+1;
    std::vector<int> nnz(size, 0);

    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            nnz[compute_side_index(b(), ghost_box, axis)] = std::min(size,U_stencil_sz);
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
        nnz[compute_cell_index(b(), ghost_box)] = std::min(size,P_stencil_sz);
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
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (Box<NDIM>::Iterator b(side_boxes[axis]); b; b++)
        {
            Index<NDIM> i = b();
            const int mat_row = compute_side_index(i, ghost_box, axis);

            std::vector<int> mat_cols(U_stencil_sz,-1);
            std::vector<double> mat_vals(U_stencil_sz,0.0);

            mat_cols[0] = mat_row;
            mat_vals[0] = (rho/dt) + 0.5*lambda;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index<NDIM> shift = 0;
                shift(d) = 1;
                const Index<NDIM> u_left = i - shift;
                const Index<NDIM> u_rght = i + shift;
                mat_cols[2*d+1] = compute_side_index(u_left, ghost_box, axis);
                mat_cols[2*d+2] = compute_side_index(u_rght, ghost_box, axis);

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

        for (unsigned int axis = 0; axis < NDIM; ++axis)
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
    for (unsigned int axis = 0; axis < NDIM; ++axis)
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

void
modifyRhsForBcs(
    Vec& v,
    const SideData<NDIM,double>& U_data,
    const CellData<NDIM,double>& P_data,
    const INSCoefs& problem_coefs,
    const double /*dt*/,
    const Box<NDIM>& box,
    const Box<NDIM>& ghost_box,
    const double* const dx)
{
    int ierr;

    const double mu = problem_coefs.getMu();
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            Index<NDIM> i = b();
            const int idx = compute_side_index(i, ghost_box, axis);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                Index<NDIM> shift = 0;
                shift(d) = 1;
                const Index<NDIM> u_left = i - shift;
                const Index<NDIM> u_rght = i + shift;
                if (!side_box.contains(u_left))
                {
                    ierr = VecSetValue(v, idx, +0.5*mu*U_data(SideIndex<NDIM>(u_left, axis, SideIndex<NDIM>::Lower))/(dx[d]*dx[d]), ADD_VALUES); IBTK_CHKERRQ(ierr);
                }
                if (!side_box.contains(u_rght))
                {
                    ierr = VecSetValue(v, idx, +0.5*mu*U_data(SideIndex<NDIM>(u_rght, axis, SideIndex<NDIM>::Lower))/(dx[d]*dx[d]), ADD_VALUES); IBTK_CHKERRQ(ierr);
                }
            }

            Index<NDIM> shift = 0;
            shift(axis) = 1;
            const Index<NDIM> p_left = i - shift;
            const Index<NDIM> p_rght = i;
            if (!box.contains(p_left))
            {
                ierr = VecSetValue(v, idx, +P_data(p_left)/dx[axis], ADD_VALUES); IBTK_CHKERRQ(ierr);
            }
            if (!box.contains(p_rght))
            {
                ierr = VecSetValue(v, idx, -P_data(p_rght)/dx[axis], ADD_VALUES); IBTK_CHKERRQ(ierr);
            }
        }
    }

    ierr = VecAssemblyBegin(v);  IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(v);  IBTK_CHKERRQ(ierr);
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

    for (unsigned int axis = 0; axis < NDIM; ++axis)
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

    const double omega = 0.65;

    double U;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const Index<NDIM>& i = b();
            const SideIndex<NDIM> s_i(i, axis, SideIndex<NDIM>::Lower);
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
    const Pointer<Database>& input_db)
    : d_object_name(object_name),
      d_is_initialized(false),
      d_solution(NULL),
      d_rhs(NULL),
      d_gcw(GHOSTS),
      d_patch_side_bc_box_overlap(),
      d_patch_cell_bc_box_overlap(),
      d_patch_side_smoother_bc_boxes(),
      d_patch_cell_smoother_bc_boxes(),
      d_hierarchy(),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_hier_bdry_fill_ops(),
      d_hier_math_ops(),
      d_in_initialize_operator_state(false),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_problem_coefs(),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_smoother_choice("additive"),
      d_U_prolongation_method("CONSTANT_REFINE"),
      d_P_prolongation_method("LINEAR_REFINE"),
      d_U_restriction_method("CONSERVATIVE_COARSEN"),
      d_P_restriction_method("CONSERVATIVE_COARSEN"),
      d_preconditioner(NULL),
      d_coarse_solver_choice("block_jacobi"),
      d_coarse_solver_tol(1.0e-6),
      d_coarse_solver_max_its(10),
      d_context(NULL),
      d_side_scratch_idx(-1),
      d_cell_scratch_idx(-1),
      d_U_bc_op(NULL),
      d_default_U_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_U_bc_coef", Pointer<Database>(NULL))),
      d_U_bc_coefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_U_bc_coef)),
      d_P_bc_op(NULL),
      d_default_P_bc_coef(new LocationIndexRobinBcCoefs<NDIM>(d_object_name+"::default_P_bc_coef", Pointer<Database>(NULL))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_current_time(0.0),
      d_new_time(0.0),
      d_U_cf_bdry_op(),
      d_P_cf_bdry_op(),
      d_U_op_stencil_fill_pattern(),
      d_P_op_stencil_fill_pattern(),
      d_U_synch_fill_pattern(),
      d_U_prolongation_refine_operator(),
      d_P_prolongation_refine_operator(),
      d_prolongation_refine_patch_strategy(),
      d_prolongation_refine_algorithm(),
      d_prolongation_refine_schedules(),
      d_U_restriction_coarsen_operator(),
      d_P_restriction_coarsen_operator(),
      d_restriction_coarsen_algorithm(),
      d_restriction_coarsen_schedules(),
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
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_U_bc_coef->setBoundaryValue(2*d  ,0.0);
        d_default_U_bc_coef->setBoundaryValue(2*d+1,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d  ,0.0);
        d_default_P_bc_coef->setBoundarySlope(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(d_default_U_bc_coef),d_default_P_bc_coef);

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
    IBAMR_DO_ONCE(
        t_restrict_residual         = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::restrictResidual()");
        t_prolong_error             = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::prolongError()");
        t_prolong_error_and_correct = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::prolongErrorAndCorrect()");
        t_smooth_error              = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::smoothError()");
        t_solve_coarsest_level      = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel()");
        t_compute_residual          = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::computeResidual()");
        t_initialize_operator_state = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::initializeOperatorState()");
        t_deallocate_operator_state = TimerManager::getManager()->getTimer("INSStaggeredBoxRelaxationFACOperator::deallocateOperatorState()");
                  );
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
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs,
    RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (U_bc_coefs[d] != NULL)
        {
            d_U_bc_coefs[d] = U_bc_coefs[d];
        }
        else
        {
            d_U_bc_coefs[d] = d_default_U_bc_coef;
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
INSStaggeredBoxRelaxationFACOperator::setTimeInterval(
    const double current_time,
    const double new_time)
{
    d_current_time = current_time;
    d_new_time = new_time;
    return;
}// setTimeInterval

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

///
///  The following routines:
///
///      setFACPreconditioner(),
///      restrictResidual(),
///      prolongError(),
///      prolongErrorAndCorrect(),
///      smoothError(),
///      solveCoarsestLevel(),
///      computeResidual(),
///      initializeOperatorState(),
///      deallocateOperatorState()
///
///  are concrete implementations of functions declared in the
///  FACPreconditionerStrategy abstract base class.
///

void
INSStaggeredBoxRelaxationFACOperator::setFACPreconditioner(
    ConstPointer<FACPreconditioner> preconditioner)
{
    if (d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::setFACPreconditioner()\n"
                   << "  cannot be called while operator state is initialized" << std::endl);
    }
    d_preconditioner = preconditioner;
    sanityCheck();
    return;
}// setFACPreconditioner

void
INSStaggeredBoxRelaxationFACOperator::restrictResidual(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBAMR_TIMER_START(t_restrict_residual);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    if (U_src_idx != U_dst_idx)
    {
        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        hier_sc_data_ops.copyData(U_dst_idx, U_src_idx, interior_only);
    }

    if (P_src_idx != P_dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, dst_ln, dst_ln);
        static const bool interior_only = false;
        hier_sc_data_ops.copyData(P_dst_idx, P_src_idx, interior_only);
    }

    xeqScheduleRestriction(dst_idxs, src_idxs, dst_ln);

    IBAMR_TIMER_STOP(t_restrict_residual);
    return;
}// restrictResidual

void
INSStaggeredBoxRelaxationFACOperator::prolongError(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBAMR_TIMER_START(t_prolong_error);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    // Refine the correction from the coarse level src data directly into the
    // fine level error.
    xeqScheduleProlongation(dst_idxs, src_idxs, dst_ln);

    IBAMR_TIMER_STOP(t_prolong_error);
    return;
}// prolongError

void
INSStaggeredBoxRelaxationFACOperator::prolongErrorAndCorrect(
    const SAMRAIVectorReal<NDIM,double>& src,
    SAMRAIVectorReal<NDIM,double>& dst,
    int dst_ln)
{
    IBAMR_TIMER_START(t_prolong_error_and_correct);

    const int U_src_idx = src.getComponentDescriptorIndex(0);
    const int P_src_idx = src.getComponentDescriptorIndex(1);
    const std::pair<int,int> src_idxs = std::make_pair(U_src_idx,P_src_idx);

    const int U_dst_idx = dst.getComponentDescriptorIndex(0);
    const int P_dst_idx = dst.getComponentDescriptorIndex(1);
    const std::pair<int,int> dst_idxs = std::make_pair(U_dst_idx,P_dst_idx);

    const std::pair<int,int> scratch_idxs = std::make_pair(d_side_scratch_idx,d_cell_scratch_idx);

    // Prolong the correction from the coarse level src data into the fine level
    // scratch data and then correct the fine level dst data.
    static const bool interior_only = false;
    if (U_src_idx != U_dst_idx)
    {
        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops_coarse(d_hierarchy, dst_ln-1, dst_ln-1);
        hier_sc_data_ops_coarse.add(U_dst_idx, U_dst_idx, U_src_idx, interior_only);
    }
    if (P_src_idx != P_dst_idx)
    {
        HierarchyCellDataOpsReal<NDIM,double> hier_sc_data_ops_coarse(d_hierarchy, dst_ln-1, dst_ln-1);
        hier_sc_data_ops_coarse.add(P_dst_idx, P_dst_idx, P_src_idx, interior_only);
    }
    xeqScheduleProlongation(scratch_idxs, src_idxs, dst_ln);
    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    hier_sc_data_ops_fine.add(U_dst_idx, U_dst_idx, d_side_scratch_idx, interior_only);
    HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops_fine(d_hierarchy, dst_ln, dst_ln);
    hier_cc_data_ops_fine.add(P_dst_idx, P_dst_idx, d_cell_scratch_idx, interior_only);

    IBAMR_TIMER_STOP(t_prolong_error_and_correct);
    return;
}// prolongErrorAndCorrect

void
INSStaggeredBoxRelaxationFACOperator::smoothError(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int level_num,
    int num_sweeps,
    bool /*performing_pre_sweeps*/,
    bool /*performing_post_sweeps*/)
{
    if (num_sweeps == 0) return;

    IBAMR_TIMER_START(t_smooth_error);

    Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(level_num);
    const int U_error_idx = error.getComponentDescriptorIndex(0);
    const int P_error_idx = error.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_side_scratch_idx;
    const int P_scratch_idx = d_cell_scratch_idx;

    // Cache coarse-fine interface ghost cell values in the "scratch" data.
    if (level_num > d_coarsest_ln && num_sweeps > 1)
    {
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            Pointer<SideData<NDIM,double> >   U_error_data = error.getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
            TBOX_ASSERT(  U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                U_scratch_data->getArrayData(axis).copy(
                    U_error_data->getArrayData(axis),
                    d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                    IntVector<NDIM>(0));
            }

            Pointer<CellData<NDIM,double> >   P_error_data = error.getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM,double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
            TBOX_ASSERT(  P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
            P_scratch_data->getArrayData().copy(
                P_error_data->getArrayData(),
                d_patch_cell_bc_box_overlap[level_num][patch_counter],
                IntVector<NDIM>(0));
        }
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
                int patch_counter = 0;
                for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());

                    Pointer<SideData<NDIM,double> >   U_error_data = error.getComponentPatchData(0, *patch);
                    Pointer<SideData<NDIM,double> > U_scratch_data = patch->getPatchData(U_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
                    TBOX_ASSERT(U_ghost_box == U_scratch_data->getGhostBox());
                    TBOX_ASSERT(  U_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(U_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        U_error_data->getArrayData(axis).copy(
                            U_scratch_data->getArrayData(axis),
                            d_patch_side_bc_box_overlap[level_num][patch_counter][axis],
                            IntVector<NDIM>(0));
                    }

                    Pointer<CellData<NDIM,double> >   P_error_data = error.getComponentPatchData(1, *patch);
                    Pointer<CellData<NDIM,double> > P_scratch_data = patch->getPatchData(P_scratch_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                    const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
                    TBOX_ASSERT(P_ghost_box == P_scratch_data->getGhostBox());
                    TBOX_ASSERT(  P_error_data->getGhostCellWidth() == d_gcw);
                    TBOX_ASSERT(P_scratch_data->getGhostCellWidth() == d_gcw);
#endif
                    P_error_data->getArrayData().copy(
                        P_scratch_data->getArrayData(),
                        d_patch_cell_bc_box_overlap[level_num][patch_counter],
                        IntVector<NDIM>(0));
                }

                // Fill the non-coarse-fine interface ghost cell values.
                const std::pair<int,int> error_idxs = std::make_pair(U_error_idx,P_error_idx);
                xeqScheduleGhostFillNoCoarse(error_idxs, level_num);
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
            const std::pair<int,int> error_idxs = std::make_pair(U_error_idx, P_error_idx);
            xeqScheduleGhostFillNoCoarse(error_idxs, level_num);
        }

        // Smooth the error on the patches.
        Vec& e = d_box_e[level_num];
        Vec& r = d_box_r[level_num];
        KSP& ksp = d_box_ksp[level_num];
        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<SideData<NDIM,double> >    U_error_data = error   .getComponentPatchData(0, *patch);
            Pointer<SideData<NDIM,double> > U_residual_data = residual.getComponentPatchData(0, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& U_ghost_box = U_error_data->getGhostBox();
            TBOX_ASSERT(U_ghost_box == U_residual_data->getGhostBox());
            TBOX_ASSERT(   U_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(U_residual_data->getGhostCellWidth() == d_gcw);
#endif
            Pointer<CellData<NDIM,double> >    P_error_data = error   .getComponentPatchData(1, *patch);
            Pointer<CellData<NDIM,double> > P_residual_data = residual.getComponentPatchData(1, *patch);
#ifdef DEBUG_CHECK_ASSERTIONS
            const Box<NDIM>& P_ghost_box = P_error_data->getGhostBox();
            TBOX_ASSERT(P_ghost_box == P_residual_data->getGhostBox());
            TBOX_ASSERT(   P_error_data->getGhostCellWidth() == d_gcw);
            TBOX_ASSERT(P_residual_data->getGhostCellWidth() == d_gcw);
#endif
            // Copy updated values from other local patches.
            if (d_smoother_choice == "multiplicative")
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    const std::map<int,Box<NDIM> > side_smoother_bc_boxes = d_patch_side_smoother_bc_boxes[level_num][patch_counter][axis];
                    for (std::map<int,Box<NDIM> >::const_iterator cit = side_smoother_bc_boxes.begin();
                         cit != side_smoother_bc_boxes.end(); ++cit)
                    {
                        const int src_patch_num = cit->first;
                        const Box<NDIM>& overlap = cit->second;
                        Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                        Pointer<SideData<NDIM,double> > src_U_error_data = error.getComponentPatchData(0, *src_patch);
                        U_error_data->getArrayData(axis).copy(src_U_error_data->getArrayData(axis), overlap, IntVector<NDIM>(0));
                    }
                }

                const std::map<int,Box<NDIM> > cell_smoother_bc_boxes = d_patch_cell_smoother_bc_boxes[level_num][patch_counter];
                for (std::map<int,Box<NDIM> >::const_iterator cit = cell_smoother_bc_boxes.begin();
                     cit != cell_smoother_bc_boxes.end(); ++cit)
                {
                    const int src_patch_num = cit->first;
                    const Box<NDIM>& overlap = cit->second;
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(src_patch_num);
                    Pointer<CellData<NDIM,double> > src_P_error_data = error.getComponentPatchData(1, *src_patch);
                    P_error_data->getArrayData().copy(src_P_error_data->getArrayData(), overlap, IntVector<NDIM>(0));
                }
            }

            // Smooth the error on the patch.
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            if (BOX_GHOSTS == 0)
            {
                // In this case, we must modify the righ-hand side to include
                // boundary condition values.
                for (Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const Box<NDIM> box(i,i);
                    copyToVec(e, *U_error_data, *P_error_data, box, box);
                    copyToVec(r, *U_residual_data, *P_residual_data, box, box);
                    modifyRhsForBcs(r, *U_error_data, *P_error_data, d_problem_coefs, d_dt, box, box, dx);
                    int ierr = KSPSolve(ksp, r, e);  IBTK_CHKERRQ(ierr);
                    copyFromVec(e, *U_error_data, *P_error_data, box, box);
                }
            }
            else
            {
                // In this case, boundary condition values are handled
                // implicitly by setting ghost values in the RHS vector.
                for (Box<NDIM>::Iterator b(patch_box); b; b++)
                {
                    const Index<NDIM>& i = b();
                    const Box<NDIM> box(i,i);
                    const Box<NDIM> ghost_box = hier::Box<NDIM>::grow(box,BOX_GHOSTS);
                    copyToVec(e, *U_error_data, *P_error_data, ghost_box, ghost_box);
                    copyToVec(r, *U_error_data, *P_error_data, ghost_box, ghost_box);
                    copyToVec(r, *U_residual_data, *P_residual_data, box, ghost_box);
                    int ierr = KSPSolve(ksp, r, e);  IBTK_CHKERRQ(ierr);
                    copyFromVec(e, *U_error_data, *P_error_data, box, ghost_box);
                }
            }
        }
    }

    // Synchronize data along patch boundaries.
    xeqScheduleSideDataSynch(U_error_idx, level_num);

    IBAMR_TIMER_STOP(t_smooth_error);
    return;
}// smoothError

bool
INSStaggeredBoxRelaxationFACOperator::solveCoarsestLevel(
    SAMRAIVectorReal<NDIM,double>& error,
    const SAMRAIVectorReal<NDIM,double>& residual,
    int coarsest_ln)
{
    IBAMR_TIMER_START(t_solve_coarsest_level);

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(coarsest_ln == d_coarsest_ln);
#endif

    smoothError(error, residual, coarsest_ln, d_coarse_solver_max_its, false, false);

    IBAMR_TIMER_STOP(t_solve_coarsest_level);
    return true;
}// solveCoarsestLevel

void
INSStaggeredBoxRelaxationFACOperator::computeResidual(
    SAMRAIVectorReal<NDIM,double>& residual,
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs,
    int level_num)
{
    IBAMR_TIMER_START(t_compute_residual);

    if (!d_preconditioner.isNull() && d_preconditioner->getNumPreSmoothingSweeps() == 0)
    {
        // Compute the residual, r = f - A*u = f - A*0.
        residual.copyVector(Pointer<SAMRAIVectorReal<NDIM,double> >(const_cast<SAMRAIVectorReal<NDIM,double>*>(&rhs),false), false);
    }
    else
    {
        // Compute the residual, r = f - A*u.
        const int U_res_idx = residual.getComponentDescriptorIndex(0);
        const int U_sol_idx = solution.getComponentDescriptorIndex(0);
        const int U_rhs_idx = rhs.getComponentDescriptorIndex(0);

        const Pointer<SideVariable<NDIM,double> > U_res_sc_var = residual.getComponentVariable(0);
        const Pointer<SideVariable<NDIM,double> > U_sol_sc_var = solution.getComponentVariable(0);
        const Pointer<SideVariable<NDIM,double> > U_rhs_sc_var = rhs.getComponentVariable(0);

        const int P_res_idx = residual.getComponentDescriptorIndex(1);
        const int P_sol_idx = solution.getComponentDescriptorIndex(1);
        const int P_rhs_idx = rhs.getComponentDescriptorIndex(1);

        const Pointer<CellVariable<NDIM,double> > P_res_cc_var = residual.getComponentVariable(1);
        const Pointer<CellVariable<NDIM,double> > P_sol_cc_var = solution.getComponentVariable(1);
        const Pointer<CellVariable<NDIM,double> > P_rhs_cc_var = rhs.getComponentVariable(1);

        // NOTE: Here, we assume that the residual is to be computed only during
        // pre-sweeps and only for zero initial guesses, so that we need to
        // compute A*u ONLY on levels level_num and level_num-1.

        // Fill ghost-cell values.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        Pointer<VariableFillPattern<NDIM> > sc_fill_pattern = new SideNoCornersFillPattern(GHOSTS, false, false, true);
        Pointer<VariableFillPattern<NDIM> > cc_fill_pattern = new CellNoCornersFillPattern(GHOSTS, false, false, true);
        InterpolationTransactionComponent U_scratch_component(U_sol_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_U_bc_coefs, sc_fill_pattern);
        InterpolationTransactionComponent P_scratch_component(P_sol_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef , cc_fill_pattern);
        std::vector<InterpolationTransactionComponent> U_P_components(2);
        U_P_components[0] = U_scratch_component;
        U_P_components[1] = P_scratch_component;
        if (d_hier_bdry_fill_ops[level_num].isNull())
        {
            d_hier_bdry_fill_ops[level_num] = new HierarchyGhostCellInterpolation();
            d_hier_bdry_fill_ops[level_num]->initializeOperatorState(U_P_components, d_hierarchy, level_num, level_num);
        }
        else
        {
            d_hier_bdry_fill_ops[level_num]->resetTransactionComponents(U_P_components);
        }
        d_hier_bdry_fill_ops[level_num]->setHomogeneousBc(true);
        d_hier_bdry_fill_ops[level_num]->fillData(d_new_time);
        if (level_num > d_coarsest_ln) xeqScheduleGhostFillNoCoarse(std::make_pair(U_sol_idx, P_sol_idx), level_num-1);

        // Compute the residual, r = f - A*u.  We assume that u=0 for all levels
        // coarser than level_num, and therefore that A*u = 0 on all levels
        // coarser than level_num-1.  (A*u may be non-zero on level_num-1
        // because of coarse-grid corrections at coarse-fine interfaces.)
        if (d_hier_math_ops[level_num].isNull())
        {
            std::ostringstream stream;
            stream << d_object_name << "::hier_math_ops_" << level_num;
            d_hier_math_ops[level_num] = new HierarchyMathOps(stream.str(), d_hierarchy, std::max(d_coarsest_ln,level_num-1), level_num);
        }

        const double rho = d_problem_coefs.getRho();
        const double mu = d_problem_coefs.getMu();
        const double lambda = d_problem_coefs.getLambda();
        PoissonSpecifications helmholtz_spec("");
        helmholtz_spec.setCConstant((rho/d_dt)+0.5*lambda);
        helmholtz_spec.setDConstant(          -0.5*mu    );

        static const bool cf_bdry_synch = true;
        d_hier_math_ops[level_num]->grad(
            U_res_idx, U_res_sc_var,
            cf_bdry_synch,
            1.0, P_sol_idx, P_sol_cc_var, NULL, d_new_time);
        d_hier_math_ops[level_num]->laplace(
            U_res_idx, U_res_sc_var,
            helmholtz_spec,
            U_sol_idx, U_sol_sc_var, NULL, d_new_time,
            1.0,
            U_res_idx, U_res_sc_var);

        d_hier_math_ops[level_num]->div(
            P_res_idx, P_res_cc_var,
            -1.0, U_sol_idx, U_sol_sc_var, NULL, d_new_time,
            cf_bdry_synch);

        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(d_hierarchy, std::max(d_coarsest_ln,level_num-1), level_num);
        hier_sc_data_ops.axpy(U_res_idx, -1.0, U_res_idx, U_rhs_idx, false);
        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(d_hierarchy, std::max(d_coarsest_ln,level_num-1), level_num);
        hier_cc_data_ops.axpy(P_res_idx, -1.0, P_res_idx, P_rhs_idx, false);
    }

    IBAMR_TIMER_STOP(t_compute_residual);
    return;
}// computeResidual

void
INSStaggeredBoxRelaxationFACOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& solution,
    const SAMRAIVectorReal<NDIM,double>& rhs)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

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

    // Setup level operators.
    d_hier_bdry_fill_ops.resize(d_finest_ln+1, NULL);
    d_hier_math_ops.resize(d_finest_ln+1, NULL);
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        d_hier_bdry_fill_ops[ln].setNull();
        d_hier_math_ops[ln].setNull();
    }

    // Allocate scratch data.
    for (int ln = std::max(d_coarsest_ln, coarsest_reset_ln); ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_side_scratch_idx)) level->allocatePatchData(d_side_scratch_idx);
        if (!level->checkAllocated(d_cell_scratch_idx)) level->allocatePatchData(d_cell_scratch_idx);
    }

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = d_hierarchy->getGridGeometry();
    IBAMR_DO_ONCE(
        geometry->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
        geometry->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
                  );

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
    d_U_restriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_U_restriction_method);
    var_db->mapIndexToVariable(d_cell_scratch_idx, var);
    d_P_restriction_coarsen_operator = geometry->lookupCoarsenOperator(var, d_P_restriction_method);
    d_U_ghostfill_nocoarse_refine_operator = NULL;
    d_P_ghostfill_nocoarse_refine_operator = NULL;
    d_U_synch_refine_operator = NULL;

    // Make space for saving communication schedules.  There is no need to
    // delete the old schedules first because we have deallocated the solver
    // state above.
    d_U_bc_op = new CartSideRobinPhysBdryOp(d_side_scratch_idx, d_U_bc_coefs, false);
    d_P_bc_op = new CartCellRobinPhysBdryOp(d_cell_scratch_idx, d_P_bc_coef , false);

    d_U_op_stencil_fill_pattern = new SideNoCornersFillPattern(GHOSTS, false, false, false);
    d_P_op_stencil_fill_pattern = new CellNoCornersFillPattern(GHOSTS, false, false, false);
    d_U_synch_fill_pattern = new SideSynchCopyFillPattern();

    std::vector<RefinePatchStrategy<NDIM>*> prolongation_refine_patch_strategies;
    prolongation_refine_patch_strategies.push_back(d_U_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_P_cf_bdry_op);
    prolongation_refine_patch_strategies.push_back(d_U_bc_op);
    prolongation_refine_patch_strategies.push_back(d_P_bc_op);
    d_prolongation_refine_patch_strategy = new RefinePatchStrategySet(
        prolongation_refine_patch_strategies.begin(), prolongation_refine_patch_strategies.end(), false);

    d_prolongation_refine_schedules.resize(d_finest_ln+1);
    d_restriction_coarsen_schedules.resize(d_finest_ln+1);
    d_ghostfill_nocoarse_refine_schedules.resize(d_finest_ln+1);
    d_U_synch_refine_schedules.resize(d_finest_ln+1);

    d_prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    d_restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
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

    d_restriction_coarsen_algorithm->registerCoarsen(
        d_side_scratch_idx,
        rhs.getComponentDescriptorIndex(0),
        d_U_restriction_coarsen_operator);
    d_restriction_coarsen_algorithm->registerCoarsen(
        d_cell_scratch_idx,
        rhs.getComponentDescriptorIndex(1),
        d_P_restriction_coarsen_operator);

    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        solution.getComponentDescriptorIndex(0),
        d_U_ghostfill_nocoarse_refine_operator,
        d_U_op_stencil_fill_pattern);
    d_ghostfill_nocoarse_refine_algorithm->registerRefine(
        solution.getComponentDescriptorIndex(1),
        solution.getComponentDescriptorIndex(1),
        solution.getComponentDescriptorIndex(1),
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

        d_ghostfill_nocoarse_refine_schedules[dst_ln] =
            d_ghostfill_nocoarse_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln), d_U_P_bc_op);

        d_U_synch_refine_schedules[dst_ln] =
            d_U_synch_refine_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln));
    }

    d_ghostfill_nocoarse_refine_schedules[d_coarsest_ln] =
        d_ghostfill_nocoarse_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln), d_U_P_bc_op);

    d_U_synch_refine_schedules[d_coarsest_ln] =
        d_U_synch_refine_algorithm->createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_ln));

    for (int dst_ln = d_coarsest_ln; dst_ln < d_finest_ln; ++dst_ln)
    {
        d_restriction_coarsen_schedules[dst_ln] =
            d_restriction_coarsen_algorithm->createSchedule(
                d_hierarchy->getPatchLevel(dst_ln  ),
                d_hierarchy->getPatchLevel(dst_ln+1));
    }

    // Initialize the box relaxation data on each level of the patch hierarchy.
    d_box_op.resize(d_finest_ln+1);
    d_box_e.resize(d_finest_ln+1);
    d_box_r.resize(d_finest_ln+1);
    d_box_ksp.resize(d_finest_ln+1);

    const Box<NDIM> box(Index<NDIM>(0),Index<NDIM>(0));
    const Box<NDIM> ghost_box = hier::Box<NDIM>::grow(box,BOX_GHOSTS);

    const double* const dx_coarsest = geometry->getDx();
    blitz::TinyVector<double,NDIM> dx;

    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        const IntVector<NDIM>& ratio = d_hierarchy->getPatchLevel(ln)->getRatio();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dx[d] = dx_coarsest[d]/static_cast<double>(ratio(d));
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
        ierr = KSPSetUp(d_box_ksp[ln]);  IBTK_CHKERRQ(ierr);
    }

    // Get overlap information for setting patch boundary conditions.
    d_patch_side_bc_box_overlap.resize(d_finest_ln+1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_side_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box,axis);
                const Box<NDIM> side_ghost_box = Box<NDIM>::grow(side_box, 1);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis] = BoxList<NDIM>(side_ghost_box);
                d_patch_side_bc_box_overlap[ln][patch_counter][axis].removeIntersections(side_box);
            }
        }
    }

    d_patch_cell_bc_box_overlap.resize(d_finest_ln+1);
    for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
        d_patch_cell_bc_box_overlap[ln].resize(num_local_patches);

        int patch_counter = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++patch_counter)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Box<NDIM>& ghost_box = Box<NDIM>::grow(patch_box, 1);

            d_patch_cell_bc_box_overlap[ln][patch_counter] = BoxList<NDIM>(ghost_box);
            d_patch_cell_bc_box_overlap[ln][patch_counter].removeIntersections(patch_box);
        }
    }

    // Get overlap information for re-setting patch boundary conditions during
    // multiplicative smoothing.
    if (d_smoother_choice == "multiplicative")
    {
        d_patch_side_smoother_bc_boxes.resize(d_finest_ln+1);
        for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
            d_patch_side_smoother_bc_boxes[ln].resize(num_local_patches);

            int patch_counter1 = 0;
            for (PatchLevel<NDIM>::Iterator p1(level); p1; p1++, ++patch_counter1)
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    d_patch_side_smoother_bc_boxes[ln][patch_counter1][axis].clear();
                }

                Pointer<Patch<NDIM> > dst_patch = level->getPatch(p1());
                const Box<NDIM>& dst_patch_box = dst_patch->getBox();
                const Box<NDIM>& dst_ghost_box = Box<NDIM>::grow(dst_patch_box, 1);

                int patch_counter2 = 0;
                for (PatchLevel<NDIM>::Iterator p2(level); patch_counter2 < patch_counter1; p2++, ++patch_counter2)
                {
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(p2());
                    const Box<NDIM>& src_patch_box = src_patch->getBox();

                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        const Box<NDIM> overlap =
                            SideGeometry<NDIM>::toSideBox(dst_ghost_box,axis) *
                            SideGeometry<NDIM>::toSideBox(src_patch_box,axis);
                        if (!overlap.empty())
                        {
                            d_patch_side_smoother_bc_boxes[ln][patch_counter1][axis][p2()] = overlap;
                        }
                    }
                }
            }
        }

        d_patch_cell_smoother_bc_boxes.resize(d_finest_ln+1);
        for (int ln = coarsest_reset_ln; ln <= finest_reset_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

            const int num_local_patches = level->getProcessorMapping().getLocalIndices().getSize();
            d_patch_cell_smoother_bc_boxes[ln].resize(num_local_patches);

            int patch_counter1 = 0;
            for (PatchLevel<NDIM>::Iterator p1(level); p1; p1++, ++patch_counter1)
            {
                d_patch_cell_smoother_bc_boxes[ln][patch_counter1].clear();

                Pointer<Patch<NDIM> > dst_patch = level->getPatch(p1());
                const Box<NDIM>& dst_patch_box = dst_patch->getBox();
                const Box<NDIM>& dst_ghost_box = Box<NDIM>::grow(dst_patch_box, 1);

                int patch_counter2 = 0;
                for (PatchLevel<NDIM>::Iterator p2(level); patch_counter2 < patch_counter1; p2++, ++patch_counter2)
                {
                    Pointer<Patch<NDIM> > src_patch = level->getPatch(p2());
                    const Box<NDIM>& src_patch_box = src_patch->getBox();
                    const Box<NDIM> overlap = dst_ghost_box * src_patch_box;
                    if (!overlap.empty())
                    {
                        d_patch_cell_smoother_bc_boxes[ln][patch_counter1][p2()] = overlap;
                    }
                }
            }
        }
    }
    else
    {
        d_patch_side_smoother_bc_boxes.clear();
        d_patch_cell_smoother_bc_boxes.clear();
    }

    // Indicate that the operator is initialized.
    d_is_initialized = true;
    d_in_initialize_operator_state = false;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
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

    IBAMR_TIMER_START(t_deallocate_operator_state);

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
            d_patch_side_smoother_bc_boxes.resize(0);
            d_patch_cell_smoother_bc_boxes.resize(0);

            d_hierarchy.setNull();
            d_coarsest_ln = -1;
            d_finest_ln   = -1;

            d_hier_bdry_fill_ops.clear();
            d_hier_math_ops.clear();

            d_U_prolongation_refine_operator    .setNull();
            d_U_cf_bdry_op                      .setNull();
            d_prolongation_refine_patch_strategy.setNull();
            d_prolongation_refine_algorithm     .setNull();
            d_prolongation_refine_schedules     .resize(0);

            d_U_restriction_coarsen_operator.setNull();
            d_restriction_coarsen_algorithm .setNull();
            d_restriction_coarsen_schedules .resize(0);

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

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
}// deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleProlongation(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(true);
    d_U_cf_bdry_op->setPatchDataIndex(U_dst_idx);

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setPhysicalBcCoef(d_P_bc_coef);
    d_P_bc_op->setHomogeneousBc(true);
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
INSStaggeredBoxRelaxationFACOperator::xeqScheduleRestriction(
    const std::pair<int,int>& dst_idxs,
    const std::pair<int,int>& src_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    const int U_src_idx = src_idxs.first;

    const int P_dst_idx = dst_idxs.second;
    const int P_src_idx = src_idxs.second;

    CoarsenAlgorithm<NDIM> coarsener;
    coarsener.registerCoarsen(U_dst_idx, U_src_idx, d_U_restriction_coarsen_operator);
    coarsener.registerCoarsen(P_dst_idx, P_src_idx, d_P_restriction_coarsen_operator);
    coarsener.resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    d_restriction_coarsen_schedules[dst_ln]->coarsenData();
    d_restriction_coarsen_algorithm->resetSchedule(d_restriction_coarsen_schedules[dst_ln]);
    return;
}// xeqScheduleRestriction

void
INSStaggeredBoxRelaxationFACOperator::xeqScheduleGhostFillNoCoarse(
    const std::pair<int,int>& dst_idxs,
    const int dst_ln)
{
    const int U_dst_idx = dst_idxs.first;
    d_U_bc_op->setPatchDataIndex(U_dst_idx);
    d_U_bc_op->setPhysicalBcCoefs(d_U_bc_coefs);
    d_U_bc_op->setHomogeneousBc(true);

    const int P_dst_idx = dst_idxs.second;
    d_P_bc_op->setPatchDataIndex(P_dst_idx);
    d_P_bc_op->setPhysicalBcCoef(d_P_bc_coef);
    d_P_bc_op->setHomogeneousBc(true);

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

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d_U_bc_coefs[d] == NULL)
        {
            TBOX_ERROR(d_object_name << ":\n"
                       << "  invalid velocity physical bc object at depth = " << d << std::endl);
        }
    }

    if (d_P_bc_coef == NULL)
    {
        TBOX_ERROR(d_object_name << ":\n"
                   << "  invalid pressure physical bc object" << std::endl);
    }
    return;
}// sanityCheck

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSStaggeredBoxRelaxationFACOperator>;

//////////////////////////////////////////////////////////////////////////////
