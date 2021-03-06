// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#include "ibtk/HierarchyMathOps.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/LinearSolver.h"

#include "tbox/Pointer.h"

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
KrylovLinearSolver::setHierarchyMathOps(Pointer<HierarchyMathOps> hier_math_ops)
{
    LinearSolver::setHierarchyMathOps(hier_math_ops);
    if (d_A) d_A->setHierarchyMathOps(d_hier_math_ops);
    if (d_pc_solver) d_pc_solver->setHierarchyMathOps(d_hier_math_ops);
    return;
} // setHierarchyMathOps

void
KrylovLinearSolver::setHomogeneousBc(bool homogeneous_bc)
{
    LinearSolver::setHomogeneousBc(homogeneous_bc);
    if (d_A) d_A->setHomogeneousBc(homogeneous_bc);
    return;
} // setHomogeneousBc

void
KrylovLinearSolver::setSolutionTime(const double solution_time)
{
    LinearSolver::setSolutionTime(solution_time);
    if (d_A) d_A->setSolutionTime(solution_time);
    if (d_pc_solver) d_pc_solver->setSolutionTime(solution_time);
    return;
} // setSolutionTime

void
KrylovLinearSolver::setTimeInterval(const double current_time, const double new_time)
{
    LinearSolver::setTimeInterval(current_time, new_time);
    if (d_A) d_A->setTimeInterval(current_time, new_time);
    if (d_pc_solver) d_pc_solver->setTimeInterval(current_time, new_time);
    return;
} // setTimeInterval

void
KrylovLinearSolver::setOperator(Pointer<LinearOperator> A)
{
    Pointer<LinearOperator> A_old = d_A;
    d_A = A;
    if (d_A)
    {
        d_A->setHomogeneousBc(d_homogeneous_bc);
        d_A->setSolutionTime(d_solution_time);
        d_A->setTimeInterval(d_current_time, d_new_time);
        if (d_is_initialized && (d_A != A_old) && d_A)
        {
            d_A->initializeOperatorState(*d_x, *d_b);
        }
    }
    return;
} // setOperator

Pointer<LinearOperator>
KrylovLinearSolver::getOperator() const
{
    return d_A;
} // getOperator

void
KrylovLinearSolver::setPreconditioner(Pointer<LinearSolver> pc_solver)
{
    Pointer<LinearSolver> pc_solver_old = d_pc_solver;
    d_pc_solver = pc_solver;
    if (d_pc_solver)
    {
        d_pc_solver->setHomogeneousBc(true);
        d_pc_solver->setSolutionTime(d_solution_time);
        d_pc_solver->setTimeInterval(d_current_time, d_new_time);
        if (d_is_initialized && (d_pc_solver != pc_solver_old) && d_pc_solver)
        {
            d_pc_solver->initializeSolverState(*d_b, *d_b);
        }
    }
    return;
} // setPreconditioner

Pointer<LinearSolver>
KrylovLinearSolver::getPreconditioner() const
{
    return d_pc_solver;
} // getPreconditioner

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
