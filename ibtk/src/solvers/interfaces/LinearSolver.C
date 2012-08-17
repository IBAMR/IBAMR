// Filename: LinearSolver.C
// Created on 09 Sep 2003 by Boyce Griffith
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

#include "LinearSolver.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LinearSolver::LinearSolver(
    const std::string& object_name,
    bool homogeneous_bc)
    : GeneralSolver(object_name, homogeneous_bc),
      d_nullspace_contains_constant_vector(false),
      d_nullspace_basis_vecs(),
      d_initial_guess_nonzero(true),
      d_rel_residual_tol(1.0e-5),
      d_abs_residual_tol(1.0e-50),
      d_max_iterations(10000),
      d_current_its(0),
      d_current_residual_norm(std::numeric_limits<double>::quiet_NaN())
{
    // intentionally blank
    return;
}// LinearSolver()

LinearSolver::~LinearSolver()
{
    // intentionally blank
    return;
}// ~LinearSolver()

void
LinearSolver::setNullspace(
    const bool nullspace_containsconstant_vector,
    const std::vector<Pointer<SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs)
{
    d_nullspace_contains_constant_vector = nullspace_containsconstant_vector;
    d_nullspace_basis_vecs = nullspace_basis_vecs;
    return;
}// setNullspace

void
LinearSolver::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    d_initial_guess_nonzero = initial_guess_nonzero;
    return;
}// setInitialGuessNonzero

bool
LinearSolver::getNullspaceContainsConstantVector() const
{
    return d_nullspace_contains_constant_vector;
}// getNullspaceContainsConstantVector

const std::vector<Pointer<SAMRAIVectorReal<NDIM,double> > >&
LinearSolver::getNullspaceBasisVectors() const
{
    return d_nullspace_basis_vecs;
}// getNullspaceBasisVectors

bool
LinearSolver::getInitialGuessNonzero() const
{
    return d_initial_guess_nonzero;
}// getInitialGuessNonzero

void
LinearSolver::setMaxIterations(
    int max_iterations)
{
    d_max_iterations = max_iterations;
    return;
}// setMaxIterations

int
LinearSolver::getMaxIterations() const
{
    return d_max_iterations;
}// getMaxIterations

void
LinearSolver::setAbsoluteTolerance(
    double abs_residual_tol)
{
    d_abs_residual_tol = abs_residual_tol;
    return;
}// setAbsoluteTolerance

double
LinearSolver::getAbsoluteTolerance() const
{
    return d_abs_residual_tol;
}// getAbsoluteTolerance

void
LinearSolver::setRelativeTolerance(
    double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    return;
}// setRelativeTolerance

double
LinearSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
}// getRelativeTolerance

int
LinearSolver::getNumIterations() const
{
    return d_current_its;
}// getNumIterations

double
LinearSolver::getResidualNorm() const
{
    return d_current_residual_norm;
}// getResidualNorm

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
