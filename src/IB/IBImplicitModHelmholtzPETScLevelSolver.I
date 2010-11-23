// Filename: IBImplicitModHelmholtzPETScLevelSolver.I
// Created on 27 Sep 2010 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline void
IBImplicitModHelmholtzPETScLevelSolver::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    d_initial_guess_nonzero = initial_guess_nonzero;
    return;
}// setInitialGuessNonzero

inline bool
IBImplicitModHelmholtzPETScLevelSolver::getInitialGuessNonzero() const
{
    return d_initial_guess_nonzero;
}// getInitialGuessNonzero

inline void
IBImplicitModHelmholtzPETScLevelSolver::setMaxIterations(
    int max_iterations)
{
    d_max_iterations = max_iterations;
    return;
}// setMaxIterations

inline int
IBImplicitModHelmholtzPETScLevelSolver::getMaxIterations() const
{
    return d_max_iterations;
}// getMaxIterations

inline void
IBImplicitModHelmholtzPETScLevelSolver::setAbsoluteTolerance(
    double abs_residual_tol)
{
    d_abs_residual_tol = abs_residual_tol;
    return;
}//setAbsoluteTolerance

inline double
IBImplicitModHelmholtzPETScLevelSolver::getAbsoluteTolerance() const
{
    return d_abs_residual_tol;
}// getAbsoluteTolerance

inline void
IBImplicitModHelmholtzPETScLevelSolver::setRelativeTolerance(
    double rel_residual_tol)
{
    d_rel_residual_tol = rel_residual_tol;
    return;
}//setRelativeTolerance

inline double
IBImplicitModHelmholtzPETScLevelSolver::getRelativeTolerance() const
{
    return d_rel_residual_tol;
}// getRelativeTolerance

inline int
IBImplicitModHelmholtzPETScLevelSolver::getNumIterations() const
{
    return d_current_its;
}// getNumIterations

inline double
IBImplicitModHelmholtzPETScLevelSolver::getResidualNorm() const
{
    return d_current_residual_norm;
}// getResidualNorm

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
