// Filename: LinearSolver.cpp
// Created on 09 Sep 2003 by Boyce Griffith
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

#include <ostream>
#include <vector>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/GeneralSolver.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LinearSolver::LinearSolver()
    : d_initial_guess_nonzero(true), d_nullspace_contains_constant_vec(false), d_nullspace_basis_vecs()
{
    d_max_iterations = 10000;
    d_rel_residual_tol = 1.0e-5;
    d_abs_residual_tol = 1.0e-50;
    return;
} // LinearSolver()

LinearSolver::~LinearSolver()
{
    // intentionally blank
    return;
} // ~LinearSolver()

void
LinearSolver::setNullspace(const bool nullspace_containsconstant_vec,
                           const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >& nullspace_basis_vecs)
{
    d_nullspace_contains_constant_vec = nullspace_containsconstant_vec;
    d_nullspace_basis_vecs = nullspace_basis_vecs;
    return;
} // setNullspace

void
LinearSolver::setInitialGuessNonzero(bool initial_guess_nonzero)
{
    d_initial_guess_nonzero = initial_guess_nonzero;
    return;
} // setInitialGuessNonzero

bool
LinearSolver::getNullspaceContainsConstantVector() const
{
    return d_nullspace_contains_constant_vec;
} // getNullspaceContainsConstantVector

const std::vector<Pointer<SAMRAIVectorReal<NDIM, double> > >&
LinearSolver::getNullspaceBasisVectors() const
{
    return d_nullspace_basis_vecs;
} // getNullspaceBasisVectors

bool
LinearSolver::getInitialGuessNonzero() const
{
    return d_initial_guess_nonzero;
} // getInitialGuessNonzero

void
LinearSolver::printClassData(std::ostream& stream)
{
    GeneralSolver::printClassData(stream);
    stream << "initial_guess_nonzero = " << d_initial_guess_nonzero << "\n"
           << "nullspace_contains_constant_vec = " << d_nullspace_contains_constant_vec << "\n"
           << "nullspace_basis_vecs.size() = " << d_nullspace_basis_vecs.size() << "\n";
    return;
} // printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
