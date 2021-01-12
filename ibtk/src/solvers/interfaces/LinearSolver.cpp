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

#include "ibtk/LinearSolver.h"

#include "Box.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include <ostream>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LinearSolver::LinearSolver()
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
