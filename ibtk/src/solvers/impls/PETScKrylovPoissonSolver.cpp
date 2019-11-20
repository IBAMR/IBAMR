// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
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

#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

#include "tbox/Database.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScKrylovPoissonSolver::PETScKrylovPoissonSolver(std::string object_name,
                                                   Pointer<Database> input_db,
                                                   std::string default_options_prefix)
    : PETScKrylovLinearSolver(std::move(object_name), input_db, std::move(default_options_prefix))
{
    // intentionally blank
    return;
} // PETScKrylovPoissonSolver()

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
