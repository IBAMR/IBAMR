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

#include "ibamr/PETScKrylovStaggeredStokesSolver.h"

#include "tbox/Database.h"

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PETScKrylovStaggeredStokesSolver::PETScKrylovStaggeredStokesSolver(const std::string& object_name,
                                                                   Pointer<Database> input_db,
                                                                   const std::string& default_options_prefix)
    : PETScKrylovLinearSolver(object_name, input_db, default_options_prefix)
{
    // intentionally blank
    return;
} // PETScKrylovStaggeredStokesSolver()

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
