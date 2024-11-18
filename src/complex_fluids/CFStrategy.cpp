// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "ibamr/CFStrategy.h"

#include "ibtk/ibtk_utilities.h"

#include "tbox/Utilities.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <unsupported/Eigen/MatrixFunctions>
IBTK_ENABLE_EXTRA_WARNINGS

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

// Namespace
namespace IBAMR
{
CFStrategy::CFStrategy(std::string object_name) : d_object_name(std::move(object_name))
{
    return;
} // Constructor

} // namespace IBAMR
