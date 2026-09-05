// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/IBKernel.h>

namespace
{
constexpr IBTK::IBKernel separate_translation_unit_kernel = IBTK::IBKernel::IB_4;
}

bool
ib_kernel_static_initialization_valid()
{
    return separate_translation_unit_kernel == IBTK::IBKernel("IB_4");
}
