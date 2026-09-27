// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_Experimental_CartesianCoupling
#define included_IBTK_Experimental_CartesianCoupling
#include <ibtk/private/CartesianCoupling.h>
namespace IBTK::Experimental
{
using IBTK::detail::TensorProductMode;
template <DataCentering C>
using CartesianCoupling = IBTK::detail::CartesianCoupling<C>;
} // namespace IBTK::Experimental
#endif
