// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_IBTK_Experimental_SideCoupling
#define included_IBTK_Experimental_SideCoupling

#include <ibtk/config.h>

#include <CartesianCoupling.h>

namespace IBTK::Experimental
{
/*! \brief Patch-local, matrix-free coupling of side-centered field components. */
using SideCoupling = CartesianCoupling<DataCentering::SIDE>;
} // namespace IBTK::Experimental

#endif
