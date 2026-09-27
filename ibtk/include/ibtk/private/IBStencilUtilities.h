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

#ifndef included_IBTK_IBStencilUtilities
#define included_IBTK_IBStencilUtilities
#include <ibtk/config.h>
namespace IBTK::detail
{
/*! \brief Place a stencil using a cell index and the physical lattice origin. */
inline int ib_stencil_lower(double x, double x_lower, double dx, int domain_lower, int cell, int width, double offset);
/*! \brief Measure a point from the stencil's first grid location in mesh units. */
inline double ib_stencil_coordinate(double x, double x_lower, double dx, int domain_lower, int lower, double offset);
} // namespace IBTK::detail
#include <ibtk/private/IBStencilUtilities-inl.h>
#endif
