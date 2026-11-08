// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_MultiphaseEx3_InitialConditions
#define included_MultiphaseEx3_InitialConditions

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_utilities.h>

#include <CellVariable.h>

#include <string>
#include <vector>

namespace MultiphaseEx3
{
SAMRAI::tbox::Pointer<IBTK::CartGridFunction>
make_sphere_initial_condition(const std::string& object_name,
                              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> var,
                              const IBTK::VectorNd& center,
                              double radius);

SAMRAI::tbox::Pointer<IBTK::CartGridFunction>
make_velocity_initial_condition(const std::string& object_name,
                                const IBTK::VectorNd& center,
                                double radius,
                                double num_interface_cells,
                                const std::vector<double>& inside_velocity,
                                const std::vector<double>& outside_velocity);
} // namespace MultiphaseEx3

#endif
