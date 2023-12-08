// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_falling_sphere_generate_sphere
#define included_falling_sphere_generate_sphere

#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <vector>

// Adapted from the sphereGen3d generator by Namu Patel; retain its point ordering.
inline void
generate_sphere(const unsigned int& /*structure_number*/,
                const int& /*level_number*/,
                int& num_vertices,
                std::vector<IBTK::Point>& vertex_posn,
                void* ctx)
{
    auto* input_db = static_cast<SAMRAI::tbox::Database*>(ctx);
    const double radius = input_db->getDouble("R");
    const int num_points = input_db->getInteger("num_points_per_diameter");
    if (radius <= 0.0)
    {
        TBOX_ERROR("The sphere radius must be positive.\n");
    }
    if (num_points < 3)
    {
        TBOX_ERROR("num_points_per_diameter must be at least 3 to generate a nondegenerate sphere.");
    }
    const double ds = 2.0 * radius / num_points;
    vertex_posn.clear();
    for (int k = 0; k < num_points; ++k)
    {
        const double z = -radius + k * ds;
        for (int i = 0; i < num_points; ++i)
        {
            const double x = -radius + i * ds;
            for (int j = 0; j < num_points; ++j)
            {
                const double y = -radius + j * ds;
                if (x * x + y * y + z * z < radius * radius)
                {
                    vertex_posn.emplace_back(x, y, z);
                }
            }
        }
    }
    num_vertices = static_cast<int>(vertex_posn.size());
}

#endif
