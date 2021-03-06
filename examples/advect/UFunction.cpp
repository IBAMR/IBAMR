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

#include "UFunction.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

UFunction::UFunction(const string& object_name, Pointer<GridGeometry<NDIM> > grid_geom, Pointer<Database> input_db)
    : CartGridFunction(object_name),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_X(),
      d_init_type("UNIFORM"),
      d_uniform_u()
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(grid_geom);
#endif

    // Default initial values.
    const double* const x_upper = d_grid_geom->getXUpper();
    const double* const x_lower = d_grid_geom->getXLower();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_X[d] = x_lower[d] + 0.5 * (x_upper[d] - x_lower[d]);
        d_uniform_u[d] = 1.0;
    }

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    return;
} // UFunction

UFunction::~UFunction()
{
    // intentionally blank
    return;
} // ~UFunction

void
UFunction::setDataOnPatch(const int data_idx,
                          Pointer<Variable<NDIM> > /*var*/,
                          Pointer<Patch<NDIM> > patch,
                          const double /*data_time*/,
                          const bool /*initial_time*/,
                          Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<FaceData<NDIM, double> > u_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(u_data);
#endif

    if (d_init_type == "UNIFORM")
    {
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            u_data->getArrayData(axis).fillAll(d_uniform_u[axis]);
        }
    }
    else if (d_init_type == "VORTEX")
    {
        const Box<NDIM>& patch_box = patch->getBox();
        const hier::Index<NDIM>& patch_lower = patch_box.lower();
        Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();

        VectorNd X;

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (FaceIterator<NDIM> it(patch_box, axis); it; it++)
            {
                const FaceIndex<NDIM>& i = it();
                const hier::Index<NDIM>& cell_idx = i.toCell(1);

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] +
                           dx[d] * (static_cast<double>(cell_idx(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                }

                // 2D vortex
                if (axis == 0)
                {
                    (*u_data)(i) = (X[1] - d_X[axis]);
                }
                else if (axis == 1)
                {
                    (*u_data)(i) = (d_X[axis] - X[0]);
                }
                else
                {
                    (*u_data)(i) = 0.0;
                }
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::setDataOnPatch()\n"
                                 << "  invalid initialization type " << d_init_type << "\n");
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
UFunction::getFromInput(Pointer<Database> db)
{
    if (db)
    {
        if (db->keyExists("X"))
        {
            db->getDoubleArray("X", d_X.data(), NDIM);
        }

        d_init_type = db->getStringWithDefault("init_type", d_init_type);

        if (d_init_type == "UNIFORM")
        {
            if (db->keyExists("uniform_u"))
            {
                db->getDoubleArray("uniform_u", d_uniform_u.data(), NDIM);
            }
        }
        else if (d_init_type == "VORTEX")
        {
            // intentionally blank
        }
        else
        {
            TBOX_ERROR(d_object_name << "::getFromInput()\n"
                                     << "  invalid initialization type " << d_init_type << "\n");
        }
    }
    return;
} // getFromInput

//////////////////////////////////////////////////////////////////////////////
