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

#include "QInit.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <SAMRAI_config.h>

#include <array>

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

QInit::QInit(const string& object_name, Pointer<GridGeometry<NDIM> > grid_geom, Pointer<Database> input_db)
    : CartGridFunction(object_name),
      d_object_name(object_name),
      d_grid_geom(grid_geom),
      d_X(),
      d_init_type("GAUSSIAN"),
      d_gaussian_kappa(0.01),
      d_zalesak_r(0.15),
      d_zalesak_slot_w(0.025),
      d_zalesak_slot_l(0.1)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(grid_geom);
#endif
    d_object_name = object_name;
    d_grid_geom = grid_geom;
#if !defined(NDEBUG)
    TBOX_ASSERT(d_grid_geom);
#endif

    // Default initial values.
    const double* const x_upper = d_grid_geom->getXUpper();
    const double* const x_lower = d_grid_geom->getXLower();

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_X[d] = x_lower[d] + 0.5 * (x_upper[d] - x_lower[d]);
    }

    d_init_type = "GAUSSIAN";

    d_gaussian_kappa = 0.01;

    d_zalesak_r = 0.15;
    d_zalesak_slot_w = 0.025;
    d_zalesak_slot_l = 0.1;

    // Initialize object with data read from the input database.
    getFromInput(input_db);

    return;
} // QInit

QInit::~QInit()
{
    // intentionally blank
    return;
} // ~QInit

void
QInit::setDataOnPatch(const int data_idx,
                      Pointer<Variable<NDIM> > /*var*/,
                      Pointer<Patch<NDIM> > patch,
                      const double data_time,
                      const bool /*initial_time*/,
                      Pointer<PatchLevel<NDIM> > /*level*/)
{
    Pointer<CellData<NDIM, double> > Q_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(Q_data);
#endif
    const Box<NDIM>& patch_box = patch->getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    double r_squared;
    VectorNd X;
    const double t = data_time;

    Q_data->fillAll(0.0);

    if (d_init_type == "GAUSSIAN")
    {
        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            // NOTE: This assumes the lattice of Gaussians are being advected
            // and diffused in the unit square.
            std::array<int, NDIM> offset;
            for (offset[0] = -2; offset[0] <= 2; ++(offset[0]))
            {
                for (offset[1] = -2; offset[1] <= 2; ++(offset[1]))
                {
#if (NDIM > 2)
                    for (offset[2] = -2; offset[2] <= 2; ++(offset[2]))
                    {
#endif
                        r_squared = 0.0;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                            r_squared += pow(X[d] - (d_X[d] + static_cast<double>(offset[d])), 2.0);
                        }
                        (*Q_data)(i) += exp(-r_squared / (4.0 * d_gaussian_kappa)) /
                                        pow(4.0 * M_PI * d_gaussian_kappa * (1.0 + t), 0.5 * static_cast<double>(NDIM));
#if (NDIM > 2)
                    }
#endif
                }
            }
        }
    }
    else if (d_init_type == "ZALESAK")
    {
        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            r_squared = 0.0;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                r_squared += pow((X[d] - d_X[d]), 2.0);
            }
            if ((sqrt(r_squared) > d_zalesak_r) ||
                ((abs(X[0] - d_X[0]) < d_zalesak_slot_w) && (X[1] - d_X[1]) < d_zalesak_slot_l))
            {
                (*Q_data)(i) = 0.0;
            }
            else
            {
                (*Q_data)(i) = 1.0;
            }
        }
    }
    else
    {
        TBOX_ERROR(d_object_name << "::initializeDataOnPatch()\n"
                                 << "  invalid initialization type " << d_init_type << "\n");
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
QInit::getFromInput(Pointer<Database> db)
{
    if (db)
    {
        if (db->keyExists("X"))
        {
            db->getDoubleArray("X", d_X.data(), NDIM);
        }

        d_init_type = db->getStringWithDefault("init_type", d_init_type);

        if (d_init_type == "GAUSSIAN")
        {
            d_gaussian_kappa = db->getDoubleWithDefault("kappa", d_gaussian_kappa);
        }
        else if (d_init_type == "ZALESAK")
        {
            d_zalesak_r = db->getDoubleWithDefault("zalesak_r", d_zalesak_r);
            d_zalesak_slot_w = db->getDoubleWithDefault("zalesak_slot_w", d_zalesak_slot_w);
            d_zalesak_slot_l = db->getDoubleWithDefault("zalesak_slot_l", d_zalesak_slot_l);
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
