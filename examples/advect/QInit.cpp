// Filename: QInit.cpp
// Created on 19 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "QInit.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <IBAMR_config.h>
#include <SAMRAI_config.h>

#include "boost/array.hpp"

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
    const Index<NDIM>& patch_lower = patch_box.lower();
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
            const Index<NDIM>& i = ic();
            // NOTE: This assumes the lattice of Gaussians are being advected
            // and diffused in the unit square.
            boost::array<int, NDIM> offset;
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
            const Index<NDIM>& i = ic();
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
                                 << "  invalid initialization type "
                                 << d_init_type
                                 << "\n");
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
                                     << "  invalid initialization type "
                                     << d_init_type
                                     << "\n");
        }
    }
    return;
} // getFromInput

//////////////////////////////////////////////////////////////////////////////
