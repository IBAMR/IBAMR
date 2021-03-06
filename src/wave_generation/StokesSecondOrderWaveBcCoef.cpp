// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/StokesSecondOrderWaveBcCoef.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Database.h"

#include <cmath>
#include <utility>

#include "ibamr/namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int EXTENSIONS_FILLABLE = 128;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesSecondOrderWaveBcCoef::StokesSecondOrderWaveBcCoef(std::string object_name,
                                                         const int comp_idx,
                                                         Pointer<Database> input_db,
                                                         Pointer<CartesianGridGeometry<NDIM> > grid_geom)
    : d_object_name(std::move(object_name)),
      d_comp_idx(comp_idx),
      d_muparser_bcs(d_object_name + "::muParser", input_db, grid_geom),
      d_grid_geom(grid_geom)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
    TBOX_ASSERT(input_db);
#endif
    // Get wave parameters.
    getFromInput(input_db);

    return;
} // StokesSecondOrderWaveBcCoef

StokesSecondOrderWaveBcCoef::~StokesSecondOrderWaveBcCoef()
{
    // intentionally left-blank.
    return;
} // ~StokesSecondOrderWaveBcCoef

void
StokesSecondOrderWaveBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                        Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                        Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                        const Pointer<Variable<NDIM> >& variable,
                                        const Patch<NDIM>& patch,
                                        const BoundaryBox<NDIM>& bdry_box,
                                        double fill_time) const
{
    // Get pgeom info.
    const Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    // Compute a representative grid spacing
    double vol_cell = 1.0;
    for (int d = 0; d < NDIM; ++d) vol_cell *= dx[d];
    auto alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

    const int dir = (NDIM == 2) ? 1 : (NDIM == 3) ? 2 : -1;

    // We take location index = 0, i.e., x_lower to be the wave inlet.
    const unsigned int location_index = bdry_box.getLocationIndex();
    if (location_index != 0)
    {
        d_muparser_bcs.setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
    }
    else
    {
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM>& bc_coef_box = (acoef_data ? acoef_data->getBox() :
                                        bcoef_data ? bcoef_data->getBox() :
                                        gcoef_data ? gcoef_data->getBox() :
                                                     Box<NDIM>());
#if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data || bc_coef_box == acoef_data->getBox());
        TBOX_ASSERT(!bcoef_data || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(!gcoef_data || bc_coef_box == gcoef_data->getBox());
#endif

        const double H = 2.0 * d_amplitude;
        const double fac1 = (0.5 * H * d_gravity * d_wave_number) / d_omega;
        const double fac2 = (3.0 * H * H * d_omega * d_wave_number) / 16.0;
        double dof_posn[NDIM];
        for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const hier::Index<NDIM>& i = b();
            if (acoef_data) (*acoef_data)(i, 0) = 1.0;
            if (bcoef_data) (*bcoef_data)(i, 0) = 0.0;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d != bdry_normal_axis)
                {
                    dof_posn[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                }
                else
                {
                    dof_posn[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                }
            }

            // Compute a numerical heaviside at the boundary from the analytical wave
            // elevation
            const double theta = d_wave_number * dof_posn[0] - d_omega * fill_time;
            const double kd = d_wave_number * d_depth;
            const double z_plus_d = dof_posn[dir];
            const double eta = 0.5 * H * cos(theta) + H * H * d_wave_number / 16.0 * cosh(kd) * (2.0 + cosh(2.0 * kd)) /
                                                          std::pow(sinh(kd), 3.0) * cos(2.0 * theta);
            const double phi = -eta + (z_plus_d - d_depth);
            double h_phi;
            if (phi < -alpha)
                h_phi = 1.0;
            else if (std::abs(phi) <= alpha)
                h_phi = 1.0 - (0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha));
            else
                h_phi = 0.0;

            if (d_comp_idx == 0)
            {
                if (gcoef_data)
                {
                    (*gcoef_data)(i, 0) = h_phi * (fac1 * cosh(d_wave_number * (z_plus_d)) * cos(theta) / cosh(kd) +
                                                   fac2 * cosh(2.0 * d_wave_number * z_plus_d) * cos(2.0 * theta) /
                                                       std::pow(sinh(kd), 4.0));
                }
            }

            if (d_comp_idx == 1)
            {
#if (NDIM == 2)
                if (gcoef_data)
                {
                    (*gcoef_data)(i, 0) = h_phi * (fac1 * sinh(d_wave_number * z_plus_d) * sin(theta) / cosh(kd) +
                                                   fac2 * sinh(2.0 * d_wave_number * z_plus_d) * sin(2.0 * theta) /
                                                       std::pow(sinh(kd), 4.0));
                }
#elif (NDIM == 3)
                if (gcoef_data) (*gcoef_data)(i, 0) = 0.0;
#endif
            }

#if (NDIM == 3)
            if (d_comp_idx == 2)
            {
                if (gcoef_data)
                {
                    (*gcoef_data)(i, 0) = h_phi * (fac1 * sinh(d_wave_number * z_plus_d) * sin(theta) / cosh(kd) +
                                                   fac2 * sinh(2.0 * d_wave_number * z_plus_d) * sin(2.0 * theta) /
                                                       std::pow(sinh(kd), 4.0));
                }
            }
#endif
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
StokesSecondOrderWaveBcCoef::numberOfExtensionsFillable() const
{
    return EXTENSIONS_FILLABLE;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
void
StokesSecondOrderWaveBcCoef::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> wave_db = input_db->getDatabase("wave_parameters_db");
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db->isDatabase("wave_parameters_db"));
#endif

    d_depth = wave_db->getDouble("depth");
    d_omega = wave_db->getDouble("omega");
    d_gravity = wave_db->getDouble("gravitational_constant");
    d_wave_number = wave_db->getDouble("wave_number");
    d_amplitude = wave_db->getDouble("amplitude");
    d_num_interface_cells = wave_db->getDouble("num_interface_cells");

    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR
