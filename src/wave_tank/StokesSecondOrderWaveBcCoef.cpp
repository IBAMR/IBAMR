// Filename: StokesSecondOrderWaveBcCoef.cpp
// Created on 1 Sep 2018 by Amneet Bhalla
//
// Copyright (c) 2002-2018, Amneet Bhalla
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
#include "ibamr/namespaces.h"
#include "ibtk/muParserRobinBcCoefs.h"
#include <ibamr/AdvDiffHierarchyIntegrator.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int EXTENSIONS_FILLABLE = 128;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesSecondOrderWaveBcCoef::StokesSecondOrderWaveBcCoef(const std::string& object_name,
                                                         const int comp_idx,
                                                         Pointer<Database> input_db,
                                                         Pointer<CartesianGridGeometry<NDIM> > grid_geom,
                                                         Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
                                                         Pointer<CellVariable<NDIM, double> > ls_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(adv_diff_solver);
    TBOX_ASSERT(ls_var);
    TBOX_ASSERT(input_db);
#endif

    d_object_name = object_name;
    d_comp_idx = comp_idx;
    d_grid_geom = grid_geom;
    d_adv_diff_solver = adv_diff_solver;
    d_ls_var = ls_var;

    // Create muParser object for boundaries except wave inlet.
    d_muparser_bcs = new muParserRobinBcCoefs(object_name + "::muParser", input_db, grid_geom);

    // Get wave parameters.
    getFromInput(input_db);

    return;
} // StokesSecondOrderWaveBcCoef

StokesSecondOrderWaveBcCoef::~StokesSecondOrderWaveBcCoef()
{
    delete d_muparser_bcs;
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
    // Get level set patch data object.
    // Note that we are using current context here, because the new context gives a different value of
    // level set field with cycle number. Therefore, the boundary condition changes with cycle number for
    // new context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int ls_idx = var_db->mapVariableAndContextToIndex(d_ls_var, d_adv_diff_solver->getCurrentContext());
    Pointer<CellData<NDIM, double> > ls_data = patch.getPatchData(ls_idx);

    // Get pgeom info.
    const Box<NDIM>& patch_box = patch.getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    // Get physical domain box extents on patch box level.
    const SAMRAI::hier::Box<NDIM> domain_box =
        SAMRAI::hier::Box<NDIM>::refine(d_grid_geom->getPhysicalDomain()[0], pgeom->getRatio());
    const IntVector<NDIM>& domain_box_lower = domain_box.lower();
    const double* const domain_x_lower = d_grid_geom->getXLower();
    const int dir = (NDIM == 2) ? 1 : (NDIM == 3) ? 2 : -1;

    // We take location index = 0, i.e., x_lower to be the wave inlet.
    const unsigned int location_index = bdry_box.getLocationIndex();
    if (location_index != 0)
    {
        d_muparser_bcs->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
    }
    else
    {
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM>& bc_coef_box =
            (acoef_data ? acoef_data->getBox() :
                          bcoef_data ? bcoef_data->getBox() : gcoef_data ? gcoef_data->getBox() : Box<NDIM>());
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
            const Index<NDIM>& i = b();
            // const double phi = (*ls_data)(i,0);
            const double dz = dx[dir];

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
            const double theta = d_wave_number * dof_posn[0] - d_omega * fill_time;
            const double kd = d_wave_number * d_depth;
            const double z_cell_surface =
                (domain_x_lower[dir] + dz * (static_cast<double>(i(dir) - domain_box_lower(dir)) + 0.5)) - d_depth;
	    const double eta = 0.5 * H * cos(theta) + H * H * d_wave_number/16.0 * cosh(kd) * (2.0 + cosh(2.0*kd))/ std::pow(sinh(kd), 3.0) * cos(2.0*theta);
            const double phi = -eta + z_cell_surface;
            double h_phi;
            const double dw = 2.0*dz;
            if (phi < -dz) h_phi = 1.0;
            else if (std::abs(phi) <= dz) h_phi = 1.0 - (0.5 + 0.5 * phi / dz + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / dz));
            else h_phi = 0.0;
            // const double vol_frac = (std::max(std::min(eta - z_cell_surface, dw / 2.0), -dw / 2.0) + dw / 2.0) / dw;
            const double vol_frac = h_phi;
            // const double vol_frac = 1.0;

            if (d_comp_idx == 0)
            {
                if (gcoef_data)
                {
                    const double z_plus_d = z_cell_surface + d_depth;
                    (*gcoef_data)(i, 0) =
                        vol_frac * (fac1 * cosh(d_wave_number * (z_plus_d)) * cos(theta) / cosh(kd) +
                                    fac2 * cosh(2.0 * d_wave_number * z_plus_d) * cos(2.0 * theta) / std::pow(sinh(kd), 4.0));
                }
            }

            if (d_comp_idx == 1)
            {
#if (NDIM == 2)
                if (gcoef_data)
                {
                    const double z_plus_d = z_cell_surface - dz / 2.0 + d_depth;
                    (*gcoef_data)(i, 0) =
                        vol_frac * (fac1 * sinh(d_wave_number * z_plus_d) * sin(theta) / cosh(kd) +
                                    fac2 * sinh(2.0 * d_wave_number * z_plus_d) * sin(2.0 * theta) / std::pow(sinh(kd), 4.0));
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
                    const double z_plus_d = z_cell_surface - dz / 2.0 + d_depth;
                    (*gcoef_data)(i, 0) =
                        vol_frac * (fac1 * sinh(d_wave_number * z_plus_d) * sin(theta) / cosh(kd) +
                                    fac2 * sinh(2.0 * kd * z_plus_d) * sin(2.0 * theta) / std::pow(sinh(kd), 4.0));
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

    return;
} // getFromInput

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR
