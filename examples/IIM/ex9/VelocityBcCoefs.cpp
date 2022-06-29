// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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

#include "VelocityBcCoefs.h"

#include <CartesianPatchGeometry.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////
namespace
{
// static double U1 = 6.3e1;
// static double U2 = 6.15e1;
//
// static double t_load = 0.5;
//
// static std::string inflow_data;
//
// static double wall = 0.24;
// static double d_in = 2.19;
// static double d_out = 7.62;
// static double h1_in = 2.72;
// static double h2_in = -2.72;
// static double h_out = 0.0;
// static double z_min = -2.95;
// static double z_max = 8.95;

}

/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

// VelocityBcCoefs::VelocityBcCoefs(const CirculationModel* circ_model, const int comp_idx)
//    : d_circ_model(circ_model), d_comp_idx(comp_idx) //set d_circ = circ_model, d_comp_idx = comp_idx

VelocityBcCoefs::VelocityBcCoefs(const INSHierarchyIntegrator* fluid_solver, const BcData bc_data, const int comp_idx)
    : d_fluid_solver(fluid_solver), d_comp_idx(comp_idx), d_bc_data(bc_data)

//  d_U1(bc_data.U1),
//  d_U2(d_U2),
//  d_wall(bc_data.wall),
//  d_d_in(bc_data.d_in),
//  d_d_out(bc_data.d_out),
//  d_h1_in(bc_data.h1_in),
//  d_h2_in(bc_data.h2_in),
//  d_h_out(bc_data.h_out),
//  d_z_min(bc_data.z_min),
//  d_z_max(bc_data.z_max)

// set d_circ = circ_model, d_comp_idx = comp_idx
{
    // this class extends Robin BC class
    // d_comp_idx is the component of velocity we are concerned with

    // intentionally blank
    return;
} // VelocityBcCoefs

VelocityBcCoefs::~VelocityBcCoefs()
{
    // intentionally blank
    return;
} // ~VelocityBcCoefs

void
VelocityBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                            Pointer<ArrayData<NDIM, double> >& bcoef_data,
                            Pointer<ArrayData<NDIM, double> >& gcoef_data,
                            const Pointer<Variable<NDIM> >& /*variable*/,
                            const Patch<NDIM>& patch,
                            const BoundaryBox<NDIM>& bdry_box,
                            double fill_time) const
{
    // bdry_box is for filling ghost data
    // getLocationIndex says where boundary is relative to given patch
    //    face (codimension 1):
    //    x_lo: 0
    //    x_hi: 1
    //    y_lo: 2
    //    y_hi: 3
    //    z_lo: 4
    //    z_hi: 5

    const int location_index = bdry_box.getLocationIndex(); // refers to faces of domain
    const int axis = location_index / 2;                    // returns 0 for x direction, 1 for y, 2 for z
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data.isNull());
#endif
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif

    const Box<NDIM>& patch_box = patch.getBox();
    const hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower(); // pointer to lower coordinates of patch

    //    U1 = d_input_db->getDoubleWithDefault("U1", U1);
    //    U2 = d_input_db->getDoubleWithDefault("U2", U2);
    //    wall = d_input_db->getDoubleWithDefault("WALL", wall);
    //    d_in = d_input_db->getDoubleWithDefault("D_IN", d_in);
    //    d_out = d_input_db->getDoubleWithDefault("D_OUT", d_out);
    //    h1_in = d_input_db->getDoubleWithDefault("H1_IN", h1_in);
    //    h2_in = d_input_db->getDoubleWithDefault("H2_IN", h2_in);
    //    h_out = d_input_db->getDoubleWithDefault("H_OUT", h_out);
    //    z_min = d_input_db->getDoubleWithDefault("Z_MIN", z_min);
    //    z_max = d_input_db->getDoubleWithDefault("Z_MAX", z_max);

    // Loops through cells of bc_coef_box, which is box over which we need BC's
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const hier::Index<NDIM>& i = bc();
        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);

        // set velocity to zero on x and y faces... for components d_comp_idx
        if (axis == 0 || axis == 1)
        {
            a = 1.0;
            b = 0.0;
            g = 0.0;
        }

        // z_lo face
        else if (location_index == 4)
        {
            if (d_comp_idx != axis)
            {
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }
            else
            {
                const double r_in = 0.5 * d_bc_data.d_in - d_bc_data.wall;         // units cm
                const Point& posn1 = Point(0.0, d_bc_data.h1_in, d_bc_data.z_min); // d_circ_model->d_posn[side];
                const Point& posn2 = Point(0.0, d_bc_data.h2_in, d_bc_data.z_min);

                //                const double r_in = 0.5 * d_in - wall; //units cm
                //                const Point& posn1 = Point(0.0,h1_in,z_min); //d_circ_model->d_posn[side];
                //                const Point& posn2 = Point(0.0,h2_in,z_min);
                double X[NDIM];
                double r_sq1 = 0.0;
                double r_sq2 = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                    if (d != axis) r_sq1 += pow(X[d] - posn1[d], 2.0); // sq of dist from posn1 on plane of boundary
                    if (d != axis) r_sq2 += pow(X[d] - posn2[d], 2.0); // sq of dist from posn2 on plane of boundary
                }
                const double r1 = sqrt(r_sq1);
                const double r2 = sqrt(r_sq2);

                const double u_in = parabolic_flow(
                    fill_time, X[1], X[1] > 0 ? r1 / r_in : r2 / r_in); // d_circ_model->d_psrc[side]; //pressure source

                // if inside circle of radius rsrc, then bc is u = h(t), if outside then bc is u = 0
                if (fill_time <= d_bc_data.t_load)
                {
                    a = (r1 <= r_in || r2 <= r_in) ? 0.0 : 1.0;
                    b = (r1 <= r_in || r2 <= r_in) ? 1.0 : 0.0;
                    g = 0.0;
                }
                else
                {
                    a = 1.0;
                    b = 0.0;
                    g = (r1 <= r_in || r2 <= r_in) ? u_in : 0.0;
                }
            }
        }

        // z_hi face
        else if (location_index == 5)
        {
            if (d_comp_idx != axis)
            {
                a = 1.0;
                b = 0.0;
                g = 0.0;
            }
            else
            {
                const double r_out = 0.5 * d_bc_data.d_out - d_bc_data.wall;      // units cm
                const Point& posn = Point(0.0, d_bc_data.h_out, d_bc_data.z_max); // d_circ_model->d_posn[side];

                //                const double r_out = 0.5 * d_out - wall; //units cm
                //                const Point& posn = Point(0.0,h_out,z_max); //d_circ_model->d_posn[side];
                double X[NDIM];
                double r_sq = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_lower(d)) + (d == axis ? 0.0 : 0.5));
                    if (d != axis) r_sq += pow(X[d] - posn[d], 2.0); // sq of dist from posn1 on plane of boundary
                }
                const double r = sqrt(r_sq);

                a = (r <= r_out) ? 0.0 : 1.0;
                b = (r <= r_out) ? 1.0 : 0.0;
                g = 0.0;
            }
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
VelocityBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable

double
VelocityBcCoefs::parabolic_flow(double t, double y, double r) const
{
    //    double t_load = 0.5;
    double c = y > 0.0 ? d_bc_data.U1 : d_bc_data.U2;

    //    double c = y > 0.0 ? U1 : U2;
    double U = c * (1.0 - r * r);
    return time_ramp(t) * U;
}

double
VelocityBcCoefs::time_ramp(double t) const
{
    return (t < d_bc_data.t_load) ? 0.0 : 1.0; //? -16.0*t*t*t + 12.0*t*t : 1.0;

    //    return (t < t_load) ? -16.0*t*t*t + 12.0*t*t : 1.0;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
