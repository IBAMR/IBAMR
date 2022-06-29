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

#include "FeedbackForcer.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <tbox/Utilities.h>

#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <SideData.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// static std::string inflow_data;

// static double wall = 0.24;
// static double d_in = 2.19;
// static double d_out = 7.62;
// static double h1_in = 2.72;
// static double h2_in = -2.72;
// static double h_out = 0.0;
// static double z_min = -2.95;
// static double z_max = 8.95;

inline double
smooth_kernel(const double r)
{
    return std::abs(r) < 1.0 ? 0.5 * (cos(M_PI * r) + 1.0) : 0.0;
} // smooth_kernel
} // namespace

////////////////////////////// PUBLIC ///////////////////////////////////////

FeedbackForcer::FeedbackForcer(const INSHierarchyIntegrator* fluid_solver,
                               const Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                               const BcData& bc_data)
    : d_fluid_solver(fluid_solver), d_patch_hierarchy(patch_hierarchy), d_bc_data(bc_data)
{
    // d_bc_coefs is a vector of pointers to RobinBcCeofStrategy objects, which are themselves NDIM-d
    // intentionally blank
    return;
} // FeedbackForcer

FeedbackForcer::~FeedbackForcer()
{
    // intentionally blank
    return;
} // ~FeedbackForcer

bool
FeedbackForcer::isTimeDependent() const
{
    return true;
} // isTimeDependent

void
FeedbackForcer::setDataOnPatch(const int data_idx,
                               Pointer<Variable<NDIM> > /*var*/,
                               Pointer<Patch<NDIM> > patch,
                               const double /*data_time*/,
                               const bool initial_time,
                               Pointer<PatchLevel<NDIM> > /*patch_level*/)
{
    //    wall = d_input_db->getDoubleWithDefault("WALL", wall);
    //    d_in = d_input_db->getDoubleWithDefault("D_IN", d_in);
    //    d_out = d_input_db->getDoubleWithDefault("D_OUT", d_out);
    //    h1_in = d_input_db->getDoubleWithDefault("H1_IN", h1_in);
    //    h2_in = d_input_db->getDoubleWithDefault("H2_IN", h2_in);
    //    h_out = d_input_db->getDoubleWithDefault("H_OUT", h_out);
    //    z_min = d_input_db->getDoubleWithDefault("Z_MIN", z_min);
    //    z_max = d_input_db->getDoubleWithDefault("Z_MAX", z_max);

    Pointer<SideData<NDIM, double> > F_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(F_data);
#endif
    F_data->fillAll(0.0);
    if (initial_time) return;
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    const double dt = d_fluid_solver->getCurrentTimeStepSize();
    const double rho = d_fluid_solver->getStokesSpecifications()->getRho();
    const double kappa = cycle_num >= 0 ? 0.25 * rho / dt : 0.0;
    Pointer<SideData<NDIM, double> > U_current_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    Pointer<SideData<NDIM, double> > U_new_data =
        patch->getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(U_current_data);
#endif
    const Box<NDIM>& patch_box = patch->getBox();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    const double* const x_lower = pgeom->getXLower();
    const double* const x_upper = pgeom->getXUpper();
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = d_patch_hierarchy->getGridGeometry();
    const double* const dx_coarsest = grid_geometry->getDx();
    double dx_finest[NDIM];
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    const IntVector<NDIM>& finest_ratio = d_patch_hierarchy->getPatchLevel(finest_ln)->getRatio();
    for (int d = 0; d < NDIM; ++d)
    {
        dx_finest[d] = dx_coarsest[d] / static_cast<double>(finest_ratio(d));
    }
    const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], pgeom->getRatio());

    // Clamp the velocity along the top/bottom of the domain (but not in the
    // interior of the structure).
    static const int axis = 2; // concerned with z_lo and z_hi faces
    const double L = 4.0 * dx_finest[axis];
    const int offset = static_cast<int>(L / dx[axis]);
    for (int side = 0; side <= 1; ++side)
    {
        const bool is_lower = side == 0;
        if (pgeom->getTouchesRegularBoundary(axis, side))
        {
            Box<NDIM> bdry_box = domain_box;
            if (side == 0)
            {
                bdry_box.upper(axis) = domain_box.lower(axis) + offset;
            }
            else
            {
                bdry_box.lower(axis) = domain_box.upper(axis) - offset;
            }
            bdry_box = bdry_box * patch_box;
            for (int component = 0; component < NDIM; ++component)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(bdry_box, component)); b; b++)
                {
                    double r_in = 0.5 * d_bc_data.d_in - d_bc_data.wall;   // d_bc_data.wall;
                    double r_out = 0.5 * d_bc_data.d_out - d_bc_data.wall; // d_bc_data.wall;
                    const Point& posn1 = Point(0.0, d_bc_data.h1_in, d_bc_data.z_min);
                    const Point& posn2 = Point(0.0, d_bc_data.h2_in, d_bc_data.z_min);
                    const Point& posn3 = Point(0.0, d_bc_data.h_out, d_bc_data.z_max);

                    //                    double r_in = 0.5 * d_in  - wall; //d_bc_data.wall;
                    //                    double r_out = 0.5 * d_out - wall; //d_bc_data.wall;
                    //                    const Point& posn1 = Point(0.0, h1_in, z_min);
                    //                    const Point& posn2 = Point(0.0, h2_in, z_min);
                    //                    const Point& posn3 = Point(0.0, h_out, z_max);
                    double r_sq1 = 0.0, r_sq2 = 0.0, r_sq3 = 0.0;

                    const hier::Index<NDIM>& i = b();
                    const SideIndex<NDIM> i_s(i, component, SideIndex<NDIM>::Lower);
                    const double U_current = U_current_data ? (*U_current_data)(i_s) : 0.0;
                    const double U_new = U_new_data ? (*U_new_data)(i_s) : 0.0;
                    const double U = (cycle_num > 0) ? 0.5 * (U_new + U_current) : U_current;

                    double X[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == component ? 0.0 : 0.5));
                        if (d != axis) r_sq1 += pow(X[d] - posn1[d], 2.0);
                        if (d != axis) r_sq2 += pow(X[d] - posn2[d], 2.0);
                        if (d != axis) r_sq3 += pow(X[d] - posn3[d], 2.0);
                    }
                    const double r1 = sqrt(r_sq1);
                    const double r2 = sqrt(r_sq2);
                    const double r3 = sqrt(r_sq3);

                    double fac;
                    if (component == axis)
                    {
                        if (side == 0)
                        {
                            fac = (r1 <= r_in || r2 <= r_in) ? 0.0 : 1.0;
                        }
                        // side 1
                        else
                        {
                            fac = (r3 <= r_out) ? 0.0 : 1.0;
                        }
                    }
                    // x, y component
                    else
                    {
                        fac = 1.0;
                    }
                    const double x_bdry = (is_lower ? x_lower[axis] : x_upper[axis]);
                    fac *= smooth_kernel((X[axis] - x_bdry) / L);
                    (*F_data)(i_s) = -fac * kappa * U;
                }
            }
        }
    }
    return;
} // setDataOnPatch

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
