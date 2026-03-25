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

#include <ibamr/INSSGSKinematics.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSTurbulenceModel.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int TENSOR_DEPTH = NDIM * (NDIM + 1) / 2;
constexpr double PI = 3.141592653589793238462643383279502884;

using Vector = std::array<double, NDIM>;
using SymTensor = std::array<double, TENSOR_DEPTH>;

enum class ManufacturedCase
{
    UNIFORM,
    RIGID_ROTATION,
    SIMPLE_SHEAR,
    TRIGONOMETRIC
};

struct SmagorinskyParameters
{
    double rho = 1.0;
    double smagorinsky_constant = 0.17;
    double filter_width_scale = 1.0;
    double max_turbulent_viscosity = std::numeric_limits<double>::max();
};

void
allocate_patch_data(const int idx, const double data_time, const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, data_time);
    }
}

void
deallocate_patch_data(const int idx, const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(idx)) level->deallocatePatchData(idx);
    }
}

hier::IntVector<NDIM>
axis_shift(const int axis, const int offset)
{
    hier::IntVector<NDIM> shift(0);
    shift(axis) = offset;
    return shift;
}

ManufacturedCase
parse_case(const std::string& case_name)
{
    if (case_name == "UNIFORM") return ManufacturedCase::UNIFORM;
    if (case_name == "RIGID_ROTATION") return ManufacturedCase::RIGID_ROTATION;
    if (case_name == "SIMPLE_SHEAR") return ManufacturedCase::SIMPLE_SHEAR;
    if (case_name == "TRIGONOMETRIC") return ManufacturedCase::TRIGONOMETRIC;
    TBOX_ERROR("Unsupported manufactured SGS case " << case_name << "\n");
    return ManufacturedCase::UNIFORM;
}

Vector
evaluate_velocity(const ManufacturedCase manufactured_case, const Vector& X)
{
    Vector U = {};
    switch (manufactured_case)
    {
    case ManufacturedCase::UNIFORM:
#if (NDIM == 2)
        U = { 0.35, -0.20 };
#endif
#if (NDIM == 3)
        U = { 0.35, -0.20, 0.15 };
#endif
        break;
    case ManufacturedCase::RIGID_ROTATION:
    {
        const double omega = 1.1;
        const double x = X[0] - 0.5;
        const double y = X[1] - 0.5;
        U[0] = -omega * y;
        U[1] = omega * x;
#if (NDIM == 3)
        U[2] = 0.0;
#endif
        break;
    }
    case ManufacturedCase::SIMPLE_SHEAR:
    {
        const double gamma = 0.8;
        U[0] = gamma * (X[1] - 0.5);
        U[1] = 0.0;
#if (NDIM == 3)
        U[2] = 0.0;
#endif
        break;
    }
    case ManufacturedCase::TRIGONOMETRIC:
    {
        const double qx = 2.0 * PI * X[0];
        const double qy = 2.0 * PI * X[1];
#if (NDIM == 2)
        U[0] = std::sin(qx) * std::cos(qy);
        U[1] = std::cos(qx) * std::sin(qy);
#endif
#if (NDIM == 3)
        const double qz = 2.0 * PI * X[2];
        U[0] = std::sin(qx) * std::cos(qy) * std::cos(qz);
        U[1] = std::cos(qx) * std::sin(qy) * std::cos(qz);
        U[2] = std::cos(qx) * std::cos(qy) * std::sin(qz);
#endif
        break;
    }
    }
    return U;
}

template <class IndexType>
Vector
compute_point(const IndexType& idx,
              const int axis,
              const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
              const hier::Index<NDIM>& patch_lower)
{
    Vector X = {};
    const double* const x_lower = patch_geom->getXLower();
    const double* const dx = patch_geom->getDx();
    for (int d = 0; d < NDIM; ++d)
    {
        double offset = 0.5;
        if constexpr (std::is_same<IndexType, CellIndex<NDIM>>::value)
        {
            offset = 0.5;
        }
        else if constexpr (std::is_same<IndexType, NodeIndex<NDIM>>::value)
        {
            offset = 0.0;
        }
        else
        {
            offset = d == axis ? 0.0 : 0.5;
        }
        X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_lower(d)) + offset);
    }
    return X;
}

double
evaluate_side_velocity(const ManufacturedCase manufactured_case,
                       const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                       const hier::Index<NDIM>& patch_lower,
                       const SideIndex<NDIM>& side_idx)
{
    const Vector X = compute_point(side_idx, side_idx.getAxis(), patch_geom, patch_lower);
    const Vector U = evaluate_velocity(manufactured_case, X);
    return U[side_idx.getAxis()];
}

SymTensor
compute_exact_cell_strain(const ManufacturedCase manufactured_case,
                          const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                          const hier::Index<NDIM>& patch_lower,
                          const CellIndex<NDIM>& idx)
{
    const double* const dx = patch_geom->getDx();
    SymTensor S = {};

#if (NDIM == 2)
    const CellIndex<NDIM> idx_x_minus(idx - axis_shift(0, 1));
    const CellIndex<NDIM> idx_x_plus(idx + axis_shift(0, 1));
    const CellIndex<NDIM> idx_y_minus(idx - axis_shift(1, 1));
    const CellIndex<NDIM> idx_y_plus(idx + axis_shift(1, 1));

    const double du0_dx0 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 0, SideIndex<NDIM>::Lower))) /
        dx[0];
    const double du1_dx1 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 1, SideIndex<NDIM>::Lower))) /
        dx[1];
    const double du0_dx1 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_plus, 0, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_plus, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_minus, 0, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_minus, 0, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[1]);
    const double du1_dx0 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_plus, 1, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_plus, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_minus, 1, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_minus, 1, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[0]);

    S[0] = du0_dx0;
    S[1] = du1_dx1;
    S[2] = 0.5 * (du0_dx1 + du1_dx0);
#endif

#if (NDIM == 3)
    const CellIndex<NDIM> idx_x_minus(idx - axis_shift(0, 1));
    const CellIndex<NDIM> idx_x_plus(idx + axis_shift(0, 1));
    const CellIndex<NDIM> idx_y_minus(idx - axis_shift(1, 1));
    const CellIndex<NDIM> idx_y_plus(idx + axis_shift(1, 1));
    const CellIndex<NDIM> idx_z_minus(idx - axis_shift(2, 1));
    const CellIndex<NDIM> idx_z_plus(idx + axis_shift(2, 1));

    const double du0_dx0 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 0, SideIndex<NDIM>::Lower))) /
        dx[0];
    const double du1_dx1 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 1, SideIndex<NDIM>::Lower))) /
        dx[1];
    const double du2_dx2 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 2, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx, 2, SideIndex<NDIM>::Lower))) /
        dx[2];
    const double du1_dx2 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_plus, 1, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_plus, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_minus, 1, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_minus, 1, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[2]);
    const double du2_dx1 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_plus, 2, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_plus, 2, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_minus, 2, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_minus, 2, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[1]);
    const double du0_dx2 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_plus, 0, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_plus, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_minus, 0, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_z_minus, 0, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[2]);
    const double du2_dx0 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_plus, 2, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_plus, 2, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_minus, 2, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_minus, 2, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[0]);
    const double du0_dx1 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_plus, 0, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_plus, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_minus, 0, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_y_minus, 0, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[1]);
    const double du1_dx0 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_plus, 1, SideIndex<NDIM>::Lower)) +
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_plus, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_minus, 1, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(idx_x_minus, 1, SideIndex<NDIM>::Upper))) /
        (4.0 * dx[0]);

    S[0] = du0_dx0;
    S[1] = du1_dx1;
    S[2] = du2_dx2;
    S[3] = 0.5 * (du1_dx2 + du2_dx1);
    S[4] = 0.5 * (du0_dx2 + du2_dx0);
    S[5] = 0.5 * (du0_dx1 + du1_dx0);
#endif
    return S;
}

#if (NDIM == 2)
SymTensor
compute_exact_node_strain(const ManufacturedCase manufactured_case,
                          const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                          const hier::Index<NDIM>& patch_lower,
                          const NodeIndex<NDIM>& idx)
{
    const double* const dx = patch_geom->getDx();
    const CellIndex<NDIM> c(idx);
    const CellIndex<NDIM> c_x_minus(idx - axis_shift(0, 1));
    const CellIndex<NDIM> c_y_minus(idx - axis_shift(1, 1));
    SymTensor S = {};

    const double du0_dx0_upper =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower))) /
        dx[0];
    const double du0_dx0_lower =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_y_minus, 0, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_y_minus, 0, SideIndex<NDIM>::Lower))) /
        dx[0];
    const double du1_dx1_right =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower))) /
        dx[1];
    const double du1_dx1_left =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_x_minus, 1, SideIndex<NDIM>::Upper)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_x_minus, 1, SideIndex<NDIM>::Lower))) /
        dx[1];
    const double du0_dx1 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_y_minus, 0, SideIndex<NDIM>::Lower))) /
        dx[1];
    const double du1_dx0 =
        (evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
         evaluate_side_velocity(
             manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_x_minus, 1, SideIndex<NDIM>::Lower))) /
        dx[0];

    S[0] = 0.5 * (du0_dx0_upper + du0_dx0_lower);
    S[1] = 0.5 * (du1_dx1_right + du1_dx1_left);
    S[2] = 0.5 * (du0_dx1 + du1_dx0);
    return S;
}
#endif

#if (NDIM == 3)
SymTensor
compute_exact_edge_strain(const ManufacturedCase manufactured_case,
                          const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                          const hier::Index<NDIM>& patch_lower,
                          const int axis,
                          const EdgeIndex<NDIM>& e_idx)
{
    const double* const dx = patch_geom->getDx();
    const hier::Index<NDIM>& i = e_idx;
    const CellIndex<NDIM> c(i);
    const CellIndex<NDIM> c_xm(i - axis_shift(0, 1));
    const CellIndex<NDIM> c_xp(i + axis_shift(0, 1));
    const CellIndex<NDIM> c_ym(i - axis_shift(1, 1));
    const CellIndex<NDIM> c_yp(i + axis_shift(1, 1));
    const CellIndex<NDIM> c_zm(i - axis_shift(2, 1));
    const CellIndex<NDIM> c_zp(i + axis_shift(2, 1));
    const CellIndex<NDIM> c_xym(i - axis_shift(0, 1) - axis_shift(1, 1));
    const CellIndex<NDIM> c_xzm(i - axis_shift(0, 1) - axis_shift(2, 1));
    const CellIndex<NDIM> c_yzm(i - axis_shift(1, 1) - axis_shift(2, 1));
    SymTensor S = {};

    if (axis == 0)
    {
        const double du0_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[0]);
        const double du0_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[1]);
        const double du0_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[2]);
        const double du1_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xp, 1, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(manufactured_case,
                                    patch_geom,
                                    patch_lower,
                                    SideIndex<NDIM>(c_xp - axis_shift(2, 1), 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[0]);
        const double du1_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[1]);
        const double du1_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower))) /
            dx[2];
        const double du2_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xp, 2, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(manufactured_case,
                                    patch_geom,
                                    patch_lower,
                                    SideIndex<NDIM>(c_xp - axis_shift(1, 1), 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[0]);
        const double du2_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower))) /
            dx[1];
        const double du2_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[2]);

        S[0] = du0_dx0;
        S[1] = du1_dx1;
        S[2] = du2_dx2;
        S[3] = 0.5 * (du1_dx2 + du2_dx1);
        S[4] = 0.5 * (du0_dx2 + du2_dx0);
        S[5] = 0.5 * (du0_dx1 + du1_dx0);
    }
    else if (axis == 1)
    {
        const double du0_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[0]);
        const double du0_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yp, 0, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(manufactured_case,
                                    patch_geom,
                                    patch_lower,
                                    SideIndex<NDIM>(c_yp - axis_shift(2, 1), 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[1]);
        const double du0_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower))) /
            dx[2];
        const double du1_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[0]);
        const double du1_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xp, 1, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(manufactured_case,
                                    patch_geom,
                                    patch_lower,
                                    SideIndex<NDIM>(c_xp - axis_shift(2, 1), 1, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[1]);
        const double du1_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[2]);
        const double du2_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xp, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[0]);
        const double du2_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[1]);
        const double du2_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[2]);

        S[0] = du0_dx0;
        S[1] = du1_dx1;
        S[2] = du2_dx2;
        S[3] = 0.5 * (du1_dx2 + du2_dx1);
        S[4] = 0.5 * (du0_dx2 + du2_dx0);
        S[5] = 0.5 * (du0_dx1 + du1_dx0);
    }
    else
    {
        const double du0_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[0]);
        const double du0_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 0, SideIndex<NDIM>::Lower))) /
            dx[1];
        const double du0_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zp, 0, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(manufactured_case,
                                    patch_geom,
                                    patch_lower,
                                    SideIndex<NDIM>(c_zp - axis_shift(1, 1), 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 0, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_yzm, 0, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[2]);
        const double du1_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower))) /
            dx[0];
        const double du1_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 1, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[1]);
        const double du1_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zp, 1, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(manufactured_case,
                                    patch_geom,
                                    patch_lower,
                                    SideIndex<NDIM>(c_zp - axis_shift(0, 1), 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_zm, 1, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xzm, 1, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[2]);
        const double du2_dx0 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[0]);
        const double du2_dx1 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
            (2.0 * dx[1]);
        const double du2_dx2 =
            (evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Upper)) +
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Upper)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xm, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_ym, 2, SideIndex<NDIM>::Lower)) -
             evaluate_side_velocity(
                 manufactured_case, patch_geom, patch_lower, SideIndex<NDIM>(c_xym, 2, SideIndex<NDIM>::Lower))) /
            (4.0 * dx[2]);

        S[0] = du0_dx0;
        S[1] = du1_dx1;
        S[2] = du2_dx2;
        S[3] = 0.5 * (du1_dx2 + du2_dx1);
        S[4] = 0.5 * (du0_dx2 + du2_dx0);
        S[5] = 0.5 * (du0_dx1 + du1_dx0);
    }

    return S;
}
#endif

double
compute_strain_magnitude(const SymTensor& S)
{
#if (NDIM == 2)
    return std::sqrt(2.0 * (S[0] * S[0] + S[1] * S[1] + 2.0 * S[2] * S[2]));
#endif
#if (NDIM == 3)
    return std::sqrt(2.0 * (S[0] * S[0] + S[1] * S[1] + S[2] * S[2] + 2.0 * (S[3] * S[3] + S[4] * S[4] + S[5] * S[5])));
#endif
}

double
compute_turbulent_viscosity(const SymTensor& S,
                            const SmagorinskyParameters& params,
                            const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom)
{
    const double* const dx = patch_geom->getDx();
    double delta_vol = 1.0;
    for (int d = 0; d < NDIM; ++d) delta_vol *= dx[d];
    const double delta = params.filter_width_scale * std::pow(delta_vol, 1.0 / static_cast<double>(NDIM));
    const double cs_delta_sq = std::pow(params.smagorinsky_constant * delta, 2.0);
    return params.rho * std::min(cs_delta_sq * compute_strain_magnitude(S), params.max_turbulent_viscosity);
}

double
compute_exact_tau_diag(const ManufacturedCase manufactured_case,
                       const SmagorinskyParameters& params,
                       const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                       const hier::Index<NDIM>& patch_lower,
                       const CellIndex<NDIM>& idx,
                       const int comp)
{
    const SymTensor S = compute_exact_cell_strain(manufactured_case, patch_geom, patch_lower, idx);
    return 2.0 * compute_turbulent_viscosity(S, params, patch_geom) * S[comp];
}

#if (NDIM == 2)
double
compute_exact_tau_shear(const ManufacturedCase manufactured_case,
                        const SmagorinskyParameters& params,
                        const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                        const hier::Index<NDIM>& patch_lower,
                        const NodeIndex<NDIM>& idx)
{
    const SymTensor S = compute_exact_node_strain(manufactured_case, patch_geom, patch_lower, idx);
    return 2.0 * compute_turbulent_viscosity(S, params, patch_geom) * S[2];
}
#endif

#if (NDIM == 3)
double
compute_exact_tau_shear(const ManufacturedCase manufactured_case,
                        const SmagorinskyParameters& params,
                        const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                        const hier::Index<NDIM>& patch_lower,
                        const int axis,
                        const EdgeIndex<NDIM>& idx)
{
    const SymTensor S = compute_exact_edge_strain(manufactured_case, patch_geom, patch_lower, axis, idx);
    const int shear_comp = axis == 0 ? 3 : (axis == 1 ? 4 : 5);
    return 2.0 * compute_turbulent_viscosity(S, params, patch_geom) * S[shear_comp];
}
#endif

double
compute_exact_force(const ManufacturedCase manufactured_case,
                    const SmagorinskyParameters& params,
                    const Pointer<CartesianPatchGeometry<NDIM>>& patch_geom,
                    const hier::Index<NDIM>& patch_lower,
                    const SideIndex<NDIM>& s_idx)
{
    const double* const dx = patch_geom->getDx();

#if (NDIM == 2)
    if (s_idx.getAxis() == 0)
    {
        const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
        const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
        const NodeIndex<NDIM> n_l(c_u, NodeIndex<NDIM>::LowerLeft);
        const NodeIndex<NDIM> n_u(c_u, NodeIndex<NDIM>::UpperLeft);
        return (compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_u, 0) -
                compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_l, 0)) /
                   dx[0] +
               (compute_exact_tau_shear(manufactured_case, params, patch_geom, patch_lower, n_u) -
                compute_exact_tau_shear(manufactured_case, params, patch_geom, patch_lower, n_l)) /
                   dx[1];
    }
    else
    {
        const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
        const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
        const NodeIndex<NDIM> n_l(c_u, NodeIndex<NDIM>::LowerLeft);
        const NodeIndex<NDIM> n_u(c_u, NodeIndex<NDIM>::LowerRight);
        return (compute_exact_tau_shear(manufactured_case, params, patch_geom, patch_lower, n_u) -
                compute_exact_tau_shear(manufactured_case, params, patch_geom, patch_lower, n_l)) /
                   dx[0] +
               (compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_u, 1) -
                compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_l, 1)) /
                   dx[1];
    }
#endif

#if (NDIM == 3)
    if (s_idx.getAxis() == 0)
    {
        const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
        const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
        const hier::Index<NDIM>& i = c_u;
        return (compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_u, 0) -
                compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_l, 0)) /
                   dx[0] +
               (compute_exact_tau_shear(manufactured_case,
                                        params,
                                        patch_geom,
                                        patch_lower,
                                        2,
                                        EdgeIndex<NDIM>(i + axis_shift(1, 1), 2, 0)) -
                compute_exact_tau_shear(
                    manufactured_case, params, patch_geom, patch_lower, 2, EdgeIndex<NDIM>(i, 2, 0))) /
                   dx[1] +
               (compute_exact_tau_shear(manufactured_case,
                                        params,
                                        patch_geom,
                                        patch_lower,
                                        1,
                                        EdgeIndex<NDIM>(i + axis_shift(2, 1), 1, 0)) -
                compute_exact_tau_shear(
                    manufactured_case, params, patch_geom, patch_lower, 1, EdgeIndex<NDIM>(i, 1, 0))) /
                   dx[2];
    }
    else if (s_idx.getAxis() == 1)
    {
        const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
        const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
        const hier::Index<NDIM>& i = c_u;
        return (compute_exact_tau_shear(manufactured_case,
                                        params,
                                        patch_geom,
                                        patch_lower,
                                        2,
                                        EdgeIndex<NDIM>(i + axis_shift(0, 1), 2, 0)) -
                compute_exact_tau_shear(
                    manufactured_case, params, patch_geom, patch_lower, 2, EdgeIndex<NDIM>(i, 2, 0))) /
                   dx[0] +
               (compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_u, 1) -
                compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_l, 1)) /
                   dx[1] +
               (compute_exact_tau_shear(manufactured_case,
                                        params,
                                        patch_geom,
                                        patch_lower,
                                        0,
                                        EdgeIndex<NDIM>(i + axis_shift(2, 1), 0, 0)) -
                compute_exact_tau_shear(
                    manufactured_case, params, patch_geom, patch_lower, 0, EdgeIndex<NDIM>(i, 0, 0))) /
                   dx[2];
    }
    else
    {
        const CellIndex<NDIM> c_l = s_idx.toCell(SideIndex<NDIM>::Lower);
        const CellIndex<NDIM> c_u = s_idx.toCell(SideIndex<NDIM>::Upper);
        const hier::Index<NDIM>& i = c_u;
        return (compute_exact_tau_shear(manufactured_case,
                                        params,
                                        patch_geom,
                                        patch_lower,
                                        1,
                                        EdgeIndex<NDIM>(i + axis_shift(0, 1), 1, 0)) -
                compute_exact_tau_shear(
                    manufactured_case, params, patch_geom, patch_lower, 1, EdgeIndex<NDIM>(i, 1, 0))) /
                   dx[0] +
               (compute_exact_tau_shear(manufactured_case,
                                        params,
                                        patch_geom,
                                        patch_lower,
                                        0,
                                        EdgeIndex<NDIM>(i + axis_shift(1, 1), 0, 0)) -
                compute_exact_tau_shear(
                    manufactured_case, params, patch_geom, patch_lower, 0, EdgeIndex<NDIM>(i, 0, 0))) /
                   dx[1] +
               (compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_u, 2) -
                compute_exact_tau_diag(manufactured_case, params, patch_geom, patch_lower, c_l, 2)) /
                   dx[2];
    }
#endif
}

template <class IndexType>
bool
is_interior_with_width(const hier::Box<NDIM>& patch_box, const IndexType& idx, const int width)
{
    for (int d = 0; d < NDIM; ++d)
    {
        if (idx(d) <= patch_box.lower()(d) + width - 1 || idx(d) >= patch_box.upper()(d) - width + 1) return false;
    }
    return true;
}

void
fill_side_velocity(const int U_idx,
                   const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                   const ManufacturedCase manufactured_case,
                   const double data_time)
{
    allocate_patch_data(U_idx, data_time, patch_hierarchy);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(U_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM> s_idx = si();
                    (*U_data)(s_idx) = evaluate_side_velocity(manufactured_case, patch_geom, patch_lower, s_idx);
                }
            }
        }
    }
}

double
verify_cell_strain(const ManufacturedCase manufactured_case,
                   const Pointer<INSSGSKinematics>& kinematics,
                   const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy)
{
    double max_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> S_data =
                patch->getPatchData(kinematics->getCellCenteredStrainRatePatchDataIndex());
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();
            const Box<NDIM> S_box = S_data->getGhostBox();
            for (CellIterator<NDIM> ci(S_box); ci; ci++)
            {
                const CellIndex<NDIM> idx = ci();
                if (!is_interior_with_width(patch->getBox(), idx, 1)) continue;
                const SymTensor S_exact = compute_exact_cell_strain(manufactured_case, patch_geom, patch_lower, idx);
                for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                {
                    max_error = std::max(max_error, std::abs((*S_data)(idx, comp) - S_exact[comp]));
                }
            }
        }
    }
    pout << "cell strain error: " << max_error << "\n";
    return max_error;
}

#if (NDIM == 2)
double
verify_native_shear_strain(const ManufacturedCase manufactured_case,
                           const Pointer<INSSGSKinematics>& kinematics,
                           const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy)
{
    double max_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> S_data =
                patch->getPatchData(kinematics->getNodeCenteredStrainRatePatchDataIndex());
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();
            const Box<NDIM> S_box = S_data->getGhostBox();
            for (NodeIterator<NDIM> ni(S_box); ni; ni++)
            {
                const NodeIndex<NDIM> idx = ni();
                if (!is_interior_with_width(patch->getBox(), idx, 1)) continue;
                const SymTensor S_exact = compute_exact_node_strain(manufactured_case, patch_geom, patch_lower, idx);
                for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                {
                    max_error = std::max(max_error, std::abs((*S_data)(idx, comp) - S_exact[comp]));
                }
            }
        }
    }
    pout << "node strain error: " << max_error << "\n";
    return max_error;
}
#endif

#if (NDIM == 3)
double
verify_native_shear_strain(const ManufacturedCase manufactured_case,
                           const Pointer<INSSGSKinematics>& kinematics,
                           const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy)
{
    double max_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<EdgeData<NDIM, double>> S_data =
                patch->getPatchData(kinematics->getEdgeCenteredStrainRatePatchDataIndex());
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const ArrayData<NDIM, double>& S_array = S_data->getArrayData(axis);
                const Box<NDIM> S_box = S_array.getBox();
                for (EdgeIterator<NDIM> ei(S_box, axis); ei; ei++)
                {
                    const EdgeIndex<NDIM> idx = ei();
                    if (!is_interior_with_width(patch->getBox(), idx, 1)) continue;
                    const SymTensor S_exact =
                        compute_exact_edge_strain(manufactured_case, patch_geom, patch_lower, axis, idx);
                    for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                    {
                        max_error = std::max(max_error, std::abs(S_array(idx, comp) - S_exact[comp]));
                    }
                }
            }
        }
    }
    pout << "edge strain error: " << max_error << "\n";
    return max_error;
}
#endif

double
verify_turbulence_force(const ManufacturedCase manufactured_case,
                        const SmagorinskyParameters& params,
                        const int F_idx,
                        const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy)
{
    double max_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> F_data = patch->getPatchData(F_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM> s_idx = si();
                    if (!is_interior_with_width(patch->getBox(), s_idx, 2)) continue;
                    const double F_exact =
                        compute_exact_force(manufactured_case, params, patch_geom, patch_lower, s_idx);
                    max_error = std::max(max_error, std::abs((*F_data)(s_idx)-F_exact));
                }
            }
        }
    }
    pout << "SGS force error: " << max_error << "\n";
    return max_error;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> verification_db = app_initializer->getComponentDatabase("ManufacturedSGS");

    const ManufacturedCase manufactured_case = parse_case(verification_db->getString("case_name"));
    const bool compact_output = verification_db->getBoolWithDefault("compact_output", false);
    const double verification_tol = verification_db->getDoubleWithDefault("verification_tol", 1.0e-10);
    const double data_time = verification_db->getDoubleWithDefault("data_time", 0.0);

    Pointer<INSStaggeredHierarchyIntegrator> time_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
        "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
    time_integrator->registerVelocityInitialConditions(u_init);
    Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
        "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
    time_integrator->registerPressureInitialConditions(p_init);

    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
    std::vector<Pointer<muParserRobinBcCoefs>> u_bc_coef_storage(NDIM);
    if (periodic_shift.min() <= 0)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            const std::string bc_name = "u_bc_coefs_" + std::to_string(d);
            const std::string db_name = "VelocityBcCoefs_" + std::to_string(d);
            u_bc_coef_storage[d] =
                new muParserRobinBcCoefs(bc_name, app_initializer->getComponentDatabase(db_name), grid_geometry);
            u_bc_coefs[d] = u_bc_coef_storage[d];
        }
        time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
    }

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    Pointer<INSTurbulenceModel> turbulence_model = time_integrator->getTurbulenceModel();
    if (turbulence_model.isNull())
    {
        TBOX_ERROR("Expected INSStaggeredHierarchyIntegrator to construct a SMAGORINSKY turbulence model.\n");
    }

    const Pointer<Database> model_db =
        app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator")->getDatabase("TurbulenceModel");
    SmagorinskyParameters params;
    params.rho = time_integrator->getStokesSpecifications()->getRho();
    params.smagorinsky_constant = model_db->getDoubleWithDefault("smagorinsky_constant", params.smagorinsky_constant);
    params.filter_width_scale = model_db->getDoubleWithDefault("filter_width_scale", params.filter_width_scale);
    params.max_turbulent_viscosity =
        model_db->getDoubleWithDefault("max_turbulent_viscosity", params.max_turbulent_viscosity);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> U_var = time_integrator->getVelocityVariable();
    const int U_idx = var_db->mapVariableAndContextToIndex(U_var, time_integrator->getCurrentContext());
    fill_side_velocity(U_idx, patch_hierarchy, manufactured_case, data_time);

    Pointer<INSSGSKinematics> kinematics = new INSSGSKinematics("EX8::Kinematics");
    kinematics->fillGhostedVelocity(
        U_idx, time_integrator->getVelocityBoundaryConditions(), patch_hierarchy, data_time);
    kinematics->computeCellCenteredStrainRate(patch_hierarchy, data_time);
#if (NDIM == 2)
    kinematics->computeNodeCenteredStrainRate(patch_hierarchy, data_time);
#endif
#if (NDIM == 3)
    kinematics->computeEdgeCenteredStrainRate(patch_hierarchy, data_time);
#endif

    Pointer<VariableContext> ctx = var_db->getContext("EX8::ForceContext");
    Pointer<SideVariable<NDIM, double>> F_var = new SideVariable<NDIM, double>("F_sgs");
    const int F_idx = var_db->registerVariableAndContext(F_var, ctx, IntVector<NDIM>(0));
    allocate_patch_data(F_idx, data_time, patch_hierarchy);
    turbulence_model->computeTurbulenceForce(F_idx,
                                             F_var,
                                             U_idx,
                                             U_var,
                                             time_integrator->getVelocityBoundaryConditions(),
                                             patch_hierarchy,
                                             data_time,
                                             *time_integrator->getStokesSpecifications());

    const double cell_strain_error = verify_cell_strain(manufactured_case, kinematics, patch_hierarchy);
    const double native_strain_error = verify_native_shear_strain(manufactured_case, kinematics, patch_hierarchy);
    const double force_error = verify_turbulence_force(manufactured_case, params, F_idx, patch_hierarchy);
    const double max_error = std::max(cell_strain_error, std::max(native_strain_error, force_error));
    const bool success = max_error <= verification_tol;

    if (compact_output)
    {
        pout << "case = " << verification_db->getString("case_name") << "\n";
        pout << "verification = " << (success ? "PASS" : "FAIL") << "\n";
    }
    else
    {
        pout << "Maximum SGS verification error = " << max_error << "\n";
        pout << (success ? "Manufactured SGS verification passed.\n" : "Manufactured SGS verification failed.\n");
    }

    deallocate_patch_data(F_idx, patch_hierarchy);
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
