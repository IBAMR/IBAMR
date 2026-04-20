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

#ifndef included_tests_navier_stokes_stokes_petsc_mat_utilities_side_transfer_test_helpers
#define included_tests_navier_stokes_stokes_petsc_mat_utilities_side_transfer_test_helpers

#include <ibtk/IndexUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <Patch.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>

#include <array>
#include <cmath>
#include <string>

#include <ibtk/app_namespaces.h>

namespace StokesPETScMatUtilitiesSideTransferTests
{
enum class ProfileType
{
    AFFINE,
    PIECEWISE_RT0,
    NONLINEAR,
    UNKNOWN
};

inline ProfileType
string_to_profile_type(const std::string& profile)
{
    ProfileType profile_type = ProfileType::UNKNOWN;
    if (profile == "affine") profile_type = ProfileType::AFFINE;
    if (profile == "piecewise_rt0") profile_type = ProfileType::PIECEWISE_RT0;
    if (profile == "nonlinear") profile_type = ProfileType::NONLINEAR;
    if (profile_type == ProfileType::UNKNOWN) TBOX_ERROR("Unknown profile_type = " << profile << "\n");
    return profile_type;
}

inline int
compute_piecewise_region(const double x)
{
    if (x < 1.0 / 3.0) return 0;
    if (x < 2.0 / 3.0) return 1;
    return 2;
}

inline double
compute_piecewise_rt0_value(const VectorNd& X, const int axis)
{
    int transverse_sum = 0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis) continue;
        transverse_sum += compute_piecewise_region(X[d]);
    }
    const int mod3 = transverse_sum % 3;
    const int mod5 = (2 * transverse_sum + axis) % 5;
    const double c0 = 0.2 * static_cast<double>(axis + 1) + 0.04 * static_cast<double>(mod3);
    const double c1 = 0.15 + 0.02 * static_cast<double>(mod5);
    return c0 + c1 * X[axis];
}

inline double
compute_nonlinear_value(const VectorNd& X, const int axis)
{
    const double pi = 3.14159265358979323846;
    double transverse_sum = 0.0;
    double transverse_prod = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (d == axis) continue;
        transverse_sum += X[d];
        transverse_prod *= (X[d] + 0.15 * static_cast<double>(d + 1));
    }
    return std::sin(2.0 * pi * X[axis]) + 0.25 * std::cos(pi * transverse_sum) + 0.12 * X[axis] * X[axis] +
           0.04 * transverse_prod;
}

inline void
set_affine_side_field(Pointer<SideData<NDIM, double>> u_data,
                      Pointer<Patch<NDIM>> patch,
                      const std::array<double, NDIM>& coeffs)
{
    const Box<NDIM>& patch_box = patch->getBox();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const VectorNd X = IBTK::IndexUtilities::getSideCenter(*patch, i_s);
            (*u_data)(i_s) = coeffs[axis] * X[axis];
        }
    }
}

inline void
set_piecewise_rt0_side_field(Pointer<SideData<NDIM, double>> u_data, Pointer<Patch<NDIM>> patch)
{
    const Box<NDIM>& patch_box = patch->getBox();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const VectorNd X = IBTK::IndexUtilities::getSideCenter(*patch, i_s);
            (*u_data)(i_s) = compute_piecewise_rt0_value(X, axis);
        }
    }
}

inline void
set_nonlinear_side_field(Pointer<SideData<NDIM, double>> u_data, Pointer<Patch<NDIM>> patch)
{
    const Box<NDIM>& patch_box = patch->getBox();
    for (int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
        for (Box<NDIM>::Iterator b(side_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
            const VectorNd X = IBTK::IndexUtilities::getSideCenter(*patch, i_s);
            (*u_data)(i_s) = compute_nonlinear_value(X, axis);
        }
    }
}

inline void
set_test_profile_side_field(Pointer<SideData<NDIM, double>> u_data,
                            Pointer<Patch<NDIM>> patch,
                            const ProfileType profile_type,
                            const std::array<double, NDIM>& coeffs)
{
    switch (profile_type)
    {
    case ProfileType::AFFINE:
        set_affine_side_field(u_data, patch, coeffs);
        break;
    case ProfileType::PIECEWISE_RT0:
        set_piecewise_rt0_side_field(u_data, patch);
        break;
    case ProfileType::NONLINEAR:
        set_nonlinear_side_field(u_data, patch);
        break;
    case ProfileType::UNKNOWN:
    default:
        TBOX_ERROR("Unknown ProfileType encountered.\n");
    }
}

inline void
check_nontrivial(const std::string& label,
                 const double samrai_max_norm,
                 const double petsc_max_norm,
                 const double tol,
                 int& test_failures)
{
    if (IBTK::abs_equal_eps(samrai_max_norm, 0.0, tol) || IBTK::abs_equal_eps(petsc_max_norm, 0.0, tol))
    {
        ++test_failures;
        pout << "nontriviality check failed for " << label << ": SAMRAI max norm = " << samrai_max_norm
             << ", PETSc max norm = " << petsc_max_norm << ", tolerance = " << tol << "\n";
    }
}

inline void
check_max_norm_consistency(const std::string& label,
                           const double samrai_max_norm,
                           const double petsc_max_norm,
                           const double tol,
                           int& test_failures)
{
    if (!IBTK::rel_equal_eps(samrai_max_norm, petsc_max_norm, tol))
    {
        ++test_failures;
        pout << "max-norm consistency check failed for " << label << ": SAMRAI max norm = " << samrai_max_norm
             << ", PETSc max norm = " << petsc_max_norm << ", tolerance = " << tol << "\n";
    }
}
} // namespace StokesPETScMatUtilitiesSideTransferTests

#endif
