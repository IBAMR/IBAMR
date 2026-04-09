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

// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSAveragingTurbulenceStatistics.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/snapshot_utilities.h>

// Set up application namespace declarations
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <set>
#include <string>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int TENSOR_DEPTH = NDIM * (NDIM + 1) / 2;
constexpr double PI = 3.141592653589793238462643383279502884;

using Vector = std::array<double, NDIM>;
using SymTensor = std::array<double, TENSOR_DEPTH>;
using LinearVectorCoefficients = std::array<std::array<double, NDIM + 1>, NDIM>;

struct ManufacturedVelocityCoefficients
{
    LinearVectorCoefficients mean;
    LinearVectorCoefficients cosine;
    LinearVectorCoefficients sine;
};

ManufacturedVelocityCoefficients
get_manufactured_velocity_coefficients()
{
#if (NDIM == 2)
    return { { { { 0.45, 0.0, 0.0 }, { -0.20, 0.0, 0.0 } } },
             { { { 0.70, 0.0, 0.0 }, { -0.35, 0.0, 0.0 } } },
             { { { -0.25, 0.0, 0.0 }, { 0.55, 0.0, 0.0 } } } };
#endif
#if (NDIM == 3)
    return { { { { 0.45, 0.0, 0.0, 0.0 }, { -0.20, 0.0, 0.0, 0.0 }, { 0.15, 0.0, 0.0, 0.0 } } },
             { { { 0.70, 0.0, 0.0, 0.0 }, { -0.35, 0.0, 0.0, 0.0 }, { 0.28, 0.0, 0.0, 0.0 } } },
             { { { -0.25, 0.0, 0.0, 0.0 }, { 0.55, 0.0, 0.0, 0.0 }, { -0.18, 0.0, 0.0, 0.0 } } } };
#endif
}

Vector
evaluate_linear_vector(const LinearVectorCoefficients& coeffs, const Vector& X)
{
    Vector values = {};
    for (int comp = 0; comp < NDIM; ++comp)
    {
        values[comp] = coeffs[comp][0];
        for (int d = 0; d < NDIM; ++d)
        {
            values[comp] += coeffs[comp][d + 1] * X[d];
        }
    }
    return values;
}

Vector
evaluate_mean_velocity(const ManufacturedVelocityCoefficients& coeffs, const Vector& X)
{
    return evaluate_linear_vector(coeffs.mean, X);
}

Vector
evaluate_velocity(const ManufacturedVelocityCoefficients& coeffs,
                  const Vector& X,
                  const double time,
                  const double period)
{
    const double theta = 2.0 * PI * time / period;
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const Vector mean = evaluate_linear_vector(coeffs.mean, X);
    const Vector cosine = evaluate_linear_vector(coeffs.cosine, X);
    const Vector sine = evaluate_linear_vector(coeffs.sine, X);

    Vector velocity = {};
    for (int comp = 0; comp < NDIM; ++comp)
    {
        velocity[comp] = mean[comp] + cosine[comp] * c + sine[comp] * s;
    }
    return velocity;
}

SymTensor
evaluate_reynolds_stress(const ManufacturedVelocityCoefficients& coeffs, const Vector& X)
{
    const Vector cosine = evaluate_linear_vector(coeffs.cosine, X);
    const Vector sine = evaluate_linear_vector(coeffs.sine, X);
    SymTensor R = {};

#if (NDIM == 2)
    R[0] = 0.5 * (cosine[0] * cosine[0] + sine[0] * sine[0]);
    R[1] = 0.5 * (cosine[1] * cosine[1] + sine[1] * sine[1]);
    R[2] = 0.5 * (cosine[0] * cosine[1] + sine[0] * sine[1]);
#endif

#if (NDIM == 3)
    R[0] = 0.5 * (cosine[0] * cosine[0] + sine[0] * sine[0]);
    R[1] = 0.5 * (cosine[1] * cosine[1] + sine[1] * sine[1]);
    R[2] = 0.5 * (cosine[2] * cosine[2] + sine[2] * sine[2]);
    R[3] = 0.5 * (cosine[1] * cosine[2] + sine[1] * sine[2]);
    R[4] = 0.5 * (cosine[0] * cosine[2] + sine[0] * sine[2]);
    R[5] = 0.5 * (cosine[0] * cosine[1] + sine[0] * sine[1]);
#endif

    return R;
}

double
evaluate_tke(const ManufacturedVelocityCoefficients& coeffs, const Vector& X)
{
    const SymTensor R = evaluate_reynolds_stress(coeffs, X);
#if (NDIM == 2)
    return 0.5 * (R[0] + R[1]);
#endif
#if (NDIM == 3)
    return 0.5 * (R[0] + R[1] + R[2]);
#endif
}

void
fill_side_velocity(const int U_idx,
                   const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                   const ManufacturedVelocityCoefficients& coeffs,
                   const double data_time,
                   const double period)
{
    IBTK::allocate_patch_data(U_idx, data_time, patch_hierarchy);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> U_data = patch->getPatchData(U_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const dx = patch_geom->getDx();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    const SideIndex<NDIM> side_idx = si();
                    Vector X = {};
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] +
                               dx[d] * (static_cast<double>(side_idx(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5));
                    }
                    const Vector velocity = evaluate_velocity(coeffs, X, data_time, period);
                    (*U_data)(side_idx) = velocity[axis];
                }
            }
        }
    }
}

double
verify_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                       const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                       const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                       const ManufacturedVelocityCoefficients& coeffs,
                       const double snapshot_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::VerificationContext");
    Pointer<CellVariable<NDIM, double>> U_mean_var = new CellVariable<NDIM, double>("U_mean_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> R_var = new CellVariable<NDIM, double>("R_cc", TENSOR_DEPTH);
    const int U_mean_idx = var_db->registerVariableAndContext(U_mean_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));

    IBTK::allocate_patch_data({ U_mean_idx, R_idx }, snapshot_time, patch_hierarchy);

    IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                     U_mean_idx,
                                     snapshot_time,
                                     patch_hierarchy,
                                     "CONSERVATIVE_LINEAR_REFINE");
    statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);

    double max_mean_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> U_mean_data = patch->getPatchData(U_mean_idx);
            Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const dx = patch_geom->getDx();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                Vector X = {};
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_lower(d)) + 0.5);
                }

                const Vector U_exact = evaluate_mean_velocity(coeffs, X);
                const SymTensor R_exact = evaluate_reynolds_stress(coeffs, X);
                const double k_exact = evaluate_tke(coeffs, X);

                for (int comp = 0; comp < NDIM; ++comp)
                {
                    max_mean_error = std::max(max_mean_error, std::abs((*U_mean_data)(idx, comp) - U_exact[comp]));
                }
                for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                {
                    max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp) - R_exact[comp]));
                }

#if (NDIM == 2)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                max_tke_error = std::max(max_tke_error, std::abs(k_num - k_exact));
            }
        }
    }

    IBTK::deallocate_patch_data({ U_mean_idx, R_idx }, patch_hierarchy);

    pout << "CELL statistics errors: |<U>-exact|_max = " << max_mean_error << ", |R-exact|_max = " << max_reynolds_error
         << ", |k-exact|_max = " << max_tke_error << "\n";
    return std::max(max_mean_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                       const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                       const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                       const ManufacturedVelocityCoefficients& coeffs,
                       const double snapshot_time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::VerificationContext");
    Pointer<NodeVariable<NDIM, double>> U_mean_var = new NodeVariable<NDIM, double>("U_mean_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> R_var = new NodeVariable<NDIM, double>("R_nc", TENSOR_DEPTH, false);
    const int U_mean_idx = var_db->registerVariableAndContext(U_mean_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));

    IBTK::allocate_patch_data({ U_mean_idx, R_idx }, snapshot_time, patch_hierarchy);

    IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                     U_mean_idx,
                                     snapshot_time,
                                     patch_hierarchy,
                                     "LINEAR_REFINE");
    statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);

    double max_mean_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> U_mean_data = patch->getPatchData(U_mean_idx);
            Pointer<NodeData<NDIM, double>> R_data = patch->getPatchData(R_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const x_lower = patch_geom->getXLower();
            const double* const dx = patch_geom->getDx();
            const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

            for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
            {
                const NodeIndex<NDIM>& idx = ni();
                Vector X = {};
                for (int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * static_cast<double>(idx(d) - patch_lower(d));
                }

                const Vector U_exact = evaluate_mean_velocity(coeffs, X);
                const SymTensor R_exact = evaluate_reynolds_stress(coeffs, X);
                const double k_exact = evaluate_tke(coeffs, X);

                for (int comp = 0; comp < NDIM; ++comp)
                {
                    max_mean_error = std::max(max_mean_error, std::abs((*U_mean_data)(idx, comp) - U_exact[comp]));
                }
                for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                {
                    max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp) - R_exact[comp]));
                }

#if (NDIM == 2)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                max_tke_error = std::max(max_tke_error, std::abs(k_num - k_exact));
            }
        }
    }

    IBTK::deallocate_patch_data({ U_mean_idx, R_idx }, patch_hierarchy);

    pout << "NODE statistics errors: |<U>-exact|_max = " << max_mean_error << ", |R-exact|_max = " << max_reynolds_error
         << ", |k-exact|_max = " << max_tke_error << "\n";
    return std::max(max_mean_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_periodic_phase_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                      const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                      const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                      const ManufacturedVelocityCoefficients& coeffs,
                                      const double sample_period)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::PeriodicVerificationContext");
    Pointer<CellVariable<NDIM, double>> U_phase_var = new CellVariable<NDIM, double>("U_phase_cc", NDIM);
    Pointer<CellVariable<NDIM, double>> R_var = new CellVariable<NDIM, double>("R_phase_cc", TENSOR_DEPTH);
    const int U_phase_idx = var_db->registerVariableAndContext(U_phase_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));

    double max_phase_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (const double snapshot_time : statistics->getAveragedVelocityManager().getSnapshotTimePoints())
    {
        IBTK::allocate_patch_data({ U_phase_idx, R_idx }, snapshot_time, patch_hierarchy);

        IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                         U_phase_idx,
                                         snapshot_time,
                                         patch_hierarchy,
                                         "CONSERVATIVE_LINEAR_REFINE");
        statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> U_phase_data = patch->getPatchData(U_phase_idx);
                Pointer<CellData<NDIM, double>> R_data = patch->getPatchData(R_idx);
                Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                const double* const x_lower = patch_geom->getXLower();
                const double* const dx = patch_geom->getDx();
                const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    Vector X = {};
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * (static_cast<double>(idx(d) - patch_lower(d)) + 0.5);
                    }

                    const Vector U_exact = evaluate_velocity(coeffs, X, snapshot_time, sample_period);
                    for (int comp = 0; comp < NDIM; ++comp)
                    {
                        max_phase_error =
                            std::max(max_phase_error, std::abs((*U_phase_data)(idx, comp) - U_exact[comp]));
                    }
                    for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                    {
                        max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp)));
                    }
#if (NDIM == 2)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                    max_tke_error = std::max(max_tke_error, std::abs(k_num));
                }
            }
        }

        IBTK::deallocate_patch_data({ U_phase_idx, R_idx }, patch_hierarchy);
    }

    pout << "CELL periodic-phase errors: |U_phase-exact|_max = " << max_phase_error
         << ", |R|_max = " << max_reynolds_error << ", |k|_max = " << max_tke_error << "\n";
    return std::max(max_phase_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_periodic_phase_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                      const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                      const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                      const ManufacturedVelocityCoefficients& coeffs,
                                      const double sample_period)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("EX7::PeriodicVerificationContext");
    Pointer<NodeVariable<NDIM, double>> U_phase_var = new NodeVariable<NDIM, double>("U_phase_nc", NDIM, false);
    Pointer<NodeVariable<NDIM, double>> R_var = new NodeVariable<NDIM, double>("R_phase_nc", TENSOR_DEPTH, false);
    const int U_phase_idx = var_db->registerVariableAndContext(U_phase_var, ctx, IntVector<NDIM>(1));
    const int R_idx = var_db->registerVariableAndContext(R_var, ctx, IntVector<NDIM>(0));

    double max_phase_error = 0.0;
    double max_reynolds_error = 0.0;
    double max_tke_error = 0.0;
    for (const double snapshot_time : statistics->getAveragedVelocityManager().getSnapshotTimePoints())
    {
        IBTK::allocate_patch_data({ U_phase_idx, R_idx }, snapshot_time, patch_hierarchy);

        IBTK::fill_snapshot_on_hierarchy(statistics->getAveragedVelocityManager().getSnapshotCache(),
                                         U_phase_idx,
                                         snapshot_time,
                                         patch_hierarchy,
                                         "LINEAR_REFINE");
        statistics->fillReynoldsStressSnapshot(R_idx, R_var, snapshot_time, patch_hierarchy, hier_math_ops);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<NodeData<NDIM, double>> U_phase_data = patch->getPatchData(U_phase_idx);
                Pointer<NodeData<NDIM, double>> R_data = patch->getPatchData(R_idx);
                Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                const double* const x_lower = patch_geom->getXLower();
                const double* const dx = patch_geom->getDx();
                const hier::Index<NDIM>& patch_lower = patch->getBox().lower();

                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = ni();
                    Vector X = {};
                    for (int d = 0; d < NDIM; ++d)
                    {
                        X[d] = x_lower[d] + dx[d] * static_cast<double>(idx(d) - patch_lower(d));
                    }

                    const Vector U_exact = evaluate_velocity(coeffs, X, snapshot_time, sample_period);
                    for (int comp = 0; comp < NDIM; ++comp)
                    {
                        max_phase_error =
                            std::max(max_phase_error, std::abs((*U_phase_data)(idx, comp) - U_exact[comp]));
                    }
                    for (int comp = 0; comp < TENSOR_DEPTH; ++comp)
                    {
                        max_reynolds_error = std::max(max_reynolds_error, std::abs((*R_data)(idx, comp)));
                    }
#if (NDIM == 2)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1));
#endif
#if (NDIM == 3)
                    const double k_num = 0.5 * ((*R_data)(idx, 0) + (*R_data)(idx, 1) + (*R_data)(idx, 2));
#endif
                    max_tke_error = std::max(max_tke_error, std::abs(k_num));
                }
            }
        }

        IBTK::deallocate_patch_data({ U_phase_idx, R_idx }, patch_hierarchy);
    }

    pout << "NODE periodic-phase errors: |U_phase-exact|_max = " << max_phase_error
         << ", |R|_max = " << max_reynolds_error << ", |k|_max = " << max_tke_error << "\n";
    return std::max(max_phase_error, std::max(max_reynolds_error, max_tke_error));
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

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

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    Pointer<INSAveragingTurbulenceStatistics> statistics = time_integrator->getTurbulenceStatistics();
    if (statistics.isNull())
    {
        TBOX_ERROR(
            "Expected INSStaggeredHierarchyIntegrator to construct an AVERAGED_VELOCITY turbulence statistics "
            "object.\n");
    }

    Pointer<Database> verification_db = app_initializer->getComponentDatabase("ManufacturedStatistics");
    const double sample_start_time = verification_db->getDoubleWithDefault("sample_start_time", 0.0);
    const double sample_period = verification_db->getDouble("sample_period");
    const int num_samples = verification_db->getInteger("num_samples");
    const int num_periods = verification_db->getIntegerWithDefault("num_periods", 1);
    const bool periodic_mode = verification_db->getBoolWithDefault("periodic_mode", false);
    const bool compact_output = verification_db->getBoolWithDefault("compact_output", false);
    const double verification_tol = verification_db->getDoubleWithDefault("verification_tol", 1.0e-10);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> U_var = time_integrator->getVelocityVariable();
    const int U_idx = var_db->mapVariableAndContextToIndex(U_var, time_integrator->getCurrentContext());

    Pointer<IBTK::HierarchyMathOps> hier_math_ops = new IBTK::HierarchyMathOps("HierarchyMathOps", patch_hierarchy);
    hier_math_ops->resetLevels(0, patch_hierarchy->getFinestLevelNumber());

    const ManufacturedVelocityCoefficients coeffs = get_manufactured_velocity_coefficients();
    for (int period = 0; period < num_periods; ++period)
    {
        for (int sample = 0; sample < num_samples; ++sample)
        {
            const double phase_time =
                sample_start_time + sample_period * static_cast<double>(sample) / static_cast<double>(num_samples);
            const double data_time =
                periodic_mode ? phase_time + static_cast<double>(period) * sample_period : phase_time;
            fill_side_velocity(U_idx, patch_hierarchy, coeffs, data_time, sample_period);
            statistics->updateStatistics(U_idx,
                                         U_var,
                                         time_integrator->getVelocityBoundaryConditions(),
                                         data_time,
                                         patch_hierarchy,
                                         hier_math_ops);
        }
    }

    const DataCentering data_centering = statistics->getDataCentering();
    double max_error = std::numeric_limits<double>::quiet_NaN();
    if (periodic_mode)
    {
        switch (data_centering)
        {
        case DataCentering::CELL:
            max_error = verify_periodic_phase_cell_statistics(
                statistics, patch_hierarchy, hier_math_ops, coeffs, sample_period);
            break;
        case DataCentering::NODE:
            max_error = verify_periodic_phase_node_statistics(
                statistics, patch_hierarchy, hier_math_ops, coeffs, sample_period);
            break;
        default:
            TBOX_ERROR("Unsupported analysis centering " << enum_to_string<DataCentering>(data_centering) << "\n");
        }
    }
    else
    {
        const double snapshot_time = *(statistics->getAveragedVelocityManager().getSnapshotTimePoints().begin());
        switch (data_centering)
        {
        case DataCentering::CELL:
            max_error = verify_cell_statistics(statistics, patch_hierarchy, hier_math_ops, coeffs, snapshot_time);
            break;
        case DataCentering::NODE:
            max_error = verify_node_statistics(statistics, patch_hierarchy, hier_math_ops, coeffs, snapshot_time);
            break;
        default:
            TBOX_ERROR("Unsupported analysis centering " << enum_to_string<DataCentering>(data_centering) << "\n");
        }
    }

    const bool success = max_error <= verification_tol;
    if (compact_output)
    {
        pout << "mode = " << (periodic_mode ? "PERIODIC_PHASE" : "RUNNING_MEAN") << "\n";
        pout << "analysis_centering = " << enum_to_string<DataCentering>(data_centering) << "\n";
        pout << "verification = " << (success ? "PASS" : "FAIL") << "\n";
    }
    else
    {
        pout << "Maximum verification error = " << max_error << "\n";
        pout << (success ? "Manufactured turbulence statistics verification passed.\n" :
                           "Manufactured turbulence statistics verification failed.\n");
    }
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
