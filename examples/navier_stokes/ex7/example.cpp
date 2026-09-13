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
#include <ibtk/CartGridFunction.h>
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
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int TENSOR_DEPTH = NDIM * (NDIM + 1) / 2;
constexpr double PI = 3.141592653589793238462643383279502884;

using IntegerVector = std::array<int, NDIM>;
using Vector = std::array<double, NDIM>;
using SymTensor = std::array<double, TENSOR_DEPTH>;

struct HarmonicMode
{
    Vector amplitude = {};
    Vector wave_number = {};
    double omega = 0.0;
    double phase = 0.0;
};

struct PeriodicFlowModel
{
    double rho = std::numeric_limits<double>::quiet_NaN();
    double mu = std::numeric_limits<double>::quiet_NaN();
    double period = std::numeric_limits<double>::quiet_NaN();
    Vector x_lower = {};
    Vector x_upper = {};
    Vector mean_velocity = {};
    std::vector<HarmonicMode> steady_modes;
    std::vector<HarmonicMode> fluctuating_modes;
};

double
dot(const Vector& lhs, const Vector& rhs)
{
    double result = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        result += lhs[d] * rhs[d];
    }
    return result;
}

double
squared_norm(const Vector& X)
{
    return dot(X, X);
}

void
scale(Vector& X, const double factor)
{
    for (int d = 0; d < NDIM; ++d)
    {
        X[d] *= factor;
    }
}

void
add_scaled(Vector& dst, const Vector& src, const double factor)
{
    for (int d = 0; d < NDIM; ++d)
    {
        dst[d] += factor * src[d];
    }
}

#if (NDIM == 3)
Vector
cross(const Vector& lhs, const Vector& rhs)
{
    Vector result = {};
    result[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
    result[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
    result[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
    return result;
}
#endif

double
phase_angle(const HarmonicMode& mode, const Vector& X, const double time)
{
    return dot(mode.wave_number, X) + mode.omega * time + mode.phase;
}

double
map_to_period(const double period, double time)
{
    if (IBTK::abs_equal_eps(period, 0.0)) return 0.0;
    time = std::fmod(time, period);
    if (time < 0.0) time += period;
    if (IBTK::abs_equal_eps(time, period)) time = 0.0;
    return time;
}

bool
is_periodic_sample_time(const std::set<double>& snapshot_time_points,
                        const double period,
                        const double time,
                        const double tol = 1.0e-8)
{
    const double phase_time = map_to_period(period, time);
    auto it = snapshot_time_points.lower_bound(phase_time);
    if (it != snapshot_time_points.end() && IBTK::abs_equal_eps(*it, phase_time, tol)) return true;
    if (it != snapshot_time_points.begin())
    {
        --it;
        if (IBTK::abs_equal_eps(*it, phase_time, tol)) return true;
    }
    return false;
}

Vector
evaluate_mean_velocity(const PeriodicFlowModel& model, const Vector& X)
{
    Vector velocity = model.mean_velocity;
    for (const HarmonicMode& mode : model.steady_modes)
    {
        add_scaled(velocity, mode.amplitude, std::cos(phase_angle(mode, X, 0.0)));
    }
    return velocity;
}

Vector
evaluate_velocity(const PeriodicFlowModel& model, const Vector& X, const double time)
{
    Vector velocity = evaluate_mean_velocity(model, X);
    for (const HarmonicMode& mode : model.fluctuating_modes)
    {
        add_scaled(velocity, mode.amplitude, std::cos(phase_angle(mode, X, time)));
    }
    return velocity;
}

Vector
evaluate_time_derivative(const PeriodicFlowModel& model, const Vector& X, const double time)
{
    Vector time_derivative = {};
    for (const HarmonicMode& mode : model.fluctuating_modes)
    {
        add_scaled(time_derivative, mode.amplitude, -mode.omega * std::sin(phase_angle(mode, X, time)));
    }
    return time_derivative;
}

Vector
evaluate_laplacian(const PeriodicFlowModel& model, const Vector& X, const double time)
{
    Vector laplacian = {};
    for (const HarmonicMode& mode : model.steady_modes)
    {
        add_scaled(laplacian, mode.amplitude, -squared_norm(mode.wave_number) * std::cos(phase_angle(mode, X, 0.0)));
    }
    for (const HarmonicMode& mode : model.fluctuating_modes)
    {
        add_scaled(laplacian, mode.amplitude, -squared_norm(mode.wave_number) * std::cos(phase_angle(mode, X, time)));
    }
    return laplacian;
}

Vector
evaluate_convective_term(const PeriodicFlowModel& model, const Vector& X, const double time)
{
    const Vector velocity = evaluate_velocity(model, X, time);
    Vector convective_term = {};
    for (const HarmonicMode& mode : model.steady_modes)
    {
        add_scaled(
            convective_term, mode.amplitude, -std::sin(phase_angle(mode, X, 0.0)) * dot(velocity, mode.wave_number));
    }
    for (const HarmonicMode& mode : model.fluctuating_modes)
    {
        add_scaled(
            convective_term, mode.amplitude, -std::sin(phase_angle(mode, X, time)) * dot(velocity, mode.wave_number));
    }
    return convective_term;
}

Vector
evaluate_body_force(const PeriodicFlowModel& model, const Vector& X, const double time)
{
    const Vector time_derivative = evaluate_time_derivative(model, X, time);
    const Vector convective_term = evaluate_convective_term(model, X, time);
    const Vector laplacian = evaluate_laplacian(model, X, time);

    Vector body_force = {};
    for (int d = 0; d < NDIM; ++d)
    {
        body_force[d] = model.rho * (time_derivative[d] + convective_term[d]) - model.mu * laplacian[d];
    }
    return body_force;
}

SymTensor
evaluate_reynolds_stress(const PeriodicFlowModel& model)
{
    SymTensor reynolds_stress = {};
    for (const HarmonicMode& mode : model.fluctuating_modes)
    {
#if (NDIM == 2)
        reynolds_stress[0] += 0.5 * mode.amplitude[0] * mode.amplitude[0];
        reynolds_stress[1] += 0.5 * mode.amplitude[1] * mode.amplitude[1];
        reynolds_stress[2] += 0.5 * mode.amplitude[0] * mode.amplitude[1];
#endif
#if (NDIM == 3)
        reynolds_stress[0] += 0.5 * mode.amplitude[0] * mode.amplitude[0];
        reynolds_stress[1] += 0.5 * mode.amplitude[1] * mode.amplitude[1];
        reynolds_stress[2] += 0.5 * mode.amplitude[2] * mode.amplitude[2];
        reynolds_stress[3] += 0.5 * mode.amplitude[1] * mode.amplitude[2];
        reynolds_stress[4] += 0.5 * mode.amplitude[0] * mode.amplitude[2];
        reynolds_stress[5] += 0.5 * mode.amplitude[0] * mode.amplitude[1];
#endif
    }
    return reynolds_stress;
}

double
evaluate_tke(const PeriodicFlowModel& model)
{
    const SymTensor reynolds_stress = evaluate_reynolds_stress(model);
#if (NDIM == 2)
    return 0.5 * (reynolds_stress[0] + reynolds_stress[1]);
#endif
#if (NDIM == 3)
    return 0.5 * (reynolds_stress[0] + reynolds_stress[1] + reynolds_stress[2]);
#endif
}

template <class VectorEvaluator>
void
fill_side_data(const int data_idx,
               const Pointer<Patch<NDIM>>& patch,
               VectorEvaluator&& evaluator,
               const std::string& object_name)
{
    Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(data_idx);
    if (side_data.isNull())
    {
        TBOX_ERROR(object_name << ": expected side-centered destination data.\n");
    }

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
                X[d] =
                    x_lower[d] + dx[d] * (static_cast<double>(side_idx(d) - patch_lower(d)) + (axis == d ? 0.0 : 0.5));
            }
            const Vector values = evaluator(X);
            (*side_data)(side_idx) = values[axis];
        }
    }
}

class ExactPeriodicVelocityField : public CartGridFunction
{
public:
    ExactPeriodicVelocityField(std::string object_name, const PeriodicFlowModel& model)
        : d_object_name(std::move(object_name)), d_model(model)
    {
        return;
    }

    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatch(int data_idx,
                        Pointer<Variable<NDIM>> /*var*/,
                        Pointer<Patch<NDIM>> patch,
                        const double data_time,
                        const bool /*initial_time*/,
                        Pointer<PatchLevel<NDIM>> /*patch_level*/) override
    {
        fill_side_data(
            data_idx,
            patch,
            [this, data_time](const Vector& X) { return evaluate_velocity(d_model, X, data_time); },
            d_object_name);
    }

private:
    std::string d_object_name;
    PeriodicFlowModel d_model;
};

class ExactPeriodicBodyForce : public CartGridFunction
{
public:
    ExactPeriodicBodyForce(std::string object_name, const PeriodicFlowModel& model)
        : d_object_name(std::move(object_name)), d_model(model)
    {
        return;
    }

    bool isTimeDependent() const override
    {
        return true;
    }

    void setDataOnPatch(int data_idx,
                        Pointer<Variable<NDIM>> /*var*/,
                        Pointer<Patch<NDIM>> patch,
                        const double data_time,
                        const bool /*initial_time*/,
                        Pointer<PatchLevel<NDIM>> /*patch_level*/) override
    {
        fill_side_data(
            data_idx,
            patch,
            [this, data_time](const Vector& X) { return evaluate_body_force(d_model, X, data_time); },
            d_object_name);
    }

private:
    std::string d_object_name;
    PeriodicFlowModel d_model;
};

template <class DatabaseValueType>
DatabaseValueType get_value_with_default(const Pointer<Database>& input_db,
                                         const std::string& key,
                                         const DatabaseValueType& default_value);

template <>
double
get_value_with_default(const Pointer<Database>& input_db, const std::string& key, const double& default_value)
{
    return input_db->getDoubleWithDefault(key, default_value);
}

template <>
int
get_value_with_default(const Pointer<Database>& input_db, const std::string& key, const int& default_value)
{
    return input_db->getIntegerWithDefault(key, default_value);
}

template <class DefaultValueType>
DefaultValueType
read_component_value(const Pointer<Database>& input_db,
                     const std::string& key_prefix,
                     const int component,
                     const DefaultValueType& default_value)
{
    return get_value_with_default(input_db, key_prefix + "_" + std::to_string(component), default_value);
}

Vector
make_orthogonal_direction(const Vector& wave_number, const double amplitude, std::mt19937* generator = nullptr)
{
    Vector direction = {};
#if (NDIM == 2)
    direction[0] = -wave_number[1];
    direction[1] = wave_number[0];
#endif
#if (NDIM == 3)
    Vector trial = {};
    if (generator)
    {
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        for (int d = 0; d < NDIM; ++d)
        {
            trial[d] = dist(*generator);
        }
    }
    else
    {
        trial[0] = 1.0;
        trial[1] = 0.0;
        trial[2] = 0.0;
    }

    const double projection = dot(trial, wave_number) / squared_norm(wave_number);
    direction = trial;
    add_scaled(direction, wave_number, -projection);
    if (squared_norm(direction) < 1.0e-12)
    {
        trial = { { 0.0, 1.0, 0.0 } };
        const double fallback_projection = dot(trial, wave_number) / squared_norm(wave_number);
        direction = trial;
        add_scaled(direction, wave_number, -fallback_projection);
    }
#endif

    scale(direction, amplitude / std::sqrt(squared_norm(direction)));
    return direction;
}

HarmonicMode
make_mode(const PeriodicFlowModel& model,
          const IntegerVector& wave_numbers,
          const int temporal_harmonic,
          const double amplitude,
          const double phase,
          std::mt19937* generator = nullptr)
{
    HarmonicMode mode;
    for (int d = 0; d < NDIM; ++d)
    {
        const double domain_length = model.x_upper[d] - model.x_lower[d];
        mode.wave_number[d] = 2.0 * PI * static_cast<double>(wave_numbers[d]) / domain_length;
    }
    mode.omega = 2.0 * PI * static_cast<double>(temporal_harmonic) / model.period;
    mode.phase = phase;
    mode.amplitude = make_orthogonal_direction(mode.wave_number, amplitude, generator);
    return mode;
}

PeriodicFlowModel
build_periodic_flow_model(const Pointer<CartesianGridGeometry<NDIM>>& grid_geometry,
                          const Pointer<Database>& input_db,
                          const double rho,
                          const double mu)
{
    PeriodicFlowModel model;
    model.rho = rho;
    model.mu = mu;
    model.period = input_db->getDouble("period");
    for (int d = 0; d < NDIM; ++d)
    {
        model.x_lower[d] = grid_geometry->getXLower()[d];
        model.x_upper[d] = grid_geometry->getXUpper()[d];
    }

#if (NDIM == 2)
    const Vector default_mean_velocity = { { 0.30, -0.15 } };
    const double default_steady_amplitude = 0.12;
    const double default_fluctuation_amplitude = 0.10;
    const double default_noise_amplitude = 0.03;
    const int default_noise_num_modes = 2;
#endif
#if (NDIM == 3)
    const Vector default_mean_velocity = { { 0.25, -0.12, 0.08 } };
    const double default_steady_amplitude = 0.10;
    const double default_fluctuation_amplitude = 0.08;
    const double default_noise_amplitude = 0.02;
    const int default_noise_num_modes = 3;
#endif

    for (int d = 0; d < NDIM; ++d)
    {
        model.mean_velocity[d] = read_component_value(input_db, "mean_velocity", d, default_mean_velocity[d]);
    }

    const double steady_amplitude = input_db->getDoubleWithDefault("steady_mode_amplitude", default_steady_amplitude);
    const double fluctuation_amplitude =
        input_db->getDoubleWithDefault("fluctuation_amplitude", default_fluctuation_amplitude);
    const double noise_amplitude = input_db->getDoubleWithDefault("noise_amplitude", default_noise_amplitude);
    const int noise_num_modes = input_db->getIntegerWithDefault("noise_num_modes", default_noise_num_modes);
    const int noise_seed = input_db->getIntegerWithDefault("noise_seed", 12345);
    const int max_spatial_harmonic = input_db->getIntegerWithDefault("max_spatial_harmonic", 3);

#if (NDIM == 2)
    model.steady_modes.push_back(make_mode(model, { { 1, 1 } }, 0, steady_amplitude, 0.30));
    model.fluctuating_modes.push_back(make_mode(model, { { 2, 1 } }, 1, fluctuation_amplitude, -0.40));
    model.fluctuating_modes.push_back(make_mode(model, { { 1, 2 } }, 2, 0.8 * fluctuation_amplitude, 0.55));
#endif
#if (NDIM == 3)
    model.steady_modes.push_back(make_mode(model, { { 1, 1, 1 } }, 0, steady_amplitude, 0.30));
    model.fluctuating_modes.push_back(make_mode(model, { { 2, 1, 1 } }, 1, fluctuation_amplitude, -0.35));
    model.fluctuating_modes.push_back(make_mode(model, { { 1, 2, 1 } }, 2, 0.85 * fluctuation_amplitude, 0.45));
    model.fluctuating_modes.push_back(make_mode(model, { { 1, 1, 2 } }, 3, 0.70 * fluctuation_amplitude, -0.25));
#endif

    std::mt19937 generator(noise_seed);
    std::uniform_int_distribution<int> wave_number_dist(1, max_spatial_harmonic);
    std::uniform_real_distribution<double> phase_dist(-PI, PI);
    const double noise_mode_amplitude =
        noise_num_modes > 0 ? noise_amplitude / std::sqrt(static_cast<double>(noise_num_modes)) : 0.0;
    const int first_noise_harmonic = static_cast<int>(model.fluctuating_modes.size()) + 1;
    for (int mode = 0; mode < noise_num_modes; ++mode)
    {
        IntegerVector wave_numbers = {};
        for (int d = 0; d < NDIM; ++d)
        {
            wave_numbers[d] = wave_number_dist(generator);
        }
        model.fluctuating_modes.push_back(make_mode(
            model, wave_numbers, first_noise_harmonic + mode, noise_mode_amplitude, phase_dist(generator), &generator));
    }

    return model;
}

double
verify_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                       const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                       const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                       const PeriodicFlowModel& model,
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

    const SymTensor R_exact = evaluate_reynolds_stress(model);
    const double k_exact = evaluate_tke(model);
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

                const Vector U_exact = evaluate_mean_velocity(model, X);

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

    pout << "CELL running-mean errors: |<U>-exact|_max = " << max_mean_error
         << ", |R-exact|_max = " << max_reynolds_error << ", |k-exact|_max = " << max_tke_error << "\n";
    return std::max(max_mean_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                       const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                       const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                       const PeriodicFlowModel& model,
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

    const SymTensor R_exact = evaluate_reynolds_stress(model);
    const double k_exact = evaluate_tke(model);
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

                const Vector U_exact = evaluate_mean_velocity(model, X);

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

    pout << "NODE running-mean errors: |<U>-exact|_max = " << max_mean_error
         << ", |R-exact|_max = " << max_reynolds_error << ", |k-exact|_max = " << max_tke_error << "\n";
    return std::max(max_mean_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_periodic_phase_cell_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                      const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                      const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                      const PeriodicFlowModel& model)
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

                    const Vector U_exact = evaluate_velocity(model, X, snapshot_time);
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

    pout << "CELL phase-average errors: |U_phase-exact|_max = " << max_phase_error
         << ", |R|_max = " << max_reynolds_error << ", |k|_max = " << max_tke_error << "\n";
    return std::max(max_phase_error, std::max(max_reynolds_error, max_tke_error));
}

double
verify_periodic_phase_node_statistics(const Pointer<INSAveragingTurbulenceStatistics>& statistics,
                                      const Pointer<PatchHierarchy<NDIM>>& patch_hierarchy,
                                      const Pointer<IBTK::HierarchyMathOps>& hier_math_ops,
                                      const PeriodicFlowModel& model)
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

                    const Vector U_exact = evaluate_velocity(model, X, snapshot_time);
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

    pout << "NODE phase-average errors: |U_phase-exact|_max = " << max_phase_error
         << ", |R|_max = " << max_reynolds_error << ", |k|_max = " << max_tke_error << "\n";
    return std::max(max_phase_error, std::max(max_reynolds_error, max_tke_error));
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "turbulence_stats.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    const bool dump_viz_data = app_initializer->dumpVizData();
    const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
    const int viz_dump_interval = app_initializer->getVizDumpInterval();
    const bool dump_restart_data = app_initializer->dumpRestartData();
    const int restart_dump_interval = app_initializer->getRestartDumpInterval();
    const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
    const bool dump_timer_data = app_initializer->dumpTimerData();
    const int timer_dump_interval = app_initializer->getTimerDumpInterval();

    Pointer<INSStaggeredHierarchyIntegrator> time_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    if (grid_geometry->getPeriodicShift().min() <= 0)
    {
        TBOX_ERROR("examples/navier_stokes/ex7 requires a fully periodic domain.\n");
    }
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

    Pointer<Database> ins_db = app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator");
    const PeriodicFlowModel model = build_periodic_flow_model(grid_geometry,
                                                              app_initializer->getComponentDatabase("PeriodicFlow"),
                                                              ins_db->getDouble("rho"),
                                                              ins_db->getDouble("mu"));

    time_integrator->registerVelocityInitialConditions(new ExactPeriodicVelocityField("u_init", model));
    if (input_db->keyExists("PressureInitialConditions"))
    {
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);
    }
    time_integrator->registerBodyForceFunction(new ExactPeriodicBodyForce("f_fcn", model));

    Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
    if (uses_visit)
    {
        time_integrator->registerVisItDataWriter(visit_data_writer);
    }

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    Pointer<INSAveragingTurbulenceStatistics> statistics = time_integrator->getTurbulenceStatistics();
    if (statistics.isNull())
    {
        TBOX_ERROR(
            "Expected INSStaggeredHierarchyIntegrator to construct an AVERAGED_VELOCITY turbulence statistics "
            "object.\n");
    }
    const std::set<double>& snapshot_time_points = statistics->getAveragedVelocityManager().getSnapshotTimePoints();
    const bool periodic_mode = snapshot_time_points.size() > 1;

    Pointer<IBTK::HierarchyMathOps> hier_math_ops = new IBTK::HierarchyMathOps("HierarchyMathOps", patch_hierarchy);
    hier_math_ops->resetLevels(0, patch_hierarchy->getFinestLevelNumber());

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM, double>> U_var = time_integrator->getVelocityVariable();
    const int U_idx = var_db->mapVariableAndContextToIndex(U_var, time_integrator->getCurrentContext());

    // Record the exact initialized state so phase-averaged runs that start at t = 0
    // have a corresponding statistics snapshot.
    statistics->updateStatistics(U_idx,
                                 U_var,
                                 time_integrator->getVelocityBoundaryConditions(),
                                 time_integrator->getIntegratorTime(),
                                 patch_hierarchy,
                                 hier_math_ops);

    int iteration_num = time_integrator->getIntegratorStep();
    double loop_time = time_integrator->getIntegratorTime();
    if (dump_viz_data && uses_visit)
    {
        time_integrator->setupPlotData();
        visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
    }

    const double loop_time_end = time_integrator->getEndTime();
    double dt = 0.0;
    while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
    {
        iteration_num = time_integrator->getIntegratorStep();
        loop_time = time_integrator->getIntegratorTime();

        pout << "\n";
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        pout << "At beginning of timestep # " << iteration_num << "\n";
        pout << "Simulation time is " << loop_time << "\n";

        dt = time_integrator->getMaximumTimeStepSize();
        time_integrator->advanceHierarchy(dt);
        loop_time = time_integrator->getIntegratorTime();

        if (!periodic_mode || is_periodic_sample_time(snapshot_time_points, model.period, loop_time))
        {
            const int U_current_idx = var_db->mapVariableAndContextToIndex(U_var, time_integrator->getCurrentContext());
            statistics->updateStatistics(U_current_idx,
                                         U_var,
                                         time_integrator->getVelocityBoundaryConditions(),
                                         loop_time,
                                         patch_hierarchy,
                                         hier_math_ops);
        }

        pout << "\n";
        pout << "At end       of timestep # " << iteration_num << "\n";
        pout << "Simulation time is " << loop_time << "\n";
        pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        pout << "\n";

        iteration_num += 1;
        const bool last_step = !time_integrator->stepsRemaining();
        if (dump_viz_data && uses_visit && viz_dump_interval > 0 &&
            (iteration_num % viz_dump_interval == 0 || last_step))
        {
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }
        if (dump_restart_data && restart_dump_interval > 0 && (iteration_num % restart_dump_interval == 0 || last_step))
        {
            RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
        }
        if (dump_timer_data && timer_dump_interval > 0 && (iteration_num % timer_dump_interval == 0 || last_step))
        {
            TimerManager::getManager()->print(plog);
        }
    }

    const DataCentering data_centering = statistics->getDataCentering();
    double max_error = std::numeric_limits<double>::quiet_NaN();
    if (periodic_mode)
    {
        switch (data_centering)
        {
        case DataCentering::CELL:
            max_error = verify_periodic_phase_cell_statistics(statistics, patch_hierarchy, hier_math_ops, model);
            break;
        case DataCentering::NODE:
            max_error = verify_periodic_phase_node_statistics(statistics, patch_hierarchy, hier_math_ops, model);
            break;
        default:
            TBOX_ERROR("Unsupported analysis centering " << IBAMR::enum_to_string<IBAMR::DataCentering>(data_centering)
                                                         << "\n");
        }
        pout << "Periodic steady state flag = " << (statistics->isAtSteadyState() ? "TRUE" : "FALSE") << "\n";
    }
    else
    {
        const double snapshot_time = *(statistics->getAveragedVelocityManager().getSnapshotTimePoints().begin());
        switch (data_centering)
        {
        case DataCentering::CELL:
            max_error = verify_cell_statistics(statistics, patch_hierarchy, hier_math_ops, model, snapshot_time);
            break;
        case DataCentering::NODE:
            max_error = verify_node_statistics(statistics, patch_hierarchy, hier_math_ops, model, snapshot_time);
            break;
        default:
            TBOX_ERROR("Unsupported analysis centering " << IBAMR::enum_to_string<IBAMR::DataCentering>(data_centering)
                                                         << "\n");
        }
    }

    bool enforce_verification = false;
    double verification_tol = std::numeric_limits<double>::infinity();
    if (input_db->keyExists("StatisticsDiagnostics"))
    {
        Pointer<Database> diagnostics_db = app_initializer->getComponentDatabase("StatisticsDiagnostics");
        enforce_verification = diagnostics_db->getBoolWithDefault("enforce_verification", false);
        verification_tol = diagnostics_db->getDoubleWithDefault("verification_tol", verification_tol);
    }

    pout << "Maximum diagnostic error = " << max_error << "\n";
    if (!enforce_verification) return EXIT_SUCCESS;

    const bool success = max_error <= verification_tol;
    pout << (success ? "Turbulence-statistics diagnostics passed.\n" : "Turbulence-statistics diagnostics failed.\n");
    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
