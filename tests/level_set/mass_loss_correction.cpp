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

#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/MemoryDatabase.h>
#include <tbox/RestartManager.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <tuple>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
enum class Geometry
{
    PLANE,
    CURVE,
    CORNER_CURVE,
    STEEP,
    NEGATIVE,
    POSITIVE,
    CUTOFF,
    NEAR_CUTOFF
};

double
field_value(const Geometry field, const double x, const double y, const double alpha, const double length)
{
    switch (field)
    {
    case Geometry::PLANE:
        return x - 0.5 * length;
    case Geometry::CURVE:
        return std::hypot(x - 0.5 * length, y - 0.5 * length) - 0.2 * length;
    case Geometry::CORNER_CURVE:
        return std::hypot(x - 0.25 * length, y - 0.25 * length) - 0.125 * length;
    case Geometry::STEEP:
        return 4.0 * (x - 0.5 * length);
    case Geometry::NEGATIVE:
        return -length;
    case Geometry::POSITIVE:
        return length;
    case Geometry::CUTOFF:
        return -alpha;
    case Geometry::NEAR_CUTOFF:
        return -alpha * (1.0 - 1.0e-6);
    default:
        TBOX_ERROR("Unknown correction-test geometry\n");
    }
    return 0.0;
}

void
fill_correction_fields(Pointer<PatchHierarchy<NDIM>> hierarchy,
                       const std::array<int, 3>& indices,
                       const Geometry field,
                       const Geometry solid,
                       const double length)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> phi = patch->getPatchData(indices[0]);
            Pointer<CellData<NDIM, double>> psi = patch->getPatchData(indices[1]);
            Pointer<CellData<NDIM, double>> other = patch->getPatchData(indices[2]);
            Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
            const Box<NDIM>& box = patch->getBox();
            double dv = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                dv *= geom->getDx()[d];
            }
            const double alpha = std::pow(dv, 1.0 / NDIM);
            for (Box<NDIM>::Iterator i(box); i; i++)
            {
                const double x = geom->getXLower()[0] + (i()(0) - box.lower(0) + 0.5) * geom->getDx()[0];
                const double y = geom->getXLower()[1] + (i()(1) - box.lower(1) + 0.5) * geom->getDx()[1];
                (*phi)(i()) = field_value(field, x, y, alpha, length);
                (*psi)(i()) = field_value(solid, y, x, alpha, length);
                (*other)(i()) = 17.0;
            }
        }
    }
}

void
check_algebraic_corrections(Pointer<AdvDiffHierarchyIntegrator> integrator,
                            const std::vector<Pointer<CellVariable<NDIM, double>>>& variables,
                            Pointer<Database> input,
                            const double length)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::array<int, 3> current, next, original;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, finest_ln);
    for (int k = 0; k < 3; ++k)
    {
        current[k] = var_db->mapVariableAndContextToIndex(variables[k], integrator->getCurrentContext());
        next[k] = var_db->mapVariableAndContextToIndex(variables[k], integrator->getNewContext());
        original[k] = var_db->registerClonedPatchDataIndex(variables[k], current[k]);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(next[k], 0.1);
            hierarchy->getPatchLevel(ln)->allocatePatchData(original[k], 0.0);
        }
    }
    const double domain_volume = std::pow(length, NDIM);
    const double eps = std::numeric_limits<double>::epsilon();
    int cells = 0;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            cells += level->getPatch(p())->getBox().size();
        }
    }
    const double accumulation_eps = (2.0 * IBTK_MPI::sumReduction(cells) + 8.0) * eps;
    const double identity_tolerance = accumulation_eps / (1.0 - accumulation_eps) * domain_volume;
    const std::string failure = input->getStringWithDefault("failure", "");
    const bool exact_initial_success = input->getBoolWithDefault("exact_initial_success", false);
    // The third coordinate is extruded for these algebraic cases.
    using Case = std::tuple<Geometry, Geometry, double>;
    const Case cases[] = {
        { Geometry::NEGATIVE, Geometry::POSITIVE, 1.0 },   { Geometry::POSITIVE, Geometry::POSITIVE, 0.0 },
        { Geometry::NEGATIVE, Geometry::POSITIVE, 0.0 },   { Geometry::POSITIVE, Geometry::POSITIVE, 1.0 },
        { Geometry::PLANE, Geometry::NEGATIVE, 0.0 },      { Geometry::PLANE, Geometry::CUTOFF, 0.0 },
        { Geometry::PLANE, Geometry::NEAR_CUTOFF, 0.5 },   { Geometry::STEEP, Geometry::POSITIVE, 0.5 },
        { Geometry::STEEP, Geometry::POSITIVE, 0.6 },      { Geometry::PLANE, Geometry::POSITIVE, 0.4 },
        { Geometry::PLANE, Geometry::POSITIVE, 0.6 },      { Geometry::CURVE, Geometry::POSITIVE, 0.8 },
        { Geometry::CURVE, Geometry::POSITIVE, 0.0001 },   { Geometry::CURVE, Geometry::PLANE, 0.4 },
        { Geometry::CURVE, Geometry::PLANE, 0.6 },         { Geometry::CORNER_CURVE, Geometry::POSITIVE, 0.04 },
        { Geometry::PLANE, Geometry::POSITIVE, -1.0e-15 }, { Geometry::PLANE, Geometry::POSITIVE, 1.0 + 1.0e-15 }
    };
    int case_number = 0;
    for (const bool three_phase : { false, true })
    {
        std::vector<Pointer<CellVariable<NDIM, double>>> fields{ variables[0] };
        if (three_phase)
        {
            fields.push_back(variables[1]);
        }
        for (const Case& test : cases)
        {
            const Geometry field = failure.empty() ? std::get<0>(test) : Geometry::CURVE;
            const Geometry solid = failure.empty() ? std::get<1>(test) : Geometry::POSITIVE;
            fill_correction_fields(hierarchy, current, field, solid, length);
            for (int k = 0; k < 3; ++k)
            {
                data_ops.copyData(next[k], current[k]);
                data_ops.copyData(original[k], current[k]);
            }
            LevelSetUtilities::LevelSetContainer container(integrator, fields);
            const std::vector<double> before = three_phase ?
                                                   LevelSetUtilities::computeHeavisideIntegrals3PhaseFlows(container) :
                                                   LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(container);
            const int phase = three_phase ? 1 : 0;
            const double capacity = before[0] + before[1];
            double target = std::get<2>(test) * capacity;
            if (!three_phase && (solid == Geometry::NEGATIVE || solid == Geometry::CUTOFF))
            {
                target = before[phase];
            }
            Pointer<Database> controls = new MemoryDatabase("FixerInput");
            if (exact_initial_success)
            {
                controls->putInteger("max_its", 1);
                controls->putDouble("rel_tol", 0.0);
                controls->putDouble("abs_tol", 0.0);
                target = before[phase];
            }
            if (input->keyExists("abs_tol"))
            {
                controls->putDouble("abs_tol", input->getDouble("abs_tol"));
            }
            if (failure == "interval_zero" || failure == "interval_negative")
            {
                controls->putInteger("correction_interval", failure == "interval_zero" ? 0 : -1);
            }
            if (failure == "width_input")
            {
                controls->putDouble("half_width", 0.0);
            }
            if (failure == "relative_tolerance")
            {
                controls->putDouble("rel_tol", 1.0);
            }
            if (failure == "absolute_tolerance")
            {
                controls->putDouble("abs_tol", -1.0);
            }
            if (failure == "budget_zero" || failure == "budget_exhausted")
            {
                controls->putInteger("max_its", failure == "budget_zero" ? 0 : 1);
                target = 0.8 * capacity;
            }
            if (failure == "width_restart" || failure == "target_restart")
            {
                {
                    LevelSetUtilities::LevelSetMassLossFixer writer("Fixer", integrator, fields, controls);
                    writer.setInitialVolume(before[phase]);
                    if (failure == "width_restart")
                    {
                        writer.getLevelSetContainer().setInterfaceHalfWidth(-1.0);
                    }
                    else
                    {
                        writer.setTargetVolume(std::numeric_limits<double>::quiet_NaN());
                    }
                    RestartManager::getManager()->writeRestartFile("invalid_restart", 1);
                }
                if (!RestartManager::getManager()->openRestartFile("invalid_restart", 1, IBTK_MPI::getNodes()))
                {
                    TBOX_ERROR("Could not open invalid-state restart fixture\n");
                }
            }
            LevelSetUtilities::LevelSetMassLossFixer fixer("Fixer", integrator, fields, controls, false);
            fixer.setInitialVolume(before[phase]);
            if (failure == "target_nonfinite")
            {
                target = std::numeric_limits<double>::quiet_NaN();
            }
            else if (failure == "target_infeasible")
            {
                target = 2.0 * capacity;
            }
            if (failure != "target_restart")
            {
                fixer.setTargetVolume(target);
            }
            if (failure == "width_setter")
            {
                fixer.getLevelSetContainer().setInterfaceHalfWidth(-1.0);
            }
            if (failure == "rank_nonfinite" && IBTK_MPI::getRank() == IBTK_MPI::getNodes() - 1)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
                PatchLevel<NDIM>::Iterator p(level);
                TBOX_ASSERT(p);
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> phi = patch->getPatchData(next[0]);
                (*phi)(patch->getBox().lower()) = std::numeric_limits<double>::quiet_NaN();
            }
            if (three_phase)
            {
                LevelSetUtilities::fixMassLoss3PhaseFlows(0.0, 0.1, false, 1, &fixer);
            }
            else
            {
                LevelSetUtilities::fixMassLoss2PhaseFlows(0.0, 0.1, false, 1, &fixer);
            }
            if (!failure.empty())
            {
                // Unexpected continuation must be reported as success to the
                // test-local parent, which requires the real fatal diagnostic.
                break;
            }
            if (fixer.getTargetVolume() != target)
            {
                TBOX_ERROR("Correction regression: stored application target changed\n");
            }
            const double q = fixer.getLagrangeMultiplier();
            double shift_error = 0.0;
            int local_support = 0;
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
                    double cell_volume = 1.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        cell_volume *= geom->getDx()[d];
                    }
                    const double alpha = std::pow(cell_volume, 1.0 / NDIM);
                    for (int k = 0; k < 3; ++k)
                    {
                        Pointer<CellData<NDIM, double>> old_data = patch->getPatchData(original[k]);
                        Pointer<CellData<NDIM, double>> current_data = patch->getPatchData(current[k]);
                        Pointer<CellData<NDIM, double>> new_data = patch->getPatchData(next[k]);
                        for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                        {
                            const double value = (*new_data)(i());
                            if ((*current_data)(i()) != (*old_data)(i()) || (k != 0 && value != (*old_data)(i())) ||
                                !std::isfinite(value))
                            {
                                TBOX_ERROR("Correction regression: changed CURRENT, solid, or unrelated data\n");
                            }
                            if (k == 0)
                            {
                                local_support += smooth_delta((*old_data)(i()), alpha) > 0.0;
                                shift_error = std::max(shift_error, std::abs((value - (*old_data)(i())) - q));
                            }
                        }
                    }
                }
            }
            if (NDIM == 2 && IBTK_MPI::getNodes() == 2 && finest_ln == 0 && field == Geometry::CORNER_CURVE)
            {
                const int minimum_support = IBTK_MPI::minReduction(local_support);
                const int maximum_support = IBTK_MPI::maxReduction(local_support);
                if (minimum_support != 0 || maximum_support == 0)
                {
                    TBOX_ERROR("Correction regression: missing rank without local interface support\n");
                }
            }
            shift_error = IBTK_MPI::maxReduction(shift_error);
            if (!(shift_error <= 8.0 * eps * std::max(length, std::abs(q))) || fixer.getTime() != 0.1)
            {
                TBOX_ERROR("Correction regression: nonuniform shift or incorrect correction time\n");
            }
            data_ops.copyData(current[0], next[0]);
            const std::vector<double> after = three_phase ?
                                                  LevelSetUtilities::computeHeavisideIntegrals3PhaseFlows(container) :
                                                  LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(container);
            if (!(std::abs(after[0] + after[1] - capacity) <= identity_tolerance) ||
                (three_phase && !(std::abs(after[0] + after[1] + after[2] - domain_volume) <= identity_tolerance)))
            {
                TBOX_ERROR("Correction regression: composite phase-sum identity\n");
            }
            const double tolerance =
                input->getDoubleWithDefault("abs_tol", 64.0 * eps * capacity) + 1.0e-12 * std::abs(target);
            if (!(std::abs(after[phase] - target) <= (exact_initial_success ? 0.0 : tolerance)))
            {
                TBOX_ERROR("Correction regression: phase volume missed the specified tolerance\n");
            }
            if (std::abs(before[phase] - target) <= tolerance && q != 0.0)
            {
                TBOX_ERROR("Correction regression: initial success changed the field\n");
            }
            if (std::abs(before[phase] - target) > tolerance &&
                !((three_phase ? q : -q) * (target - before[phase]) > 0.0))
            {
                TBOX_ERROR("Correction regression: wrong shift sign\n");
            }
            if (field == Geometry::PLANE && solid == Geometry::POSITIVE &&
                (std::get<2>(test) == 0.4 || std::get<2>(test) == 0.6) && finest_ln == 0)
            {
                // For the symmetric plane, the continuum location agrees with
                // the discrete regularization to within one cell width.
                const double expected_q = (three_phase ? 1.0 : -1.0) * (target / capacity - 0.5) * length;
                if (!(std::abs(q - expected_q) <= length / 16.0))
                {
                    TBOX_ERROR("Correction regression: planar interface location\n");
                }
            }
            plog << "case = " << case_number++ << "; phases = " << (three_phase ? 3 : 2)
                 << "; target, volume, shift: " << target / domain_volume << ' ' << after[phase] / domain_volume << ' '
                 << q / length << '\n';
        }
        if (!failure.empty())
        {
            break;
        }
    }
    for (int k = 0; k < 3; ++k)
    {
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(next[k]);
            hierarchy->getPatchLevel(ln)->deallocatePatchData(original[k]);
        }
        var_db->removePatchDataIndex(original[k]);
    }
}
void
check_lifecycle_callback(const double current_time,
                         const double new_time,
                         const bool skip_synchronize,
                         const int num_cycles,
                         void* ctx)
{
    std::pair<LevelSetUtilities::LevelSetMassLossFixer*, bool>* callback =
        static_cast<std::pair<LevelSetUtilities::LevelSetMassLossFixer*, bool>*>(ctx);
    Pointer<AdvDiffHierarchyIntegrator> integrator =
        callback->first->getLevelSetContainer().getAdvDiffHierarchyIntegrator();
    Pointer<PatchHierarchy<NDIM>> hierarchy = integrator->getPatchHierarchy();
    const int idx = VariableDatabase<NDIM>::getDatabase()->mapVariableAndContextToIndex(
        callback->first->getLevelSetContainer().getLevelSetVariable(), integrator->getNewContext());
    const bool skipped = integrator->getIntegratorStep() % callback->first->getCorrectionInterval() != 0;
    std::vector<double> before;
    if (skipped)
    {
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> data = patch->getPatchData(idx);
                for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                {
                    before.push_back((*data)(i()));
                }
            }
        }
    }
    if (callback->second)
    {
        LevelSetUtilities::fixMassLoss3PhaseFlows(
            current_time, new_time, skip_synchronize, num_cycles, callback->first);
    }
    else
    {
        LevelSetUtilities::fixMassLoss2PhaseFlows(
            current_time, new_time, skip_synchronize, num_cycles, callback->first);
    }
    if (skipped)
    {
        std::size_t j = 0;
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> data = patch->getPatchData(idx);
                for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                {
                    if ((*data)(i()) != before[j++])
                    {
                        TBOX_ERROR("Correction regression: skipped callback changed NEW data\n");
                    }
                }
            }
        }
    }
}

void
check_correction_lifecycle(Pointer<AdvDiffHierarchyIntegrator> integrator,
                           const std::vector<Pointer<CellVariable<NDIM, double>>>& variables,
                           Pointer<Database> input,
                           const double length)
{
    const bool three_phase = input->getInteger("phases") == 3;
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    const bool test_restart = input->getBoolWithDefault("test_restart", false);
    const double domain_volume = std::pow(length, NDIM);
    Pointer<PatchHierarchy<NDIM>> hierarchy = integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    std::array<int, 3> current;
    for (int k = 0; k < 3; ++k)
    {
        current[k] = var_db->mapVariableAndContextToIndex(variables[k], integrator->getCurrentContext());
    }
    if (!from_restart)
    {
        fill_correction_fields(hierarchy, current, Geometry::PLANE, Geometry::PLANE, length);
    }
    std::vector<Pointer<CellVariable<NDIM, double>>> fields{ variables[0] };
    if (three_phase)
    {
        fields.push_back(variables[1]);
    }
    LevelSetUtilities::LevelSetContainer container(integrator, fields);
    const auto volumes = [&]()
    {
        return three_phase ? LevelSetUtilities::computeHeavisideIntegrals3PhaseFlows(container) :
                             LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(container);
    };
    const int phase = three_phase ? 1 : 0;
    const double capacity = (three_phase ? 0.5 : 1.0) * domain_volume;
    Pointer<Database> controls = new MemoryDatabase("LifecycleInput");
    controls->putInteger("correction_interval", 2);
    controls->putInteger("max_its", from_restart ? 96 : 64);
    controls->putDouble("rel_tol", from_restart ? 5.0e-13 : 1.0e-12);
    controls->putDouble("half_width", from_restart ? 2.0 : 1.0);
    LevelSetUtilities::LevelSetMassLossFixer fixer("LifecycleFixer", integrator, fields, controls);
    fixer.setInitialVolume(volumes()[phase]);
    if (from_restart)
    {
        if (!std::isfinite(fixer.getInitialVolume()) || fixer.getInitialVolume() == fixer.getTargetVolume() ||
            fixer.getTargetVolume() != 0.4 * capacity || fixer.getMaxIterations() != 96 ||
            fixer.getErrorRelTolerance() != 5.0e-13 || fixer.getLevelSetContainer().getInterfaceHalfWidth() != 1.0)
        {
            TBOX_ERROR("Correction regression: restart/input precedence\n");
        }
    }
    std::pair<LevelSetUtilities::LevelSetMassLossFixer*, bool> callback(&fixer, three_phase);
    integrator->registerPostprocessIntegrateHierarchyCallback(check_lifecycle_callback, &callback);
    const int finest_ln = hierarchy->getFinestLevelNumber();
    while (integrator->getIntegratorStep() < 3)
    {
        const int step = integrator->getIntegratorStep();
        if (step == 0 || step == 2)
        {
            fixer.setTargetVolume((step == 0 ? 0.4 : 0.6) * capacity);
        }
        const double previous_q = fixer.getLagrangeMultiplier();
        const double previous_time = fixer.getTime();
        integrator->advanceHierarchy(0.01);
        const double volume = volumes()[phase];
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * capacity +
                                 fixer.getErrorRelTolerance() * std::abs(fixer.getTargetVolume());
        if (!(std::abs(volume - fixer.getTargetVolume()) <= tolerance))
        {
            TBOX_ERROR("Correction regression: CURRENT volume after synchronization and context swap\n");
        }
        if (step == 1)
        {
            if (fixer.getTime() != previous_time ||
                !(fixer.getLagrangeMultiplier() == previous_q ||
                  (std::isnan(previous_q) && std::isnan(fixer.getLagrangeMultiplier()))))
            {
                TBOX_ERROR("Correction regression: skipped event changed logging state\n");
            }
        }
        else if (fixer.getTime() != integrator->getIntegratorTime())
        {
            TBOX_ERROR("Correction regression: incorrect pre-increment correction schedule\n");
        }
        plog << "completed step = " << step + 1
             << "; normalized target, CURRENT volume: " << fixer.getTargetVolume() / domain_volume << ' '
             << volume / domain_volume << '\n';
        if (test_restart && !from_restart && step == 0)
        {
            RestartManager::getManager()->writeRestartFile("restart", 1);
        }
    }
    if (test_restart)
    {
        std::vector<double> values{ fixer.getInitialVolume() };
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int idx : current)
                {
                    Pointer<CellData<NDIM, double>> data = patch->getPatchData(idx);
                    for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                    {
                        values.push_back((*data)(i()));
                    }
                }
            }
        }
        const std::string filename = "uninterrupted-" + std::to_string(IBTK_MPI::getRank()) + ".dat";
        if (!from_restart)
        {
            std::ofstream reference(filename, std::ios::binary);
            reference.write(reinterpret_cast<const char*>(values.data()), values.size() * sizeof(double));
            if (!reference)
            {
                TBOX_ERROR("Cannot write uninterrupted correction reference\n");
            }
        }
        else
        {
            std::vector<double> reference(values.size());
            std::ifstream stream(filename, std::ios::binary);
            stream.read(reinterpret_cast<char*>(reference.data()), reference.size() * sizeof(double));
            if (!stream || stream.peek() != std::ifstream::traits_type::eof())
            {
                TBOX_ERROR("Invalid uninterrupted correction reference\n");
            }
            double error = 0.0;
            for (std::size_t i = 0; i < values.size(); ++i)
            {
                if (!std::isfinite(values[i]) || !std::isfinite(reference[i]))
                {
                    TBOX_ERROR("Correction regression: nonfinite restart comparison data\n");
                }
                error = std::max(error, std::abs(values[i] - reference[i]));
            }
            error = IBTK_MPI::maxReduction(error);
            if (!(error <= 128.0 * std::numeric_limits<double>::epsilon() * length))
            {
                TBOX_ERROR("Correction regression: restart differs from uninterrupted fields\n");
            }
            plog << "normalized restart field difference: " << error / length << '\n';
        }
    }
}

void
copy_transport_distance(int idx, Pointer<HierarchyMathOps> ops, double, bool, void* ctx)
{
    std::pair<int, int>* source = static_cast<std::pair<int, int>*>(ctx);
    ++source->second;
    Pointer<PatchHierarchy<NDIM>> hierarchy = ops->getPatchHierarchy();
    HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, hierarchy->getFinestLevelNumber());
    data_ops.copyData(idx, source->first);
}

void
reinitialize_transport(int idx,
                       Pointer<HierarchyMathOps> ops,
                       const int step,
                       const double time,
                       const bool initial_time,
                       const bool regrid_time,
                       void* ctx)
{
    if (regrid_time)
    {
        return;
    }
    RelaxationLSMethod* relaxation = static_cast<RelaxationLSMethod*>(ctx);
    relaxation->setReinitializeLSData(true);
    relaxation->initializeLSData(idx, ops, step, time, initial_time);
}

void
check_corrected_transport(Pointer<AdvDiffHierarchyIntegrator> integrator,
                          const std::vector<Pointer<CellVariable<NDIM, double>>>& variables,
                          LocationIndexRobinBcCoefs<NDIM>& bc)
{
    Pointer<PatchHierarchy<NDIM>> hierarchy = integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int finest_ln = hierarchy->getFinestLevelNumber();
    std::array<int, 3> current;
    for (int k = 0; k < 3; ++k)
    {
        current[k] = var_db->mapVariableAndContextToIndex(variables[k], integrator->getCurrentContext());
    }
    fill_correction_fields(hierarchy, current, Geometry::CURVE, Geometry::POSITIVE, 1.0);
    HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, finest_ln);
    data_ops.copyData(current[2], current[0]);
    std::array<std::pair<int, int>, 2> sources{ std::make_pair(current[0], 0), std::make_pair(current[2], 0) };
    std::array<Pointer<RelaxationLSMethod>, 2> relaxations;
    for (int k = 0; k < 2; ++k)
    {
        Pointer<Database> db = new MemoryDatabase("TransportRelaxationInput");
        db->putInteger("max_iterations", 2);
        db->putString("order", "THIRD_ORDER_ENO");
        db->putString("time_stepping_scheme", "TVD_RK2");
        relaxations[k] = new RelaxationLSMethod("TransportRelaxation" + std::to_string(k), db, false);
        relaxations[k]->registerInterfaceNeighborhoodLocatingFcn(copy_transport_distance, &sources[k]);
        relaxations[k]->registerPhysicalBoundaryCondition(&bc);
        integrator->registerResetFunction(variables[2 * k], reinitialize_transport, relaxations[k].getPointer());
    }
    LevelSetUtilities::LevelSetContainer corrected(integrator, variables[0]);
    LevelSetUtilities::LevelSetContainer control(integrator, variables[2]);
    const double target = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(corrected)[0];
    LevelSetUtilities::LevelSetMassLossFixer fixer("TransportFixer", integrator, { variables[0] }, nullptr, false);
    fixer.setInitialVolume(target);
    integrator->registerPostprocessIntegrateHierarchyCallback(LevelSetUtilities::fixMassLoss2PhaseFlows, &fixer);
    for (int step = 0; step < 3; ++step)
    {
        integrator->advanceHierarchy(0.005);
        const double volume = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(corrected)[0];
        if (!(std::abs(volume - target) <= 64.0 * std::numeric_limits<double>::epsilon() + 1.0e-12 * target))
        {
            TBOX_ERROR("Correction regression: transported/reinitialized phase volume\n");
        }
    }
    double distance_error = 0.0, control_error = 0.0, change = 0.0, difference = 0.0, mesh_width = 0.0;
    const int weight_idx = integrator->getHierarchyMathOps()->getCellWeightPatchDescriptorIndex();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> phi = patch->getPatchData(current[0]);
            Pointer<CellData<NDIM, double>> other = patch->getPatchData(current[2]);
            Pointer<CellData<NDIM, double>> weight = patch->getPatchData(weight_idx);
            Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
            const Box<NDIM>& box = patch->getBox();
            mesh_width = std::max(mesh_width, geom->getDx()[0]);
            for (Box<NDIM>::Iterator i(box); i; i++)
            {
                const double x = geom->getXLower()[0] + (i()(0) - box.lower(0) + 0.5) * geom->getDx()[0];
                const double y = geom->getXLower()[1] + (i()(1) - box.lower(1) + 0.5) * geom->getDx()[1];
                const double distance = std::hypot(x - 0.5, y - 0.5) - 0.2;
                const double value = (*phi)(i());
                const double other_value = (*other)(i());
                if (!std::isfinite(value) || !std::isfinite(other_value))
                {
                    TBOX_ERROR("Correction regression: nonfinite transported field\n");
                }
                if (std::abs(distance) < 0.15)
                {
                    distance_error += std::abs(value - distance) * (*weight)(i());
                    control_error += std::abs(other_value - distance) * (*weight)(i());
                }
                change = std::max(change, std::abs(other_value - distance));
                difference = std::max(difference, std::abs(value - other_value));
            }
        }
    }
    distance_error = IBTK_MPI::sumReduction(distance_error);
    control_error = IBTK_MPI::sumReduction(control_error);
    change = IBTK_MPI::maxReduction(change);
    difference = IBTK_MPI::maxReduction(difference);
    mesh_width = IBTK_MPI::maxReduction(mesh_width);
    // The exact transported distance changes by at most max|u|*time.
    // Allow two cells of spatial error around the original circular interface.
    const double geometry_bound = std::sqrt(2.0) * integrator->getIntegratorTime() + 2.0 * mesh_width;
    if (!(distance_error <= geometry_bound && control_error <= geometry_bound && difference <= 2.0 * mesh_width &&
          change > 128.0 * std::numeric_limits<double>::epsilon()) ||
        sources[0].second != 3 || sources[1].second != 3)
    {
        TBOX_ERROR("Correction regression: transported interface geometry or inactive reinitialization: "
                   << distance_error << ", " << control_error << ", " << difference << ", " << change
                   << ", calls = " << sources[0].second << ", " << sources[1].second << '\n');
    }
    const double volume = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(corrected)[0];
    const double control_volume = LevelSetUtilities::computeHeavisideIntegrals2PhaseFlows(control)[0];
    plog << "transported target, corrected volume, control volume: " << target << ' ' << volume << ' ' << control_volume
         << '\n';
    plog << "near-interface distance L1, control L1, field difference: " << distance_error << ' ' << control_error
         << ' ' << difference << '\n';
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
        plog << std::setprecision(16);
        Pointer<AdvDiffHierarchyIntegrator> integrator =
            new AdvDiffSemiImplicitHierarchyIntegrator("AdvDiff", app->getComponentDatabase("AdvDiff"));
        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("Hierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tagger = new StandardTagAndInitialize<NDIM>(
            "Tagger", integrator, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load = new LoadBalancer<NDIM>("Load", app->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "Gridding", app->getComponentDatabase("GriddingAlgorithm"), tagger, boxes, load);
        Pointer<FaceVariable<NDIM, double>> velocity = new FaceVariable<NDIM, double>("velocity");
        integrator->registerAdvectionVelocity(velocity);
        integrator->setAdvectionVelocityFunction(
            velocity,
            new muParserCartGridFunction("velocity_function", app->getComponentDatabase("Velocity"), geometry));
        LocationIndexRobinBcCoefs<NDIM> bc;
        for (int face = 0; face < 2 * NDIM; ++face)
        {
            bc.setBoundarySlope(face, 0.0);
        }
        std::vector<Pointer<CellVariable<NDIM, double>>> variables;
        for (int k = 0; k < 3; ++k)
        {
            Pointer<CellVariable<NDIM, double>> var = new CellVariable<NDIM, double>("phi" + std::to_string(k));
            integrator->registerTransportedQuantity(var);
            integrator->setDiffusionCoefficient(var, 0.0);
            integrator->setConvectiveDifferencingType(var, ADVECTIVE);
            integrator->setAdvectionVelocity(var, velocity);
            integrator->setPhysicalBcCoef(var, &bc);
            integrator->setInitialConditions(var,
                                             new muParserCartGridFunction("initial" + std::to_string(k),
                                                                          app->getComponentDatabase("Initial"),
                                                                          geometry));
            variables.push_back(var);
        }
        integrator->initializePatchHierarchy(hierarchy, gridding);
        const double length = geometry->getXUpper()[0] - geometry->getXLower()[0];
        if (app->getInputDatabase()->getStringWithDefault("mode", "algebraic") == "lifecycle")
        {
            check_correction_lifecycle(integrator, variables, app->getInputDatabase(), length);
        }
        else if (app->getInputDatabase()->getStringWithDefault("mode", "algebraic") == "transport")
        {
            check_corrected_transport(integrator, variables, bc);
        }
        else
        {
            check_algebraic_corrections(integrator, variables, app->getInputDatabase(), length);
        }
    }
}
