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

#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonPointRelaxationFACOperator.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PoissonFACPreconditioner.h>
#include <ibtk/SAMRAIScopedVectorCopy.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
using HierarchyVector = SAMRAIVectorReal<NDIM, double>;

double
active_range_difference(const HierarchyVector& first, const HierarchyVector& second, int lower, int upper)
{
    double maximum = 0.0;
    int nonfinite = 0;
    Pointer<PatchHierarchy<NDIM>> hierarchy = first.getPatchHierarchy();
    for (int ln = lower; ln <= upper; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> x = patch->getPatchData(first.getComponentDescriptorIndex(0)),
                                            y = patch->getPatchData(second.getComponentDescriptorIndex(0)),
                                            weight = patch->getPatchData(first.getControlVolumeIndex(0));
            for (CellIterator<NDIM> c(patch->getBox()); c; c++)
            {
                if ((*weight)(c()) == 0.0)
                {
                    continue;
                }
                const double difference = std::abs((*x)(c()) - (*y)(c()));
                nonfinite += !std::isfinite(difference);
                maximum = std::max(maximum, difference);
            }
        }
    }
    if (IBTK_MPI::sumReduction(nonfinite) != 0)
    {
        TBOX_ERROR("FAC prefix test: nonfinite active comparison data\n");
    }
    return IBTK_MPI::maxReduction(maximum);
}

// Count the actual traversal without substituting any numerical operation.
class CountingFACOperator : public CCPoissonPointRelaxationFACOperator
{
public:
    CountingFACOperator(Pointer<Database> db, int levels)
        : CCPoissonPointRelaxationFACOperator("fac_operator", db, ""),
          d_visits(levels, 0),
          d_check_prefix(db->getBoolWithDefault("check_prefix", false))
    {
    }

    void prolongError(const HierarchyVector& source, HierarchyVector& destination, int level) override
    {
        ++d_visits[level];
        CCPoissonPointRelaxationFACOperator::prolongError(source, destination, level);
    }

    void prolongErrorAndCorrect(const HierarchyVector& source, HierarchyVector& destination, int level) override
    {
        ++d_visits[level];
        CCPoissonPointRelaxationFACOperator::prolongErrorAndCorrect(source, destination, level);
    }

    bool solveCoarsestLevel(HierarchyVector& error, const HierarchyVector& residual, int level) override
    {
        ++d_visits[level];
        return CCPoissonPointRelaxationFACOperator::solveCoarsestLevel(error, residual, level);
    }

    // Compare real residuals starting from identical, independently owned
    // below-range backing. Retain the expected child prefix until restriction.
    void computeResidual(HierarchyVector& residual,
                         const HierarchyVector& solution,
                         const HierarchyVector& rhs,
                         int lower,
                         int upper) override
    {
        if (!d_check_prefix || lower == 0 || solution.getCoarsestLevelNumber() == 0)
        {
            CCPoissonPointRelaxationFACOperator::computeResidual(residual, solution, rhs, lower, upper);
            return;
        }
        Pointer<PatchHierarchy<NDIM>> hierarchy = solution.getPatchHierarchy();
        HierarchyVector full("prefix_full", hierarchy, 0, hierarchy->getFinestLevelNumber());
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> full_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, 0, hierarchy->getFinestLevelNumber());
        full.addComponent(solution.getComponentVariable(0),
                          solution.getComponentDescriptorIndex(0),
                          solution.getControlVolumeIndex(0),
                          full_ops);
        SAMRAIScopedVectorCopy<double> actual_storage(full), zero_storage(full);
        SAMRAIScopedVectorDuplicate<double> actual_residual_storage(full), zero_residual_storage(full);
        HierarchyVector& actual = actual_storage;
        HierarchyVector& zero = zero_storage;
        HierarchyVector& actual_residual = actual_residual_storage;
        HierarchyVector& zero_residual = zero_residual_storage;
        HierarchyCellDataOpsReal<NDIM, double> range_ops(hierarchy, lower, upper);
        range_ops.setToScalar(zero.getComponentDescriptorIndex(0), 0.0, false);
        CCPoissonPointRelaxationFACOperator::computeResidual(actual_residual, actual, rhs, lower, upper);
        CCPoissonPointRelaxationFACOperator::computeResidual(zero_residual, zero, rhs, lower, upper);
        // The child starts after the actual parent's synchronization. Compare
        // its zero-state contribution with the parent's original zero state.
        if (upper > lower + 1)
        {
            d_parent_residual = std::make_unique<SAMRAIScopedVectorCopy<double>>(actual_residual);
            d_child_zero = std::make_unique<SAMRAIScopedVectorCopy<double>>(actual);
            HierarchyVector& child_zero = *d_child_zero;
            range_ops.setToScalar(child_zero.getComponentDescriptorIndex(0), 0.0, false);
            SAMRAIScopedVectorCopy<double> child_storage(child_zero);
            SAMRAIScopedVectorDuplicate<double> child_residual_storage(full);
            HierarchyVector& child = child_storage;
            HierarchyVector& child_residual = child_residual_storage;
            CCPoissonPointRelaxationFACOperator::computeResidual(child_residual, child, rhs, lower, upper - 1);
            const double discrepancy = active_range_difference(zero_residual, child_residual, lower, upper - 2);
            if (discrepancy > 1.0e-10)
            {
                TBOX_ERROR("FAC prefix test: parent-child zero-state contributions differ\n");
            }
            plog << "parent-child boundary discrepancy " << lower << ' ' << upper << ' ' << discrepancy << '\n';
        }
        d_expected_prefix = std::make_unique<SAMRAIScopedVectorCopy<double>>(actual_residual);
        HierarchyVector& expected = *d_expected_prefix;
        d_prefix_lower = lower;
        d_prefix_upper = upper - 2;
        for (int ln = 0; ln <= upper; ++ln)
        {
            double response = 0.0, discarded = 0.0, boundary = 0.0, actual_change = 0.0, zero_change = 0.0,
                   backing_difference = 0.0;
            int nonfinite = 0;
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> u = patch->getPatchData(solution.getComponentDescriptorIndex(0)),
                                                a = patch->getPatchData(actual.getComponentDescriptorIndex(0)),
                                                z = patch->getPatchData(zero.getComponentDescriptorIndex(0)),
                                                w = patch->getPatchData(solution.getControlVolumeIndex(0));
                for (CellIterator<NDIM> c(u->getGhostBox()); c; c++)
                {
                    nonfinite += !std::isfinite((*a)(c())) || !std::isfinite((*z)(c())) || !std::isfinite((*u)(c()));
                    backing_difference = std::max(backing_difference, std::abs((*a)(c()) - (*z)(c())));
                    actual_change = std::max(actual_change, std::abs((*a)(c()) - (*u)(c())));
                    zero_change = std::max(zero_change, std::abs((*z)(c()) - (ln < lower ? (*u)(c()) : 0.0)));
                }
                if (ln >= lower)
                {
                    Pointer<CellData<NDIM, double>> ra = patch->getPatchData(
                                                        actual_residual.getComponentDescriptorIndex(0)),
                                                    rz = patch->getPatchData(
                                                        zero_residual.getComponentDescriptorIndex(0)),
                                                    f = patch->getPatchData(rhs.getComponentDescriptorIndex(0)),
                                                    e = patch->getPatchData(expected.getComponentDescriptorIndex(0));
                    for (CellIterator<NDIM> c(patch->getBox()); c; c++)
                    {
                        (*e)(c()) = (*f)(c()) + ((*ra)(c()) - (*rz)(c()));
                        if ((*w)(c()) == 0.0)
                        {
                            continue;
                        }
                        nonfinite +=
                            !std::isfinite((*ra)(c())) || !std::isfinite((*rz)(c())) || !std::isfinite((*f)(c()));
                        response = std::max(response, std::abs((*ra)(c()) - (*rz)(c())));
                        discarded = std::max(discarded, std::abs((*ra)(c()) - (*f)(c())));
                        boundary = std::max(boundary, std::abs((*rz)(c()) - (*f)(c())));
                    }
                }
            }
            if (IBTK_MPI::sumReduction(nonfinite) != 0)
            {
                TBOX_ERROR("FAC prefix test: nonfinite evaluation data\n");
            }
            plog << "range " << lower << " " << upper << " level " << ln << " response "
                 << IBTK_MPI::maxReduction(response) << " discarded " << IBTK_MPI::maxReduction(discarded)
                 << " boundary " << IBTK_MPI::maxReduction(boundary) << " actual_sync "
                 << IBTK_MPI::maxReduction(actual_change) << " zero_sync " << IBTK_MPI::maxReduction(zero_change)
                 << " synchronized_difference " << IBTK_MPI::maxReduction(backing_difference) << '\n';
        }
        CCPoissonPointRelaxationFACOperator::computeResidual(residual, solution, rhs, lower, upper);
        HierarchyCellDataOpsReal<NDIM, double> check_ops(hierarchy, lower, upper);
        check_ops.subtract(actual_residual.getComponentDescriptorIndex(0),
                           actual_residual.getComponentDescriptorIndex(0),
                           residual.getComponentDescriptorIndex(0),
                           true);
        const double discrepancy = check_ops.L1Norm(actual_residual.getComponentDescriptorIndex(0));
        if (!std::isfinite(discrepancy) || discrepancy != 0.0)
        {
            TBOX_ERROR("FAC prefix test: paired evaluations changed the actual residual\n");
        }
    }

    void restrictResidual(const HierarchyVector& source, HierarchyVector& destination, int level) override
    {
        if (d_expected_prefix && d_prefix_upper >= d_prefix_lower)
        {
            HierarchyVector& expected = *d_expected_prefix;
            const double maximum = active_range_difference(source, expected, d_prefix_lower, d_prefix_upper);
            if (maximum > 1.0e-10)
            {
                TBOX_ERROR("FAC prefix test: child RHS lost the live correction response\n");
            }
            // Independently evaluate the actual child equation at zero with
            // its inherited backing and the RHS supplied by the production cycle.
            HierarchyVector& child_zero = *d_child_zero;
            SAMRAIScopedVectorDuplicate<double> child_residual_storage(child_zero);
            HierarchyVector& child_residual = child_residual_storage;
            CCPoissonPointRelaxationFACOperator::computeResidual(
                child_residual, child_zero, source, d_prefix_lower, level);
            const double handoff =
                active_range_difference(child_residual, *d_parent_residual, d_prefix_lower, d_prefix_upper);
            if (handoff > 1.0e-10)
            {
                TBOX_ERROR("FAC prefix test: zero child residual does not preserve the parent residual\n");
            }
            plog << "zero-child residual discrepancy " << handoff << '\n';
            plog << "child prefix discrepancy " << maximum << '\n';
        }
        d_expected_prefix.reset();
        d_parent_residual.reset();
        d_child_zero.reset();
        CCPoissonPointRelaxationFACOperator::restrictResidual(source, destination, level);
    }

    void resetVisits()
    {
        std::fill(d_visits.begin(), d_visits.end(), 0);
    }

    const std::vector<int>& getVisits() const
    {
        return d_visits;
    }

private:
    std::vector<int> d_visits;
    bool d_check_prefix;
    int d_prefix_lower = 0, d_prefix_upper = -1;
    std::unique_ptr<SAMRAIScopedVectorCopy<double>> d_expected_prefix, d_parent_residual, d_child_zero;
};

double
data_difference(Pointer<PatchHierarchy<NDIM>> hierarchy,
                int first,
                int second,
                int weight,
                bool active_only,
                bool interiors_only = false)
{
    double maximum = 0.0;
    int nonfinite = 0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> x = patch->getPatchData(first), y = patch->getPatchData(second),
                                            w = patch->getPatchData(weight);
            for (CellIterator<NDIM> c(active_only || interiors_only ? patch->getBox() : x->getGhostBox()); c; c++)
            {
                if (active_only && (*w)(c()) == 0.0)
                {
                    continue;
                }
                const double value = std::abs((*x)(c()) - (*y)(c()));
                nonfinite += !std::isfinite(value);
                maximum = std::max(maximum, value);
            }
        }
    }
    if (IBTK_MPI::sumReduction(nonfinite) != 0)
    {
        TBOX_ERROR("FAC cycle test: nonfinite comparison data\n");
    }
    return IBTK_MPI::maxReduction(maximum);
}
// Exercise vector views separately from their caller-owned backing allocation.
void
check_ranges(Pointer<PatchHierarchy<NDIM>> hierarchy,
             Pointer<CellVariable<NDIM, double>> variable,
             int weight,
             Pointer<Database> input,
             const PoissonSpecifications& specification,
             RobinBcCoefStrategy<NDIM>* boundary)
{
    const int finest = hierarchy->getFinestLevelNumber();
    VariableDatabase<NDIM>* variable_db = VariableDatabase<NDIM>::getDatabase();
    std::array<int, 4> indices;
    for (int k = 0; k < 4; ++k)
    {
        indices[k] = variable_db->registerVariableAndContext(
            variable, variable_db->getContext("range_" + std::to_string(k)), IntVector<NDIM>(1));
        for (int ln = 0; ln <= finest; ++ln)
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(indices[k]);
        }
    }
    Pointer<HierarchyCellDataOpsReal<NDIM, double>> caller_ops =
        new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, 0, finest);
    const auto reset = [&]()
    {
        for (int ln = 0; ln <= finest; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<CellData<NDIM, double>> x = level->getPatch(p())->getPatchData(indices[0]),
                                                b = level->getPatch(p())->getPatchData(indices[1]);
                x->fillAll(7.0);
                for (CellIterator<NDIM> c(b->getGhostBox()); c; c++)
                {
                    (*b)(c()) = 1.0 + 0.01 * (c()(0) + 2 * c()(1));
                }
            }
        }
        caller_ops->copyData(indices[2], indices[1], false);
    };
    const auto view = [&](int index, int lower, int upper)
    {
        Pointer<HierarchyVector> vector = new HierarchyVector("range", hierarchy, lower, upper);
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> range_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, lower, upper);
        vector->addComponent(variable, index, weight, range_ops);
        return vector;
    };
    Pointer<CountingFACOperator> strategy = new CountingFACOperator(input->getDatabase("fac_db"), finest + 1);
    PoissonFACPreconditioner fac("range_fac", strategy, input->getDatabase("fac_db"), "");
    fac.setPoissonSpecifications(specification);
    fac.setPhysicalBcCoef(boundary);
    fac.setHomogeneousBc(true);
    fac.setNumPostSmoothingSweeps(2);
    if (input->getDatabase("fac_db")->getBoolWithDefault("check_prefix", false))
    {
        reset();
        fac.setMGCycleType(V_CYCLE);
        fac.setNumPreSmoothingSweeps(2);
        Pointer<HierarchyVector> x = view(indices[0], 1, finest), b = view(indices[1], 1, finest);
        fac.initializeSolverState(*x, *b);
        for (int repeat = 0; repeat < 2; ++repeat)
        {
            reset();
            plog << "initialized solve " << repeat << '\n';
            fac.solveSystem(*x, *b);
            if (repeat == 0)
            {
                caller_ops->copyData(indices[3], indices[0], false);
            }
            else if (data_difference(hierarchy, indices[0], indices[3], weight, false) != 0.0)
            {
                TBOX_ERROR("FAC prefix test: initialized solve changed on reuse\n");
            }
            // Freeze the actual end-of-cycle backing in independent copies.
            // This measures f - A(x; g_end), not an externally prescribed solve.
            Pointer<HierarchyVector> full = view(indices[0], 0, finest);
            SAMRAIScopedVectorCopy<double> evaluation_storage(full), zero_storage(full);
            SAMRAIScopedVectorDuplicate<double> residual_storage(full), zero_residual_storage(full);
            HierarchyVector& evaluation = evaluation_storage;
            HierarchyVector& zero = zero_storage;
            HierarchyVector& residual = residual_storage;
            HierarchyVector& zero_residual = zero_residual_storage;
            HierarchyCellDataOpsReal<NDIM, double> range_ops(hierarchy, 1, finest);
            range_ops.setToScalar(zero.getComponentDescriptorIndex(0), 0.0, false);
            strategy->CCPoissonPointRelaxationFACOperator::computeResidual(residual, evaluation, *b, 1, finest);
            strategy->CCPoissonPointRelaxationFACOperator::computeResidual(zero_residual, zero, *b, 1, finest);
            const double residual_norm = range_ops.L2Norm(residual.getComponentDescriptorIndex(0), weight),
                         zero_norm = range_ops.L2Norm(zero_residual.getComponentDescriptorIndex(0), weight);
            if (!std::isfinite(residual_norm) || !std::isfinite(zero_norm))
            {
                TBOX_ERROR("FAC prefix test: nonfinite completed-cycle residual\n");
            }
            plog << "completed-cycle residual " << residual_norm << " zero-state residual " << zero_norm << '\n';
        }
        fac.deallocateSolverState();
    }
    else if (input->keyExists("reject_cycle"))
    {
        reset();
        fac.setMGCycleType(string_to_enum<MGCycleType>(input->getString("reject_cycle")));
        fac.setNumPreSmoothingSweeps(input->getIntegerWithDefault("pre_sweeps", 2));
        fac.setMGCycleMultiplicity(input->getIntegerWithDefault("multiplicity", 2));
        const int lower = input->getIntegerWithDefault("range_lower", 1);
        Pointer<HierarchyVector> x = view(indices[0], lower, finest), b = view(indices[1], lower, finest);
        const bool initial_setup = input->getBoolWithDefault("initial_setup", false);
        if (!initial_setup)
        {
            // Recheck inputs after a supported initialization and configuration change.
            const MGCycleType cycle = fac.getMGCycleType();
            const int pre = fac.getNumPreSmoothingSweeps();
            fac.setMGCycleType(V_CYCLE);
            fac.setNumPreSmoothingSweeps(0);
            fac.initializeSolverState(*x, *b);
            fac.setMGCycleType(cycle);
            fac.setNumPreSmoothingSweeps(pre);
        }
        if (input->getBoolWithDefault("missing_backing", false))
        {
            hierarchy->getPatchLevel(0)->deallocatePatchData(indices[0]);
        }
        if (initial_setup)
        {
            fac.initializeSolverState(*x, *b);
        }
        else
        {
            fac.solveSystem(*x, *b);
        }
        fac.deallocateSolverState();
    }
    else
    {
        std::vector<std::pair<int, int>> ranges = { { 1, 1 }, { finest, finest }, { 1, finest } };
        if (finest > 2)
        {
            ranges.emplace_back(2, finest);
        }
        for (const auto& range : ranges)
        {
            const int lower = range.first, upper = range.second;
            Pointer<HierarchyVector> x = view(indices[0], lower, upper), b = view(indices[1], lower, upper);
            const bool single = lower == upper;
            const std::vector<std::pair<MGCycleType, int>> cycles = { { V_CYCLE, 1 },  { W_CYCLE, 2 },
                                                                      { F_CYCLE, 2 },  { FMG_CYCLE, 1 },
                                                                      { MU_CYCLE, 1 }, { MU_CYCLE, 2 },
                                                                      { MU_CYCLE, 3 } };
            for (const auto& cycle_case : cycles)
            {
                const MGCycleType cycle = cycle_case.first;
                const int multiplicity = cycle_case.second;
                fac.setMGCycleMultiplicity(multiplicity);
                if (!single && (cycle == W_CYCLE || cycle == F_CYCLE || cycle == FMG_CYCLE ||
                                (cycle == MU_CYCLE && multiplicity > 1)))
                {
                    continue;
                }
                for (int pre : { 0, 2 })
                {
                    fac.setMGCycleType(cycle);
                    fac.setNumPreSmoothingSweeps(0);
                    reset();
                    fac.initializeSolverState(*x, *b);
                    fac.setNumPreSmoothingSweeps(pre);
                    for (int repeat = 0; repeat < 2; ++repeat)
                    {
                        reset();
                        strategy->resetVisits();
                        fac.solveSystem(*x, *b);
                        // Compare every RHS entry, including ghosts and levels
                        // outside the view, except the legacy in-place V path.
                        if (!(cycle == V_CYCLE && pre == 0) &&
                            data_difference(hierarchy, indices[1], indices[2], weight, false) != 0.0)
                        {
                            TBOX_ERROR("FAC range test: RHS changed\n");
                        }
                        int changed_below = 0;
                        for (int ln = 0; ln <= finest; ++ln)
                        {
                            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                            {
                                Pointer<CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(indices[0]);
                                for (CellIterator<NDIM> c(data->getGhostBox()); c; c++)
                                {
                                    if (!std::isfinite((*data)(c())))
                                    {
                                        TBOX_ERROR("FAC range test: nonfinite solution\n");
                                    }
                                    if (ln < lower && (*data)(c()) != 7.0)
                                    {
                                        ++changed_below;
                                    }
                                    if ((ln > upper || ln < lower - 1) && (*data)(c()) != 7.0)
                                    {
                                        TBOX_ERROR("FAC range test: unexpected outside-range mutation\n");
                                    }
                                }
                            }
                        }
                        changed_below = IBTK_MPI::sumReduction(changed_below);
                        const bool below_access = !single && pre > 0;
                        if ((changed_below != 0) != below_access)
                        {
                            TBOX_ERROR("FAC range test: unexpected below-range mutation\n");
                        }
                        if (single && strategy->getVisits()[lower] != 1)
                        {
                            TBOX_ERROR("FAC range test: single-level cycle did not coarse-solve once\n");
                        }
                        if (repeat == 0)
                        {
                            caller_ops->copyData(indices[3], indices[0], false);
                        }
                        else if (data_difference(hierarchy, indices[0], indices[3], weight, false) != 0.0)
                        {
                            TBOX_ERROR("FAC range test: reused solve changed interiors or ghosts\n");
                        }
                    }
                    fac.deallocateSolverState();
                    // A fresh state with another component index must reproduce
                    // every entry, including the permitted below-range writes.
                    reset();
                    caller_ops->setToScalar(indices[2], 7.0, false);
                    Pointer<HierarchyVector> alternate = view(indices[2], lower, upper);
                    fac.initializeSolverState(*x, *b);
                    fac.solveSystem(*alternate, *b);
                    fac.deallocateSolverState();
                    if (data_difference(hierarchy, indices[2], indices[3], weight, false) != 0.0)
                    {
                        TBOX_ERROR("FAC range test: fresh solve changed interiors or ghosts\n");
                    }
                    if (single)
                    {
                        reset();
                        for (int ln = 0; ln < lower; ++ln)
                        {
                            hierarchy->getPatchLevel(ln)->deallocatePatchData(indices[0]);
                            hierarchy->getPatchLevel(ln)->deallocatePatchData(indices[1]);
                        }
                        fac.solveSystem(*x, *b);
                        for (int ln = 0; ln < lower; ++ln)
                        {
                            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                            level->allocatePatchData(indices[0]);
                            level->allocatePatchData(indices[1]);
                            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                            {
                                Pointer<CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(indices[0]);
                                data->fillAll(7.0);
                            }
                        }
                        if (data_difference(hierarchy, indices[0], indices[3], weight, false) != 0.0)
                        {
                            TBOX_ERROR("FAC range test: single-level solve depends on backing data\n");
                        }
                    }
                    if (cycle == MU_CYCLE && multiplicity == 1)
                    {
                        reset();
                        fac.setMGCycleType(V_CYCLE);
                        fac.solveSystem(*x, *b);
                        // The default V path does not promise identical ghosts.
                        // Compare every interior, including covered cells.
                        if (data_difference(hierarchy, indices[0], indices[3], weight, false, pre == 0) != 0.0)
                        {
                            TBOX_ERROR("FAC range test: multiplicity one differs from V\n");
                        }
                    }
                    plog << "range = " << lower << ' ' << upper << "; cycle = " << enum_to_string(cycle);
                    if (cycle == MU_CYCLE)
                    {
                        plog << "; mu = " << multiplicity;
                    }
                    plog << "; pre = " << pre << "; correction L2 = " << view(indices[3], lower, upper)->L2Norm()
                         << '\n';
                }
            }
            // The default V path and single-level solves need no backing data
            // below the view. Check with both caller components absent there.
            reset();
            for (int ln = 0; ln < lower; ++ln)
            {
                hierarchy->getPatchLevel(ln)->deallocatePatchData(indices[0]);
                hierarchy->getPatchLevel(ln)->deallocatePatchData(indices[1]);
            }
            for (MGCycleType cycle : { V_CYCLE, MU_CYCLE })
            {
                fac.setMGCycleType(cycle);
                fac.setMGCycleMultiplicity(1);
                fac.setNumPreSmoothingSweeps(single ? 2 : 0);
                fac.solveSystem(*x, *b);
            }
            for (int ln = 0; ln < lower; ++ln)
            {
                hierarchy->getPatchLevel(ln)->allocatePatchData(indices[0]);
                hierarchy->getPatchLevel(ln)->allocatePatchData(indices[1]);
            }
        }
    }
    for (int ln = 0; ln <= finest; ++ln)
    {
        for (int index : indices)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(index);
        }
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "fac_cycles.log");
        Pointer<Database> input = app->getInputDatabase();
        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tags = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> balance =
            new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), tags, boxes, balance);
        gridding->makeCoarsestLevel(hierarchy, 0.0);
        for (int ln = 0; gridding->levelCanBeRefined(ln); ++ln)
        {
            gridding->makeFinerLevel(hierarchy, 0.0, 0.0, 1);
            if (!hierarchy->finerLevelExists(ln))
            {
                break;
            }
        }
        const int finest = hierarchy->getFinestLevelNumber();
        TBOX_ASSERT(finest + 1 == input->getInteger("levels"));

        enum Field
        {
            EXACT,
            RHS,
            SECOND_RHS,
            COMBINED_RHS,
            EVALUATION,
            RESIDUAL,
            CORRECTION,
            RESULT,
            SECOND_RESULT,
            THIRD_RESULT,
            SNAPSHOT,
            ITERATE,
            ERROR,
            V_RESULT,
            W_RESULT,
            OPS_PROBE,
            FIELD_COUNT
        };
        HierarchyMathOps math_ops("math_ops", hierarchy);
        const int weight = math_ops.getCellWeightPatchDescriptorIndex();
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> shared_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, 0, finest);
        HierarchyCellDataOpsReal<NDIM, double>& ops = *shared_ops;
        std::array<int, FIELD_COUNT> indices;
        std::array<std::unique_ptr<HierarchyVector>, FIELD_COUNT> vectors;
        Pointer<CellVariable<NDIM, double>> variable = new CellVariable<NDIM, double>("u");
        VariableDatabase<NDIM>* variable_db = VariableDatabase<NDIM>::getDatabase();
        LocationIndexRobinBcCoefs<NDIM> boundary("boundary", nullptr);
        for (int d = 0; d < NDIM; ++d)
        {
            boundary.setBoundaryValue(2 * d, 0.0);
            boundary.setBoundaryValue(2 * d + 1, 0.0);
        }
        PoissonSpecifications specification("specification");
        specification.setCConstant(1.0);
        specification.setDConstant(-1.0);
        if (input->getBoolWithDefault("range_test", false))
        {
            Logger::getInstance()->setAbortAppender(new TestAppender());
            check_ranges(hierarchy, variable, weight, input, specification, &boundary);
            return 0;
        }
        for (int f = 0; f < FIELD_COUNT; ++f)
        {
            const std::string name = "field_" + std::to_string(f);
            indices[f] =
                variable_db->registerVariableAndContext(variable, variable_db->getContext(name), IntVector<NDIM>(1));
            for (int ln = 0; ln <= finest; ++ln)
            {
                hierarchy->getPatchLevel(ln)->allocatePatchData(indices[f], 0.0);
            }
            vectors[f] = std::make_unique<HierarchyVector>(name, hierarchy, 0, finest);
            vectors[f]->addComponent(variable, indices[f], weight, shared_ops);
            ops.setToScalar(indices[f], 0.0, false);
        }
        const auto copy = [&](Field destination, Field source)
        { ops.copyData(indices[destination], indices[source], false); };
        const auto difference = [&](Field first, Field second, bool active_only = true)
        { return data_difference(hierarchy, indices[first], indices[second], weight, active_only); };
        const auto check_data_ops = [&]()
        {
            // Caller-owned operations must still span the complete hierarchy,
            // even when FAC creates vectors for truncated correction equations.
            for (int ln = 0; ln <= finest; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(indices[OPS_PROBE]);
                    data->fillAll(-1.0);
                }
            }
            ops.setToScalar(indices[OPS_PROBE], 1.0, false);
            int unchanged = 0;
            for (int ln = 0; ln <= finest; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(indices[OPS_PROBE]);
                    for (CellIterator<NDIM> c(data->getGhostBox()); c; c++)
                    {
                        unchanged += (*data)(c()) != 1.0;
                    }
                }
            }
            if (IBTK_MPI::sumReduction(unchanged) != 0)
            {
                TBOX_ERROR("FAC cycle test: FAC changed the caller's data-operation level range\n");
            }
        };

        const double pi = std::acos(-1.0);
        for (int ln = 0; ln <= finest; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            const double n = input->getInteger("N") * level->getRatio()(0);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> exact = patch->getPatchData(indices[EXACT]),
                                                second = patch->getPatchData(indices[EVALUATION]);
                for (CellIterator<NDIM> c(patch->getBox()); c; c++)
                {
                    const double x = (c()(0) + 0.5) / n, y = (c()(1) + 0.5) / n;
                    (*exact)(c()) =
                        std::sin(pi * x) * std::sin(pi * y) + 0.2 * std::sin(3.0 * pi * x) * std::sin(2.0 * pi * y);
                    (*second)(c()) = std::sin(2.0 * pi * x) * std::sin(3.0 * pi * y);
                }
            }
        }
        CCLaplaceOperator laplace("laplace", input->getDatabase("fac_db"));
        laplace.setPoissonSpecifications(specification);
        laplace.setPhysicalBcCoef(&boundary);
        laplace.setHomogeneousBc(true);
        laplace.initializeOperatorState(*vectors[EVALUATION], *vectors[RESIDUAL]);
        laplace.apply(*vectors[EVALUATION], *vectors[SECOND_RHS]);
        copy(EVALUATION, EXACT);
        laplace.apply(*vectors[EVALUATION], *vectors[RHS]);
        ops.linearSum(indices[COMBINED_RHS], 0.31, indices[RHS], -0.7, indices[SECOND_RHS], false);

        Pointer<CountingFACOperator> strategy = new CountingFACOperator(input->getDatabase("fac_db"), finest + 1);
        PoissonFACPreconditioner fac("fac", strategy, input->getDatabase("fac_db"), "");
        fac.setPoissonSpecifications(specification);
        fac.setPhysicalBcCoef(&boundary);
        fac.setHomogeneousBc(true);
        const auto apply = [&](Field destination, Field source)
        {
            check_data_ops();
            // The legacy no-presmoothing V path may modify covered RHS storage.
            copy(SNAPSHOT, source);
            ops.setToScalar(indices[destination], 7.0, false);
            strategy->resetVisits();
            fac.solveSystem(*vectors[destination], *vectors[source]);
            check_data_ops();
            const bool legacy_v = fac.getMGCycleType() == V_CYCLE && fac.getNumPreSmoothingSweeps() == 0;
            if (difference(source, SNAPSHOT, legacy_v) != 0.0)
            {
                TBOX_ERROR("FAC cycle test: RHS data was modified\n");
            }
            copy(source, SNAPSHOT);
        };
        const auto residual = [&]()
        {
            copy(EVALUATION, ITERATE);
            laplace.apply(*vectors[EVALUATION], *vectors[RESIDUAL]);
            ops.subtract(indices[RESIDUAL], indices[RHS], indices[RESIDUAL], false);
        };
        TBOX_ASSERT(string_to_enum<MGCycleType>("MU_CYCLE") == MU_CYCLE);
        TBOX_ASSERT(string_to_enum<MGCycleType>("mu_cycle") == UNKNOWN_MG_CYCLE_TYPE);
        TBOX_ASSERT(string_to_enum<MGCycleType>("Mu_Cycle") == UNKNOWN_MG_CYCLE_TYPE);
        TBOX_ASSERT(string_to_enum<MGCycleType>("MU") == UNKNOWN_MG_CYCLE_TYPE);
        TBOX_ASSERT(string_to_enum<MGCycleType>("MU-CYCLE") == UNKNOWN_MG_CYCLE_TYPE);
        TBOX_ASSERT(enum_to_string(MU_CYCLE) == "MU_CYCLE");
        struct CycleCase
        {
            MGCycleType type;
            const char* name;
            int multiplicity;
        };
        const std::array<CycleCase, 5> cycles = { { { V_CYCLE, "V", 1 },
                                                    { W_CYCLE, "W", 2 },
                                                    { MU_CYCLE, "MU3", 3 },
                                                    { F_CYCLE, "F", 2 },
                                                    { FMG_CYCLE, "FMG", 1 } } };
        plog << "levels = " << finest + 1 << '\n';
        for (int pre : { 0, 2 })
        {
            fac.setNumPreSmoothingSweeps(pre);
            fac.setNumPostSmoothingSweeps(2);
            fac.initializeSolverState(*vectors[CORRECTION], *vectors[RHS]);
            check_data_ops();
            plog << "pre_sweeps = " << pre << '\n';
            for (const CycleCase& cycle : cycles)
            {
                fac.setMGCycleType(cycle.type);
                fac.setMGCycleMultiplicity(cycle.multiplicity);
                TBOX_ASSERT(fac.getMGCycleMultiplicity() == cycle.multiplicity);
                apply(RESULT, RHS);
                plog << cycle.name << " visits =";
                for (int ln = finest; ln >= 0; --ln)
                {
                    const int depth = finest - ln;
                    int expected = cycle.type == F_CYCLE ? depth + 1 : 1;
                    if (cycle.type == FMG_CYCLE)
                    {
                        // Each nested V-cycle visits this level once. Noncoarse levels
                        // also receive the initial FMG prolongation.
                        expected = depth + (ln > 0 ? 2 : 1);
                    }
                    else if (cycle.type != F_CYCLE)
                    {
                        for (int i = 0; i < depth; ++i)
                        {
                            expected *= cycle.multiplicity;
                        }
                    }
                    TBOX_ASSERT(strategy->getVisits()[ln] == expected);
                    plog << ' ' << strategy->getVisits()[ln];
                }
                plog << '\n';
                if (cycle.type == V_CYCLE)
                {
                    copy(V_RESULT, RESULT);
                }
                if (cycle.type == W_CYCLE)
                {
                    copy(W_RESULT, RESULT);
                }
                copy(ITERATE, RESULT);
                residual();
                plog << cycle.name << " first correction L2 = " << vectors[RESULT]->L2Norm()
                     << "; residual L2 = " << vectors[RESIDUAL]->L2Norm() << '\n';
                if (cycle.type == FMG_CYCLE && !(vectors[RESIDUAL]->L2Norm() < vectors[RHS]->L2Norm()))
                {
                    TBOX_ERROR("FAC cycle test: FMG did not reduce the residual\n");
                }

                apply(SECOND_RESULT, SECOND_RHS);
                apply(THIRD_RESULT, COMBINED_RHS);
                ops.linearSum(indices[ERROR], 0.31, indices[RESULT], -0.7, indices[SECOND_RESULT], false);
                if (!(difference(THIRD_RESULT, ERROR) < 1.0e-11))
                {
                    TBOX_ERROR("FAC cycle test: FAC map is not linear\n");
                }

                // Change both vector indices after initialization, then rebuild solver state.
                copy(COMBINED_RHS, RHS);
                apply(THIRD_RESULT, COMBINED_RHS);
                if (!(difference(THIRD_RESULT, RESULT) < 1.0e-12))
                {
                    TBOX_ERROR("FAC cycle test: solve depends on vector data indices\n");
                }
                fac.deallocateSolverState();
                fac.initializeSolverState(*vectors[SECOND_RESULT], *vectors[SECOND_RHS]);
                check_data_ops();
                apply(THIRD_RESULT, RHS);
                if (!(difference(THIRD_RESULT, RESULT) < 1.0e-12))
                {
                    TBOX_ERROR("FAC cycle test: reinitialization changed the correction\n");
                }
                ops.linearSum(indices[COMBINED_RHS], 0.31, indices[RHS], -0.7, indices[SECOND_RHS], false);

                ops.setToScalar(indices[ERROR], 0.0, false);
                apply(CORRECTION, ERROR);
                if (difference(CORRECTION, ERROR) != 0.0)
                {
                    TBOX_ERROR("FAC cycle test: zero RHS produced a nonzero correction\n");
                }
                copy(ITERATE, EXACT);
                for (int ln = 0; ln < finest; ++ln)
                {
                    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM>> patch = level->getPatch(p());
                        Pointer<CellData<NDIM, double>> x = patch->getPatchData(indices[ITERATE]),
                                                        w = patch->getPatchData(weight);
                        for (CellIterator<NDIM> c(patch->getBox()); c; c++)
                        {
                            if ((*w)(c()) == 0.0)
                            {
                                (*x)(c()) += 0.125;
                            }
                        }
                    }
                }
                residual();
                if (!(vectors[RESIDUAL]->L2Norm() < 1.0e-11))
                {
                    TBOX_ERROR("FAC cycle test: covered storage changed the exact residual\n");
                }
                apply(CORRECTION, RESIDUAL);
                ops.add(indices[ITERATE], indices[ITERATE], indices[CORRECTION], false);
                if (!(difference(ITERATE, EXACT) < 1.0e-11))
                {
                    TBOX_ERROR("FAC cycle test: exact solution is not a fixed point\n");
                }

                for (double starting_scale : { 0.0, 0.25 })
                {
                    ops.scale(indices[ITERATE], starting_scale, indices[EXACT], false);
                    residual();
                    const double initial_residual = vectors[RESIDUAL]->L2Norm();
                    ops.subtract(indices[ERROR], indices[ITERATE], indices[EXACT], false);
                    const double initial_error = vectors[ERROR]->L2Norm();
                    for (int iteration = 0; iteration < 15; ++iteration)
                    {
                        apply(CORRECTION, RESIDUAL);
                        ops.add(indices[ITERATE], indices[ITERATE], indices[CORRECTION], false);
                        residual();
                    }
                    ops.subtract(indices[ERROR], indices[ITERATE], indices[EXACT], false);
                    const double final_residual = vectors[RESIDUAL]->L2Norm(), final_error = vectors[ERROR]->L2Norm();
                    if (!(final_residual < initial_residual))
                    {
                        TBOX_ERROR("FAC cycle test: composite residual did not decrease\n");
                    }
                    if (!(final_error < initial_error))
                    {
                        TBOX_ERROR("FAC cycle test: solution error did not decrease\n");
                    }
                    plog << cycle.name << " start = " << (starting_scale == 0.0 ? "zero" : "nonzero")
                         << "; final residual ratio = " << final_residual / initial_residual
                         << "; final error ratio = " << final_error / initial_error << '\n';
                }
            }
            fac.setMGCycleType(MU_CYCLE);
            fac.setMGCycleMultiplicity(1);
            apply(CORRECTION, RHS);
            if (!(difference(CORRECTION, V_RESULT) < 1.0e-12))
            {
                TBOX_ERROR("FAC cycle test: mu=1 differs from V\n");
            }
            fac.setMGCycleMultiplicity(2);
            apply(CORRECTION, RHS);
            if (!(difference(CORRECTION, W_RESULT) < 1.0e-12))
            {
                TBOX_ERROR("FAC cycle test: mu=2 differs from W\n");
            }
            fac.deallocateSolverState();
        }

        // Start with the V path without general-cycle scratch, then grow and
        // reuse its scratch without rebuilding state when smoothing changes.
        fac.setMGCycleType(V_CYCLE);
        fac.setNumPreSmoothingSweeps(0);
        fac.initializeSolverState(*vectors[CORRECTION], *vectors[RHS]);
        check_data_ops();
        apply(RESULT, RHS);
        fac.setNumPreSmoothingSweeps(2);
        apply(SECOND_RESULT, RHS);
        fac.setNumPreSmoothingSweeps(0);
        apply(CORRECTION, RHS);
        if (!(difference(CORRECTION, RESULT) < 1.0e-12))
        {
            TBOX_ERROR("FAC cycle test: disabling presmoothing changed the V correction\n");
        }
        fac.deallocateSolverState();
        fac.setNumPreSmoothingSweeps(2);
        fac.initializeSolverState(*vectors[SECOND_RESULT], *vectors[SECOND_RHS]);
        check_data_ops();
        apply(CORRECTION, RHS);
        if (!(difference(CORRECTION, SECOND_RESULT) < 1.0e-12))
        {
            TBOX_ERROR("FAC cycle test: changing presmoothing differs from fresh initialization\n");
        }
        fac.deallocateSolverState();
        plog << "Presmoothing changes match fresh initialization\n";
        plog << "Caller data operations retain the full hierarchy range\n";
        laplace.deallocateOperatorState();
        for (int ln = 0; ln <= finest; ++ln)
        {
            for (int index : indices)
            {
                hierarchy->getPatchLevel(ln)->deallocatePatchData(index);
            }
        }
    }
    return 0;
}
