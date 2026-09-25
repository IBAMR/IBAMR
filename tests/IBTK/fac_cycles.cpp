// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
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

#include <ibtk/app_namespaces.h>

namespace
{
using HierarchyVector = SAMRAIVectorReal<NDIM, double>;

// Count the actual traversal without substituting any numerical operation.
class CountingFACOperator : public CCPoissonPointRelaxationFACOperator
{
public:
    CountingFACOperator(Pointer<Database> db, int levels)
        : CCPoissonPointRelaxationFACOperator("fac_operator", db, ""), d_visits(levels, 0)
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
};

void
require(bool condition, const std::string& message)
{
    if (!condition) TBOX_ERROR("FAC cycle test: " << message << '\n');
}

double
data_difference(Pointer<PatchHierarchy<NDIM>> hierarchy, int first, int second, int weight, bool active_only)
{
    double maximum = 0.0;
    int nonfinite = 0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            auto patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> x = patch->getPatchData(first), y = patch->getPatchData(second),
                                            w = patch->getPatchData(weight);
            for (CellIterator<NDIM> c(active_only ? patch->getBox() : x->getGhostBox()); c; c++)
            {
                if (active_only && (*w)(c()) == 0.0) continue;
                const double value = std::abs((*x)(c()) - (*y)(c()));
                nonfinite += !std::isfinite(value);
                maximum = std::max(maximum, value);
            }
        }
    }
    require(IBTK_MPI::sumReduction(nonfinite) == 0, "nonfinite comparison data");
    return IBTK_MPI::maxReduction(maximum);
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
            if (!hierarchy->finerLevelExists(ln)) break;
        }
        const int finest = hierarchy->getFinestLevelNumber();
        require(finest + 1 == input->getInteger("levels"), "incomplete test hierarchy");

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
        auto& ops = *shared_ops;
        std::array<int, FIELD_COUNT> indices;
        std::array<std::unique_ptr<HierarchyVector>, FIELD_COUNT> vectors;
        Pointer<CellVariable<NDIM, double>> variable = new CellVariable<NDIM, double>("u");
        auto* variable_db = VariableDatabase<NDIM>::getDatabase();
        for (int f = 0; f < FIELD_COUNT; ++f)
        {
            const std::string name = "field_" + std::to_string(f);
            indices[f] =
                variable_db->registerVariableAndContext(variable, variable_db->getContext(name), IntVector<NDIM>(1));
            for (int ln = 0; ln <= finest; ++ln) hierarchy->getPatchLevel(ln)->allocatePatchData(indices[f], 0.0);
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
                    for (CellIterator<NDIM> c(data->getGhostBox()); c; c++) unchanged += (*data)(c()) != 1.0;
                }
            }
            require(IBTK_MPI::sumReduction(unchanged) == 0, "FAC changed the caller's data-operation level range");
        };

        const double pi = std::acos(-1.0);
        for (int ln = 0; ln <= finest; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            const double n = input->getInteger("N") * level->getRatio()(0);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                auto patch = level->getPatch(p());
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
        LocationIndexRobinBcCoefs<NDIM> boundary("boundary", nullptr);
        for (int d = 0; d < NDIM; ++d)
        {
            boundary.setBoundaryValue(2 * d, 0.0);
            boundary.setBoundaryValue(2 * d + 1, 0.0);
        }
        PoissonSpecifications specification("specification");
        specification.setCConstant(1.0);
        specification.setDConstant(-1.0);
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
            require(difference(source, SNAPSHOT, legacy_v) == 0.0, "RHS data was modified");
            copy(source, SNAPSHOT);
        };
        const auto residual = [&]()
        {
            copy(EVALUATION, ITERATE);
            laplace.apply(*vectors[EVALUATION], *vectors[RESIDUAL]);
            ops.subtract(indices[RESIDUAL], indices[RHS], indices[RESIDUAL], false);
        };
        struct CycleCase
        {
            MGCycleType type;
            const char* name;
            int multiplicity;
        };
        const std::array<CycleCase, 4> cycles = {
            { { V_CYCLE, "V", 1 }, { W_CYCLE, "W", 2 }, { MU_CYCLE, "MU3", 3 }, { F_CYCLE, "F", 2 } }
        };
        plog << "levels = " << finest + 1 << '\n';
        for (int pre : { 0, 2 })
        {
            fac.setNumPreSmoothingSweeps(pre);
            fac.setNumPostSmoothingSweeps(2);
            fac.initializeSolverState(*vectors[CORRECTION], *vectors[RHS]);
            check_data_ops();
            plog << "pre_sweeps = " << pre << '\n';
            for (const auto& cycle : cycles)
            {
                fac.setMGCycleType(cycle.type);
                fac.setMGCycleMultiplicity(cycle.multiplicity);
                require(fac.getMGCycleMultiplicity() == cycle.multiplicity, "multiplicity accessor mismatch");
                apply(RESULT, RHS);
                plog << cycle.name << " visits =";
                for (int ln = finest; ln >= 0; --ln)
                {
                    const int depth = finest - ln;
                    int expected = cycle.type == F_CYCLE ? depth + 1 : 1;
                    if (cycle.type != F_CYCLE)
                        for (int i = 0; i < depth; ++i) expected *= cycle.multiplicity;
                    require(strategy->getVisits()[ln] == expected, "incorrect cycle traversal");
                    plog << ' ' << strategy->getVisits()[ln];
                }
                plog << '\n';
                if (cycle.type == V_CYCLE) copy(V_RESULT, RESULT);
                if (cycle.type == W_CYCLE) copy(W_RESULT, RESULT);
                copy(ITERATE, RESULT);
                residual();
                plog << cycle.name << " first correction L2 = " << vectors[RESULT]->L2Norm()
                     << "; residual L2 = " << vectors[RESIDUAL]->L2Norm() << '\n';

                apply(SECOND_RESULT, SECOND_RHS);
                apply(THIRD_RESULT, COMBINED_RHS);
                ops.linearSum(indices[ERROR], 0.31, indices[RESULT], -0.7, indices[SECOND_RESULT], false);
                require(difference(THIRD_RESULT, ERROR) < 1.0e-11, "FAC map is not linear");

                // Change both vector indices after initialization, then rebuild solver state.
                copy(COMBINED_RHS, RHS);
                apply(THIRD_RESULT, COMBINED_RHS);
                require(difference(THIRD_RESULT, RESULT) < 1.0e-12, "solve depends on vector data indices");
                fac.deallocateSolverState();
                fac.initializeSolverState(*vectors[SECOND_RESULT], *vectors[SECOND_RHS]);
                check_data_ops();
                apply(THIRD_RESULT, RHS);
                require(difference(THIRD_RESULT, RESULT) < 1.0e-12, "reinitialization changed the correction");
                ops.linearSum(indices[COMBINED_RHS], 0.31, indices[RHS], -0.7, indices[SECOND_RHS], false);

                ops.setToScalar(indices[ERROR], 0.0, false);
                apply(CORRECTION, ERROR);
                require(difference(CORRECTION, ERROR) == 0.0, "zero RHS produced a nonzero correction");
                copy(ITERATE, EXACT);
                for (int ln = 0; ln < finest; ++ln)
                {
                    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        auto patch = level->getPatch(p());
                        Pointer<CellData<NDIM, double>> x = patch->getPatchData(indices[ITERATE]),
                                                        w = patch->getPatchData(weight);
                        for (CellIterator<NDIM> c(patch->getBox()); c; c++)
                            if ((*w)(c()) == 0.0) (*x)(c()) += 0.125;
                    }
                }
                residual();
                require(vectors[RESIDUAL]->L2Norm() < 1.0e-11, "covered storage changed the exact residual");
                apply(CORRECTION, RESIDUAL);
                ops.add(indices[ITERATE], indices[ITERATE], indices[CORRECTION], false);
                require(difference(ITERATE, EXACT) < 1.0e-11, "exact solution is not a fixed point");

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
                    require(final_residual < initial_residual, "composite residual did not decrease");
                    require(final_error < initial_error, "solution error did not decrease");
                    plog << cycle.name << " start = " << (starting_scale == 0.0 ? "zero" : "nonzero")
                         << "; final residual ratio = " << final_residual / initial_residual
                         << "; final error ratio = " << final_error / initial_error << '\n';
                }
                plog << cycle.name << " invariants, linearity, reuse, and convergence passed\n";
            }
            fac.setMGCycleType(MU_CYCLE);
            fac.setMGCycleMultiplicity(1);
            apply(CORRECTION, RHS);
            require(difference(CORRECTION, V_RESULT) < 1.0e-12, "mu=1 differs from V");
            fac.setMGCycleMultiplicity(2);
            apply(CORRECTION, RHS);
            require(difference(CORRECTION, W_RESULT) < 1.0e-12, "mu=2 differs from W");
            plog << "mu=1 matches V; mu=2 matches W\n";
            fac.setMGCycleType(FMG_CYCLE);
            apply(RESULT, RHS);
            copy(ITERATE, RESULT);
            residual();
            require(vectors[RESIDUAL]->L2Norm() < vectors[RHS]->L2Norm(), "FMG did not reduce the residual");
            ops.setToScalar(indices[ERROR], 0.0, false);
            apply(CORRECTION, ERROR);
            require(difference(CORRECTION, ERROR) == 0.0, "FMG zero RHS produced a nonzero correction");
            plog << "FMG zero-input and nested initialization passed\n";
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
        require(difference(CORRECTION, RESULT) < 1.0e-12, "disabling presmoothing changed the V correction");
        fac.deallocateSolverState();
        fac.setNumPreSmoothingSweeps(2);
        fac.initializeSolverState(*vectors[SECOND_RESULT], *vectors[SECOND_RHS]);
        check_data_ops();
        apply(CORRECTION, RHS);
        require(difference(CORRECTION, SECOND_RESULT) < 1.0e-12,
                "changing presmoothing differs from fresh initialization");
        fac.deallocateSolverState();
        plog << "Presmoothing changes match fresh initialization\n";
        plog << "Caller data operations retain the full hierarchy range\n";
        laplace.deallocateOperatorState();
        for (int ln = 0; ln <= finest; ++ln)
            for (int index : indices) hierarchy->getPatchLevel(ln)->deallocatePatchData(index);
    }
    return 0;
}
