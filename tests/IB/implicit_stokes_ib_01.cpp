// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBKernelEvaluators.h>
#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/IBOperatorRegistry.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <map>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int NUM_POINTS = 16;
int evaluator_calls = 0;

// An application-compiled evaluator with the supplied four-point kernel's weights.
struct CustomKernel
{
    std::array<double, 4> operator()(double r) const;
};

std::array<double, 4>
CustomKernel::operator()(const double r) const
{
    ++evaluator_calls;
    return IBKernelEvaluatorIB4{}(r);
}

void
count_integration_cycles(double /*current_time*/, double /*new_time*/, const int cycle_num, void* ctx)
{
    int* callback_count = static_cast<int*>(ctx);
    if (cycle_num != *callback_count)
    {
        TBOX_ERROR("Integration callback did not run exactly once per cycle.\n");
    }
    ++*callback_count;
}

struct RegridState
{
    int count = 0;
    int expected_levels = 1;
};

void
count_regrids(Pointer<BasePatchHierarchy<NDIM>> hierarchy, double /*data_time*/, bool /*initial_time*/, void* ctx)
{
    auto* state = static_cast<RegridState*>(ctx);
    if (hierarchy->getNumberOfLevels() != state->expected_levels)
    {
        TBOX_ERROR("Regrid did not retain the assigned hierarchy levels.\n");
    }
    ++state->count;
}

void
generate_structure(const unsigned int& structure,
                   const int& level,
                   int& num_vertices,
                   std::vector<IBTK::Point>& positions,
                   void* ctx)
{
    const int finest_level = *static_cast<const int*>(ctx);
    num_vertices = (structure == 0 && level == finest_level) ? NUM_POINTS : 0;
    positions.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * k / NUM_POINTS;
        positions[k](0) = 0.5 + 0.18 * std::cos(theta);
        positions[k](1) = 0.5 + 0.25 * std::sin(theta);
    }
}

void
generate_springs(
    const unsigned int& structure,
    const int& level,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>& specs,
    void* ctx)
{
    const int finest_level = *static_cast<const int*>(ctx);
    if (structure != 0 || level != finest_level)
    {
        return;
    }
    for (int k = 0; k < NUM_POINTS; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % NUM_POINTS };
        if (edge.first > edge.second)
        {
            std::swap(edge.first, edge.second);
        }
        spring_map.emplace(edge.first, edge);
        IBRedundantInitializer::SpringSpec spec;
        spec.force_fcn_idx = 0;
        spec.parameters = { 5.0, 0.0 };
        specs.emplace(edge, spec);
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    // Keep optional-library warnings and build paths out of compared output.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
        Pointer<Database> input = app->getInputDatabase();
        const bool custom_kernel = input->getBoolWithDefault("custom_kernel", false);
        if (custom_kernel)
        {
            IBOperatorRegistry::register_interpolation_matrix_sc(IBKernel("APP_FOUR_POINT"),
                                                                 IBKernelTensorProductEvaluator{ CustomKernel{} });
        }
        Pointer<INSStaggeredHierarchyIntegrator> ins = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator", app->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
        Pointer<IBImplicitStaggeredHierarchyIntegrator> integrator = new IBImplicitStaggeredHierarchyIntegrator(
            "IBHierarchyIntegrator", app->getComponentDatabase("IBHierarchyIntegrator"), method, ins);
        int finest_level = input->getIntegerWithDefault("MAX_LEVELS", 1) - 1;
        const bool verify_regrid = input->getBoolWithDefault("VERIFY_REGRID", false);
        RegridState regrid_state;
        regrid_state.expected_levels = finest_level + 1;
        if (verify_regrid)
        {
            integrator->registerRegridHierarchyCallback(count_regrids, &regrid_state);
        }
        int callback_count = 0;
        integrator->registerIntegrateHierarchyCallback(count_integration_cycles, &callback_count);
        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tagging = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", integrator, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), tagging, boxes, load_balancer);
        Pointer<IBRedundantInitializer> initializer =
            new IBRedundantInitializer("IBRedundantInitializer", app->getComponentDatabase("IBRedundantInitializer"));
        initializer->setStructureNamesOnLevel(finest_level, { "ellipse" });
        initializer->registerInitStructureFunction(generate_structure, &finest_level);
        initializer->registerInitSpringDataFunction(generate_springs, &finest_level);
        method->registerLInitStrategy(initializer);
        method->registerIBLagrangianForceFunction(new IBStandardForceGen());
        ins->registerVelocityInitialConditions(new muParserCartGridFunction(
            "initial_velocity", app->getComponentDatabase("VelocityInitialConditions"), geometry));
        integrator->initializePatchHierarchy(hierarchy, gridding);
        method->freeLInitStrategy();
        initializer.setNull();

        VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
        const int u = variables->mapVariableAndContextToIndex(ins->getVelocityVariable(), ins->getCurrentContext());
        const int p = variables->mapVariableAndContextToIndex(ins->getPressureVariable(), ins->getCurrentContext());
        HierarchySideDataOpsReal<NDIM, double> velocity_ops(hierarchy, 0, finest_level);
        HierarchyCellDataOpsReal<NDIM, double> pressure_ops(hierarchy, 0, finest_level);
        Pointer<HierarchyMathOps> math_ops = ins->getHierarchyMathOps();
        plog << std::scientific << std::setprecision(8);
        for (int step = 0; step < 3; ++step)
        {
            const int previous_calls = evaluator_calls;
            callback_count = 0;
            // A shorter final step also exercises timestep-dependent matrix setup.
            integrator->advanceHierarchy(step < 2 ? 0.001 : 0.0005);
            if (callback_count != integrator->getNumberOfCycles())
            {
                TBOX_ERROR("Integration callback count does not match the number of cycles.\n");
            }
            // Geometric norms are independent of redistribution's Lagrangian vector ordering.
            finest_level = hierarchy->getFinestLevelNumber();
            velocity_ops.resetLevels(0, finest_level);
            pressure_ops.resetLevels(0, finest_level);
            Vec positions = method->getLDataManager()->getLData("X", finest_level)->getVec();
            Vec centered_positions = nullptr;
            PetscErrorCode ierr = VecDuplicate(positions, &centered_positions);
            IBTK_CHKERRQ(ierr);
            ierr = VecCopy(positions, centered_positions);
            IBTK_CHKERRQ(ierr);
            ierr = VecShift(centered_positions, -0.5);
            IBTK_CHKERRQ(ierr);
            double radius_norm = 0.0;
            ierr = VecNorm(centered_positions, NORM_2, &radius_norm);
            IBTK_CHKERRQ(ierr);
            const double velocity_norm = velocity_ops.L2Norm(u, math_ops->getSideWeightPatchDescriptorIndex());
            const double pressure_norm = pressure_ops.L2Norm(p, math_ops->getCellWeightPatchDescriptorIndex());
            if (!std::isfinite(radius_norm) || radius_norm <= 0.0 || !std::isfinite(velocity_norm) ||
                !std::isfinite(pressure_norm))
            {
                TBOX_ERROR("Non-finite or trivial state after hierarchy advancement.\n");
            }
            ierr = VecDestroy(&centered_positions);
            IBTK_CHKERRQ(ierr);
            plog << "step " << step + 1 << " time " << integrator->getIntegratorTime() << " velocity_L2 "
                 << velocity_norm << " pressure_L2 " << pressure_norm << " radius_L2 " << radius_norm << '\n';
            if (custom_kernel && evaluator_calls == previous_calls)
            {
                TBOX_ERROR("The selected application evaluator was not used during advancement.\n");
            }
        }
        if (verify_regrid)
        {
            if (regrid_state.count != 3)
            {
                TBOX_ERROR("Expected the initial automatic regrid and two subsequent regrids.\n");
            }
            plog << "regrids " << regrid_state.count << " levels " << hierarchy->getNumberOfLevels() << '\n';
        }
    }
    return 0;
}
