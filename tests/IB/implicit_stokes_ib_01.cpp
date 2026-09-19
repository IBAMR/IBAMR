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

#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/ib_kernel_evaluators.h>
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
#include <memory>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int NUM_POINTS = 16;

// One fixed level: regridding and multilevel hierarchies are not exercised.

// An application-compiled evaluator with the supplied four-point kernel's weights.
struct CustomKernel
{
    CustomKernel(double scale, std::shared_ptr<const int> lifetime);
    static constexpr std::size_t get_stencil_width();

    template <typename Output, std::floating_point Input>
    Output evaluate(const Input& r) const requires requires(const Input& x)
    {
        IBKernelEvaluators::IB4{}.template evaluate<Output>(x);
    };
    std::unique_ptr<const double> scale;
    std::shared_ptr<const int> lifetime;
};

CustomKernel::CustomKernel(double value, std::shared_ptr<const int> token)
    : scale(std::make_unique<const double>(value)), lifetime(std::move(token))
{
}

struct CycleData
{
    int count = 0;
    IBMethod* method = nullptr;
    // The converged "ordinary" (Newton-solved) new-time coupling positions from the end of the previous
    // cycle. reinitializeOperatorsAndSolvers() reruns updateFixedLEOperators() at the start of every cycle
    // (not just once per timestep), which freezes whatever the previous cycle left as the new-time position
    // guess; so a later cycle's fixed positions must equal this snapshot exactly.
    Vec previous_ordinary = nullptr;
};

constexpr std::size_t
CustomKernel::get_stencil_width()
{
    return 4;
}

template <typename Output, std::floating_point Input>
Output
CustomKernel::evaluate(const Input& r) const requires requires(const Input& x)
{
    IBKernelEvaluators::IB4{}.template evaluate<Output>(x);
}
{
    Output values = IBKernelEvaluators::IB4{}.template evaluate<Output>(r);
    using Coefficient = ib_kernel_weights_value_t<Output>;
    const Coefficient factor = *scale;
    for (std::size_t i = 0; i < get_stencil_width(); ++i)
    {
        values[i] *= factor;
    }
    return values;
}

void
count_integration_cycles(double /*current_time*/, double new_time, const int cycle_num, void* ctx)
{
    CycleData* data = static_cast<CycleData*>(ctx);
    int* callback_count = &data->count;
    if (cycle_num != *callback_count)
    {
        TBOX_ERROR("Integration callback did not run exactly once per cycle.\n");
    }
    ++*callback_count;
    std::vector<Pointer<LData>>* ordinary = nullptr;
    bool* needs_fill = nullptr;
    data->method->getPositionData(&ordinary, &needs_fill, IBTK::TimePoint::NEW_TIME);
    Vec fixed = data->method->getFinestLevelLECouplingPositions(new_time);
    Vec ordinary_vec = (*ordinary)[0]->getVec();
    if (cycle_num > 0)
    {
        Vec drift = nullptr;
        IBTK_CHKERRQ(VecDuplicate(fixed, &drift));
        IBTK_CHKERRQ(VecWAXPY(drift, -1.0, fixed, data->previous_ordinary));
        double drift_norm = 0.0;
        IBTK_CHKERRQ(VecNorm(drift, NORM_2, &drift_norm));
        IBTK_CHKERRQ(VecDestroy(&drift));
        if (!(drift_norm < 1.0e-10))
        {
            TBOX_ERROR("Fixed coupling positions did not freeze the previous cycle's converged positions.\n");
        }
    }
    // Secondary sanity guard: confirm the fixed and ordinary positions are actually independent quantities
    // (not, say, both accidentally reading the same underlying data).
    Vec difference = nullptr;
    IBTK_CHKERRQ(VecDuplicate(fixed, &difference));
    IBTK_CHKERRQ(VecWAXPY(difference, -1.0, fixed, ordinary_vec));
    double norm = 0.0;
    IBTK_CHKERRQ(VecNorm(difference, NORM_2, &norm));
    IBTK_CHKERRQ(VecDestroy(&difference));
    if (!(norm > 1.0e-10))
    {
        TBOX_ERROR("Fixed coupling positions were not distinguished from updated positions.\n");
    }
    // Record this cycle's converged positions as the value the next cycle's fixed positions must reproduce.
    if (data->previous_ordinary) IBTK_CHKERRQ(VecDestroy(&data->previous_ordinary));
    IBTK_CHKERRQ(VecDuplicate(ordinary_vec, &data->previous_ordinary));
    IBTK_CHKERRQ(VecCopy(ordinary_vec, data->previous_ordinary));
}

void
generate_structure(const unsigned int& structure,
                   const int& level,
                   int& num_vertices,
                   std::vector<IBTK::Point>& positions,
                   void*)
{
    num_vertices = (structure == 0 && level == 0) ? NUM_POINTS : 0;
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
    void*)
{
    if (structure != 0 || level != 0)
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
        spec.parameters = { 5.0, 0.05 };
        specs.emplace(edge, spec);
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without Silo.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
        Pointer<Database> input = app->getInputDatabase();
        const bool custom_kernel = input->getBoolWithDefault("custom_kernel", false);
        Pointer<INSStaggeredHierarchyIntegrator> ins = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator", app->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
        Pointer<IBImplicitStaggeredHierarchyIntegrator> integrator = new IBImplicitStaggeredHierarchyIntegrator(
            "IBHierarchyIntegrator", app->getComponentDatabase("IBHierarchyIntegrator"), method, ins);
        std::weak_ptr<const int> evaluator_lifetime;
        if (custom_kernel)
        {
            auto lifetime = std::make_shared<const int>(1);
            evaluator_lifetime = lifetime;
            integrator->setJacobianOperatorBuilder(IBOperatorBuilder(
                IBKernelEvaluatorTensorProduct{ CustomKernel{ 1.0, lifetime }, CustomKernel{ 1.0, lifetime } }));
        }
        CycleData cycles{ 0, method.getPointer() };
        integrator->registerIntegrateHierarchyCallback(count_integration_cycles, &cycles);
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
        initializer->setStructureNamesOnLevel(0, { "ellipse" });
        initializer->registerInitStructureFunction(generate_structure);
        initializer->registerInitSpringDataFunction(generate_springs);
        method->registerLInitStrategy(initializer);
        method->registerIBLagrangianForceFunction(new IBStandardForceGen());
        ins->registerVelocityInitialConditions(new muParserCartGridFunction(
            "initial_velocity", app->getComponentDatabase("VelocityInitialConditions"), geometry));
        integrator->initializePatchHierarchy(hierarchy, gridding);
        if (input->getBoolWithDefault("late_kernel", false))
        {
            Pointer<Logger::Appender> abort_appender = new TestAppender();
            Logger::getInstance()->setAbortAppender(abort_appender);
            integrator->setJacobianOperatorBuilder(
                IBOperatorBuilder(IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB4{} }));
            return 0;
        }
        method->freeLInitStrategy();
        initializer.setNull();

        VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
        const int u = variables->mapVariableAndContextToIndex(ins->getVelocityVariable(), ins->getCurrentContext());
        const int p = variables->mapVariableAndContextToIndex(ins->getPressureVariable(), ins->getCurrentContext());
        HierarchySideDataOpsReal<NDIM, double> velocity_ops(hierarchy, 0, 0);
        HierarchyCellDataOpsReal<NDIM, double> pressure_ops(hierarchy, 0, 0);
        Pointer<HierarchyMathOps> math_ops = ins->getHierarchyMathOps();
        Vec centered_positions = nullptr;
        PetscErrorCode ierr = VecDuplicate(method->getLDataManager()->getLData("X", 0)->getVec(), &centered_positions);
        IBTK_CHKERRQ(ierr);
        // The Jacobian kernel enters only the convergence of the nonlinear solve, so report what the kernel requires
        // and the input key that selected it, since fixtures with the same ghost width would otherwise be
        // indistinguishable.
        plog << "jacobian kernel = " << input->getStringWithDefault("JAC_KERNEL", "IB_4")
             << ", ghost width = " << integrator->getJacobianOperatorBuilder().getMinimumGhostWidth() << '\n';
        plog << std::scientific << std::setprecision(8);
        const double dt = input->getDoubleWithDefault("dt", 0.001);
        for (int step = 0; step < 3; ++step)
        {
            cycles.count = 0;
            // A shorter final step also exercises timestep-dependent matrix setup.
            integrator->advanceHierarchy(step < 2 ? dt : 0.5 * dt);
            if (cycles.count != integrator->getNumberOfCycles())
            {
                TBOX_ERROR("Integration callback count does not match the number of cycles.\n");
            }
            // Geometric norms are independent of redistribution's Lagrangian vector ordering.
            ierr = VecCopy(method->getLDataManager()->getLData("X", 0)->getVec(), centered_positions);
            IBTK_CHKERRQ(ierr);
            ierr = VecShift(centered_positions, -0.5);
            IBTK_CHKERRQ(ierr);
            double radius_norm = 0.0;
            ierr = VecNorm(centered_positions, NORM_2, &radius_norm);
            IBTK_CHKERRQ(ierr);
            plog << "step " << step + 1 << " time " << integrator->getIntegratorTime() << " velocity_L2 "
                 << velocity_ops.L2Norm(u, math_ops->getSideWeightPatchDescriptorIndex()) << " pressure_L2 "
                 << pressure_ops.L2Norm(p, math_ops->getCellWeightPatchDescriptorIndex()) << " radius_L2 "
                 << radius_norm << '\n';
            if (custom_kernel && evaluator_lifetime.expired())
            {
                TBOX_ERROR("The integrator did not retain the application evaluator.\n");
            }
        }
        ierr = VecDestroy(&centered_positions);
        IBTK_CHKERRQ(ierr);
        if (cycles.previous_ordinary) IBTK_CHKERRQ(VecDestroy(&cycles.previous_ordinary));
    }
    return 0;
}
