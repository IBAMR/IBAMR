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
#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ib_kernel_evaluators.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Array.h>
#include <tbox/Logger.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <BoundaryBox.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideIndex.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int NUM_POINTS = 16;

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
    int finest_level = 0;
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
    Vec ordinary_vec = (*ordinary)[data->finest_level]->getVec();
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

// The level and the center and radii of the elliptical structure.
struct StructureSpec
{
    int finest_level = 0;
    std::array<double, 2> center = { 0.5, 0.5 }, radii = { 0.18, 0.25 };
};

void
generate_structure(const unsigned int& structure,
                   const int& level,
                   int& num_vertices,
                   std::vector<IBTK::Point>& positions,
                   void* ctx)
{
    const StructureSpec& spec = *static_cast<const StructureSpec*>(ctx);
    num_vertices = (structure == 0 && level == spec.finest_level) ? NUM_POINTS : 0;
    positions.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * k / NUM_POINTS;
        positions[k](0) = spec.center[0] + spec.radii[0] * std::cos(theta);
        positions[k](1) = spec.center[1] + spec.radii[1] * std::sin(theta);
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
    if (structure != 0 || level != static_cast<const StructureSpec*>(ctx)->finest_level)
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

// Return the largest magnitude of the velocity components on the physical boundary normal to it, after requiring that
// they equal the values that bc_coefs prescribe at data_time.
double
check_boundary_normal_velocity(const int u_idx,
                               Pointer<PatchHierarchy<NDIM>> hierarchy,
                               const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                               const double data_time)
{
    double max_velocity = 0.0, max_error = 0.0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            const SAMRAI::tbox::Array<BoundaryBox<NDIM>> boundary_boxes =
                PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            for (int n = 0; n < boundary_boxes.size(); ++n)
            {
                const int axis = boundary_boxes[n].getLocationIndex() / 2;
                const BoundaryBox<NDIM> trimmed_box =
                    PhysicalBoundaryUtilities::trimBoundaryCodim1Box(boundary_boxes[n], *patch);
                const Box<NDIM> coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_box);
                Pointer<ArrayData<NDIM, double>> acoef = new ArrayData<NDIM, double>(coef_box, 1);
                Pointer<ArrayData<NDIM, double>> bcoef = new ArrayData<NDIM, double>(coef_box, 1);
                Pointer<ArrayData<NDIM, double>> gcoef = new ArrayData<NDIM, double>(coef_box, 1);
                bc_coefs[axis]->setBcCoefs(acoef, bcoef, gcoef, nullptr, *patch, trimmed_box, data_time);
                for (Box<NDIM>::Iterator b(coef_box); b; b++)
                {
                    const double velocity = (*u_data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
                    max_velocity = std::max(max_velocity, std::abs(velocity));
                    max_error = std::max(max_error, std::abs(velocity - (*gcoef)(b(), 0) / (*acoef)(b(), 0)));
                }
            }
        }
    }
    if (!std::isfinite(max_error) || max_error > 1.0e-12 * std::max(1.0, max_velocity))
    {
        TBOX_ERROR("The boundary normal velocity differs from its prescribed value by " << max_error << ".\n");
    }
    return max_velocity;
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
            integrator->registerJacobianOperatorBuilder(
                IBKernelTensorProduct("APP_CUSTOM"),
                IBOperatorBuilder(
                    IBKernelEvaluatorTensorProduct{ CustomKernel{ 1.0, lifetime }, CustomKernel{ 1.0, lifetime } }));
        }
        int finest_level = input->getIntegerWithDefault("MAX_LEVELS", 1) - 1;
        StructureSpec structure_spec;
        structure_spec.finest_level = finest_level;
        if (input->keyExists("STRUCTURE_CENTER"))
        {
            input->getDoubleArray("STRUCTURE_CENTER", structure_spec.center.data(), 2);
        }
        const bool verify_regrid = input->getBoolWithDefault("VERIFY_REGRID", false);
        RegridState regrid_state;
        regrid_state.expected_levels = finest_level + 1;
        if (verify_regrid)
        {
            integrator->registerRegridHierarchyCallback(count_regrids, &regrid_state);
        }
        CycleData cycles{ 0, method.getPointer(), finest_level };
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
        initializer->setStructureNamesOnLevel(finest_level, { "ellipse" });
        initializer->registerInitStructureFunction(generate_structure, &structure_spec);
        initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
        method->registerLInitStrategy(initializer);
        method->registerIBLagrangianForceFunction(new IBStandardForceGen());
        ins->registerVelocityInitialConditions(new muParserCartGridFunction(
            "initial_velocity", app->getComponentDatabase("VelocityInitialConditions"), geometry));
        // A nonperiodic domain prescribes the velocity on its physical boundary.
        const bool physical_boundary = geometry->getPeriodicShift().min() == 0;
        std::vector<std::unique_ptr<muParserRobinBcCoefs>> u_bc_coef_objects;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        if (physical_boundary)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                u_bc_coef_objects.push_back(std::make_unique<muParserRobinBcCoefs>(
                    "u_bc_coefs_" + std::to_string(d),
                    app->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                    geometry));
                u_bc_coefs[d] = u_bc_coef_objects.back().get();
            }
            ins->registerPhysicalBoundaryConditions(u_bc_coefs);
        }
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
        HierarchySideDataOpsReal<NDIM, double> velocity_ops(hierarchy, 0, finest_level);
        HierarchyCellDataOpsReal<NDIM, double> pressure_ops(hierarchy, 0, finest_level);
        Pointer<HierarchyMathOps> math_ops = ins->getHierarchyMathOps();
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
            if (custom_kernel && evaluator_lifetime.expired())
            {
                TBOX_ERROR("The integrator did not retain the application evaluator.\n");
            }
        }
        if (physical_boundary)
        {
            plog << "max_boundary_normal_velocity "
                 << check_boundary_normal_velocity(u, hierarchy, u_bc_coefs, integrator->getIntegratorTime()) << '\n';
        }
        if (verify_regrid)
        {
            if (regrid_state.count != 3)
            {
                TBOX_ERROR("Expected the initial automatic regrid and two subsequent regrids.\n");
            }
            plog << "regrids " << regrid_state.count << " levels " << hierarchy->getNumberOfLevels() << '\n';
        }
        if (cycles.previous_ordinary) IBTK_CHKERRQ(VecDestroy(&cycles.previous_ordinary));
    }
    return 0;
}
