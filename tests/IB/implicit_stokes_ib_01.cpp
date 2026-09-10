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
        spec.parameters = { 5.0, 0.0 };
        specs.emplace(edge, spec);
    }
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
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
        initializer->setStructureNamesOnLevel(0, { "ellipse" });
        initializer->registerInitStructureFunction(generate_structure);
        initializer->registerInitSpringDataFunction(generate_springs);
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
        HierarchySideDataOpsReal<NDIM, double> velocity_ops(hierarchy, 0, 0);
        HierarchyCellDataOpsReal<NDIM, double> pressure_ops(hierarchy, 0, 0);
        Pointer<HierarchyMathOps> math_ops = ins->getHierarchyMathOps();
        Vec centered_positions = nullptr;
        PetscErrorCode ierr = VecDuplicate(method->getLDataManager()->getLData("X", 0)->getVec(), &centered_positions);
        IBTK_CHKERRQ(ierr);
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
            if (custom_kernel && evaluator_calls == previous_calls)
            {
                TBOX_ERROR("The selected application evaluator was not used during advancement.\n");
            }
        }
        ierr = VecDestroy(&centered_positions);
        IBTK_CHKERRQ(ierr);
    }
    return 0;
}
