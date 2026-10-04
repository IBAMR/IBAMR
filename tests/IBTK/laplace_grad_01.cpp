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

#include <ibtk/AppInitializer.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SideNoCornersFillPattern.h>
#include <ibtk/muParserCartGridFunction.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideIterator.h>
#include <StandardTagAndInitialize.h>

#include <cmath>
#include <vector>

#include <ibtk/app_namespaces.h>

// A test program to check that HierarchyMathOps::laplace_grad() agrees with HierarchyMathOps::grad() followed by
// HierarchyMathOps::laplace() and yields the expected order of accuracy.

namespace
{
int
count_unequal_values(Pointer<PatchHierarchy<NDIM>> hierarchy, const int a_idx, const int b_idx)
{
    int count = 0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> a_data = patch->getPatchData(a_idx);
            Pointer<SideData<NDIM, double>> b_data = patch->getPatchData(b_idx);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> s(patch->getBox(), axis); s; s++)
                {
                    if ((*a_data)(s()) != (*b_data)(s())) ++count;
                }
            }
        }
    }
    return IBTK_MPI::sumReduction(count);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "laplace_grad.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
        Pointer<SideVariable<NDIM, double>> v_var = new SideVariable<NDIM, double>("v");
        // r holds the result of laplace_grad(), s the result of grad() followed by laplace(), f the exact result,
        // and e a difference.
        Pointer<SideVariable<NDIM, double>> r_var = new SideVariable<NDIM, double>("r");
        Pointer<SideVariable<NDIM, double>> s_var = new SideVariable<NDIM, double>("s");
        Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f");
        Pointer<SideVariable<NDIM, double>> e_var = new SideVariable<NDIM, double>("e");

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
        const int v_idx = var_db->registerVariableAndContext(v_var, ctx, IntVector<NDIM>(0));
        const int r_idx = var_db->registerVariableAndContext(r_var, ctx, IntVector<NDIM>(0));
        const int s_idx = var_db->registerVariableAndContext(s_var, ctx, IntVector<NDIM>(0));
        const int f_idx = var_db->registerVariableAndContext(f_var, ctx, IntVector<NDIM>(0));
        const int e_idx = var_db->registerVariableAndContext(e_var, ctx, IntVector<NDIM>(0));

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int level_number = 0;
        bool done = false;
        while (!done && gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, /*tag_buffer*/ 1);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int idx : { u_idx, p_idx, v_idx, r_idx, s_idx, f_idx, e_idx })
            {
                level->allocatePatchData(idx, 0.0);
            }
        }

        const auto set_data = [&](const std::string& name, const int idx, Pointer<Variable<NDIM>> var)
        {
            muParserCartGridFunction fcn(name, app_initializer->getComponentDatabase(name), grid_geometry);
            fcn.setDataOnPatchHierarchy(idx, var, patch_hierarchy, 0.0);
        };
        set_data("u", u_idx, u_var);
        set_data("p", p_idx, p_var);
        set_data("v", v_idx, v_var);
        set_data("f", f_idx, f_var);

        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        const std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        const std::vector<InterpolationTransactionComponent> transaction_comps = {
            { u_idx,
              "NONE",
              true,
              "CUBIC_COARSEN",
              "LINEAR",
              false,
              u_bc_coefs,
              new SideNoCornersFillPattern(1, /*overwrite_interior*/ false) },
            { p_idx,
              "NONE",
              true,
              "CUBIC_COARSEN",
              "LINEAR",
              false,
              nullptr,
              new CellNoCornersFillPattern(1, /*overwrite_interior*/ false) }
        };
        Pointer<HierarchyGhostCellInterpolation> ghost_fill = new HierarchyGhostCellInterpolation();
        ghost_fill->initializeOperatorState(transaction_comps, patch_hierarchy, 0, finest_ln);

        const double C = input_db->getDouble("C");
        const double D = input_db->getDouble("D");
        const double gamma = input_db->getDouble("gamma");
        const double delta = input_db->getDouble("delta");
        PoissonSpecifications poisson_spec("poisson_spec");
        poisson_spec.setCConstant(C);
        poisson_spec.setDConstant(D);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.laplace_grad(
            r_idx, r_var, poisson_spec, u_idx, u_var, gamma, p_idx, p_var, ghost_fill, 0.0, delta, v_idx, v_var);

        ghost_fill->fillData(0.0);
        hier_math_ops.grad(
            s_idx, s_var, /*cf_bdry_synch*/ true, gamma, p_idx, p_var, nullptr, 0.0, delta, v_idx, v_var);
        hier_math_ops.laplace(s_idx, s_var, poisson_spec, u_idx, u_var, nullptr, 0.0, 1.0, s_idx, s_var);

        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        HierarchySideDataOpsReal<NDIM, double> sc_data_ops(patch_hierarchy, 0, finest_ln);
        sc_data_ops.subtract(e_idx, r_idx, s_idx);
        plog << "|laplace_grad - laplace(grad)|_oo = " << std::abs(sc_data_ops.maxNorm(e_idx, h_sc_idx)) << "\n";
        plog << "values that differ = " << count_unequal_values(patch_hierarchy, r_idx, s_idx) << "\n";

        sc_data_ops.subtract(e_idx, r_idx, f_idx);
        plog << "|e|_oo = " << sc_data_ops.maxNorm(e_idx, h_sc_idx) << "\n";
        plog << "|e|_2  = " << sc_data_ops.L2Norm(e_idx, h_sc_idx) << "\n";
        plog << "|e|_1  = " << sc_data_ops.L1Norm(e_idx, h_sc_idx) << "\n";
    }
}
