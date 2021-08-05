// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAIVectorReal.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// Verify that we can correctly refine a piecewise linear solution with the
// RT0 refinement class. Since the RT0 element is a vector-valued element that
// is, on Cartesian grids,
//
//     RT0_K = (P^1(x) * P^0(y), P^0(x) * P^1(y))
//
// we expect refining a vector field which is in that space to have zero error.

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // this test only works in serial
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "rt0.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");
        Pointer<SideVariable<NDIM, double> > u_sc_var = new SideVariable<NDIM, double>("u_sc");
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx);
        Pointer<SideVariable<NDIM, double> > exact_sc_var = new SideVariable<NDIM, double>("exact_sc");
        const int exact_sc_idx = var_db->registerVariableAndContext(exact_sc_var, ctx);
        // u_cc_var is only for plotting (and testing): uncomment if output is desired
// #define DO_PLOT
#ifdef DO_PLOT
        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx);
        Pointer<CellVariable<NDIM, double> > exact_cc_var = new CellVariable<NDIM, double>("exact_cc", NDIM);
        const int exact_cc_idx = var_db->registerVariableAndContext(exact_cc_var, ctx);
#endif

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(exact_sc_idx, 0.0);
#ifdef DO_PLOT
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(exact_cc_idx, 0.0);
#endif
        }

        Pointer<VisItDataWriter<NDIM> > visit_writer = app_initializer->getVisItDataWriter();

        // The rest is just book-keeping, this is the actual test:
        auto do_test = [&](const std::string& db_u_fcn_name, const int coarse_level_n)
        {
            muParserCartGridFunction u_fcn(
                db_u_fcn_name, app_initializer->getComponentDatabase(db_u_fcn_name), grid_geometry);
            u_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
            u_fcn.setDataOnPatchHierarchy(exact_sc_idx, exact_sc_var, patch_hierarchy, 0.0);

            solv::SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, 1);
            u_vec.addComponent(u_sc_var, u_sc_idx);

            const int fine_level_n = coarse_level_n + 1;
            Pointer<PatchLevel<NDIM> > level_0 = patch_hierarchy->getPatchLevel(coarse_level_n);
            Pointer<PatchLevel<NDIM> > level_1 = patch_hierarchy->getPatchLevel(fine_level_n);

            // there should only be one patch on each patch level
            Pointer<SideData<NDIM, double> > u_sc_0_data = level_0->getPatch(0)->getPatchData(u_sc_idx);
            Pointer<SideData<NDIM, double> > u_sc_1_data = level_1->getPatch(0)->getPatchData(u_sc_idx);
            const Box<NDIM> patch_box_0 = level_0->getPatch(0)->getBox();
            const Box<NDIM> patch_box_1 = level_1->getPatch(0)->getBox();

            Pointer<SideData<NDIM, double> > exact_sc_0_data = level_0->getPatch(0)->getPatchData(exact_sc_idx);
            Pointer<SideData<NDIM, double> > exact_sc_1_data = level_1->getPatch(0)->getPatchData(exact_sc_idx);

            const IntVector<NDIM> ratio = level_1->getRatioToCoarserLevel();
            IBTK::CartSideDoubleRT0Refine refine_op;
            refine_op.refine(*level_1->getPatch(0), *level_0->getPatch(0), u_sc_idx, u_sc_idx, patch_box_1, ratio);

            solv::SAMRAIVectorReal<NDIM, double> exact_vec("e", patch_hierarchy, coarse_level_n, fine_level_n);
            exact_vec.addComponent(exact_sc_var, exact_sc_idx);
            exact_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&u_vec, false),
                               Pointer<SAMRAIVectorReal<NDIM, double> >(&exact_vec, false));

            pout << "test results for " << db_u_fcn_name << '\n';
            pout << "max norm of u_sc: " << u_vec.maxNorm() << '\n';
            pout << "max norm of exact - refined: " << std::abs(exact_vec.maxNorm()) << '\n';
        };

        do_test("constant_function", 0);
        do_test("linear_x", 0);
        do_test("linear_y", 0);
        if (NDIM == 3)
        {
            do_test("linear_z", 0);
        }
        do_test("both_linear", 0);

        for (int coarse_level_n = 0; coarse_level_n < gridding_algorithm->getMaxLevels() - 1; ++coarse_level_n)
        {
            do_test("non_rt", coarse_level_n);
        }

#ifdef DO_PLOT
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        visit_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, nullptr, 0.0, false);
        visit_writer->registerPlotQuantity(exact_cc_var->getName(), "VECTOR", exact_cc_idx);
        hier_math_ops.interp(exact_cc_idx, exact_cc_var, exact_sc_idx, exact_sc_var, nullptr, 0.0, false);
#endif

#ifdef DO_PLOT
        visit_writer->writePlotData(patch_hierarchy, 0, 0.0);
#endif
    }
}
