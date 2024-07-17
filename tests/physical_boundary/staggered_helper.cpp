// ---------------------------------------------------------------------
//
// Copyright (c) 2024 - 2024 by the IBAMR developers
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
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/StaggeredPhysicalBoundaryHelper.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Array.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <ibtk/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown
        // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitializeNd> error_detector = new StandardTagAndInitializeNd(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        Pointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

        // Create boundary conditions object
        std::vector<RobinBcCoefStrategyNd*> bc_coefs(NDIM, nullptr);
        for (int d = 0; d < NDIM; ++d)
        {
            std::string bc_coef_name = "bc_" + std::to_string(d);
            bc_coefs[d] = new muParserRobinBcCoefs(bc_coef_name, input_db->getDatabase(bc_coef_name), grid_geometry);
        }

        // Create side centered variable
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CONTEXT");
        Pointer<SideVariableNd<double> > s0_var = new SideVariableNd<double>("s0_v");
        Pointer<SideVariableNd<double> > s1_var = new SideVariableNd<double>("s1_v");
        Pointer<SideVariableNd<int> > mask_var = new SideVariableNd<int>("mask");
        const int gcw = 1;
        const int s0_idx = var_db->registerVariableAndContext(s0_var, context, gcw);
        const int s1_idx = var_db->registerVariableAndContext(s1_var, context, gcw);
        const int mask_idx = var_db->registerVariableAndContext(mask_var, context);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(s0_idx);
            level->allocatePatchData(s1_idx);
            level->allocatePatchData(mask_idx);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<SideDataNd<double> > s0_data = patch->getPatchData(s0_idx);
                s0_data->fillAll(std::numeric_limits<double>::quiet_NaN());
                Pointer<SideDataNd<double> > s1_data = patch->getPatchData(s1_idx);
                s1_data->fillAll(0.0);
            }
        }

        StaggeredPhysicalBoundaryHelper phys_bc_helper;
        phys_bc_helper.cacheBcCoefData(bc_coefs, 0.0, patch_hierarchy);

        plog << "Testing masking function\n";
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<SideDataNd<double> > s0_data = patch->getPatchData(s0_idx);
                Pointer<SideDataNd<double> > s1_data = patch->getPatchData(s1_idx);
                phys_bc_helper.copyDataAtDirichletBoundaries(s0_data, s1_data, patch);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIteratorNd si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndexNd& idx = si();
                        if (!std::isnan((*s0_data)(idx)))
                            pout << "Copied value on index " << idx << " and axis " << axis << "\n";
                    }
                }
            }
        }

        plog << "Testing dirichlet boundaries\n";
        phys_bc_helper.setupMaskingFunction(mask_idx, coarsest_ln, finest_ln);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                Pointer<PatchNd> patch = level->getPatch(p());
                Pointer<SideDataNd<int> > mask_data = patch->getPatchData(mask_idx);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIteratorNd si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndexNd& idx = si();
                        if ((*mask_data)(idx) == 1)
                            pout << "Cell " << idx << " and axis " << axis << " touches boundary\n";
                    }
                }
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(mask_idx);
            level->deallocatePatchData(s1_idx);
            level->deallocatePatchData(s0_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown

    return 0;
} // main
