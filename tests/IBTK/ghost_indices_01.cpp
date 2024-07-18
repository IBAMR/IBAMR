// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
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
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <ibtk/app_namespaces.h>

// test stuff
#include "../tests.h"

// test for setting up DoF indices via PETSc.

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        SAMRAIPointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);
        SAMRAIPointer<StandardTagAndInitializeNd> error_detector = new StandardTagAndInitializeNd(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        SAMRAIPointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        // TODO: make depth, cell vs side, ghost width input arguments
        const int u_depth = input_db->getIntegerWithDefault("depth", 1);
        const IntVectorNd gcw(input_db->getIntegerWithDefault("ghost_width", 1));
        const bool use_cell = input_db->getStringWithDefault("var_type", "CELL") == "CELL";
        // disambiguate between libMesh::Variable and SAMRAI::hier::Variable
        auto u_var = use_cell ? SAMRAIPointer<hier::VariableNd>(new CellVariableNd<double>("u", u_depth)) :
                                SAMRAIPointer<hier::VariableNd>(new SideVariableNd<double>("u", u_depth));
        auto dof_var = use_cell ? SAMRAIPointer<hier::VariableNd>(new CellVariableNd<int>("dof", u_depth)) :
                                  SAMRAIPointer<hier::VariableNd>(new SideVariableNd<int>("dof", u_depth));

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, gcw);
        const int dof_idx = var_db->registerVariableAndContext(dof_var, ctx, gcw);

        // set up grid
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while ((gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }
        const int level = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= level; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_idx, 0.0);
            level->allocatePatchData(dof_idx, 0.0);
        }

        std::vector<int> n_dofs_per_proc;
        PETScVecUtilities::constructPatchLevelDOFIndices(
            n_dofs_per_proc, dof_idx, patch_hierarchy->getPatchLevel(level));

        plog << "dofs per proc:\n";
        for (const int n_dofs : n_dofs_per_proc) plog << n_dofs << '\n';

        std::ostringstream out;
        if (IBTK_MPI::getNodes() != 1)
        {
            // partitioning is only relevant when there are multiple processors
            SAMRAIPointer<PatchLevelNd> patch_level =
                patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
            const BoxArrayNd boxes = patch_level->getBoxes();
            plog << "hierarchy boxes:\n";
            for (int i = 0; i < boxes.size(); ++i) plog << boxes[i] << '\n';
            // rank is only relevant when there are multiple processors
            out << "\nrank: " << IBTK_MPI::getRank() << '\n';
        }

        SAMRAIPointer<PatchLevelNd> patch_level = patch_hierarchy->getPatchLevel(level);
        for (PatchLevelNd::Iterator p(patch_level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = patch_level->getPatch(p());
            SAMRAIPointer<CartesianPatchGeometryNd> patch_geo = patch->getPatchGeometry();
            TBOX_ASSERT(patch_geo);
            IBTK::VectorNd x_lo;
            IBTK::VectorNd x_up;
            std::copy(patch_geo->getXLower(), patch_geo->getXLower() + NDIM, x_lo.data());
            std::copy(patch_geo->getXUpper(), patch_geo->getXUpper() + NDIM, x_up.data());
            out << "x lower:\n" << x_lo << '\n';
            out << "x upper:\n" << x_up << '\n';

            BoxNd patch_box = patch->getBox();
            patch_box.grow(gcw);
            if (use_cell)
            {
                SAMRAIPointer<CellDataNd<int> > dof_data = patch->getPatchData(dof_idx);
                dof_data->print(patch_box, out);
            }
            else
            {
                SAMRAIPointer<SideDataNd<int> > dof_data = patch->getPatchData(dof_idx);
                dof_data->print(patch_box, out);
            }
        }
        IBTK_MPI::barrier();

        print_strings_on_plog_0(out.str());
    }
} // run_example
