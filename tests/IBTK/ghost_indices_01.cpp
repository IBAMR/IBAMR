// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/samrai_compatibility_names.h>
// SAMRAI INCLUDES
#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscsys.h>

#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIIntVector.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideData.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariable.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAI_config.h>

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

        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        // TODO: make depth, cell vs side, ghost width input arguments
        const int u_depth = input_db->getIntegerWithDefault("depth", 1);
        const SAMRAIIntVector gcw(input_db->getIntegerWithDefault("ghost_width", 1));
        const bool use_cell = input_db->getStringWithDefault("var_type", "CELL") == "CELL";
        // disambiguate between libMesh::Variable and SAMRAI::hier::Variable
        auto u_var = use_cell ? SAMRAIPointer<SAMRAIVariable>(new SAMRAICellVariable<double>("u", u_depth)) :
                                SAMRAIPointer<SAMRAIVariable>(new SAMRAISideVariable<double>("u", u_depth));
        auto dof_var = use_cell ? SAMRAIPointer<SAMRAIVariable>(new SAMRAICellVariable<int>("dof", u_depth)) :
                                  SAMRAIPointer<SAMRAIVariable>(new SAMRAISideVariable<int>("dof", u_depth));

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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
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
            SAMRAIPointer<SAMRAIPatchLevel> patch_level =
                patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
            const BoxArray<NDIM> boxes = patch_level->getBoxes();
            plog << "hierarchy boxes:\n";
            for (int i = 0; i < boxes.size(); ++i) plog << boxes[i] << '\n';
            // rank is only relevant when there are multiple processors
            out << "\nrank: " << IBTK_MPI::getRank() << '\n';
        }

        SAMRAIPointer<SAMRAIPatchLevel> patch_level = patch_hierarchy->getPatchLevel(level);
        for (SAMRAIPatchLevel::Iterator p(patch_level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = patch_level->getPatch(p());
            SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geo = patch->getPatchGeometry();
            TBOX_ASSERT(patch_geo);
            IBTK::VectorNd x_lo;
            IBTK::VectorNd x_up;
            std::copy(patch_geo->getXLower(), patch_geo->getXLower() + NDIM, x_lo.data());
            std::copy(patch_geo->getXUpper(), patch_geo->getXUpper() + NDIM, x_up.data());
            out << "x lower:\n" << x_lo << '\n';
            out << "x upper:\n" << x_up << '\n';

            SAMRAIBox patch_box = patch->getBox();
            patch_box.grow(gcw);
            if (use_cell)
            {
                SAMRAIPointer<SAMRAICellData<int>> dof_data = patch->getPatchData(dof_idx);
                dof_data->print(patch_box, out);
            }
            else
            {
                SAMRAIPointer<SAMRAISideData<int>> dof_data = patch->getPatchData(dof_idx);
                dof_data->print(patch_box, out);
            }
        }
        IBTK_MPI::barrier();

        print_strings_on_plog_0(out.str());
    }
} // run_example
