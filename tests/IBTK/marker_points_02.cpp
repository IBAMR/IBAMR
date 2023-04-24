// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
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
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/MarkerPatchHierarchy.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/RestartManager.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <fstream>

#include <ibtk/app_namespaces.h>

#include "../tests.h"

void
test(const MarkerPatch& marker_patch, std::ostream& output)
{
    for (unsigned int i = 0; i < marker_patch.size(); ++i)
    {
        const auto marker = marker_patch[i];
        output << "  index = " << std::get<0>(marker) << " X = ";
        for (unsigned int d = 0; d < NDIM - 1; ++d)
        {
            output << std::get<1>(marker)[d] << ", ";
        }
        output << std::get<1>(marker)[NDIM - 1] << " V = ";
        for (unsigned int d = 0; d < NDIM - 1; ++d)
        {
            output << std::get<2>(marker)[d] << ", ";
        }
        output << std::get<2>(marker)[NDIM - 1] << " contains = " << marker_patch.contains(std::get<1>(marker)) << '\n';
    }
}

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    // Skip warnings:
    Logger::getInstance()->setWarning(false);

    const auto rank = IBTK_MPI::getRank();
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db = input_db->keyExists("test") ? input_db->getDatabase("test") : nullptr;
    const bool set_velocity = test_db && test_db->getBoolWithDefault("set_velocity", false);
    const bool timestep = test_db && test_db->getBoolWithDefault("timestep", false);
    const bool collective_print = test_db && test_db->getBoolWithDefault("collective_print", false);
    const bool from_restart = RestartManager::getManager()->isFromRestart();

    int u_idx = IBTK::invalid_index;
    Pointer<hier::Variable<NDIM> > u_var;
    if (set_velocity)
    {
        // Set up the velocity on the Cartesian grid:
        u_var = new pdat::SideVariable<NDIM, double>("u_sc", 1);
        const std::string kernel = test_db->getStringWithDefault("kernel", "PIECEWISE_LINEAR");
        const int ghost_width = LEInteractor::getMinimumGhostWidth(kernel);

        auto* var_db = hier::VariableDatabase<NDIM>::getDatabase();
        tbox::Pointer<hier::VariableContext> ctx = var_db->getContext("context");
        u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(ghost_width));
    }

    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    grid_geometry->addSpatialRefineOperator(new CartesianCellDoubleLinearRefine<NDIM>());
    grid_geometry->addSpatialRefineOperator(new CartSideDoubleSpecializedLinearRefine());

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

    // Set up marker points.
    EigenAlignedVector<IBTK::Point> positions;
    EigenAlignedVector<IBTK::Vector> velocities;
    if (!from_restart)
    {
        for (unsigned int i = 0; i < 10; ++i)
        {
            for (unsigned int j = 0; j < 10; ++j)
            {
                positions.emplace_back(double(i) / 10.0, double(j) / 10.0);
                if (set_velocity)
                    velocities.emplace_back(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
                else
                    velocities.emplace_back(double(i) * 10 + j, 0.0);
            }
        }
    }
    MarkerPatchHierarchy marker_points("MarkerPoints", patch_hierarchy, positions, velocities);

    if (set_velocity)
    {
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_idx, 0.0);
        }
        HierarchySideDataOpsReal<NDIM, double> ops(patch_hierarchy);
        ops.setToScalar(u_idx, std::numeric_limits<double>::quiet_NaN(), false);
        IBTK::muParserCartGridFunction u_fcn("u", test_db->getDatabase("u"), patch_hierarchy->getGridGeometry());
        u_fcn.setDataOnPatchHierarchy(u_idx, u_var, patch_hierarchy, 0.0);

        std::vector<std::unique_ptr<muParserRobinBcCoefs> > u_bc_coefs;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coef_ptrs;

        if (test_db->keyExists("UBcCoefs_0"))
        {
            for (int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs.emplace_back(new muParserRobinBcCoefs("u_bc_coefs_" + std::to_string(d),
                                                                 test_db->getDatabase("UBcCoefs_" + std::to_string(d)),
                                                                 grid_geometry));
                u_bc_coef_ptrs.push_back(u_bc_coefs.back().get());
            }
        }

        // sync ghost data:
        using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_components(1);
        ghost_cell_components[0] = ITC(
            u_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", true, u_bc_coef_ptrs, nullptr);
        IBTK::HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        // interpolate:
        const std::string kernel = test_db->getStringWithDefault("kernel", "PIECEWISE_LINEAR");
        marker_points.setVelocities(u_idx, kernel);
    }

    if (timestep)
    {
        // Pick CFL = 1.0 to simulate the largest marker point travel possible:
        Pointer<PatchLevel<NDIM> > finest_level =
            patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
        const double dx = get_min_patch_dx(*finest_level);
        const double U_MAX = 1.0;
        const double dt = dx / U_MAX;

        const std::string kernel = test_db->getStringWithDefault("kernel", "PIECEWISE_LINEAR");
        const int num_timesteps = test_db->getIntegerWithDefault("num_timesteps", 10);
        const int starting_step = from_restart ? app_initializer->getRestartRestoreNumber() + 1 : 0;
        for (int i = starting_step; i < num_timesteps; ++i)
        {
            marker_points.midpointStep(dt, u_idx, u_idx, kernel);

            if (app_initializer->dumpRestartData() && (i % app_initializer->getRestartDumpInterval() == 0))
            {
                RestartManager::getManager()->writeRestartFile(app_initializer->getRestartDumpDirectory(), i);
            }
        }
    }

    if (collective_print)
    {
        auto all_points = marker_points.collectAllMarkers();
        for (unsigned int k = 0; k < marker_points.getNumberOfMarkers(); ++k)
        {
            plog << "X = ";
            for (unsigned int d = 0; d < NDIM - 1; ++d)
            {
                plog << all_points.first[k][d] << ", ";
            }
            plog << all_points.first[k][NDIM - 1] << " V = ";
            for (unsigned int d = 0; d < NDIM - 1; ++d)
            {
                plog << all_points.second[k][d] << ", ";
            }
            plog << all_points.second[k][NDIM - 1] << '\n';
        }
    }
    else
    {
        std::ostringstream out;
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            int local_patch_num = 0;
            Pointer<PatchLevel<NDIM> > current_level = patch_hierarchy->getPatchLevel(ln);
            for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
            {
                if (rank == current_level->getMappingForPatch(p))
                {
                    out << "level = " << ln << " patch = " << current_level->getPatch(p)->getBox() << std::endl;
                    test(marker_points.getMarkerPatch(ln, local_patch_num), out);

                    ++local_patch_num;
                }
            }
        }

        print_strings_on_plog_0(out.str());
    }
}
