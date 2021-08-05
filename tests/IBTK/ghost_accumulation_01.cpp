// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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
#include <ibtk/LEInteractor.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/SAMRAIGhostDataAccumulator.h>
#include <ibtk/muParserCartGridFunction.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <random>

#include <ibtk/app_namespaces.h>

// test stuff
#include "../tests.h"

// test for accumulating ghost data with SAMRAIGhostDataAccumulator.

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

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

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // TODO: make depth, cell vs side, ghost width input arguments
        const int u_depth = input_db->getIntegerWithDefault("depth", 1);
        const IntVector<NDIM> gcw(input_db->getIntegerWithDefault("ghost_width", 1));
        plog << "GCW: " << gcw << '\n';
        const bool use_cell = input_db->getStringWithDefault("var_type", "CELL") == "CELL";

        // disambiguate between libMesh::Variable and SAMRAI::hier::Variable
        auto u_var = use_cell ? Pointer<hier::Variable<NDIM> >(new CellVariable<NDIM, double>("u", u_depth)) :
                                Pointer<hier::Variable<NDIM> >(new SideVariable<NDIM, double>("u", u_depth));

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, gcw);

        // set up grid
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while ((gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_idx, 0.0);
        }

        // Set up values for u
        Pointer<PatchLevel<NDIM> > patch_level = patch_hierarchy->getPatchLevel(finest_ln);
        if (input_db->getStringWithDefault("fill_test", "uniform") == "uniform")
        {
            for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
                if (use_cell)
                {
                    Pointer<CellData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                    u_data->fillAll(1.0);
                }
                else
                {
                    Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                    u_data->fillAll(1.0);
                }
            }
        }
        else if (input_db->getStringWithDefault("fill_test", "uniform") == "spread")
        {
            // output values and their coordinates:
            const std::size_t n_points = 100;
            // Note that side-centered data already has NDIM components: we
            // need that many for Q
            const int Q_depth = use_cell ? 1 : NDIM;
            std::vector<double> Q_data(Q_depth * n_points);
            const int X_depth = NDIM;
            std::vector<double> X_data(X_depth * n_points);

            // populate coordinates randomly:
            std::mt19937 std_seq(42u);
            // TODO: this hard-codes the domain size
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            for (double& v : Q_data) v = distribution(std_seq);
            for (double& v : X_data) v = distribution(std_seq);

            // spread *only points in patch interiors*:
            const std::string spread_fcn = input_db->getString("IB_DELTA_FUNCTION");
            for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
                const Box<NDIM> box = patch->getBox();
                if (use_cell)
                {
                    Pointer<CellData<NDIM, double> > q_data = patch->getPatchData(u_idx);
                    q_data->fillAll(0.0);
                    LEInteractor::spread(q_data, Q_data, Q_depth, X_data, X_depth, patch, box, spread_fcn);
                }
                else
                {
                    Pointer<SideData<NDIM, double> > q_data = patch->getPatchData(u_idx);
                    q_data->fillAll(0.0);
                    LEInteractor::spread(q_data, Q_data, Q_depth, X_data, X_depth, patch, box, spread_fcn);
                }
            }
        }

        // accumulate values
        SAMRAIGhostDataAccumulator accumulator(patch_hierarchy, u_var, gcw, finest_ln, finest_ln);
        accumulator.accumulateGhostData(u_idx);

        std::ostringstream out;
        out.precision(16);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            tbox::Pointer<CartesianPatchGeometry<NDIM> > patch_geo = patch->getPatchGeometry();
            TBOX_ASSERT(patch_geo);
            IBTK::VectorNd x_lo;
            IBTK::VectorNd x_up;
            std::copy(patch_geo->getXLower(), patch_geo->getXLower() + NDIM, x_lo.data());
            std::copy(patch_geo->getXUpper(), patch_geo->getXUpper() + NDIM, x_up.data());
            out << "\nRank: " << IBTK_MPI::getRank() << '\n';
            out << "x lower:\n" << x_lo << '\n';
            out << "x upper:\n" << x_up << '\n';

            Box<NDIM> patch_box = patch->getBox();
            if (use_cell)
            {
                Pointer<CellData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                u_data->print(patch_box, out, 16);
            }
            else
            {
                Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
                u_data->print(patch_box, out, 16);
            }
        }
        IBTK_MPI::barrier();

        print_strings_on_plog_0(out.str());
    }
} // run_example
