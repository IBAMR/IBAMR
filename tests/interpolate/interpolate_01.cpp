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

// Config files

#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAI_config.h>

// Headers for basic PETSc objects
#include <petscsys.h>

// SAMRAI INCLUDES
#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianCellDoubleLinearRefine.h>
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
#include <SAMRAIRefineOperator.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAIVisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <random>

#include <ibtk/app_namespaces.h>

// test stuff
#include "../tests.h"

// This file is the main driver for interpolation tests. This verifies that
// the stencils produce consistent output; i.e., we test the stencils and not
// every overload of LEInteractor::interpolate().
//
// This test can be run in two different ways: when use_exact = TRUE is
// specified in the input file then we use an exact multilinear solution
// (hardcoded below) and assert that interpolated values match this exactly
// (all kernels except piecewise constant and discontinuous linear can
// interpolate multilinear solutions exactly). If this value is FALSE (the
// default) then we interpolate a trigonometric field and print the values to
// output.

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

        // we don't want to use a conservative refinement scheme
        SAMRAIPointer<SAMRAIRefineOperator> linear_refine = new SAMRAICartesianCellDoubleLinearRefine();
        grid_geometry->addSpatialRefineOperator(linear_refine);

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

        const int n_ghosts = LEInteractor::getMinimumGhostWidth(input_db->getString("IB_DELTA_FUNCTION"));
        SAMRAIPointer<SAMRAICellVariable<double>> u_cc_var = new SAMRAICellVariable<double>("u_cc", NDIM);
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, SAMRAIIntVector(n_ghosts));

        // Register variables for plotting.
        SAMRAIPointer<SAMRAIVisItDataWriter> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_cc_var->getName() + std::to_string(d), "SCALAR", u_cc_idx, d);
        }

        // Initialize the AMR patch hierarchy.
        const int tag_buffer = 1;
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
        }

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        u_fcn.setDataOnPatchHierarchy(u_cc_idx, u_cc_var, patch_hierarchy, 0.0);

        // Fill ghost data.
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_components(1);
        ghost_cell_components[0] = ITC(u_cc_idx,
                                       "LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // u_bc_coefs
                                       nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        // Here comes the actual test:
        {
            // output values and their coordinates:
            const std::size_t n_points = 100;
            const int Q_depth = NDIM;
            std::vector<double> Q_data(Q_depth * n_points, std::numeric_limits<double>::max());
            const int X_depth = NDIM;
            std::vector<double> X_data(X_depth * n_points);

            // required SAMRAI data:
            //
            // avoid problems with filling boundary ghost data (which isn't
            // relevant to this test) by picking the sole patch on level 1:
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(1);
            const SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(0);
            SAMRAIPointer<SAMRAICellData<double>> q_data = patch->getPatchData(u_cc_idx);
            const SAMRAIBox& interp_box = patch->getBox();

            // populate coordinates randomly:
            std::mt19937 std_seq(42u);
            const SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double width = patch_x_upper[0] - patch_x_lower[0];
            std::cout << "patch bounds:\n";
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::cout << patch_x_lower[d] << ", " << patch_x_upper[d] << '\n';
                TBOX_ASSERT(width == patch_x_upper[d] - patch_x_lower[d]);
            }
            const bool use_exact = input_db->getBoolWithDefault("use_exact", false);
            double lower = patch_x_lower[0];
            double upper = patch_x_upper[0];
            std::uniform_real_distribution<double> distribution(lower, upper);
            for (double& v : X_data) v = distribution(std_seq);

            // interpolate:
            const std::string interp_fcn = input_db->getString("IB_DELTA_FUNCTION");
            LEInteractor::interpolate(Q_data, Q_depth, X_data, X_depth, q_data, patch, interp_box, interp_fcn);

            // output:
            std::ofstream out("output");
            if (!use_exact)
            {
                out.precision(16);
                out << "x, y, z, u0, u1, u2, e0, e1, e2\n";
            }
            for (std::size_t point_n = 0; point_n < n_points; ++point_n)
            {
                if (use_exact)
                {
                    double exact = 0;
                    const double x = X_data[point_n * NDIM + 0];
                    const double y = X_data[point_n * NDIM + 1];
#if NDIM == 2
                    exact = 1 + 2 * x + 3 * y + 4 * x * y;
#else
                    {
                        const double z = X_data[point_n * NDIM + 2];
                        exact = 1 + 2 * x + 3 * y - z + 4 * x * y + 2 * x * z + 3 * x * y * z;
                    }
#endif
                    for (int d = 0; d < NDIM; ++d)
                    {
                        // the ultra-wide kernels have a lot of trouble with
                        // roundoff that is evident at different optimization
                        // settings
                        const double tol = LEInteractor::getStencilSize(interp_fcn) > 4 ? 1e-10 : 1e-12;
                        const double error = std::abs(Q_data[point_n * NDIM + d] - exact);
                        if (error > tol)
                        {
                            if (NDIM == 2)
                            {
                                std::cout << x << ", " << y << ", " << error << '\n';
                                out << x << ", " << y << ", " << error << '\n';
                            }
                            else
                            {
                                const double z = X_data[point_n * NDIM + 2];
                                std::cout << x << ", " << y << ", " << z << ", " << error << '\n';
                                out << x << ", " << y << ", " << z << ", " << error << '\n';
                            }
                        }
                    }
                }
                else
                {
                    out << X_data[point_n * NDIM + 0] << ", " << X_data[point_n * NDIM + 1] << ", "
                        << (NDIM == 3 ? X_data[point_n * NDIM + 2] : 0) << ", " << Q_data[point_n * NDIM + 0] << ", "
                        << Q_data[point_n * NDIM + 1] << ", " << (NDIM == 3 ? Q_data[point_n * NDIM + 2] : 0) << '\n';
                }
            }
            if (use_exact) out << "OK\n";
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
    } // cleanup dynamically allocated objects prior to shutdown
} // main
