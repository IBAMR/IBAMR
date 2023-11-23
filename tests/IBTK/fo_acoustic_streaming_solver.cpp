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
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/FOAcousticStreamingPETScLevelSolver.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <sstream>

#include <ibamr/namespaces.h>

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

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "fo_acoustic_streaming_solver.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
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

        Pointer<SideVariable<NDIM, double> > u_sc_var = new SideVariable<NDIM, double>("u_sc", /*depth*/ 2);
        Pointer<CellVariable<NDIM, double> > p_cc_var = new CellVariable<NDIM, double>("p_cc", /*depth*/ 2);
        Pointer<SideVariable<NDIM, double> > fu_sc_var = new SideVariable<NDIM, double>("fu_sc", /*depth*/ 2);
        Pointer<CellVariable<NDIM, double> > fp_cc_var = new CellVariable<NDIM, double>("fp_cc", /*depth*/ 2);
        Pointer<SideVariable<NDIM, double> > eu_sc_var = new SideVariable<NDIM, double>("eu_sc", /*depth*/ 2);
        Pointer<CellVariable<NDIM, double> > ep_cc_var = new CellVariable<NDIM, double>("ep_cc", /*depth*/ 2);

#if (NDIM == 2)
        Pointer<NodeVariable<NDIM, double> > mu_nc_var = new NodeVariable<NDIM, double>("mu_node");
        const int mu_nc_idx = var_db->registerVariableAndContext(mu_nc_var, ctx, IntVector<NDIM>(1));
#elif (NDIM == 3)
        Pointer<EdgeVariable<NDIM, double> > mu_ec_var = new EdgeVariable<NDIM, double>("mu_edge");
        const int mu_ec_idx = var_db->registerVariableAndContext(mu_ec_var, ctx, IntVector<NDIM>(1));
        Pointer<CellVariable<NDIM, double> > mu_cc_var = new CellVariable<NDIM, double>("mu_cc");
        const int mu_cc_idx = var_db->registerVariableAndContext(mu_cc_var, ctx, IntVector<NDIM>(0));
#endif

        Pointer<SideVariable<NDIM, double> > rho_sc_var = new SideVariable<NDIM, double>("rho_sc");
        const int rho_sc_idx = var_db->registerVariableAndContext(rho_sc_var, ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > lambda_cc_var = new CellVariable<NDIM, double>("lambda_cc");
        const int lambda_cc_idx = var_db->registerVariableAndContext(lambda_cc_var, ctx, IntVector<NDIM>(1));

        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(1));
        const int p_cc_idx = var_db->registerVariableAndContext(p_cc_var, ctx, IntVector<NDIM>(1));
        const int fu_sc_idx = var_db->registerVariableAndContext(fu_sc_var, ctx, IntVector<NDIM>(1));
        const int fp_cc_idx = var_db->registerVariableAndContext(fp_cc_var, ctx, IntVector<NDIM>(1));
        const int eu_sc_idx = var_db->registerVariableAndContext(eu_sc_var, ctx, IntVector<NDIM>(1));
        const int ep_cc_idx = var_db->registerVariableAndContext(ep_cc_var, ctx, IntVector<NDIM>(1));

        // Plotting variables
        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc", 2 * NDIM);
        Pointer<CellVariable<NDIM, double> > fu_cc_var = new CellVariable<NDIM, double>("fu_cc", 2 * NDIM);
        Pointer<CellVariable<NDIM, double> > eu_cc_var = new CellVariable<NDIM, double>("eu_cc", 2 * NDIM);

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(0));
        const int fu_cc_idx = var_db->registerVariableAndContext(fu_cc_var, ctx, IntVector<NDIM>(0));
        const int eu_cc_idx = var_db->registerVariableAndContext(eu_cc_var, ctx, IntVector<NDIM>(0));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        // visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "VECTOR", u_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("u_cc_real_" + std::to_string(d), "SCALAR", u_cc_idx, d);
            visit_data_writer->registerPlotQuantity("u_cc_imag_" + std::to_string(d), "SCALAR", u_cc_idx, NDIM + d);
        }

        // visit_data_writer->registerPlotQuantity(p_cc_var->getName(), "VECTOR", p_cc_idx);
        {
            visit_data_writer->registerPlotQuantity("p_cc_real", "SCALAR", p_cc_idx, 0);
            visit_data_writer->registerPlotQuantity("p_cc_imag", "SCALAR", p_cc_idx, 1);
        }

        // visit_data_writer->registerPlotQuantity(fu_cc_var->getName(), "VECTOR", fu_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("fu_cc_real_" + std::to_string(d), "SCALAR", fu_cc_idx, d);
            visit_data_writer->registerPlotQuantity("fu_cc_imag_" + std::to_string(d), "SCALAR", fu_cc_idx, NDIM + d);
        }

        // visit_data_writer->registerPlotQuantity(fp_cc_var->getName(), "VECTOR", fp_cc_idx);
        {
            visit_data_writer->registerPlotQuantity("fp_cc_real", "SCALAR", fp_cc_idx, 0);
            visit_data_writer->registerPlotQuantity("fp_cc_imag", "SCALAR", fp_cc_idx, 1);
        }

        // visit_data_writer->registerPlotQuantity(eu_cc_var->getName(), "VECTOR", eu_cc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("eu_cc_real_" + std::to_string(d), "SCALAR", eu_cc_idx, d);
            visit_data_writer->registerPlotQuantity("eu_cc_imag_" + std::to_string(d), "SCALAR", eu_cc_idx, NDIM + d);
        }

        // visit_data_writer->registerPlotQuantity(ep_cc_var->getName(), "VECTOR", ep_cc_idx);
        {
            visit_data_writer->registerPlotQuantity("ep_cc_real", "SCALAR", ep_cc_idx, 0);
            visit_data_writer->registerPlotQuantity("ep_cc_imag", "SCALAR", ep_cc_idx, 1);
        }

#if (NDIM == 2)
        visit_data_writer->registerPlotQuantity(mu_nc_var->getName(), "SCALAR", mu_nc_idx);
#elif (NDIM == 3)
        visit_data_writer->registerPlotQuantity(mu_cc_var->getName(), "SCALAR", mu_cc_idx);
#endif

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

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(p_cc_idx, 0.0);
            level->allocatePatchData(fu_sc_idx, 0.0);
            level->allocatePatchData(fp_cc_idx, 0.0);
            level->allocatePatchData(eu_sc_idx, 0.0);
            level->allocatePatchData(ep_cc_idx, 0.0);
#if (NDIM == 2)
            level->allocatePatchData(mu_nc_idx, 0.0);
#elif (NDIM == 3)
            level->allocatePatchData(mu_ec_idx, 0.0);
            level->allocatePatchData(mu_cc_idx, 0.0);
#endif
            level->allocatePatchData(rho_sc_idx, 0.0);
            level->allocatePatchData(lambda_cc_idx, 0.0);

            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(fu_cc_idx, 0.0);
            level->allocatePatchData(eu_cc_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> x_vec("up", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f_up", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e_up", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        x_vec.addComponent(u_sc_var, u_sc_idx, h_sc_idx);
        x_vec.addComponent(p_cc_var, p_cc_idx, h_cc_idx);
        f_vec.addComponent(fu_sc_var, fu_sc_idx, h_sc_idx);
        f_vec.addComponent(fp_cc_var, fp_cc_idx, h_cc_idx);
        e_vec.addComponent(eu_sc_var, eu_sc_idx, h_sc_idx);
        e_vec.addComponent(ep_cc_var, ep_cc_idx, h_cc_idx);

        x_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2> u_bc_coefs;
        u_bc_coefs[0].resize(NDIM);
        u_bc_coefs[1].resize(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (int comp = 0; comp < 2; ++comp)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    u_bc_coefs[comp][d] = nullptr;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "ur_bc_coefs_" << d;
                const std::string bc_coefs_name = bc_coefs_name_stream.str();

                std::ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityRealBcCoefs_" << d;
                const std::string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[0][d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                std::ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "ui_bc_coefs_" << d;
                const std::string bc_coefs_name = bc_coefs_name_stream.str();

                std::ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityImagBcCoefs_" << d;
                const std::string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[1][d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }

        std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs(NDIM, nullptr);
        if (!(periodic_shift.min() > 0))
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                rho_bc_coefs[d] = new muParserRobinBcCoefs(
                    "rho_bc_coef",
                    app_initializer->getComponentDatabase("DensityBcCoefs_" + std::to_string(d)),
                    grid_geometry);
            }
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ShearViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ShearViscosityBcCoefs"), grid_geometry);
        }

        RobinBcCoefStrategy<NDIM>* lambda_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("BulkViscosityBcCoefs"))
        {
            lambda_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("BulkViscosityBcCoefs"), grid_geometry);
        }

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction p_fcn("p", app_initializer->getComponentDatabase("p"), grid_geometry);
        muParserCartGridFunction fu_fcn("fu", app_initializer->getComponentDatabase("fu"), grid_geometry);
        muParserCartGridFunction fp_fcn("fp", app_initializer->getComponentDatabase("fp"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);
        muParserCartGridFunction rho_fcn("rho", app_initializer->getComponentDatabase("rho"), grid_geometry);
        muParserCartGridFunction lambda_fcn("lambda", app_initializer->getComponentDatabase("lambda"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(eu_sc_idx, eu_sc_var, patch_hierarchy, 0.0);
        p_fcn.setDataOnPatchHierarchy(ep_cc_idx, ep_cc_var, patch_hierarchy, 0.0);
        fu_fcn.setDataOnPatchHierarchy(fu_sc_idx, fu_sc_var, patch_hierarchy, 0.0);
        fp_fcn.setDataOnPatchHierarchy(fp_cc_idx, fp_cc_var, patch_hierarchy, 0.0);
#if (NDIM == 2)
        mu_fcn.setDataOnPatchHierarchy(mu_nc_idx, mu_nc_var, patch_hierarchy, 0.0);
#elif (NDIM == 3)
        mu_fcn.setDataOnPatchHierarchy(mu_ec_idx, mu_ec_var, patch_hierarchy, 0.0);
#endif
        rho_fcn.setDataOnPatchHierarchy(rho_sc_idx, rho_sc_var, patch_hierarchy, 0.0);
        lambda_fcn.setDataOnPatchHierarchy(lambda_cc_idx, lambda_cc_var, patch_hierarchy, 0.0);

        // Fill ghost cells of viscosity.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comp(3);
#if (NDIM == 2)
        transaction_comp[0] = InterpolationTransactionComponent(mu_nc_idx,
                                                                /*DATA_REFINE_TYPE*/ "LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ false,
                                                                /*DATA_COARSEN_TYPE*/ "CONSTANT_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                mu_bc_coef,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));
#elif (NDIM == 3)
        transaction_comp[0] = InterpolationTransactionComponent(mu_ec_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ false,
                                                                /*DATA_COARSEN_TYPE*/ "CONSERVATIVE_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                mu_bc_coef,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));
#endif

        transaction_comp[1] = InterpolationTransactionComponent(rho_sc_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ false,
                                                                /*DATA_COARSEN_TYPE*/ "CONSERVATIVE_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                rho_bc_coefs,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));

        transaction_comp[2] = InterpolationTransactionComponent(lambda_cc_idx,
                                                                /*DATA_REFINE_TYPE*/ "CONSERVATIVE_LINEAR_REFINE",
                                                                /*USE_CF_INTERPOLATION*/ false,
                                                                /*DATA_COARSEN_TYPE*/ "CONSERVATIVE_COARSEN",
                                                                /*BDRY_EXTRAP_TYPE*/ "LINEAR",
                                                                /*CONSISTENT_TYPE_2_BDRY*/ false,
                                                                lambda_bc_coef,
                                                                Pointer<VariableFillPattern<NDIM> >(NULL));

        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
        hier_bdry_fill->setHomogeneousBc(false);
        hier_bdry_fill->fillData(/*time*/ 0.0);

        // Setup the PETScLevelSolver.
        FOAcousticStreamingPETScLevelSolver petsc_solver(
            "FOAcousticStreamingPETScLevelSolver",
            app_initializer->getComponentDatabase("FOAcousticStreamingPETScLevelSolver"),
            "fo_acoustic_");
        petsc_solver.setAcousticAngularFrequency(input_db->getDouble("acoustic_angular_freq"));
        petsc_solver.setSoundSpeed(input_db->getDouble("sound_speed"));
        petsc_solver.setBoundaryConditionCoefficients(u_bc_coefs);
        petsc_solver.setMassDensityPatchDataIndex(rho_sc_idx);
#if (NDIM == 2)
        petsc_solver.setShearViscosityPatchDataIndex(mu_nc_idx);
#elif (NDIM == 3)
        petsc_solver.setShearViscosityPatchDataIndex(mu_ec_idx);
#endif
        petsc_solver.setBulkViscosityPatchDataIndex(lambda_cc_idx);
        petsc_solver.initializeSolverState(x_vec, f_vec);
        petsc_solver.solveSystem(x_vec, f_vec);

        // Compute error and print error norms.
        // e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
        //                Pointer<SAMRAIVectorReal<NDIM, double> >(&x_vec, false));
        const double e_max_norm = e_vec.maxNorm();
        const double e_l2_norm = e_vec.L2Norm();
        const double e_l1_norm = e_vec.L1Norm();
        pout << "|e|_oo = " << e_max_norm << "\n";
        pout << "|e|_2  = " << e_l2_norm << "\n";
        pout << "|e|_1  = " << e_l1_norm << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cc_idx, u_cc_var, u_sc_idx, u_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(fu_cc_idx, fu_cc_var, fu_sc_idx, fu_sc_var, NULL, 0.0, synch_cf_interface);
        hier_math_ops.interp(eu_cc_idx, eu_cc_var, eu_sc_idx, eu_sc_var, NULL, 0.0, synch_cf_interface);

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > eu_cc_data = patch->getPatchData(eu_cc_idx);
                Pointer<CellData<NDIM, double> > ep_cc_data = patch->getPatchData(ep_cc_idx);
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                    if (!intersection.empty())
                    {
                        eu_cc_data->fillAll(0.0, intersection);
                        ep_cc_data->fillAll(0.0, intersection);
                    }
                }
            }
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        // Cleanup dumb pointers
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[0][d];
            delete u_bc_coefs[1][d];

            delete rho_bc_coefs[d];
        }

        delete mu_bc_coef;
        delete lambda_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example
