// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
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

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HyperbolicLevelIntegrator.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvectorExplicitPredictorPatchOps.h>
#include <ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h>
#include <ibamr/RelaxationLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <LocationIndexRobinBcCoefs.h>
#include <TimeRefinementIntegrator.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

struct CircularInterface
{
    IBTK::Vector X0;
    double R;
};
void
circular_interface_neighborhood(int D_idx,
                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                double /*time*/,
                                bool /*initial_time*/,
                                void* /*ctx*/)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(D_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                const double* patch_X_lower = patch_geom->getXLower();
                const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                const double* const patch_dx = patch_geom->getDx();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                const double x = coord[0];
                const double y = coord[1];
                (*D_data)(ci) = 1.0 * (std::pow(x - 1.0, 2.0) + std::pow(y - 1.0, 2.0) + 0.1) *
                                (std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - 1.0);
            }
        }
    }
    return;
} // circular_interface_neighborhood

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "advect.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        // Get solver configuration options.
        bool using_refined_timestepping = false;
        if (main_db->keyExists("timestepping"))
        {
            string timestepping_method = main_db->getString("timestepping");
            if (timestepping_method == "SYNCHRONIZED")
            {
                using_refined_timestepping = false;
            }
            else
            {
                using_refined_timestepping = true;
            }
        }
        if (using_refined_timestepping)
        {
            pout << "using subcycled timestepping.\n";
        }
        else
        {
            pout << "NOT using subcycled timestepping.\n";
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor = new AdvectorExplicitPredictorPatchOps(
            "AdvectorExplicitPredictorPatchOps",
            app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<AdvectorPredictorCorrectorHyperbolicPatchOps> hyp_patch_ops =
            new AdvectorPredictorCorrectorHyperbolicPatchOps(
                "AdvectorPredictorCorrectorHyperbolicPatchOps",
                app_initializer->getComponentDatabase("AdvectorPredictorCorrectorHyperbolicPatchOps"),
                explicit_predictor,
                grid_geometry);
        Pointer<HyperbolicLevelIntegrator<NDIM> > hyp_level_integrator =
            new HyperbolicLevelIntegrator<NDIM>("HyperbolicLevelIntegrator",
                                                app_initializer->getComponentDatabase("HyperbolicLevelIntegrator"),
                                                hyp_patch_ops,
                                                true,
                                                using_refined_timestepping);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               hyp_level_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        Pointer<TimeRefinementIntegrator<NDIM> > time_integrator =
            new TimeRefinementIntegrator<NDIM>("TimeRefinementIntegrator",
                                               app_initializer->getComponentDatabase("TimeRefinementIntegrator"),
                                               patch_hierarchy,
                                               hyp_level_integrator,
                                               gridding_algorithm);

        // Setup the advection velocity.
        const bool u_is_div_free = main_db->getBoolWithDefault("u_is_div_free", false);
        if (u_is_div_free)
        {
            pout << "advection velocity u is discretely divergence free.\n";
        }
        else
        {
            pout << "advection velocity u is NOT discretely divergence free.\n";
        }
        Pointer<FaceVariable<NDIM, double> > u_var = new FaceVariable<NDIM, double>("u");
        Pointer<CartGridFunction> u_fcn = new muParserCartGridFunction(
            "u_fcn", app_initializer->getComponentDatabase("AdvectionVelocityFunction"), grid_geometry);
        hyp_patch_ops->registerAdvectionVelocity(u_var);
        hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        hyp_patch_ops->setAdvectionVelocityFunction(u_var, u_fcn);

        // Setup the advected quantity.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        pout << "solving the advection equation in "
             << IBAMR::enum_to_string<ConvectiveDifferencingType>(difference_form) << " form.\n";
        Pointer<CellVariable<NDIM, double> > Q_var = new CellVariable<NDIM, double>("Q");
        LocationIndexRobinBcCoefs<NDIM> physical_bc_coef(
            "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        hyp_patch_ops->registerTransportedQuantity(Q_var);
        hyp_patch_ops->setAdvectionVelocity(Q_var, u_var);
        hyp_patch_ops->setConvectiveDifferencingType(Q_var, difference_form);
        hyp_patch_ops->setPhysicalBcCoefs(Q_var, &physical_bc_coef);

        // Level set initial conditions
        Pointer<CartGridFunction> Q_init = new muParserCartGridFunction(
            "Q_init", app_initializer->getComponentDatabase("QInitFunction"), grid_geometry);
        hyp_patch_ops->setInitialConditions(Q_var, Q_init);

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit) hyp_patch_ops->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializeHierarchy();

        // Create inital level set
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        input_db->getDoubleArray("X0", circle.X0.data(), NDIM);

        Pointer<VariableContext> current_ctx = hyp_level_integrator->getCurrentContext();
        Pointer<VariableContext> scratch_ctx = hyp_level_integrator->getScratchContext();
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, current_ctx);
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, scratch_ctx);

        Pointer<CellVariable<NDIM, double> > E_var = new CellVariable<NDIM, double>("E");
        const int E_idx = var_db->registerVariableAndContext(E_var, scratch_ctx);

        // Heaviside
        Pointer<CellVariable<NDIM, double> > H_var = new CellVariable<NDIM, double>("H");
        const int H_idx = var_db->registerVariableAndContext(H_var, scratch_ctx);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(Q_scratch_idx))
                level->allocatePatchData(Q_scratch_idx, time_integrator->getIntegratorTime());
            if (!level->checkAllocated(E_idx)) level->allocatePatchData(E_idx, time_integrator->getIntegratorTime());
            if (!level->checkAllocated(H_idx)) level->allocatePatchData(H_idx, time_integrator->getIntegratorTime());
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        Pointer<HierarchyMathOps> hier_math_ops =
            new HierarchyMathOps("HierarchyMathOps", patch_hierarchy, coarsest_ln, finest_ln);
        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("LevelSet"));
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&circular_interface_neighborhood, (void*)&circle);
        level_set_ops->registerPhysicalBoundaryCondition(&physical_bc_coef);
        level_set_ops->initializeLSData(Q_scratch_idx,
                                        hier_math_ops,
                                        time_integrator->getIntegratorStep(),
                                        time_integrator->getIntegratorTime(),
                                        /*initial_time*/ true);

        // Compute L1 error from analytical solution
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > E_data = patch->getPatchData(E_idx);
                Pointer<CellData<NDIM, double> > Q_scratch_data = patch->getPatchData(Q_scratch_idx);
                Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(H_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());

                    // Get physical coordinates
                    IBTK::Vector coord = IBTK::Vector::Zero();
                    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    const double* patch_X_lower = patch_geom->getXLower();
                    const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
                    const double* const patch_dx = patch_geom->getDx();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord[d] =
                            patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                    }
                    const double distance =
                        std::sqrt(std::pow((coord[0] - circle.X0(0)), 2.0) + std::pow((coord[1] - circle.X0(1)), 2.0));

                    (*E_data)(ci) = distance - circle.R;

                    const double phi =
                        -(*Q_scratch_data)(ci); // This make sure phi is positive inide the circle so that H = 1.
                    (*H_data)(ci) = IBTK::discontinuous_heaviside(phi);
                }
            }
        }

        HierarchyCellDataOpsReal<NDIM, double> cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        cc_data_ops.subtract(E_idx, E_idx, Q_scratch_idx);
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
        const double EQ_domain = cc_data_ops.L1Norm(E_idx, wgt_cc_idx);

        const double Numerical_volume = cc_data_ops.integral(H_idx, wgt_cc_idx);
        const double Exact_volume = M_PI * std::pow(circle.R, 2.0);
        const double vol_error = std::abs(Numerical_volume - Exact_volume) / Exact_volume;

        double E_domain = 0.0;
        double E_interface = 0.0;
        int num_interface_pts = 0;
        // Compute L1 Norm for specific regions
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(Q_scratch_idx);
                Pointer<CellData<NDIM, double> > E_data = patch->getPatchData(E_idx);
                Pointer<CellData<NDIM, double> > W_data = patch->getPatchData(wgt_cc_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi = (*D_data)(ci);
                    const double err = (*E_data)(ci);
                    const double dV = (*W_data)(ci);
                    Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    const double* const patch_dx = patch_geom->getDx();

                    if (std::abs(phi) < 1.2 * patch_dx[0])
                    {
                        E_interface += std::abs(err) * dV;
                        num_interface_pts++;
                    }
                    if (phi > -0.8) E_domain += std::abs(err) * dV;
                }
            }
        }
        // Perform sum reduction
        num_interface_pts = IBTK_MPI::sumReduction(num_interface_pts);
        E_interface = IBTK_MPI::sumReduction(E_interface);
        E_domain = IBTK_MPI::sumReduction(E_domain);

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");

            out << "Error in Q near interface after level set re-initialization:" << std::endl
                << "L1-norm:  " << std::setprecision(10) << E_interface << std::endl;
            out << "Number of points within the interface (used to compute interface error):" << std::endl
                << num_interface_pts << std::endl;
            out << "Error in Q in entire domain (minus center) after level set re-initialization:" << std::endl
                << "L1-norm:  " << std::setprecision(10) << E_domain << std::endl;
            out << "Error in Q in entire domain (including center) after level set re-initialization:" << std::endl
                << "L1-norm:  " << std::setprecision(10) << EQ_domain << std::endl;
            out << "Volume error for the re-initialized circle using the discontinuous Heaviside function: "
                << std::setprecision(10) << vol_error << std::endl;
        }

        // Register for plotting
        visit_data_writer->registerPlotQuantity("Error", "SCALAR", E_idx);

        cc_data_ops.copyData(Q_current_idx, Q_scratch_idx);
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (level->checkAllocated(Q_scratch_idx)) level->deallocatePatchData(Q_scratch_idx);
        }

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        if (dump_viz_data && uses_visit)
        {
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

    } // cleanup dynamically allocated objects prior to shutdown
} // main
