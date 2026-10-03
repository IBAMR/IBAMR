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

// Config files

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredDivergenceFreePhysBdryOp.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredPressureBcCoef.h>
#include <ibamr/INSStaggeredVelocityBcCoef.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/RestartCleaner.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <muParser.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// A velocity given by muParser expressions function_0, function_1, ... of the
// position X_0, X_1, ....
class ParsedVelocity
{
public:
    explicit ParsedVelocity(Pointer<Database> db) : d_parsers(NDIM)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            d_parsers[d].SetExpr(db->getString("function_" + std::to_string(d)));
            for (unsigned int k = 0; k < NDIM; ++k)
            {
                d_parsers[d].DefineVar("X_" + std::to_string(k), &d_X[k]);
            }
        }
    }

    ParsedVelocity(const ParsedVelocity&) = delete;
    ParsedVelocity& operator=(const ParsedVelocity&) = delete;

    // Evaluate component axis at the position x.
    double operator()(const unsigned int axis, const std::array<double, NDIM>& x)
    {
        d_X = x;
        return d_parsers[axis].Eval();
    }

private:
    std::array<double, NDIM> d_X;
    std::vector<mu::Parser> d_parsers;
};

// Function prototypes
void check_open_boundary_ghost_values(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                      Pointer<INSHierarchyIntegrator> ins_integrator,
                                      const vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);
void check_traction_corner_accuracy(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                    Pointer<INSHierarchyIntegrator> ins_integrator,
                                    Pointer<Database> exact_velocity_db);

void check_divergence_free_extension(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                     Pointer<INSStaggeredHierarchyIntegrator> ins_integrator,
                                     Pointer<Database> exact_velocity_db,
                                     int extension_ghost_width);

void output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> ins_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

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
        // prevent a warning about timer initializations
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        // Initialize RestartCleaner if configured
        Pointer<RestartCleaner> restart_cleaner;
        if (dump_restart_data && input_db->isDatabase("RestartCleaner"))
        {
            restart_cleaner = new RestartCleaner("RestartCleaner", input_db->getDatabase("RestartCleaner"));
        }

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> time_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const bool check_ghost_values = main_db->keyExists("check_open_boundary_ghost_values") &&
                                        main_db->getBool("check_open_boundary_ghost_values");
        if (solver_type == "STAGGERED")
        {
            time_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            time_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerPressureInitialConditions(p_init);

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create body force function specification objects (when necessary).
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        const bool check_accuracy = input_db->keyExists("check_traction_corner_accuracy") &&
                                    input_db->getBool("check_traction_corner_accuracy");
        const bool check_divergence_free = input_db->keyExists("check_divergence_free_extension") &&
                                           input_db->getBool("check_divergence_free_extension");
        if (check_accuracy || check_divergence_free)
        {
            // preprocessIntegrateHierarchy() sets up the velocity boundary condition objects.
            const double current_time = time_integrator->getIntegratorTime();
            time_integrator->preprocessIntegrateHierarchy(
                current_time, current_time + time_integrator->getMaximumTimeStepSize(), 1);
            if (check_accuracy)
            {
                check_traction_corner_accuracy(
                    patch_hierarchy, time_integrator, app_initializer->getComponentDatabase("ExactVelocity"));
            }
            if (check_divergence_free)
            {
                Pointer<INSStaggeredHierarchyIntegrator> staggered_integrator = time_integrator;
                if (!staggered_integrator)
                {
                    TBOX_ERROR("check_divergence_free_extension requires the STAGGERED solver.\n");
                }
                check_divergence_free_extension(patch_hierarchy,
                                                staggered_integrator,
                                                app_initializer->getComponentDatabase("ExactVelocity"),
                                                input_db->getInteger("extension_ghost_width"));
            }
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                delete u_bc_coefs[d];
            }
            return 0;
        }

        // Deallocate initialization objects.
        app_initializer.setNull();

        // The check does not write the input database to the log file.
        if (check_ghost_values)
        {
            check_open_boundary_ghost_values(patch_hierarchy, time_integrator, u_bc_coefs);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                delete u_bc_coefs[d];
            }
            return 0;
        }

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                if (restart_cleaner)
                {
                    restart_cleaner->cleanup();
                }
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy, time_integrator, iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }

        // Report RestartCleaner status for verification
        if (restart_cleaner)
        {
            auto remaining = restart_cleaner->getAvailableRestartRestoreNumbers();
            pout << "\n"
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                 << "RestartCleaner: " << remaining.size() << " restart directories remain\n"
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<Variable<NDIM>> u_var = time_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();

        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        const Pointer<Variable<NDIM>> p_var = time_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();

        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx, loop_time);
        }

        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        Pointer<CellVariable<NDIM, double>> u_cc_var = u_var;
        if (u_cc_var)
        {
            HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            hier_cc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);

            pout << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(u_idx, wgt_cc_idx) << "\n"
                 << "  L2-norm:  " << hier_cc_data_ops.L2Norm(u_idx, wgt_cc_idx) << "\n"
                 << "  max-norm: " << hier_cc_data_ops.maxNorm(u_idx, wgt_cc_idx) << "\n";
        }

        Pointer<SideVariable<NDIM, double>> u_sc_var = u_var;
        if (u_sc_var)
        {
            HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
            pout << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n";
        }

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        pout << "Error in p at time " << loop_time - 0.5 * dt << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data && uses_visit)
        {
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main

// Apply the staggered Stokes operator to non-smooth data in two ways and compare the results.  The first is the
// library's application, which sets the normal velocity ghost values that impose TRACTION conditions and the pressure
// boundary value -g where the normal velocity is not prescribed.  The second is assembled here with divergence-free
// normal velocity ghost values and the pressure ghost values p_G = 2*p_b - p_I with p_b = 2*mu*du_n/dx_n - g.
void
check_open_boundary_ghost_values(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                 Pointer<INSHierarchyIntegrator> ins_integrator,
                                 const vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    const double fill_time = 0.5;
    const double mu = ins_integrator->getStokesSpecifications()->getMu();

    // Configure the integrator's boundary condition objects.
    const vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = ins_integrator->getVelocityBoundaryConditions();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto U_bc_coef = dynamic_cast<INSStaggeredVelocityBcCoef*>(U_bc_coefs[d]);
        U_bc_coef->setStokesSpecifications(ins_integrator->getStokesSpecifications());
        U_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
        U_bc_coef->setSolutionTime(fill_time);
    }
    RobinBcCoefStrategy<NDIM>* P_bc_coef = ins_integrator->getPressureBoundaryConditions();
    auto P_ins_bc_coef = dynamic_cast<INSStaggeredPressureBcCoef*>(P_bc_coef);
    P_ins_bc_coef->setPhysicalBcCoefs(u_bc_coefs);
    P_ins_bc_coef->setSolutionTime(fill_time);

    // Data: u and p are the operator's argument, u2 is a copy of u with divergence-free ghost values, and y1 and y2 are
    // the two results.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("open_boundary_ghost_values");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    auto register_data = [&](Pointer<Variable<NDIM>> var, const int ghosts)
    {
        const int idx = var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(ghosts));
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(idx, fill_time);
        }
        return idx;
    };
    const int u_idx = register_data(u_var, 1), p_idx = register_data(p_var, 1);
    const int u2_idx = register_data(new SideVariable<NDIM, double>("u2"), 1);
    Pointer<SideVariable<NDIM, double>> y1u_var = new SideVariable<NDIM, double>("y1u");
    Pointer<CellVariable<NDIM, double>> y1p_var = new CellVariable<NDIM, double>("y1p");
    const int y1u_idx = register_data(y1u_var, 0), y1p_idx = register_data(y1p_var, 0);
    const int y2u_idx = register_data(new SideVariable<NDIM, double>("y2u"), 0);
    const int y2p_idx = register_data(new CellVariable<NDIM, double>("y2p"), 0);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            u_data->fillAll(0.0);
            p_data->fillAll(0.0);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                {
                    (*u_data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                        std::sin(1.7 * it()(0) + 2.9 * it()(1) + 0.8 * axis);
                }
            }
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                (*p_data)(it()) = std::cos(2.3 * it()(0) + 0.9 * it()(1));
            }
        }
    }

    // Weights, vectors, and the operator.
    HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
    SAMRAIVectorReal<NDIM, double> x("x", patch_hierarchy, 0, finest_ln), y1("y1", patch_hierarchy, 0, finest_ln);
    x.addComponent(u_var, u_idx, wgt_sc_idx);
    x.addComponent(p_var, p_idx, wgt_cc_idx);
    y1.addComponent(y1u_var, y1u_idx, wgt_sc_idx);
    y1.addComponent(y1p_var, y1p_idx, wgt_cc_idx);
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    bc_helper->cacheBcCoefData(u_bc_coefs, fill_time, patch_hierarchy);
    PoissonSpecifications spec("spec");
    spec.setCConstant(1.0);
    spec.setDConstant(-mu);
    StaggeredStokesOperator op("StaggeredStokesOperator", /*homogeneous_bc*/ false);
    op.setVelocityPoissonSpecifications(spec);
    op.setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
    op.setPhysicalBoundaryHelper(bc_helper);
    op.setSolutionTime(fill_time);
    op.setTimeInterval(fill_time, fill_time);
    op.initializeOperatorState(x, y1);

    // (i) The library's application; it sets the ghost values of x.
    op.apply(x, y1);

    // (ii) Divergence-free normal velocity ghost values and the corresponding pressure ghost values (which replace
    // those of x) where the normal velocity is not prescribed, and the operator assembled from its terms.
    HierarchySideDataOpsReal<NDIM, double> sc_ops(patch_hierarchy, 0, finest_ln);
    HierarchyCellDataOpsReal<NDIM, double> cc_ops(patch_hierarchy, 0, finest_ln);
    sc_ops.copyData(u2_idx, u_idx, /*interior_only*/ false);
    bc_helper->enforceDivergenceFreeConditionAtBoundary(u2_idx);
    double max_ghost_difference = 0.0;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> u2_data = patch->getPatchData(u2_idx);
            Pointer<CellData<NDIM, double>> p_data = patch->getPatchData(p_idx);
            const tbox::Array<BoundaryBox<NDIM>>& bdry_boxes = pgeom->getCodimensionBoundaries(1);
            for (int k = 0; k < bdry_boxes.size(); ++k)
            {
                const unsigned int axis = bdry_boxes[k].getLocationIndex() / 2;
                const bool is_lower = bdry_boxes[k].getLocationIndex() % 2 == 0;
                const Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(bdry_boxes[k]);
                Pointer<ArrayData<NDIM, double>> a_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> b_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                Pointer<ArrayData<NDIM, double>> g_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                u_bc_coefs[axis]->setBcCoefs(
                    a_data, b_data, g_data, Pointer<Variable<NDIM>>(), *patch, bdry_boxes[k], fill_time);
                // Limit the boundary cells to the tangential extent of the patch.
                Box<NDIM> cell_box = patch->getBox();
                cell_box.lower(axis) -= 1;
                cell_box.upper(axis) += 1;
                for (Box<NDIM>::Iterator it(bc_coef_box * cell_box); it; it++)
                {
                    // This is the interior cell abutting a lower boundary and the ghost cell abutting an upper
                    // boundary, so its lower face is the boundary face.
                    if (!IBTK::rel_equal_eps((*b_data)(it(), 0), 1.0))
                    {
                        continue;
                    }
                    SideIndex<NDIM> s_in(it(), axis, SideIndex<NDIM>::Lower), s_out(s_in);
                    s_in(axis) += is_lower ? 1 : -1;
                    s_out(axis) += is_lower ? -1 : 1;
                    max_ghost_difference =
                        std::max(max_ghost_difference, std::abs((*u2_data)(s_out) - (*u_data)(s_out)));
                    const double du_dx =
                        (is_lower ? 1.0 : -1.0) * ((*u2_data)(s_in) - (*u2_data)(s_out)) / (2.0 * pgeom->getDx()[axis]);
                    hier::Index<NDIM> i_in(it()), i_ghost(it());
                    (is_lower ? i_ghost : i_in)(axis) -= 1;
                    (*p_data)(i_ghost) = 2.0 * (2.0 * mu * du_dx - (*g_data)(it(), 0)) - (*p_data)(i_in);
                }
            }
        }
    }
    hier_math_ops.grad(y2u_idx, u_var, /*cf_bdry_synch*/ false, 1.0, p_idx, p_var, nullptr, fill_time);
    hier_math_ops.laplace(y2u_idx, u_var, spec, u2_idx, u_var, nullptr, fill_time, 1.0, y2u_idx, u_var);
    hier_math_ops.div(y2p_idx, p_var, -1.0, u2_idx, u_var, nullptr, fill_time, /*cf_bdry_synch*/ true);
    bc_helper->copyDataAtDirichletBoundaries(y2u_idx, u2_idx);

    // Compare.
    const double velocity_norm = sc_ops.maxNorm(y1u_idx, wgt_sc_idx);
    sc_ops.subtract(y2u_idx, y1u_idx, y2u_idx);
    cc_ops.subtract(y2p_idx, y1p_idx, y2p_idx);
    const double difference = std::max(sc_ops.maxNorm(y2u_idx, wgt_sc_idx), cc_ops.maxNorm(y2p_idx, wgt_cc_idx));
    if (max_ghost_difference == 0.0 || !std::isfinite(max_ghost_difference) || !std::isfinite(velocity_norm) ||
        !std::isfinite(difference))
    {
        TBOX_ERROR("the normal velocity ghost values do not differ or a result is not finite\n");
    }
    plog << std::setprecision(10) << "Largest difference between normal velocity ghost values: " << max_ghost_difference
         << "\nLargest magnitude of the velocity component of the result: " << velocity_norm
         << "\nLargest difference between the two results: " << difference << "\n";
} // check_open_boundary_ghost_values

void
output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
            Pointer<INSHierarchyIntegrator> ins_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(format_samrai_output_filename(iteration_num, data_dump_dirname, "hier_data"));

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getVelocityVariable(),
                                                           ins_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(ins_integrator->getPressureVariable(),
                                                           ins_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data

void
check_traction_corner_accuracy(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                               Pointer<INSHierarchyIntegrator> ins_integrator,
                               Pointer<Database> exact_velocity_db)
{
    // Fill the ghost values of a velocity field that satisfies the boundary conditions and report the
    // error in the tangential ghost values outside the x boundaries: at the lower corners, at the upper corners,
    // and elsewhere.
    ParsedVelocity exact_velocity(exact_velocity_db);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const IntVector<NDIM> ghost_width(3);
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("traction_corner_accuracy"), ghost_width);
    CartSideRobinPhysBdryOp bdry_op(u_idx, ins_integrator->getVelocityBoundaryConditions(), /*homogeneous_bc*/ false);
    const double fill_time = ins_integrator->getIntegratorTime();

    double lower_end_error = 0.0;
    double upper_end_error = 0.0;
    double other_error = 0.0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            BoxArray<NDIM> domain = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
            domain.refine(pgeom->getRatio());
            const auto exact = [&](const unsigned int axis, const hier::Index<NDIM>& i)
            {
                std::array<double, NDIM> X;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[d] = x_lower[d] + dx[d] * (i(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                }
                return exact_velocity(axis, X);
            };

            // Set the values at the faces of cells in the physical domain,
            // including in the ghost cells, as a ghost cell fill from the other
            // patches would.
            u_data->fillAll(0.0);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    hier::Index<NDIM> i_lower = b();
                    i_lower(axis) -= 1;
                    if (domain.contains(b()) || domain.contains(i_lower))
                    {
                        u_data->getArrayData(axis)(b(), 0) = exact(axis, b());
                    }
                }
            }
            bdry_op.setPatchDataIndex(u_idx);
            bdry_op.setPhysicalBoundaryConditions(*patch, fill_time, ghost_width);

            // Tangential ghost values outside the x boundaries, within the
            // extent of the patch in the other directions.
            for (unsigned int axis = 1; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    const hier::Index<NDIM>& i = b();
                    bool beyond_x_boundary = i(0) < patch_box.lower(0) || i(0) > patch_box.upper(0);
                    for (unsigned int d = 1; d < NDIM; ++d)
                    {
                        beyond_x_boundary = beyond_x_boundary && side_box.lower(d) <= i(d) && i(d) <= side_box.upper(d);
                    }
                    if (!beyond_x_boundary)
                    {
                        continue;
                    }
                    // The ghost value is set using the normal velocity at the
                    // boundary faces in the rows i(axis) - 1 and i(axis). It is
                    // at a lower or upper corner of the boundary if the lower or
                    // upper row is outside the physical domain.
                    hier::Index<NDIM> i_boundary = i;
                    i_boundary(0) = i(0) < patch_box.lower(0) ? patch_box.lower(0) : patch_box.upper(0);
                    hier::Index<NDIM> i_boundary_lower = i_boundary;
                    i_boundary_lower(axis) -= 1;
                    const double error = std::abs(u_data->getArrayData(axis)(i, 0) - exact(axis, i));
                    if (!domain.contains(i_boundary_lower))
                    {
                        lower_end_error = std::max(lower_end_error, error);
                    }
                    else if (!domain.contains(i_boundary))
                    {
                        upper_end_error = std::max(upper_end_error, error);
                    }
                    else
                    {
                        other_error = std::max(other_error, error);
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
    }
    pout << std::setprecision(6) << std::scientific
         << "max tangential ghost error at the lower corners of the x boundaries = "
         << IBTK_MPI::maxReduction(lower_end_error) << '\n'
         << "max tangential ghost error at the upper corners of the x boundaries = "
         << IBTK_MPI::maxReduction(upper_end_error) << '\n'
         << "max tangential ghost error elsewhere on the x boundaries            = "
         << IBTK_MPI::maxReduction(other_error) << '\n';
    return;
} // check_traction_corner_accuracy

void
check_divergence_free_extension(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                Pointer<INSStaggeredHierarchyIntegrator> ins_integrator,
                                Pointer<Database> exact_velocity_db,
                                const int extension_ghost_width)
{
    // Fill the ghost values outside the physical boundaries of a level of
    // several patches with the divergence-free extension, using a ghost cell
    // fill that copies the values at the other patches first, and check, within
    // the width G = extension_ghost_width of each patch, that
    //  (a) the divergence in the ghost cells outside the domain vanishes for a
    //      smooth velocity and for a velocity that is not smooth,
    //  (b) every patch that stores a ghost face outside the domain stores the
    //      same value for it. The comparison is among the patches on one rank
    //      and does not distinguish levels, so it requires one rank and one
    //      level,
    //  (c) the accumulation of values outside the domain is the adjoint of the
    //      extension with homogeneous boundary conditions: for values u in the
    //      domain and y outside it, (y, E u) = (E^T y, u), where the values in
    //      the domain include those in the ghost cells inside the domain that
    //      the patches store, and it sets to zero every value outside the
    //      domain that the extension fills,
    //  (d) the extension of a smooth divergence-free velocity whose boundary
    //      data match it approximates the velocity, by the number of boundaries
    //      that the ghost face lies beyond, and
    //  (e) all ghost values are computable.
    // The velocity has ghost width G + 1 so that all ghost values within the
    // width G are computable. A direction in which the domain is periodic has
    // no boundary: ghost values there are copies, not extensions.
    if (IBTK_MPI::getNodes() > 1 || patch_hierarchy->getFinestLevelNumber() > 0)
    {
        TBOX_ERROR(
            "check_divergence_free_extension: the comparison of the values stored by different patches "
            "requires one rank and one level, but there are "
            << IBTK_MPI::getNodes() << " ranks and " << patch_hierarchy->getFinestLevelNumber() + 1 << " levels.\n");
    }
    const int G = extension_ghost_width;
    ParsedVelocity exact_velocity(exact_velocity_db);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const IntVector<NDIM> ghost_width(G + 1);
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("divergence_free_extension"), ghost_width);
    const int y_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    INSStaggeredDivergenceFreePhysBdryOp bdry_op(u_idx, ins_integrator, /*homogeneous_bc*/ false);
    INSStaggeredDivergenceFreePhysBdryOp homogeneous_bdry_op(u_idx, ins_integrator, /*homogeneous_bc*/ true);
    const double fill_time = ins_integrator->getIntegratorTime();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // A deterministic function of the face that is not smooth.
    const auto rough = [](const unsigned int axis, const hier::Index<NDIM>& i, const double scale)
    {
        static const std::array<double, 3> coefficients = { 12.9898, 78.233, 37.719 };
        double phase = 4.1 * axis;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            phase += coefficients[d] * i(d);
        }
        return std::sin(scale * phase);
    };

    double smooth_divergence = 0.0;
    double rough_divergence = 0.0;
    double patch_difference = 0.0;
    double y_dot_Eu = 0.0;
    double ETy_dot_u = 0.0;
    double max_filled_after_adjoint = 0.0;
    double num_nonfinite_divergence = 0.0;
    std::array<double, NDIM> error_by_region_size;
    error_by_region_size.fill(0.0);
    double num_nan = 0.0;
    double num_nan_outermost = 0.0;
    std::map<std::array<int, NDIM + 1>, double> face_values;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx, fill_time);
        level->allocatePatchData(y_idx, fill_time);
        BoxArray<NDIM> domain = patch_hierarchy->getGridGeometry()->getPhysicalDomain();
        domain.refine(level->getRatio());
        const Box<NDIM> domain_box = domain[0];
        const IntVector<NDIM> periodic_shift = patch_hierarchy->getGridGeometry()->getPeriodicShift(level->getRatio());

        // Whether the face i of component axis lies beyond a boundary of the
        // domain. Taking axis = NDIM gives the cell i.
        const auto is_beyond_boundary = [&](const unsigned int axis, const hier::Index<NDIM>& i)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (periodic_shift(d) == 0 &&
                    (i(d) < domain_box.lower(d) || i(d) > domain_box.upper(d) + (d == axis ? 1 : 0)))
                {
                    return true;
                }
            }
            return false;
        };

        Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
        ghost_fill_alg->registerRefine(u_idx, u_idx, u_idx, nullptr);
        Pointer<RefineSchedule<NDIM>> ghost_fill_sched = ghost_fill_alg->createSchedule(level, &bdry_op);
        Pointer<RefineSchedule<NDIM>> homogeneous_ghost_fill_sched =
            ghost_fill_alg->createSchedule(level, &homogeneous_bdry_op);

        // Call fcn(patch, axis, i) for each side index i of the patch data
        // within width of the patch that is outside the domain.
        const auto for_each_exterior_face = [&](Pointer<Patch<NDIM>> patch, const int width, const auto& fcn)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> box = SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch->getBox(), width), axis);
                for (Box<NDIM>::Iterator b(box); b; b++)
                {
                    if (is_beyond_boundary(axis, b()))
                    {
                        fcn(axis, b());
                    }
                }
            }
        };

        // Set the values in the patch interiors, and NaN elsewhere so that
        // any dependence on unset values appears in the results.
        const auto fill_ghosts = [&](const auto& value, const Pointer<RefineSchedule<NDIM>>& sched)
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
                u_data->fillAll(nan);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
                    {
                        u_data->getArrayData(axis)(b(), 0) = value(*patch, axis, b());
                    }
                }
            }
            sched->fillData(fill_time);
        };

        // The maximum divergence in the ghost cells outside the domain.
        const auto max_exterior_divergence = [&]()
        {
            double max_divergence = 0.0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
                Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
                for (Box<NDIM>::Iterator b(Box<NDIM>::grow(patch->getBox(), G)); b; b++)
                {
                    if (!is_beyond_boundary(NDIM, b()))
                    {
                        continue;
                    }
                    double divergence = 0.0;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        hier::Index<NDIM> i_upper = b();
                        i_upper(axis) += 1;
                        divergence += (u_data->getArrayData(axis)(i_upper, 0) - u_data->getArrayData(axis)(b(), 0)) /
                                      pgeom->getDx()[axis];
                    }
                    if (std::isfinite(divergence))
                    {
                        max_divergence = std::max(max_divergence, std::abs(divergence));
                    }
                    else
                    {
                        num_nonfinite_divergence += 1.0;
                    }
                }
            }
            return max_divergence;
        };

        const auto rough_value = [&](const Patch<NDIM>&, const unsigned int axis, const hier::Index<NDIM>& i)
        { return rough(axis, i, 1.0); };

        // Smooth velocity and boundary data: (a), (b), (d), (e).
        const auto exact = [&](const Patch<NDIM>& patch, const unsigned int axis, const hier::Index<NDIM>& i)
        {
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
            std::array<double, NDIM> X;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                X[d] = pgeom->getXLower()[d] +
                       pgeom->getDx()[d] * (i(d) - patch.getBox().lower(d) + (d == axis ? 0.0 : 0.5));
            }
            return exact_velocity(axis, X);
        };
        bdry_op.setPatchDataIndex(u_idx);
        fill_ghosts(exact, ghost_fill_sched);
        smooth_divergence = std::max(smooth_divergence, max_exterior_divergence());
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            for_each_exterior_face(
                patch,
                G + 1,
                [&](const unsigned int axis, const hier::Index<NDIM>& i)
                {
                    const double value = u_data->getArrayData(axis)(i, 0);
                    if (std::isnan(value))
                    {
                        const Box<NDIM> inner =
                            SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch->getBox(), G), axis);
                        if (inner.contains(i))
                        {
                            num_nan += 1.0;
                        }
                        else
                        {
                            num_nan_outermost += 1.0;
                        }
                        return;
                    }
                    const Box<NDIM> inner = SideGeometry<NDIM>::toSideBox(Box<NDIM>::grow(patch->getBox(), G), axis);
                    if (!inner.contains(i))
                    {
                        return;
                    }
                    std::array<int, NDIM + 1> key;
                    key[0] = static_cast<int>(axis);
                    int region_size = 0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        key[d + 1] = i(d);
                        if (periodic_shift(d) == 0 &&
                            (i(d) < domain_box.lower(d) || i(d) > domain_box.upper(d) + (d == axis ? 1 : 0)))
                        {
                            ++region_size;
                        }
                    }
                    const auto it = face_values.find(key);
                    if (it == face_values.end())
                    {
                        face_values[key] = value;
                    }
                    else
                    {
                        patch_difference = std::max(patch_difference, std::abs(value - it->second));
                    }
                    error_by_region_size[region_size - 1] =
                        std::max(error_by_region_size[region_size - 1], std::abs(value - exact(*patch, axis, i)));
                });
        }

        // Velocity that is not smooth: (a).
        fill_ghosts(rough_value, ghost_fill_sched);
        rough_divergence = std::max(rough_divergence, max_exterior_divergence());

        // Adjoint: (c). The values y outside the domain within width G are
        // arbitrary, and the values u in the domain are not smooth.
        homogeneous_bdry_op.setPatchDataIndex(u_idx);
        fill_ghosts(rough_value, homogeneous_ghost_fill_sched);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<SideData<NDIM, double>> y_data = patch->getPatchData(y_idx);
            y_data->fillAll(0.0);
            for_each_exterior_face(patch,
                                   G,
                                   [&](const unsigned int axis, const hier::Index<NDIM>& i)
                                   {
                                       y_data->getArrayData(axis)(i, 0) = rough(axis, i, 1.3) + 0.2;
                                       y_dot_Eu += y_data->getArrayData(axis)(i, 0) * u_data->getArrayData(axis)(i, 0);
                                   });
            homogeneous_bdry_op.setPatchDataIndex(y_idx);
            homogeneous_bdry_op.accumulateFromPhysicalBoundaryData(*patch, fill_time, ghost_width);
            homogeneous_bdry_op.setPatchDataIndex(u_idx);

            // The accumulation sets to zero every value outside the domain that
            // the extension fills, which here is every computable value within
            // the ghost width of the patch data.
            for_each_exterior_face(
                patch,
                G + 1,
                [&](const unsigned int axis, const hier::Index<NDIM>& i)
                {
                    const double value = y_data->getArrayData(axis)(i, 0);
                    if (!std::isfinite(value))
                    {
                        TBOX_ERROR("check_divergence_free_extension: the value after the accumulation of component "
                                   << axis << " at the face " << i << " outside the domain is not finite.\n");
                    }
                    max_filled_after_adjoint = std::max(max_filled_after_adjoint, std::abs(value));
                });
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u_data->getGhostBox(), axis)); b; b++)
                {
                    if (!is_beyond_boundary(axis, b()))
                    {
                        ETy_dot_u += y_data->getArrayData(axis)(b(), 0) * u_data->getArrayData(axis)(b(), 0);
                    }
                }
            }
        }
        level->deallocatePatchData(u_idx);
        level->deallocatePatchData(y_idx);
    }

    pout << std::setprecision(12)
         << "max |Div_h u| outside the domain, smooth velocity     = " << IBTK_MPI::maxReduction(smooth_divergence)
         << '\n'
         << "max |Div_h u| outside the domain, non-smooth velocity = " << IBTK_MPI::maxReduction(rough_divergence)
         << '\n'
         << "max difference between patches in a ghost value        = " << IBTK_MPI::maxReduction(patch_difference)
         << '\n'
         << "(y, E u)   = " << IBTK_MPI::sumReduction(y_dot_Eu) << '\n'
         << "(E^T y, u) = " << IBTK_MPI::sumReduction(ETy_dot_u) << '\n'
         << "max |value| outside the domain after the accumulation  = "
         << IBTK_MPI::maxReduction(max_filled_after_adjoint) << '\n';
    for (unsigned int n = 0; n < NDIM; ++n)
    {
        pout << "max error in a ghost value beyond " << n + 1 << (n == 0 ? " boundary   = " : " boundaries = ")
             << IBTK_MPI::maxReduction(error_by_region_size[n]) << '\n';
    }
    pout << "number of non-finite ghost-cell divergences within width G = "
         << IBTK_MPI::sumReduction(num_nonfinite_divergence) << '\n'
         << "number of NaN ghost values within width G        = " << IBTK_MPI::sumReduction(num_nan) << '\n'
         << "number of NaN ghost values in the outermost layer = " << IBTK_MPI::sumReduction(num_nan_outermost) << '\n';
    return;
} // check_divergence_free_extension
