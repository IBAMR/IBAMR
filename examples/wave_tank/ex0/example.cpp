// Copyright (c) 2002-2018, Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/StokesFirstOrderWaveBcCoef.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Application
#include "GravityForcing.h"
#include "LSLocateColumnInterface.h"
#include "SetDampingZone.h"
#include "SetFluidProperties.h"
#include "SetLSProperties.h"
#include "TagLSRefinementCells.h"

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

void damp_outlet_waves(double current_time,
                       double new_time,
                       bool skip_synchronize_new_state_data,
                       int num_cycles,
                       void* ctx);

struct WaveDamper
{
    Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_integrator;
    double d_x_zone_start, d_x_zone_end;
};

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
bool
run_example(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            time_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            time_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "SEMI_IMPLICIT");
        if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: SEMI_IMPLICIT");
        }
        time_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
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

        // Setup level set information
        ColumnInterface column;
        input_db->getDoubleArray("X_UR", &column.X_UR[0], NDIM);

        const string& ls_name = "level_set";
        Pointer<CellVariable<NDIM, double> > phi_var = new CellVariable<NDIM, double>(ls_name);
        adv_diff_integrator->registerTransportedQuantity(phi_var);
        adv_diff_integrator->setDiffusionCoefficient(phi_var, 0.0);
        // Set the advection velocity of the bubble.
        adv_diff_integrator->setAdvectionVelocity(phi_var, time_integrator->getAdvectionVelocityVariable());

        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateColumnInterface* ptr_LSLocateColumnInterface =
            new LSLocateColumnInterface("LSLocateColumnInterface", adv_diff_integrator, phi_var, column);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateColumnInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateColumnInterface));
        SetLSProperties* ptr_SetLSProperties = new SetLSProperties("SetLSProperties", level_set_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var, &callSetLSCallbackFunction, static_cast<void*>(ptr_SetLSProperties));

        // Setup the INS maintained material properties.
        Pointer<Variable<NDIM> > rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariable<NDIM, double>("rho");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho");
        }

        time_integrator->registerMassDensityVariable(rho_var);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_inside = input_db->getDouble("RHO_I");
        const double rho_outside = input_db->getDouble("RHO_O");
        const double mu_inside = input_db->getDouble("MU_I");
        const double mu_outside = input_db->getDouble("MU_O");
        const int ls_reinit_interval = input_db->getInteger("LS_REINIT_INTERVAL");
        const double num_interface_cells = input_db->getDouble("NUM_INTERFACE_CELLS");

        // Callback functions can either be registered with the NS integrator, or the advection-diffusion integrator
        // Note that these will set the initial conditions for density and viscosity, based on level set information
        SetFluidProperties* ptr_SetFluidProperties = new SetFluidProperties("SetFluidProperties",
                                                                            adv_diff_integrator,
                                                                            phi_var,
                                                                            rho_outside,
                                                                            rho_inside,
                                                                            mu_outside,
                                                                            mu_inside,
                                                                            ls_reinit_interval,
                                                                            num_interface_cells);
        time_integrator->registerResetFluidDensityFcn(&callSetFluidDensityCallbackFunction,
                                                      static_cast<void*>(ptr_SetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&callSetFluidViscosityCallbackFunction,
                                                        static_cast<void*>(ptr_SetFluidProperties));

        // Register callback function for tagging refined cells for level set data
        const double tag_value = input_db->getDouble("LS_TAG_VALUE");
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        TagLSRefinementCells ls_tagger;
        ls_tagger.d_ls_gas_var = phi_var;
        ls_tagger.d_tag_value = tag_value;
        ls_tagger.d_tag_abs_thresh = tag_thresh;
        ls_tagger.d_adv_diff_solver = adv_diff_integrator;
        time_integrator->registerApplyGradientDetectorCallback(&callTagLSRefinementCellsCallbackFunction,
                                                               static_cast<void*>(&ls_tagger));

        // Create Eulerian initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            time_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new StokesFirstOrderWaveBcCoef(bc_coefs_name,
                                                               /*comp*/ d,
                                                               app_initializer->getComponentDatabase(bc_coefs_db_name),
                                                               grid_geometry,
                                                               adv_diff_integrator,
                                                               phi_var);
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create a damping zone near the channel outlet to absorb water waves through momentum equation.
        // This method adds an implicit damping term in the momentum equation.
        DampingZone* ptr_DampingZone = new DampingZone(
            "OutletDampingZone", time_integrator, app_initializer->getComponentDatabase("OutletDampingZone"));
        time_integrator->registerAddDampingZoneFcn(&callDampingZoneCallbackFunction,
                                                   static_cast<void*>(ptr_DampingZone));

        // Create a damping zone near the channel outlet to absorb water waves via time-splitting approach.
        // This method uses a time splitting approach and modifies the fluid momentum in the
        // post processing step.
        WaveDamper smooth_damper;
        smooth_damper.d_ins_integrator = time_integrator;
        smooth_damper.d_x_zone_start = input_db->getDouble("X_ZONE_START");
        smooth_damper.d_x_zone_end = input_db->getDouble("X_ZONE_END");
        time_integrator->registerPostprocessIntegrateHierarchyCallback(&damp_outlet_waves,
                                                                       static_cast<void*>(&smooth_damper));

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(phi_var, phi_bc_coef);
        }
        level_set_ops->registerPhysicalBoundaryCondition(phi_bc_coef);

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("RhoBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("RhoBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("MuBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("MuBcCoefs"), grid_geometry);
            time_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        // Forcing terms
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        Pointer<CartGridFunction> grav_force = new GravityForcing("GravityForcing", time_integrator, grav_const);

        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var);

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(surface_tension_force);
        eul_forces->addFunction(grav_force);
        time_integrator->registerBodyForceFunction(eul_forces);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Remove the AppInitializer
        app_initializer.setNull();

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

        // File to write to for fluid mass data
        ofstream mass_file;

        if (!SAMRAI_MPI::getRank())
        {
            mass_file.open("mass_fluid.txt");
        }
        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
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

            // Compute the fluid mass in the domain from interpolated density
            const int rho_ins_idx = time_integrator->getLinearOperatorRhoPatchDataIndex();
#if !defined(NDEBUG)
            TBOX_ASSERT(rho_ins_idx >= 0);
#endif
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            HierarchySideDataOpsReal<NDIM, double> hier_rho_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const double mass_fluid = hier_rho_data_ops.integral(rho_ins_idx, wgt_sc_idx);

            // Write to file
            if (!SAMRAI_MPI::getRank())
            {
                mass_file << std::setprecision(13) << loop_time << "\t" << mass_fluid << std::endl;
            }

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

        // Close file
        if (!SAMRAI_MPI::getRank())
        {
            mass_file.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        // Cleanup other dumb pointers
        delete ptr_SetLSProperties;
        delete ptr_SetFluidProperties;
        delete ptr_LSLocateColumnInterface;
        delete ptr_DampingZone;

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return true;
} // run_example

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                           time_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getPressureVariable(),
                                                           time_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data

void
damp_outlet_waves(double current_time, double new_time, bool skip_synchronize_new_state_data, int num_cycles, void* ctx)
{
    pout << "Damping waves at outlet ...\n";
    WaveDamper* ptr_wave_damper = static_cast<WaveDamper*>(ctx);
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = ptr_wave_damper->d_ins_integrator->getPatchHierarchy();
    const double& x_zone_start = ptr_wave_damper->d_x_zone_start;
    const double& x_zone_end = ptr_wave_damper->d_x_zone_end;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int u_new_idx = var_db->mapVariableAndContextToIndex(ptr_wave_damper->d_ins_integrator->getVelocityVariable(),
                                                         ptr_wave_damper->d_ins_integrator->getNewContext());

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_new_idx);

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    Index<NDIM> i = it();
                    SideIndex<NDIM> i_side(i, axis, SideIndex<NDIM>::Lower);
                    const double x_posn =
                        patch_x_lower[0] +
                        patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + axis == 0 ? 0.0 : 0.5);
                    if (x_posn >= x_zone_start)
                    {
                        double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                        double gamma = 1.0 - (exp(std::pow(xtilde, 1.5)) - 1.0) / (exp(1.0) - 1.0);
                        (*u_data)(i_side, 0) *= gamma;
                    }
                }
            }
        }
    }

    return;
} // damp_outlet_waves
