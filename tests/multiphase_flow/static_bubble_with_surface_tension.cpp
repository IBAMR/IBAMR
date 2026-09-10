// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <ibamr/app_namespaces.h>

// Application
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/vc_ins_utilities.h>

#include "LSLocateCircularInterface.h"

#if (NDIM == 2)
#define SC_NORMAL_FC IBAMR_FC_FUNC(sc_normal_2d, SC_NORMAL_2D)
#else
#define SC_NORMAL_FC IBAMR_FC_FUNC(sc_normal_3d, SC_NORMAL_3D)
#endif

extern "C"
{
    void SC_NORMAL_FC(double* N00,
                      double* N01,
#if (NDIM == 3)
                      double* N02,
#endif
                      double* N10,
                      double* N11,
#if (NDIM == 3)
                      double* N12,
                      double* N20,
                      double* N21,
                      double* N22,
#endif
                      const int& N_gcw,
                      const double* U,
                      const int& U_gcw,
                      const int& ilower0,
                      const int& iupper0,
                      const int& ilower1,
                      const int& iupper1,
#if (NDIM == 3)
                      const int& ilower2,
                      const int& iupper2,
#endif
                      const double* dx);
}

// Exercise the library kernel with actual SAMRAI cell and side data, including
// the unequal ghost widths used for the Marangoni temperature gradient.
void
check_surface_tension_gradient()
{
    const Box<NDIM> box(hier::Index<NDIM>(-3), hier::Index<NDIM>(4));
    const double dx[3] = { 0.2, 0.3, 0.4 };
    for (int input_gcw = 1; input_gcw <= 3; ++input_gcw)
    {
        CellData<NDIM, double> scalar(box, 1, IntVector<NDIM>(input_gcw));
        for (Box<NDIM>::Iterator it(scalar.getGhostBox()); it; it++)
        {
            const CellIndex<NDIM> ci(it());
            double value = 0.0;
            for (int d = 0; d < NDIM; ++d)
            {
                const double x = (ci(d) + 0.5) * dx[d];
                value += (d + 1.0) * x * x + (d + 2.0) * x;
            }
            value += (ci(0) + 0.5) * dx[0] * (ci(1) + 0.5) * dx[1];
            scalar(ci) = value;
        }
        for (int output_gcw = 1; output_gcw <= 3; ++output_gcw)
        {
            SideData<NDIM, double> gradient(box, NDIM, IntVector<NDIM>(output_gcw));
            gradient.fillAll(std::numeric_limits<double>::quiet_NaN());
            SC_NORMAL_FC(gradient.getPointer(0, 0),
                         gradient.getPointer(0, 1),
#if (NDIM == 3)
                         gradient.getPointer(0, 2),
#endif
                         gradient.getPointer(1, 0),
                         gradient.getPointer(1, 1),
#if (NDIM == 3)
                         gradient.getPointer(1, 2),
                         gradient.getPointer(2, 0),
                         gradient.getPointer(2, 1),
                         gradient.getPointer(2, 2),
#endif
                         output_gcw,
                         scalar.getPointer(),
                         input_gcw,
                         box.lower(0),
                         box.upper(0),
                         box.lower(1),
                         box.upper(1),
#if (NDIM == 3)
                         box.lower(2),
                         box.upper(2),
#endif
                         dx);
            double max_error = 0.0;
            Box<NDIM> valid_box = box;
            valid_box.grow(std::min(input_gcw, output_gcw) - 1);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(valid_box, axis)); it; it++)
                {
                    const SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);
                    double x[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        x[d] = (si(d) + (d == axis ? 0.0 : 0.5)) * dx[d];
                    }
                    for (int d = 0; d < NDIM; ++d)
                    {
                        double expected = 2.0 * (d + 1.0) * x[d] + d + 2.0;
                        if (d == 0)
                        {
                            expected += x[1];
                        }
                        else if (d == 1)
                        {
                            expected += x[0];
                        }
                        const double computed = gradient(si, d);
                        if (!std::isfinite(computed))
                        {
                            TBOX_ERROR("Incorrect surface-tension gradient with input ghost width "
                                       << input_gcw << " and output ghost width " << output_gcw << "\n");
                        }
                        max_error = std::max(max_error, std::abs(computed - expected));
                    }
                }
            }
            plog << std::setprecision(13) << "Input ghost width " << input_gcw << ", output ghost width " << output_gcw
                 << ": maximum gradient error = " << max_error << '\n';
        }
    }
}

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator,
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
template <class BaseIntegrator>
class ProjectionTestIntegrator : public BaseIntegrator
{
public:
    using BaseIntegrator::BaseIntegrator;

    int getRegridProjectionCount() const
    {
        return d_regrid_projection_count;
    }

protected:
    void regridProjection(const bool initial_time) override
    {
        BaseIntegrator::regridProjection(initial_time);
        if (!initial_time)
        {
            ++d_regrid_projection_count;
        }
    }

private:
    int d_regrid_projection_count = 0;
};

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        if (input_db->getBoolWithDefault("check_gradient", false))
        {
            PIO::logOnlyNodeZero("output");
            check_surface_tension_gradient();
            return 0;
        }

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
        const bool check_regrid = input_db->getBoolWithDefault("check_regrid_projection", false);
        Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            time_integrator = new ProjectionTestIntegrator<INSVCStaggeredConservativeHierarchyIntegrator>(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            time_integrator = new ProjectionTestIntegrator<INSVCStaggeredNonConservativeHierarchyIntegrator>(
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

        // Get exact pressure jump
        const double dP_exact = input_db->getDouble("dP_exact");

        // Setup level set information
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        const string& ls_name = "level_set";
        Pointer<CellVariable<NDIM, double>> phi_var = new CellVariable<NDIM, double>(ls_name);
        adv_diff_integrator->registerTransportedQuantity(phi_var);
        adv_diff_integrator->setDiffusionCoefficient(phi_var, 0.0);
        // Set the advection velocity of the bubble.
        adv_diff_integrator->setAdvectionVelocity(phi_var, time_integrator->getAdvectionVelocityVariable());

        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateCircularInterface* ptr_LSLocateCircularInterface =
            new LSLocateCircularInterface("LSLocateCircularInterface", adv_diff_integrator, phi_var, &circle);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateCircularInterfaceCallbackFunction,
                                                                static_cast<void*>(ptr_LSLocateCircularInterface));
        IBAMR::LevelSetUtilities::SetLSProperties set_ls_properties("SetLSProperties", level_set_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var, &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(&set_ls_properties));

        // LS initial conditions
        if (input_db->keyExists("LevelSetInitialConditions"))
        {
            Pointer<CartGridFunction> phi_init = new muParserCartGridFunction(
                "phi_init", app_initializer->getComponentDatabase("LevelSetInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(phi_var, phi_init);
        }

        // Setup the INS maintained material properties.
        Pointer<Variable<NDIM>> rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariable<NDIM, double>("rho");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho");
        }
        time_integrator->registerMassDensityVariable(rho_var);

        Pointer<CellVariable<NDIM, double>> mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_inside = input_db->getDouble("RHO_I");
        const double rho_outside = input_db->getDouble("RHO_O");
        const double mu_inside = input_db->getDouble("MU_I");
        const double mu_outside = input_db->getDouble("MU_O");
        const double num_interface_cells = input_db->getDouble("NUM_INTERFACE_CELLS");

        // Callback functions can either be registered with the NS integrator, or the advection-diffusion integrator
        IBAMR::VCINSUtilities::SetFluidProperties set_fluid_properties("SetFluidProperties",
                                                                       adv_diff_integrator,
                                                                       phi_var,
                                                                       rho_outside,
                                                                       rho_inside,
                                                                       mu_outside,
                                                                       mu_inside,
                                                                       num_interface_cells);
        time_integrator->registerResetFluidDensityFcn(&IBAMR::VCINSUtilities::callSetDensityCallbackFunction,
                                                      static_cast<void*>(&set_fluid_properties));
        time_integrator->registerResetFluidViscosityFcn(&IBAMR::VCINSUtilities::callSetViscosityCallbackFunction,
                                                        static_cast<void*>(&set_fluid_properties));
        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            time_integrator->registerVelocityInitialConditions(u_init);
        }

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

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(phi_var, phi_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("RhoBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("RhoBcCoefs"), grid_geometry);
            time_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("MuBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("MuBcCoefs"), grid_geometry);
            time_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        // Set up the surface tension force
        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var);
        time_integrator->registerBodyForceFunction(surface_tension_force);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
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

        // File to write errors.
        std::ostringstream out;

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

            // Compute the fluid mass in the domain from interpolated density.
            const int rho_ins_idx = time_integrator->getLinearOperatorRhoPatchDataIndex();
#if !defined(NDEBUG)
            TBOX_ASSERT(rho_ins_idx >= 0);
#endif
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const double fluid_mass = hier_sc_data_ops.integral(rho_ins_idx, wgt_sc_idx);

            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int p_idx = var_db->mapVariableAndContextToIndex(time_integrator->getPressureVariable(),
                                                                   time_integrator->getCurrentContext());
            HierarchyCellDataOpsReal<NDIM, double> hier_p_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            const double dP = hier_p_data_ops.max(p_idx) - hier_p_data_ops.min(p_idx);

            const int u_idx = var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                                   time_integrator->getCurrentContext());
            const double Umax = hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx);
            const double UL1 = hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx);

            if (!std::isfinite(fluid_mass) || !std::isfinite(dP) || !std::isfinite(Umax) || !std::isfinite(UL1))
            {
                TBOX_ERROR("The bubble diagnostics must be finite.\n");
            }

            // Write to file.
            if (!check_regrid && !IBTK_MPI::getRank())
            {
                out << std::setprecision(16) << loop_time << "\t"
                    << "Relative pressure error = |dP - dP_exact|/dP_exact = " << std::abs(dP - dP_exact) / dP_exact
                    << "\t"
                    << "|U|_max = " << Umax << "\t"
                    << "|U|_L1 = " << UL1 << "\t"
                    << "fluid mass = " << fluid_mass << std::endl;
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

        if (check_regrid)
        {
            const int projection_count =
                conservative_form ?
                    dynamic_cast<ProjectionTestIntegrator<INSVCStaggeredConservativeHierarchyIntegrator>*>(
                        time_integrator.getPointer())
                        ->getRegridProjectionCount() :
                    dynamic_cast<ProjectionTestIntegrator<INSVCStaggeredNonConservativeHierarchyIntegrator>*>(
                        time_integrator.getPointer())
                        ->getRegridProjectionCount();
            if (projection_count == 0)
            {
                TBOX_ERROR("The test must perform a projection after regridding.\n");
            }
        }

        PIO::logOnlyNodeZero("output");
        if (check_regrid)
        {
            plog << "Regrid projection completed with finite bubble diagnostics.\n";
        }
        else
        {
            plog << out.str();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }

        // Cleanup other dumb pointers
        delete ptr_LSLocateCircularInterface;

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
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
    snprintf(temp_buf, sizeof(temp_buf), "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
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
