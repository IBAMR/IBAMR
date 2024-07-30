// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/BrinkmanPenalizationRigidBodyDynamics.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SpongeLayerForceFunction.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application specific includes.
#include "LSLocateGasInterface.h"
#include "LevelSetInitialCondition.h"

RectangularInterface rectangle;

// Struct to reset solid level set
struct SolidLevelSetResetter
{
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
    Pointer<CellVariable<NDIM, double> > ls_solid_var;
    Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd;
};

void
reset_solid_level_set_callback_fcn(double current_time, double new_time, int /*cycle_num*/, void* ctx)
{
    SolidLevelSetResetter* resetter = static_cast<SolidLevelSetResetter*>(ctx);

    // Get the new centroid of the body
    const double dt = new_time - current_time;
    Eigen::Vector3d XCOM_current = resetter->bp_rbd->getCurrentCOMPosn();
    Eigen::Vector3d XCOM_new = XCOM_current + dt * (resetter->bp_rbd->getNewCOMTransVelocity());

    IBTK::Vector c, b;
    c << XCOM_new(0), XCOM_new(1);
    b << rectangle.W / 2.0, rectangle.H / 2.0;

    // Set a large value away from the solid body.
    Pointer<PatchHierarchy<NDIM> > patch_hier = resetter->adv_diff_integrator->getPatchHierarchy();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_solid_idx =
        var_db->mapVariableAndContextToIndex(resetter->ls_solid_var, resetter->adv_diff_integrator->getNewContext());

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const hier::Index<NDIM>& ci = it();
                IBTK::Vector coord = IBTK::Vector::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                IBTK::Vector p = coord - c;
                IBTK::Vector d = p.cwiseAbs() - b;
                (*ls_solid_data)(ci) = (d.cwiseMax(0.0)).norm() + std::min(d.maxCoeff(), 0.0);
            }
        }
    }

    return;
}

void
imposed_kinematics(double /*data_time*/,
                   int /*cycle_num*/,
                   Eigen::Vector3d& U_com,
                   Eigen::Vector3d& W_com,
                   void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // imposed_kinematics

void
external_force_torque(double /*data_time*/, int /*cycle_num*/, Eigen::Vector3d& F, Eigen::Vector3d& T, void* /*ctx*/)
{
    F.setZero();
    F[1] = rectangle.rho_solid * rectangle.W * rectangle.H * rectangle.g_y;
    T.setZero();
    return;
} // imposed_kinematics

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
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IBLevelSet.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
        if (dump_restart_data && (restart_dump_interval > 0) && !restart_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(restart_dump_dirname);
        }

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Setup solid information
        rectangle.W = input_db->getDouble("WIDTH_RECTANGLE");
        rectangle.H = input_db->getDouble("HEIGHT_RECTANGLE");
        rectangle.X0[0] = input_db->getDouble("XCOM");
        rectangle.X0[1] = input_db->getDouble("YCOM");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator =
            new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        // Cartesian grid geometry and AMR algorithm objects
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               navier_stokes_integrator,
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

        // Create level sets for solid interface.
        const string& ls_name_solid = "level_set_solid";
        Pointer<CellVariable<NDIM, double> > phi_var_solid = new CellVariable<NDIM, double>(ls_name_solid);

        // Create level sets for gas/liquid interface.
        const double fluid_height = input_db->getDouble("GAS_LS_INIT");
        const string& ls_name_gas = "level_set_gas";
        Pointer<CellVariable<NDIM, double> > phi_var_gas = new CellVariable<NDIM, double>(ls_name_gas);
        Pointer<RelaxationLSMethod> level_set_gas_ops =
            new RelaxationLSMethod(ls_name_gas, app_initializer->getComponentDatabase("LevelSet_Gas"));
        LSLocateGasInterface* ptr_LSLocateGasInterface =
            new LSLocateGasInterface("LSLocateGasInterface", adv_diff_integrator, phi_var_gas, fluid_height);
        level_set_gas_ops->registerInterfaceNeighborhoodLocatingFcn(&IBAMR::callLSLocateGasInterfaceCallbackFunction,
                                                                    static_cast<void*>(ptr_LSLocateGasInterface));

        // Register the level sets with advection diffusion integrator.
        adv_diff_integrator->registerTransportedQuantity(phi_var_solid);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_solid, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_solid,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        adv_diff_integrator->registerTransportedQuantity(phi_var_gas);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_gas, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_gas,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Register the reinitialization functions for the level set variables
        LevelSetUtilities::SetLSProperties setSetLSProperties("SetLSProperties", level_set_gas_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var_gas, &LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(&setSetLSProperties));

        // Solid level set initial conditions
        Pointer<CartGridFunction> phi_solid_init = new LevelSetInitialCondition("solid_ls_init", rectangle);
        adv_diff_integrator->setInitialConditions(phi_var_solid, phi_solid_init);

        // Gas level set initial conditions
        Pointer<CartGridFunction> phi_gas_init = new muParserCartGridFunction(
            "phi_gas_init", app_initializer->getComponentDatabase("GasLevelSetInitialCondition"), grid_geometry);
        adv_diff_integrator->setInitialConditions(phi_var_gas, phi_gas_init);

        // Lagrange multiplier to conserve mass of the phases.
        std::vector<Pointer<CellVariable<NDIM, double> > > ls_vars{ phi_var_gas, phi_var_solid };
        LevelSetUtilities::LevelSetMassLossFixer level_set_fixer(
            "LevelSetMassLossFixer",
            adv_diff_integrator,
            ls_vars,
            app_initializer->getComponentDatabase("LevelSetMassFixer"),
            /*restart*/ true);
        adv_diff_integrator->registerPostprocessIntegrateHierarchyCallback(&LevelSetUtilities::fixMassLoss3PhaseFlows,
                                                                           static_cast<void*>(&level_set_fixer));

        // Reset solid geometry
        SolidLevelSetResetter solid_level_set_resetter;
        solid_level_set_resetter.adv_diff_integrator = adv_diff_integrator;
        solid_level_set_resetter.ls_solid_var = phi_var_solid;
        adv_diff_integrator->registerIntegrateHierarchyCallback(&reset_solid_level_set_callback_fcn,
                                                                static_cast<void*>(&solid_level_set_resetter));

        // Setup the advected and diffused fluid quantities.
        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<hier::Variable<NDIM> > rho_var = new SideVariable<NDIM, double>("rho");
        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_fluid = input_db->getDouble("RHO_F");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const int num_solid_interface_cells = input_db->getDouble("NUM_SOLID_INTERFACE_CELLS");
        const int num_gas_interface_cells = input_db->getDouble("NUM_GAS_INTERFACE_CELLS");
        rectangle.rho_solid = rho_solid;
        const double mu_fluid = input_db->getDouble("MU_F");
        const double mu_gas = input_db->getDouble("MU_G");
        const double mu_solid = input_db->getDoubleWithDefault("MU_S", std::numeric_limits<double>::quiet_NaN());
        const bool set_mu_solid = input_db->getBool("SET_MU_S");

        IBAMR::VCINSUtilities::SetFluidProperties setSetFluidProperties("SetFluidProperties",
                                                                        adv_diff_integrator,
                                                                        phi_var_gas,
                                                                        phi_var_solid,
                                                                        rho_fluid,
                                                                        rho_gas,
                                                                        rho_solid,
                                                                        mu_fluid,
                                                                        mu_gas,
                                                                        mu_solid,
                                                                        num_gas_interface_cells,
                                                                        num_solid_interface_cells,
                                                                        set_mu_solid);
        navier_stokes_integrator->registerResetFluidDensityFcn(&IBAMR::VCINSUtilities::callSetDensityCallbackFunction,
                                                               static_cast<void*>(&setSetFluidProperties));
        navier_stokes_integrator->registerResetFluidViscosityFcn(
            &IBAMR::VCINSUtilities::callSetViscosityCallbackFunction, static_cast<void*>(&setSetFluidProperties));

        // Register callback function for tagging refined cells for level set data
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");

        LevelSetUtilities::TagLSRefinementCells ls_gas_tagger(
            adv_diff_integrator, phi_var_gas, -tag_thresh, tag_thresh);
        LevelSetUtilities::TagLSRefinementCells ls_solid_tagger(
            adv_diff_integrator, phi_var_solid, std::numeric_limits<double>::lowest(), tag_thresh);
        navier_stokes_integrator->registerApplyGradientDetectorCallback(&LevelSetUtilities::tagLSCells,
                                                                        static_cast<void*>(&ls_gas_tagger));
        navier_stokes_integrator->registerApplyGradientDetectorCallback(&LevelSetUtilities::tagLSCells,
                                                                        static_cast<void*>(&ls_solid_tagger));

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<std::unique_ptr<RobinBcCoefStrategy<NDIM> > > u_bc_coefs(NDIM);
        if (periodic_shift.min() == 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = std::make_unique<muParserRobinBcCoefs>(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions({
                u_bc_coefs[0].get(), u_bc_coefs[1].get()
#if (NDIM == 3)
                                         ,
                    u_bc_coefs[2].get()
#endif
            });
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM> > rho_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM> > mu_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef.get());
        }

        std::unique_ptr<RobinBcCoefStrategy<NDIM> > phi_bc_coef;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = std::make_unique<muParserRobinBcCoefs>(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_gas, phi_bc_coef.get());
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef.get());

        // LS reinit boundary conditions, which is set to be the same as the BCs
        // for advection
        RobinBcCoefStrategy<NDIM>* ls_reinit_bcs = phi_bc_coef.get();
        level_set_gas_ops->registerPhysicalBoundaryCondition(ls_reinit_bcs);

        // Class GravityForcing can be utilized to apply the gravitational force $ \rho g $
        // using the density field, which includes all three phases: liquid, gas, and solid;
        // or using the flow density field, which includes only liquid and gas phases
        // and excludes the solid phase.
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        rectangle.g_y = grav_const[1];
        const string grav_type = input_db->getStringWithDefault("GRAV_TYPE", "FLOW");
        Pointer<CartGridFunction> grav_force;
        if (grav_type == "FULL")
        {
            grav_force =
                new IBAMR::VCINSUtilities::GravityForcing("FullGravityForcing", navier_stokes_integrator, grav_const);
        }
        else if (grav_type == "FLOW")
        {
            grav_force =
                new IBAMR::VCINSUtilities::GravityForcing("FlowGravityForcing",
                                                          adv_diff_integrator,
                                                          phi_var_gas,
                                                          app_initializer->getComponentDatabase("FlowGravityForcing"),
                                                          grav_const);
        }
        else
        {
            TBOX_ERROR("Unsupported GRAV_TYPE specified: " << grav_type << "\n"
                                                           << "Valid options are: FLOW, FULL");
        }

        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var_gas);

        Pointer<SpongeLayerForceFunction> sponge_fcn =
            new SpongeLayerForceFunction("SpongeLayerForceFunction",
                                         app_initializer->getComponentDatabase("SpongeLayerForceFunction"),
                                         navier_stokes_integrator,
                                         grid_geometry);

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(grav_force);
        eul_forces->addFunction(surface_tension_force);
        eul_forces->addFunction(sponge_fcn);
        navier_stokes_integrator->registerBodyForceFunction(eul_forces);

        // Configure the Brinkman penalization object to do the rigid body dynamics.
        Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd =
            new BrinkmanPenalizationRigidBodyDynamics("Brinkman Body",
                                                      phi_var_solid,
                                                      adv_diff_integrator,
                                                      navier_stokes_integrator,
                                                      app_initializer->getComponentDatabase("BrinkmanPenalization"),
                                                      /*register_for_restart*/ true);
        FreeRigidDOFVector free_dofs;
        free_dofs << 0, 1, 0;
        Eigen::Vector3d U_i = Eigen::Vector3d::Zero();
        const double mass = rectangle.rho_solid * rectangle.H * rectangle.W;
        bp_rbd->setSolveRigidBodyVelocity(free_dofs);
        bp_rbd->registerKinematicsFunction(&imposed_kinematics);
        bp_rbd->registerExternalForceTorqueFunction(&external_force_torque);
        bp_rbd->setInitialConditions(rectangle.X0, U_i, U_i, mass);
        navier_stokes_integrator->registerBrinkmanPenalizationStrategy(bp_rbd);
        solid_level_set_resetter.bp_rbd = bp_rbd;

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            navier_stokes_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        navier_stokes_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = navier_stokes_integrator->getIntegratorStep();
        double loop_time = navier_stokes_integrator->getIntegratorTime();

        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                navier_stokes_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
        }

        // Open streams to save position and velocity of the structure, as well
        // as the volume of the two phase and the Lagrange multiplier.
        ofstream rbd_stream, ft_stream, vol_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            rbd_stream.open("rbd.curve", ios_base::out | ios_base::app);
            ft_stream.open("hydro_force_torque.curve", ios_base::out | ios_base::app);

            vol_stream.open("vol.curve", ios_base::out | ios_base::app);
            vol_stream.precision(16);
            vol_stream.setf(ios::fixed, ios::floatfield);
        }

        const bool is_from_restart = RestartManager::getManager()->isFromRestart();
        if (!is_from_restart)
        {
            // Target the liquid volume in the Newton's iterations
            std::vector<double> H_integrals =
                LevelSetUtilities::computeHeavisideIntegrals3PhaseFlows(level_set_fixer.getLevelSetContainer());
            level_set_fixer.setInitialVolume(H_integrals[1]);

            // Save the initial volume of gas, liquid, and gas+solid phases, and the Lagrange multiplier in the stream.
            if (SAMRAI_MPI::getRank() == 0)
            {
                vol_stream << 0.0 << "\t" << H_integrals[0] << "\t" << H_integrals[1] << "\t" << H_integrals[2] << "\t"
                           << H_integrals[0] + H_integrals[2] << "\t" << 0.0 << std::endl;
            }
        }

        // Main time step loop.
        double loop_time_end = navier_stokes_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && navier_stokes_integrator->stepsRemaining())
        {
            iteration_num = navier_stokes_integrator->getIntegratorStep();
            loop_time = navier_stokes_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = navier_stokes_integrator->getMaximumTimeStepSize();
            pout << "Advancing hierarchy with timestep size dt = " << dt << "\n";
            navier_stokes_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            std::vector<double> H_integrals =
                LevelSetUtilities::computeHeavisideIntegrals3PhaseFlows(level_set_fixer.getLevelSetContainer());
            if (SAMRAI_MPI::getRank() == 0)
            {
                vol_stream.precision(16);
                vol_stream.setf(ios::fixed, ios::floatfield);
                vol_stream << loop_time << "\t" << H_integrals[0] << "\t" << H_integrals[1] << "\t" << H_integrals[2]
                           << "\t" << H_integrals[0] + H_integrals[2] << "\t" << level_set_fixer.getLagrangeMultiplier()
                           << std::endl;
            }

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !navier_stokes_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "Writing visualization files...\n\n";
                navier_stokes_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "Writing restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "Writing timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            if (SAMRAI_MPI::getRank() == 0)
            {
                const Eigen::Vector3d& rbd_posn = bp_rbd->getCurrentCOMPosn();
                const Eigen::Vector3d& rbd_trans_vel = bp_rbd->getCurrentCOMTransVelocity();

                Eigen::Vector3d hydro_force_pressure, hydro_force_viscous, hydro_torque_pressure, hydro_torque_viscous;
                bp_rbd->getHydrodynamicForceTorque(
                    hydro_force_pressure, hydro_force_viscous, hydro_torque_pressure, hydro_torque_viscous);

                rbd_stream.precision(12);
                rbd_stream.setf(ios::fixed, ios::floatfield);
                rbd_stream << loop_time << "\t" << rbd_posn[1] << "\t" << rbd_trans_vel[1] << std::endl;

                ft_stream.precision(12);
                ft_stream.setf(ios::fixed, ios::floatfield);
                ft_stream << loop_time << "\t" << hydro_force_pressure[0] << "\t" << hydro_force_pressure[1] << "\t"
                          << hydro_force_pressure[2] << "\t" << hydro_force_viscous[0] << "\t" << hydro_force_viscous[1]
                          << "\t" << hydro_force_viscous[2] << "\t" << hydro_torque_pressure[0] << "\t"
                          << hydro_torque_pressure[1] << "\t" << hydro_torque_pressure[2] << "\t"
                          << hydro_torque_viscous[0] << "\t" << hydro_torque_viscous[1] << "\t"
                          << hydro_torque_viscous[2] << std::endl;
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            vol_stream.close();
            rbd_stream.close();
            ft_stream.close();
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main
