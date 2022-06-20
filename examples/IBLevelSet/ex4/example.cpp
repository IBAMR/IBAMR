// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2022 by the IBAMR developers
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

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <Eigen/Geometry>

#include <ibamr/app_namespaces.h>

// Applic<ation specific includes.
#include "LevelSetInitialCondition.h"
#include "SetFluidSolidDensity.h"
#include "SetFluidSolidViscosity.h"
#include "TagLSRefinementCells.h"

// Declaring the sign of v and returning it.
inline int
sgn(double v)
{
    return ((v < 0) ? -1 : (v > 0) ? 1 : 0);
}

FoilInterface foilA;

// Struct to reset solid level set
struct SolidLevelSetResetter
{
    SolidLevelSetResetter(Pointer<AdvDiffHierarchyIntegrator> integrator,
                          Pointer<CellVariable<NDIM, double> > var,
                          Pointer<BrinkmanPenalizationRigidBodyDynamics> bp,
                          FoilInterface* foil)
        : adv_diff_integrator(integrator), ls_solid_var(var), bp_rbd(bp), ptr_foil(foil)
    {
        return;
    }

    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
    Pointer<CellVariable<NDIM, double> > ls_solid_var;
    Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd;
    FoilInterface* ptr_foil;
};

void
reset_solid_level_set_callback_fcn(double current_time, double new_time, int /*cycle_num*/, void* ctx)
{
    SolidLevelSetResetter* resetter = static_cast<SolidLevelSetResetter*>(ctx);
    const FoilInterface& foil = *(resetter->ptr_foil);

    // Get the new centroid of the body
    const double dt = new_time - current_time;
    const Eigen::Vector3d XCOM_current = resetter->bp_rbd->getCurrentCOMPosn();
    const Eigen::Vector3d XCOM_new = XCOM_current + dt * (resetter->bp_rbd->getNewCOMTransVelocity());

    // b) Rotational matrix.
    const double theta = foil.theta_0 * std::sin(2 * M_PI * foil.freq * new_time);

    const Eigen::Vector3d rot_axis(0.0, 0.0, 1.0);
    Eigen::Quaterniond q(Eigen::AngleAxisd(theta, rot_axis));
    q.normalize();

    const Eigen::Matrix3d R_mat = q.toRotationMatrix();

    const double& R = foil.R;
    const Eigen::Vector3d& X0 = foil.X0;
    const Eigen::Vector3d& X1 = foil.X1;
    const Eigen::Vector3d& X2 = foil.X2;
    const Eigen::Vector3d& X3 = foil.X3;
    Eigen::Vector3d check1, check2, check3, check4, check5, check6;

    // Rotate and translate the top, bottom and circular surface to get new coordiantes.
    const Eigen::Vector3d X1_new = R_mat * (X1 - X0) + XCOM_new;
    const Eigen::Vector3d X2_new = R_mat * (X2 - X0) + XCOM_new;
    const Eigen::Vector3d X3_new = R_mat * (X3 - X0) + XCOM_new;

    const Eigen::Vector3d X_T_new = (X1_new + X2_new + X3_new) / 3.0;

    const double slope1 = (X3_new[1] - X1_new[1]) / (X3_new[0] - X1_new[0]);
    const double slope2 = (X3_new[1] - X2_new[1]) / (X3_new[0] - X2_new[0]);
    const double slope3 = (X2_new[0] - X1_new[0]) / (X2_new[1] - X1_new[1]);

    const double y_intercept1 = X1_new[1] - slope1 * X1_new[0];
    const double y_intercept2 = X2_new[1] - slope2 * X2_new[0];
    const double x_intercept3 = X1_new[0] - slope3 * X1_new[1];

    double distance1[2], distance2[3]; // Foil has three surfaces and 1 surface for circle.

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
                Eigen::Vector3d coord = Eigen::Vector3d::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                // Distance from the semi-circle.
                distance1[0] =
                    std::sqrt(std::pow((coord[0] - XCOM_new[0]), 2.0) + std::pow((coord[1] - XCOM_new[1]), 2.0)) - R;

                // Distance from top triangle surface.
                distance2[0] = std::abs(coord[1] - slope1 * coord[0] - y_intercept1) / std::sqrt(1.0 + slope1 * slope1);

                check1 = (X1_new - X3_new).cross(X_T_new - X3_new);
                check2 = (X1_new - X3_new).cross(coord - X3_new);

                distance2[0] *= (-sgn(check1[2]) * sgn(check2[2]));

                // Distance from bottom triangle surface.
                distance2[1] = std::abs(coord[1] - slope2 * coord[0] - y_intercept2) / std::sqrt(1.0 + slope2 * slope2);

                check3 = (X2_new - X3_new).cross(X_T_new - X3_new);
                check4 = (X2_new - X3_new).cross(coord - X3_new);

                distance2[1] *= (-sgn(check3[2]) * sgn(check4[2]));

                // Distance from base of triangle.
                distance2[2] = std::abs(coord[0] - slope3 * coord[1] - x_intercept3) / std::sqrt(1.0 + slope3 * slope3);

                check5 = (X2_new - X1_new).cross(X_T_new - X1_new);
                check6 = (X2_new - X1_new).cross(coord - X1_new);

                distance2[2] *= (-sgn(check5[2]) * sgn(check6[2]));

                distance1[1] =
                    std::max({ distance2[0], distance2[1], distance2[2] }); // intersection of three traingle surfaces

                (*ls_solid_data)(ci) = std::min({ distance1[0], distance1[1] }); // union of circle and triangle
            }
        }
    }

    return;
}

void
imposed_kinematics(double data_time, int /*cycle_num*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* ctx)
{
    const FoilInterface& foil = *(static_cast<FoilInterface*>(ctx));
    const double& theta_0 = foil.theta_0;
    const double& freq = foil.freq;

    U_com.setZero();
    W_com.setZero();
    W_com[2] = theta_0 * 2 * M_PI * freq * std::cos(2 * M_PI * freq * data_time);

    return;
} // imposed_kinematics

void
external_force_torque(double /*data_time*/, int /*cycle_num*/, Eigen::Vector3d& F, Eigen::Vector3d& T, void* /*ctx*/)
{
    F.setZero();
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
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

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
        foilA.R = input_db->getDouble("R");
        foilA.mass = input_db->getDouble("MASS");

        foilA.X0[0] = input_db->getDouble("X_COM");
        foilA.X0[1] = input_db->getDouble("Y_COM");
        foilA.X1[0] = input_db->getDouble("X1");
        foilA.X1[1] = input_db->getDouble("Y1");
        foilA.X2[0] = input_db->getDouble("X2");
        foilA.X2[1] = input_db->getDouble("Y2");
        foilA.X3[0] = input_db->getDouble("X3");
        foilA.X3[1] = input_db->getDouble("Y3");

        foilA.freq = input_db->getDouble("FREQ");
        foilA.theta_0 = input_db->getDouble("THETA_0");

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

        // Create level sets for solid interfaces.
        const string& ls_foilA = "level_set_foilA";
        Pointer<CellVariable<NDIM, double> > phi_var_foilA = new CellVariable<NDIM, double>(ls_foilA);

        // Register the level sets with advection diffusion integrator.
        adv_diff_integrator->registerTransportedQuantity(phi_var_foilA);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_foilA, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_foilA,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Solid level set initial condition
        Pointer<CartGridFunction> phi_foilA_init = new LevelSetInitialCondition("phi_foilA_init", foilA);
        adv_diff_integrator->setInitialConditions(phi_var_foilA, phi_foilA_init);

        // Solid level set resetting consition
        SolidLevelSetResetter foilA_level_set_resetter(adv_diff_integrator, phi_var_foilA, /*bp_rbd*/ nullptr, &foilA);
        adv_diff_integrator->registerIntegrateHierarchyCallback(&reset_solid_level_set_callback_fcn,
                                                                static_cast<void*>(&foilA_level_set_resetter));

        // Fluid density and viscosity.
        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<hier::Variable<NDIM> > rho_var;
        rho_var = new SideVariable<NDIM, double>("rho");
        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_fluid = input_db->getDouble("RHO_F");
        SetFluidSolidDensity* ptr_setFluidSolidDensity = new SetFluidSolidDensity("SetFluidSolidDensity", rho_fluid);
        navier_stokes_integrator->registerResetFluidDensityFcn(&callSetFluidSolidDensityCallbackFunction,
                                                               static_cast<void*>(ptr_setFluidSolidDensity));

        const double mu_fluid = input_db->getDouble("MU_F");
        SetFluidSolidViscosity* ptr_setFluidSolidViscosity =
            new SetFluidSolidViscosity("SetFluidSolidViscosity", mu_fluid);
        navier_stokes_integrator->registerResetFluidViscosityFcn(&callSetFluidSolidViscosityCallbackFunction,
                                                                 static_cast<void*>(ptr_setFluidSolidViscosity));

        // Register callback function for tagging refined cells for level set data
        const double tag_value = input_db->getDouble("LS_TAG_VALUE");
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        TagLSRefinementCells ls_foilA_tagger(adv_diff_integrator, phi_var_foilA, tag_value, tag_thresh);
        navier_stokes_integrator->registerApplyGradientDetectorCallback(&callTagSolidLSRefinementCellsCallbackFunction,
                                                                        static_cast<void*>(&ls_foilA_tagger));

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

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_foilA, phi_bc_coef);

        // Configure the Brinkman penalization object to do the rigid body dynamics.
        Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd_A =
            new BrinkmanPenalizationRigidBodyDynamics("Airfoil_A",
                                                      phi_var_foilA,
                                                      adv_diff_integrator,
                                                      navier_stokes_integrator,
                                                      app_initializer->getComponentDatabase("BrinkmanPenalization"),
                                                      /*register_for_restart*/ true);

        {
            FreeRigidDOFVector free_dofs;
            free_dofs << 1, 1, 0;
            Eigen::Vector3d U_i = Eigen::Vector3d::Zero();
            Eigen::Vector3d W_i(0.0, 0.0, 2 * M_PI * foilA.freq * foilA.theta_0);
            bp_rbd_A->setSolveRigidBodyVelocity(free_dofs);
            bp_rbd_A->registerKinematicsFunction(&imposed_kinematics, &foilA);
            bp_rbd_A->registerExternalForceTorqueFunction(&external_force_torque);
            bp_rbd_A->setInitialConditions(foilA.X0, U_i, W_i, foilA.mass);
            navier_stokes_integrator->registerBrinkmanPenalizationStrategy(bp_rbd_A);
            foilA_level_set_resetter.bp_rbd = bp_rbd_A;
        }

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

        // Open streams to save position and velocity of the structure.
        ofstream rbd_foilA_stream, ft_foilA_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            rbd_foilA_stream.open("rbd.foilA", ios_base::out | ios_base::app);
            ft_foilA_stream.open("hydro_force_torque.foilA", ios_base::out | std::ios_base::app);
        }

        // Main time step loop.
        double loop_time_end = navier_stokes_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && navier_stokes_integrator->stepsRemaining())
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
                {
                    const Eigen::Vector3d& rbd_posn = bp_rbd_A->getCurrentCOMPosn();
                    const Eigen::Vector3d& rbd_trans_vel = bp_rbd_A->getCurrentCOMTransVelocity();

                    rbd_foilA_stream.precision(12);
                    rbd_foilA_stream.setf(ios::fixed, ios::floatfield);
                    rbd_foilA_stream << loop_time << "\t" << rbd_posn[0] << "\t" << rbd_posn[1] << "\t"
                                     << rbd_trans_vel[0] << "\t" << rbd_trans_vel[1] << std::endl;

                    Eigen::Vector3d hydro_force_pressure, hydro_force_viscous, hydro_torque_pressure,
                        hydro_torque_viscous;
                    bp_rbd_A->getHydrodynamicForceTorque(
                        hydro_force_pressure, hydro_force_viscous, hydro_torque_pressure, hydro_torque_viscous);
                    ft_foilA_stream.precision(12);
                    ft_foilA_stream.setf(ios::fixed, ios::floatfield);
                    ft_foilA_stream << loop_time << "\t" << hydro_force_viscous[0] << "\t" << hydro_force_viscous[1]
                                    << "\t" << hydro_force_pressure[0] << "\t" << hydro_force_pressure[1] << std::endl;
                }
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            rbd_foilA_stream.close();
            ft_foilA_stream.close();
        }

        // Delete dumb pointers.
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete ptr_setFluidSolidDensity;
        delete ptr_setFluidSolidViscosity;
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main
