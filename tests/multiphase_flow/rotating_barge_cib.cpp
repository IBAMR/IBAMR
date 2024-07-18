// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/ConstraintIBMethod.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateBargeInterface.h"
#include "LSLocateGasInterface.h"
#include "RigidBodyKinematics.h"

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

    // Several parts of the code (such as LDataManager) expect mesh files,
    // specified in the input file, to exist in the current working
    // directory. Since tests are run in temporary directories we need to
    // regenerate these input files first.
    //
    // The following is simply the bargeGen executable:
    if (IBTK_MPI::getRank() == 0)
    {
        const double Lx = 5.0;
        const double Ly = 2.5;
        const int Nx = 500 * 2 * 2;
        const int Ny = 250 * 2 * 2;

        // Dimensional parameters
        const double dx = Lx / Nx;
        const double dy = Ly / Ny;
        const double Length = 0.3;
        const double Width = Length / 3.0;
        const double Xcom = 0.0;
        const double Ycom = 0.0;
        const double theta = 15;             // in degrees
        const double t = theta * M_PI / 180; // in radians

        const int NumPtsX = std::ceil(Length / dx) + 1;
        const int NumPtsY = std::ceil(Width / dy) + 1;

        std::ofstream barge_coord_stream("barge2d.vertex");
        barge_coord_stream << NumPtsX * NumPtsY << "\n";
        // Generate Rectangular Prism
        for (int j = 1; j <= NumPtsY; ++j)
        {
            const double y = Ycom - Width / 2 + (j - 1) * dy;
            for (int i = 1; i <= NumPtsX; ++i)
            {
                const double x = Xcom - Length / 2 + (i - 1) * dy;
                const double r_x = x * std::cos(t) - y * std::sin(t);
                const double r_y = x * std::sin(t) + y * std::cos(t);
                barge_coord_stream << std::scientific << std::setprecision(7) << std::setfill('0') << r_x << "\t" << r_y
                                   << "\n";
            }
        }
    }

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "IB.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        SAMRAIPointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        if (NDIM != 2)
        {
            TBOX_ERROR("This case is presently not implemented for NDIM = 3");
        }

        // Set up the advection diffusion hierarchy integrator
        SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "PREDICTOR_CORRECTOR");
        if (adv_diff_solver_type == "PREDICTOR_CORRECTOR")
        {
            auto predictor = make_samrai_shared<AdvectorExplicitPredictorPatchOps>(
                "AdvectorExplicitPredictorPatchOps",
                app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
            adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator(
                "AdvDiffPredictorCorrectorHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"),
                predictor);
        }
        else if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: PREDICTOR_CORRECTOR, SEMI_IMPLICIT");
        }
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        const int num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        auto ib_method_ops = make_samrai_shared<ConstraintIBMethod>(
            "ConstraintIBMethod", app_initializer->getComponentDatabase("ConstraintIBMethod"), num_structures);
        SAMRAIPointer<IBHierarchyIntegrator> time_integrator = make_samrai_shared<IBExplicitHierarchyIntegrator>(
            "IBHierarchyIntegrator",
            app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
            ib_method_ops,
            navier_stokes_integrator);

        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);

        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize",
            time_integrator,
            app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
                                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                    error_detector,
                                                    box_generator,
                                                    load_balancer);

        // Configure the IB solver.
        auto ib_initializer = make_samrai_shared<IBStandardInitializer>(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        auto ib_force_fcn = make_samrai_shared<IBStandardForceGen>();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Setup level set information
        BargeInterface barge;
        barge.length = input_db->getDouble("BARGE_LENGTH");
        barge.width = input_db->getDouble("BARGE_WIDTH");

        {
            const std::vector<std::vector<double> >& structure_COM = ib_method_ops->getCurrentStructureCOM();
            barge.COM(0) = structure_COM[0][0];
            barge.COM(1) = structure_COM[0][1];
#if (NDIM == 3)
            barge.COM(2) = structure_COM[0][2];
#endif
        }
        const double fluid_height = input_db->getDouble("GAS_LS_INIT");

        const string& ls_name_solid = "level_set_solid";
        const double vol_elem = input_db->getDoubleWithDefault("VOL_ELEM", -1.0);
        SAMRAIPointer<CellVariableNd<double> > phi_var_solid =
            make_samrai_shared<CellVariableNd<double> >(ls_name_solid);
        auto level_set_solid_ops = make_samrai_shared<RelaxationLSMethod>(
            ls_name_solid, app_initializer->getComponentDatabase("LevelSet_Solid"));
        LSLocateBargeInterface setLSLocateBargeInterface("LSLocateBargeInterface",
                                                         adv_diff_integrator,
                                                         phi_var_solid,
                                                         ib_method_ops->getLDataManager(),
                                                         vol_elem,
                                                         &barge);
        level_set_solid_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateBargeInterfaceCallbackFunction,
                                                                      static_cast<void*>(&setLSLocateBargeInterface));

        const string& ls_name_gas = "level_set_gas";
        SAMRAIPointer<CellVariableNd<double> > phi_var_gas = make_samrai_shared<CellVariableNd<double> >(ls_name_gas);
        auto level_set_gas_ops =
            make_samrai_shared<RelaxationLSMethod>(ls_name_gas, app_initializer->getComponentDatabase("LevelSet_Gas"));
        LSLocateGasInterface setLSLocateGasInterface(
            "LSLocateGasInterface", adv_diff_integrator, phi_var_gas, fluid_height);
        level_set_gas_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateGasInterfaceCallbackFunction,
                                                                    static_cast<void*>(&setLSLocateGasInterface));

        adv_diff_integrator->registerTransportedQuantity(phi_var_solid);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_solid, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_solid,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        adv_diff_integrator->registerTransportedQuantity(phi_var_gas);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_gas, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_gas,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Register the reinitialization functions for the level set variables
        IBAMR::LevelSetUtilities::SetLSProperties setSetLSPropertiesSolid("SetLSPropertiesSolid", level_set_solid_ops);
        adv_diff_integrator->registerResetFunction(phi_var_solid,
                                                   &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy,
                                                   static_cast<void*>(&setSetLSPropertiesSolid));

        IBAMR::LevelSetUtilities::SetLSProperties setSetLSPropertiesGas("SetLSPropertiesGas", level_set_gas_ops);
        adv_diff_integrator->registerResetFunction(phi_var_gas,
                                                   &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy,
                                                   static_cast<void*>(&setSetLSPropertiesGas));

        // Setup the advected and diffused quantities.
        SAMRAIPointer<VariableNd> rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariableNd<double>("rho");
        }
        else
        {
            rho_var = new CellVariableNd<double>("rho");
        }

        navier_stokes_integrator->registerMassDensityVariable(rho_var);

        SAMRAIPointer<CellVariableNd<double> > mu_var = make_samrai_shared<CellVariableNd<double> >("mu");
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_fluid = input_db->getDouble("RHO_F");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const double mu_fluid = input_db->getDouble("MU_F");
        const double mu_gas = input_db->getDouble("MU_G");
        const double mu_solid = input_db->getDoubleWithDefault("MU_S", std::numeric_limits<double>::quiet_NaN());
        const bool set_mu_solid = input_db->getBool("SET_MU_S");
        const int num_solid_interface_cells = input_db->getDouble("NUM_SOLID_INTERFACE_CELLS");
        const int num_gas_interface_cells = input_db->getDouble("NUM_GAS_INTERFACE_CELLS");

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
        const double tag_min_value = -tag_thresh;
        const double tag_max_value = tag_thresh;
        IBAMR::LevelSetUtilities::TagLSRefinementCells ls_tagger(
            adv_diff_integrator, phi_var_gas, tag_min_value, tag_max_value);
        time_integrator->registerApplyGradientDetectorCallback(&IBAMR::LevelSetUtilities::tagLSCells,
                                                               static_cast<void*>(&ls_tagger));

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> u_init = make_samrai_shared<muParserCartGridFunction>(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> p_init = make_samrai_shared<muParserCartGridFunction>(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVectorNd& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategyNd*> u_bc_coefs(NDIM);
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

        RobinBcCoefStrategyNd* rho_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategyNd* mu_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategyNd* phi_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_gas, phi_bc_coef);
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef);

        // LS reinit boundary conditions, which is set to be the same as the BCs
        // for advection
        level_set_solid_ops->registerPhysicalBoundaryCondition(phi_bc_coef);
        level_set_gas_ops->registerPhysicalBoundaryCondition(phi_bc_coef);

        // Initialize objects
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        SAMRAIPointer<CartGridFunction> grav_force = make_samrai_shared<IBAMR::VCINSUtilities::GravityForcing>(
            "GravityForcing", navier_stokes_integrator, grav_const);

        auto surface_tension_force = make_samrai_shared<SurfaceTensionForceFunction>(
            "SurfaceTensionForceFunction",
            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
            adv_diff_integrator,
            phi_var_gas);

        auto eul_forces = make_samrai_shared<CartGridFunctionSet>("eulerian_forces");
        eul_forces->addFunction(grav_force);
        eul_forces->addFunction(surface_tension_force);
        time_integrator->registerBodyForceFunction(eul_forces);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Create ConstraintIBKinematics objects
        vector<SAMRAIPointer<ConstraintIBKinematics> > ibkinematics_ops_vec;
        SAMRAIPointer<ConstraintIBKinematics> ib_kinematics_op;
        // struct_0
        const string& object_name = "Barge";
        ib_kinematics_op = new RigidBodyKinematics(
            object_name,
            app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase(object_name),
            ib_method_ops->getLDataManager(),
            patch_hierarchy);
        ibkinematics_ops_vec.push_back(ib_kinematics_op);

        // Register ConstraintIBKinematics, physical boundary operators and
        // other things with ConstraintIBMethod.
        ib_method_ops->registerConstraintIBKinematics(ibkinematics_ops_vec);
        if (vol_elem > 0.0) ib_method_ops->setVolumeElement(vol_elem, 0);
        ib_method_ops->setVelocityPhysBdryOp(time_integrator->getVelocityPhysBdryOp());
        ib_method_ops->initializeHierarchyOperatorsandData();

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        // File to write for barge angle.
        ofstream output_file;
        if (!IBTK_MPI::getRank()) output_file.open("output");
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
            pout << "Advancing hierarchy with timestep size dt = " << dt << "\n";
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // Update the position of barge reference point.
            const std::vector<std::vector<double> >& structure_COM = ib_method_ops->getCurrentStructureCOM();
            barge.COM(0) = structure_COM[0][0];
            barge.COM(1) = structure_COM[0][1];
#if (NDIM == 3)
            barge.COM(2) = structure_COM[0][2];
#endif

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Write to file
            if (!IBTK_MPI::getRank())
            {
                output_file << std::setprecision(10) << loop_time << '\t' << barge.theta << std::endl;
            }

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
        }

        // Close file
        if (!IBTK_MPI::getRank()) output_file.close();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown

    return 0;
} // main
