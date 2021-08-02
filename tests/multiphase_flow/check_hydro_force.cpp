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

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBHydrodynamicSurfaceForceEvaluator.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/RelaxationLSMethod.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application specific includes.
#include "LSLocateGasInterface.h"
#include "SetFluidGasSolidDensity.h"
#include "SetFluidGasSolidViscosity.h"
#include "SetLSProperties.h"

int coarsest_ln, max_finest_ln;
double dx, ds;

// Struct to maintain the properties of the circular interface
struct CircularInterface
{
    Eigen::Vector3d X0;
    double R;
};
CircularInterface circle;

void
calculate_distance_analytically(Pointer<PatchHierarchy<NDIM> > patch_hierarchy, int E_idx)
{
    int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > E_data = patch->getPatchData(E_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                const double distance =
                    std::sqrt(std::pow((coord[0] - circle.X0(0)), 2.0) + std::pow((coord[1] - circle.X0(1)), 2.0)
#if (NDIM == 3)
                              + std::pow((coord[2] - circle.X0(2)), 2.0)
#endif
                    );

                (*E_data)(ci) = distance - circle.R;
            }
        }
    }

    return;
} // calculate_distance_analytically

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

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IBLevelSet.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Setup solid information
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator;
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

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "PREDICTOR_CORRECTOR");
        if (adv_diff_solver_type == "PREDICTOR_CORRECTOR")
        {
            Pointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps(
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
        level_set_gas_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateGasInterfaceCallbackFunction,
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
        SetLSProperties* ptr_setSetLSProperties = new SetLSProperties("SetLSProperties", nullptr, level_set_gas_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var_gas, &callSetGasLSCallbackFunction, static_cast<void*>(ptr_setSetLSProperties));

        // Setup the advected and diffused fluid quantities.
        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<hier::Variable<NDIM> > rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariable<NDIM, double>("rho");
        }
        else
        {
            rho_var = new CellVariable<NDIM, double>("rho");
        }
        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const int ls_reinit_interval = input_db->getInteger("LS_REINIT_INTERVAL");
        const double rho_fluid = input_db->getDouble("RHO_F");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const int num_solid_interface_cells = input_db->getDouble("NUM_SOLID_INTERFACE_CELLS");
        const int num_gas_interface_cells = input_db->getDouble("NUM_GAS_INTERFACE_CELLS");
        SetFluidGasSolidDensity* ptr_setFluidGasSolidDensity = new SetFluidGasSolidDensity("SetFluidGasSolidDensity",
                                                                                           adv_diff_integrator,
                                                                                           phi_var_solid,
                                                                                           phi_var_gas,
                                                                                           rho_fluid,
                                                                                           rho_gas,
                                                                                           rho_solid,
                                                                                           ls_reinit_interval,
                                                                                           num_solid_interface_cells,
                                                                                           num_gas_interface_cells);
        navier_stokes_integrator->registerResetFluidDensityFcn(&callSetFluidGasSolidDensityCallbackFunction,
                                                               static_cast<void*>(ptr_setFluidGasSolidDensity));

        const double mu_fluid = input_db->getDouble("MU_F");
        const double mu_gas = input_db->getDouble("MU_G");
        const double mu_solid = input_db->getDoubleWithDefault("MU_S", std::numeric_limits<double>::quiet_NaN());
        const bool set_mu_solid = input_db->getBool("SET_MU_S");
        SetFluidGasSolidViscosity* ptr_setFluidGasSolidViscosity =
            new SetFluidGasSolidViscosity("SetFluidGasSolidViscosity",
                                          adv_diff_integrator,
                                          phi_var_solid,
                                          phi_var_gas,
                                          mu_fluid,
                                          mu_gas,
                                          mu_solid,
                                          ls_reinit_interval,
                                          num_solid_interface_cells,
                                          num_gas_interface_cells,
                                          set_mu_solid);
        navier_stokes_integrator->registerResetFluidViscosityFcn(&callSetFluidGasSolidViscosityCallbackFunction,
                                                                 static_cast<void*>(ptr_setFluidGasSolidViscosity));

        // Create Eulerian initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

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
        adv_diff_integrator->setPhysicalBcCoef(phi_var_gas, phi_bc_coef);
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef);

        // LS reinit boundary conditions, which is set to be the same as the BCs
        // for advection
        RobinBcCoefStrategy<NDIM>* ls_reinit_bcs = phi_bc_coef;
        level_set_gas_ops->registerPhysicalBoundaryCondition(ls_reinit_bcs);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        navier_stokes_integrator->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        navier_stokes_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Compute distance function.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double> > d_var = new CellVariable<NDIM, double>("distance", 1);
        const IntVector<NDIM> no_width = 0;
        Pointer<VariableContext> main_ctx = var_db->getContext("Main");
        const int d_idx = var_db->registerVariableAndContext(d_var, main_ctx, no_width);
        int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(d_idx, 0.0);
        }
        calculate_distance_analytically(patch_hierarchy, d_idx);

        // Copy the distance function into level set of solid.
        int phi_solid_idx =
            var_db->mapVariableAndContextToIndex(phi_var_solid, adv_diff_integrator->getCurrentContext());
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, hier_finest_ln);
        hier_cc_data_ops.copyData(phi_solid_idx, d_idx);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Write out initial visualization data.
        int iteration_num = navier_stokes_integrator->getIntegratorStep();
        double loop_time = navier_stokes_integrator->getIntegratorTime();

        plog << "\n\nWriting visualization files...\n\n";
        navier_stokes_integrator->setupPlotData();
        visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);

        // Set velocity and pressure field and calculate hydrodynamic forces and torques
        // Although INS integrator already sets them the initial projection can changes it.
        int u_ins_idx = var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                             navier_stokes_integrator->getCurrentContext());
        u_init->setDataOnPatchHierarchy(
            u_ins_idx, navier_stokes_integrator->getVelocityVariable(), patch_hierarchy, 0.0);
        int p_ins_idx = var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                             navier_stokes_integrator->getCurrentContext());
        p_init->setDataOnPatchHierarchy(
            p_ins_idx, navier_stokes_integrator->getPressureVariable(), patch_hierarchy, 0.0);

        IBHydrodynamicSurfaceForceEvaluator force_evaluator("force_evaluator",
                                                            phi_var_solid,
                                                            adv_diff_integrator,
                                                            navier_stokes_integrator,
                                                            Pointer<Database>(nullptr));

        IBTK::Vector3d pressure_force, viscous_force, pressure_torque, viscous_torque;
        force_evaluator.computeHydrodynamicForceTorque(
            pressure_force, viscous_force, pressure_torque, viscous_torque, circle.X0);

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "Vertical pressure force analytical = " << M_PI * circle.R * circle.R << std::endl;
            out << "Vertical pressure force numerical  = " << pressure_force[1] << std::endl;

            out << "Horizontal pressure force analytical = " << 0.0 << std::endl;
            out << "Horizontal pressure force numerical  = " << pressure_force[0] << std::endl;
            out << "Pressure torque (should be zero): " << pressure_torque << std::endl;
        }

        // Delete dumb pointers.
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete ptr_setFluidGasSolidDensity;
        delete ptr_setFluidGasSolidViscosity;
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown

    return 0;
} // main
