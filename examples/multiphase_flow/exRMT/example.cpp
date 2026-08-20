// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config filess
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
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateColumnInterface.h"
#include "LevelSetInitialCondition.h"

/*************** */
#include <ibamr/AdvDiffPPMConvectiveOperator.h>
#include <ibamr/ConvectiveOperator.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

#include "Box.h"
#include "CellData.h"
#include "Index.h"

#include <CartesianPatchGeometry.h>
#include <FaceVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>
/**************** */

#include <ibamr/RMTForceFunction.h>
#include <ibamr/RMTMethod.h>

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

/********************************/

void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSVCStaggeredHierarchyIntegrator> time_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

// Pointer<CellVariable<NDIM, double> > w_var = new CellVariable<NDIM, double>("W_FROM_XI");
// int w_idx = IBTK::invalid_index;

// Pointer<CellVariable<NDIM, double> > phiw_err_var = new CellVariable<NDIM, double>("PHI_MINUS_W");
// int phiw_err_idx = IBTK::invalid_index;

static double
fillDispMagFromXi(Pointer<PatchHierarchy<NDIM> > hierarchy,
                  int disp_mag_idx,
                  int phi_cur_idx,
                  int xi0_cur_idx,
                  int xi1_cur_idx)
{
    double dispmax_local = 0.0;

    const int finest_ln = hierarchy->getFinestLevelNumber();

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();

            auto geom = Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
            const double* dx = geom->getDx();
            const double* x_lo = geom->getXLower();
            const hier::Index<NDIM>& ilo = box.lower();

            auto phi = Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(phi_cur_idx));
            auto xi0 = Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(xi0_cur_idx));
            auto xi1 = Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(xi1_cur_idx));
            auto disp = Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(disp_mag_idx));

            if (phi.isNull() || xi0.isNull() || xi1.isNull() || disp.isNull()) continue;

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM> ci(it());
                const double phiv = (*phi)(ci);

                if (phiv < 0.0)
                {
                    const double x = x_lo[0] + dx[0] * (double(ci(0) - ilo(0)) + 0.5);
                    const double y = x_lo[1] + dx[1] * (double(ci(1) - ilo(1)) + 0.5);

                    const double disp_x = x - (*xi0)(ci);
                    const double disp_y = y - (*xi1)(ci);
                    const double disp_mag = std::sqrt(disp_x * disp_x + disp_y * disp_y);

                    (*disp)(ci) = disp_mag;
                    dispmax_local = std::max(dispmax_local, disp_mag);
                }
                else
                {
                    (*disp)(ci) = 0.0;
                }
            }
        }
    }

    return IBTK_MPI::maxReduction(dispmax_local);
}

static double
computeVelMaxInSolid(Pointer<PatchHierarchy<NDIM> > hierarchy, int u_idx, int phi_idx)
{
    double velmax_local = 0.0;

    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();

            auto geom = Pointer<SAMRAI::geom::CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());
            const double* dx = geom->getDx();

            auto phi = Pointer<SAMRAI::pdat::CellData<NDIM, double> >(patch->getPatchData(phi_idx));
            auto u = Pointer<SAMRAI::pdat::SideData<NDIM, double> >(patch->getPatchData(u_idx));
            if (phi.isNull() || u.isNull()) continue;

            const double margin = 1.0 * std::max(dx[0], dx[1]);

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                const SAMRAI::pdat::CellIndex<NDIM> ci(it());

                if ((*phi)(ci) < -margin)
                {
                    SAMRAI::pdat::SideIndex<NDIM> sx(ci, 0, 0);
                    SAMRAI::pdat::SideIndex<NDIM> sy(ci, 1, 0);

                    const double vx = (*u)(sx, 0);
                    const double vy = (*u)(sy, 0);
                    const double vmag = std::sqrt(vx * vx + vy * vy);

                    velmax_local = std::max(velmax_local, vmag);
                }
            }
        }
    }

    return IBTK_MPI::maxReduction(velmax_local);
}

static void
computeAndWriteCentroid(Pointer<PatchHierarchy<NDIM> > hierarchy, int phi_idx, double time, const std::string& filename)
{
    double m_local = 0.0;
    double mx_local = 0.0;
    double my_local = 0.0;

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& box = patch->getBox();

            auto phi = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_idx));

            if (phi.isNull()) continue;

            auto geom = Pointer<CartesianPatchGeometry<NDIM> >(patch->getPatchGeometry());

            const double* dx = geom->getDx();
            const double* x_lo = geom->getXLower();
            const hier::Index<NDIM>& ilo = box.lower();

            const double cell_vol = dx[0] * dx[1];

            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                if ((*phi)(ci) < 0.0)
                {
                    const double x = x_lo[0] + dx[0] * (double(ci(0) - ilo(0)) + 0.5);

                    const double y = x_lo[1] + dx[1] * (double(ci(1) - ilo(1)) + 0.5);

                    m_local += cell_vol;
                    mx_local += x * cell_vol;
                    my_local += y * cell_vol;
                }
            }
        }
    }

    const double m = IBTK_MPI::sumReduction(m_local);
    const double mx = IBTK_MPI::sumReduction(mx_local);
    const double my = IBTK_MPI::sumReduction(my_local);

    if (m <= 0.0) return;

    const double xc = mx / m;
    const double yc = my / m;

    if (IBTK_MPI::getRank() == 0)
    {
        std::ofstream out(filename.c_str(), std::ios::app);
        out << std::setprecision(16) << time << " " << xc << " " << yc << "\n";
    }
}
/*********************** */

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

        // remove
        int phi_cur_idx = IBTK::invalid_index;
        int xi0_cur_idx = IBTK::invalid_index;
        int xi1_cur_idx = IBTK::invalid_index;
        int phixi_idx = IBTK::invalid_index;
        int dphi_idx = IBTK::invalid_index;

        int phi_new_idx = IBTK::invalid_index;
        int xi0_new_idx = IBTK::invalid_index;
        int xi1_new_idx = IBTK::invalid_index;

        int disp_mag_idx = IBTK::invalid_index;

        // int phi0_idx = IBTK::invalid_index;

        // Pointer<CellVariable<NDIM, double> > phi0_var = new CellVariable<NDIM, double>("PHI0");

        Pointer<CellVariable<NDIM, double> > disp_mag_var = new CellVariable<NDIM, double>("DISP_MAG");
        //

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
        // ColumnInterface column;
        // input_db->getDoubleArray("X_UR", &column.X_UR[0], NDIM);

        const string& ls_name = "level_set";
        Pointer<CellVariable<NDIM, double> > phi_var = new CellVariable<NDIM, double>(ls_name);

        /************************* */
        // IBAMR
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        // Pointer<VariableContext> phi0_ctx = var_db->getContext("PHI0_CONTEXT");
        // phi0_idx = var_db->registerVariableAndContext(phi0_var, phi0_ctx, IntVector<NDIM>(1));

        Pointer<CellVariable<NDIM, double> > xi0_var = new CellVariable<NDIM, double>("xi0");
        Pointer<CellVariable<NDIM, double> > xi1_var = new CellVariable<NDIM, double>("xi1");

        Pointer<CellVariable<NDIM, double> > phixi_var = new CellVariable<NDIM, double>("PHI_FROM_XI");

        Pointer<CellVariable<NDIM, double> > dphi_var = new CellVariable<NDIM, double>("DELTA_PHI");
        /************************* */

        adv_diff_integrator->registerTransportedQuantity(phi_var);
        adv_diff_integrator->setDiffusionCoefficient(phi_var, 0.0);

        /*************************** */
        adv_diff_integrator->registerTransportedQuantity(xi0_var);
        adv_diff_integrator->registerTransportedQuantity(xi1_var);

        const double kappa_xi = 1.0e-12;
        adv_diff_integrator->setDiffusionCoefficient(xi0_var, kappa_xi);
        adv_diff_integrator->setDiffusionCoefficient(xi1_var, kappa_xi);

        /*************************** */

        // Set the advection velocity of the bubble.
        adv_diff_integrator->setAdvectionVelocity(phi_var, time_integrator->getAdvectionVelocityVariable());

        /****************************** */
        adv_diff_integrator->setAdvectionVelocity(xi0_var, time_integrator->getAdvectionVelocityVariable());
        adv_diff_integrator->setAdvectionVelocity(xi1_var, time_integrator->getAdvectionVelocityVariable());

         /******************************* */

        ColumnInterface column;
        input_db->getDoubleArray("X_UR", &column.X_UR[0], NDIM);

        Pointer<RelaxationLSMethod> level_set_ops =
            new RelaxationLSMethod("RelaxationLSMethod", app_initializer->getComponentDatabase("RelaxationLSMethod"));
        LSLocateColumnInterface setLSLocateColumnInterface(
            "LSLocateColumnInterface", adv_diff_integrator, phi_var, column);
        level_set_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateColumnInterfaceCallbackFunction,
                                                                static_cast<void*>(&setLSLocateColumnInterface));
        IBAMR::LevelSetUtilities::SetLSProperties setSetLSProperties("SetLSProperties", level_set_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var, &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(&setSetLSProperties));

        // LS initial conditions
        // Pointer<CartGridFunction> phi_init = new LevelSetInitialCondition("ls_init", column);
        // adv_diff_integrator->setInitialConditions(phi_var, phi_init);

        /****************************** */
        Pointer<CartGridFunction> phi_init = new muParserCartGridFunction(
            "phi_init", app_initializer->getComponentDatabase("LevelSetInitialConditions"), grid_geometry);
        adv_diff_integrator->setInitialConditions(phi_var, phi_init);
        /***************************** */

        /************************ */
        Pointer<CartGridFunction> xi0_init = new muParserCartGridFunction(
            "xi0_init", app_initializer->getComponentDatabase("Xi0InitialConditions"), grid_geometry);

        Pointer<CartGridFunction> xi1_init = new muParserCartGridFunction(
            "xi1_init", app_initializer->getComponentDatabase("Xi1InitialConditions"), grid_geometry);

        adv_diff_integrator->setInitialConditions(xi0_var, xi0_init);
        adv_diff_integrator->setInitialConditions(xi1_var, xi1_init);

        /************************** */

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
        const double num_interface_cells = input_db->getDouble("NUM_INTERFACE_CELLS");

        // Callback functions can either be registered with the NS integrator, or the advection-diffusion integrator
        // Note that these will set the initial conditions for density and viscosity, based on level set information
        IBAMR::VCINSUtilities::SetFluidProperties setSetFluidProperties("SetFluidProperties",
                                                                        adv_diff_integrator,
                                                                        phi_var,
                                                                        rho_outside,
                                                                        rho_inside,
                                                                        mu_outside,
                                                                        mu_inside,
                                                                        num_interface_cells);
        time_integrator->registerResetFluidDensityFcn(&IBAMR::VCINSUtilities::callSetDensityCallbackFunction,
                                                      static_cast<void*>(&setSetFluidProperties));
        time_integrator->registerResetFluidViscosityFcn(&IBAMR::VCINSUtilities::callSetViscosityCallbackFunction,
                                                        static_cast<void*>(&setSetFluidProperties));

        /************************************ */

        /***************************** */

        // Register callback function for tagging refined cells for level set data
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        const double tag_min_value = -tag_thresh;
        const double tag_max_value = tag_thresh;
        IBAMR::LevelSetUtilities::TagLSRefinementCells ls_tagger(
            adv_diff_integrator, phi_var, tag_min_value, tag_max_value);
        time_integrator->registerApplyGradientDetectorCallback(&IBAMR::LevelSetUtilities::tagLSCells,
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
        level_set_ops->registerPhysicalBoundaryCondition(phi_bc_coef);

        /************************ */
        RobinBcCoefStrategy<NDIM>* xi0_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("Xi0BcCoefs"))
        {
            xi0_bc_coef = new muParserRobinBcCoefs(
                "xi0_bc_coef", app_initializer->getComponentDatabase("Xi0BcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(xi0_var, xi0_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* xi1_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("Xi1BcCoefs"))
        {
            xi1_bc_coef = new muParserRobinBcCoefs(
                "xi1_bc_coef", app_initializer->getComponentDatabase("Xi1BcCoefs"), grid_geometry);
            adv_diff_integrator->setPhysicalBcCoef(xi1_var, xi1_bc_coef);
        }
        /************************ */

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

        // Forcing terms
        // std::vector<double> grav_const(NDIM);
        // input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        // Pointer<CartGridFunction> grav_force =
        //     new IBAMR::VCINSUtilities::GravityForcing("GravityForcing", time_integrator, grav_const);

        Pointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var);

        Pointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");

        eul_forces->addFunction(surface_tension_force);
        // eul_forces->addFunction(grav_force);
        /******************************* */
        Pointer<Database> rmt_db = input_db->getDatabase("RMT");
        const double GGS = rmt_db->getDouble("GGS");
        const double rho_s = input_db->getDoubleWithDefault("RHO_I", 1.0);
        const double g0y = input_db->getDoubleWithDefault("GRAV_Y", 0.0);
        const double Tramp = input_db->getDoubleWithDefault("T_RAMP", 0.0);
        const double nu_s = input_db->getDoubleWithDefault("NU", 0.499);

        /******************************* */

        time_integrator->registerBodyForceFunction(eul_forces);

        // pout << ">>> RMT FORCE REGISTERED <<<\n";
        pout << ">>> RMT FORCE DEBUG MODE: NOT REGISTERED TO NS <<<\n";

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();

        // Initialize hierarchy configuration and data on all patches.

        {
            auto* vdb = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
            Pointer<SAMRAI::hier::VariableContext> ctx = adv_diff_integrator->getCurrentContext();
            const SAMRAI::hier::IntVector<NDIM> gcw(1);

            phixi_idx = vdb->registerVariableAndContext(phixi_var, ctx, gcw);
            dphi_idx = vdb->registerVariableAndContext(dphi_var, ctx, gcw);
            disp_mag_idx = vdb->registerVariableAndContext(disp_mag_var, ctx, gcw);

            // w_idx = vdb->registerVariableAndContext(w_var, ctx, gcw);

            // phiw_err_idx = var_db->registerVariableAndContext(
            //     phiw_err_var, var_db->getContext("PHI_MINUS_W_CTX"), IntVector<NDIM>(0));
        }

        // rmt_force->setPhi0AndWIndices(phi0_idx, w_idx);

        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
            // visit_data_writer->registerPlotQuantity("W_FROM_XI", "SCALAR", w_idx, 0, 1, "CELL");
        }

        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        /***************************** */

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                // if (patch->getPatchData(w_idx).isNull())
                // {
                //     patch->allocatePatchData(w_idx, 0.0);
                // }

                // auto w_data = Pointer<CellData<NDIM, double> >(patch->getPatchData(w_idx));
                // w_data->fillAll(0.0);
            }
        }

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());

                if (!patch->checkAllocated(phixi_idx)) patch->allocatePatchData(phixi_idx);
                if (!patch->checkAllocated(dphi_idx)) patch->allocatePatchData(dphi_idx);
                if (!patch->checkAllocated(disp_mag_idx)) patch->allocatePatchData(disp_mag_idx);

                // if (!patch->checkAllocated(w_idx)) patch->allocatePatchData(w_idx);
            }
        }
        /******************** */
        // remove
        int phi_scr = IBTK::invalid_index;
        int xi0_scr = IBTK::invalid_index;
        int xi1_scr = IBTK::invalid_index;

        {
            auto* vdb = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

            // CURRENT indices we will also use in time loop
            phi_cur_idx = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getCurrentContext());
            xi0_cur_idx = vdb->mapVariableAndContextToIndex(xi0_var, adv_diff_integrator->getCurrentContext());
            xi1_cur_idx = vdb->mapVariableAndContextToIndex(xi1_var, adv_diff_integrator->getCurrentContext());

            phi_new_idx = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getNewContext());
            // const int phi_scr = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getScratchContext());

            xi0_new_idx = vdb->mapVariableAndContextToIndex(xi0_var, adv_diff_integrator->getNewContext());
            // const int xi0_scr = vdb->mapVariableAndContextToIndex(xi0_var, adv_diff_integrator->getScratchContext());

            xi1_new_idx = vdb->mapVariableAndContextToIndex(xi1_var, adv_diff_integrator->getNewContext());
            // const int xi1_scr = vdb->mapVariableAndContextToIndex(xi1_var, adv_diff_integrator->getScratchContext());

            phi_scr = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getScratchContext());
            xi0_scr = vdb->mapVariableAndContextToIndex(xi0_var, adv_diff_integrator->getScratchContext());
            xi1_scr = vdb->mapVariableAndContextToIndex(xi1_var, adv_diff_integrator->getScratchContext());

            visit_data_writer->registerPlotQuantity("phi", "SCALAR", phi_cur_idx, 0, 1, "CELL");
            visit_data_writer->registerPlotQuantity("phi_from_xi", "SCALAR", phixi_idx, 0, 1, "CELL");
            visit_data_writer->registerPlotQuantity("DeltaPhi", "SCALAR", dphi_idx, 0, 1, "CELL");

            visit_data_writer->registerPlotQuantity("DISP_MAG", "SCALAR", disp_mag_idx, 0, 1, "CELL");
        }


        // Pointer<Database> rmt_db = input_db->getDatabase("RMT");

const double rmt_band_width = rmt_db->getDouble("BAND_WIDTH");
const int rmt_num_passes = rmt_db->getInteger("NUM_PASSES");
const int rmt_reinit_num_iter = rmt_db->getInteger("REINIT_NUM_ITER");
const double rmt_reinit_dtau = rmt_db->getDouble("REINIT_DTAU");

Pointer<RMTMethod> rmt_method =
    new RMTMethod("RMTMethod",
                  patch_hierarchy,
                  phi_var,
                  xi0_var,
                  xi1_var,
                  phi_cur_idx,
                  xi0_cur_idx,
                  xi1_cur_idx,
                  rmt_band_width,
                  rmt_num_passes,
                  rmt_reinit_num_iter,
                  rmt_reinit_dtau,
                  GGS,
                  nu_s,
                  rho_s,
                  g0y,
                  Tramp);

if (uses_visit)
{
    rmt_method->registerVisItDataWriter(visit_data_writer);
}

rmt_method->initializeForceFunction(phi_cur_idx,
                                    phi_new_idx,
                                    phi_scr,
                                    xi0_cur_idx,
                                    xi0_new_idx,
                                    xi0_scr,
                                    xi1_cur_idx,
                                    xi1_new_idx,
                                    xi1_scr);

rmt_method->initializeLevelData(time_integrator->getIntegratorTime());

eul_forces->addFunction(rmt_method->getForceFunction());
        /***************************** */
        /* copy phi(t=0) into phi0 ONCE */
        /***************************** */
        // TBOX_ASSERT(phi_cur_idx != IBTK::invalid_index);
        // TBOX_ASSERT(phi0_idx != IBTK::invalid_index);

        // for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        // {
        //     Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

        //     for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        //     {
        //         Pointer<Patch<NDIM> > patch = level->getPatch(p());

        //         if (patch->getPatchData(phi0_idx).isNull())
        //         {
        //             patch->allocatePatchData(phi0_idx, 0.0);
        //         }

        //         auto phi0_data = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi0_idx));
        //         auto phi_cur_data = Pointer<CellData<NDIM, double> >(patch->getPatchData(phi_cur_idx));

        //         TBOX_ASSERT(!phi0_data.isNull());
        //         TBOX_ASSERT(!phi_cur_data.isNull());

        //         phi0_data->copy(*phi_cur_data);
        //     }
        // }
        /************************* */

        /*********************************** */

        // Remove the AppInitializer
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("centroid.dat");
            out << "# time xc yc\n";
        }

        if (dump_viz_data && iteration_num % viz_dump_interval == 0)
        {
            fillDispMagFromXi(patch_hierarchy, disp_mag_idx, phi_cur_idx, xi0_cur_idx, xi1_cur_idx);

            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        // File to write to for fluid mass data
        ofstream mass_file, front_file, height_file, umax_file;

        if (!IBTK_MPI::getRank())
        {
            mass_file.open("mass_fluid.txt");
            front_file.open("front_position.txt");
            height_file.open("height_fluid.txt");

            umax_file.open("umax_vs_time.txt");
            umax_file << "# time\tumax\n";
        }
        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
      

        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())


        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            {
                pout << "\n";
                pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                pout << "At beginning of timestep # " << iteration_num << "\n";
                pout << "Simulation time is " << loop_time << "\n";
            }
            dt = time_integrator->getMaximumTimeStepSize();

            // /********************** */

            {
                auto* vdb = VariableDatabase<NDIM>::getDatabase();

                phi_cur_idx = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getCurrentContext());
                xi0_cur_idx = vdb->mapVariableAndContextToIndex(xi0_var, adv_diff_integrator->getCurrentContext());
                xi1_cur_idx = vdb->mapVariableAndContextToIndex(xi1_var, adv_diff_integrator->getCurrentContext());

                int u_cur_idx = vdb->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                                  time_integrator->getCurrentContext());

                const double* dx_global = grid_geometry->getDx();
                const double clamp_end_cut = 2.0 * dx_global[0];
            }

        
            /*************************** */

            time_integrator->advanceHierarchy(dt);

            /************************ */

            {
                auto* vdb = VariableDatabase<NDIM>::getDatabase();

                phi_cur_idx = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getCurrentContext());
            }

            computeAndWriteCentroid(patch_hierarchy, phi_cur_idx, time_integrator->getIntegratorTime(), "centroid.dat");

            // 1) get current-context indices after advanceHierarchy
            int u_cur_idx = IBTK::invalid_index;

            // 1) get current-context indices after advanceHierarchy
            {
                auto* vdb = VariableDatabase<NDIM>::getDatabase();

                phi_cur_idx = vdb->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getCurrentContext());

                xi0_cur_idx = vdb->mapVariableAndContextToIndex(xi0_var, adv_diff_integrator->getCurrentContext());

                xi1_cur_idx = vdb->mapVariableAndContextToIndex(xi1_var, adv_diff_integrator->getCurrentContext());

                // w_idx = vdb->mapVariableAndContextToIndex(w_var, adv_diff_integrator->getCurrentContext());

                u_cur_idx = vdb->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                              time_integrator->getCurrentContext());
            }

            const double* dx_global = grid_geometry->getDx();
            const double dx_max = std::max(dx_global[0], dx_global[1]);

            const double extrap_band = 4.0 * dx_max;
            const double pin_band = 2.0 * dx_max;
            const double reinit_band = 4.0 * dx_max;

            const double clamp_end_cut = 2.0 * dx_global[0];

            
            rmt_method->postprocess();


            const double disp_max =
                fillDispMagFromXi(patch_hierarchy, disp_mag_idx, phi_cur_idx, xi0_cur_idx, xi1_cur_idx);

            // pout << "[DISP_MAX] t=" << (loop_time + dt) << " disp_max=" << disp_max << "\n";

            if (!IBTK_MPI::getRank())
            {
                static std::ofstream center_file("center_disp_vs_time.txt");

                if (iteration_num == 0)
                {
                    center_file << "# time disp_max\n";
                }

                center_file << std::setprecision(12) << std::scientific << (loop_time + dt) << " " << disp_max << "\n";
            }

            /*********************** */
            loop_time += dt;

            {
                pout << "\n";
                pout << "At end       of timestep # " << iteration_num << "\n";
                pout << "Simulation time is " << loop_time << "\n";
                pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                pout << "\n";
            }
            // Compute the fluid mass in the domain from interpolated density
            const int rho_ins_idx = time_integrator->getLinearOperatorRhoPatchDataIndex();
            // #if !defined(NDEBUG)
            //             TBOX_ASSERT(rho_ins_idx >= 0);
            // #endif
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            HierarchySideDataOpsReal<NDIM, double> hier_rho_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const double mass_fluid = hier_rho_data_ops.integral(rho_ins_idx, wgt_sc_idx);

            // Compute the front position and the height of the fluid. This can be approximately done by finding
            // the maximum x and y (z in 3D) coordinate of the negative level set values.
            double fluid_front = -std::numeric_limits<double>::max();
            double fluid_height = -std::numeric_limits<double>::max();
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const int phi_idx = var_db->mapVariableAndContextToIndex(phi_var, adv_diff_integrator->getCurrentContext());
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(phi_idx);
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
                            coord[d] = patch_X_lower[d] +
                                       patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                        }

                        // Check if this coordinate is maximal for negative level set values
                        const double phi = (*D_data)(ci);
                        const int front_dim = 0;
                        const int height_dim = NDIM - 1;
                        fluid_front = (phi <= 0.0 && coord[front_dim] >= fluid_front) ? coord[front_dim] : fluid_front;
                        fluid_height =
                            (phi <= 0.0 && coord[height_dim] >= fluid_height) ? coord[height_dim] : fluid_height;
                    }
                }
            }

            // Max reduction
            fluid_front = IBTK_MPI::maxReduction(fluid_front);
            fluid_height = IBTK_MPI::maxReduction(fluid_height);

            // Write to file
            if (!IBTK_MPI::getRank())
            {
                mass_file << std::setprecision(13) << loop_time << "\t" << mass_fluid << std::endl;
                front_file << std::setprecision(13) << loop_time << "\t" << fluid_front << std::endl;
                height_file << std::setprecision(13) << loop_time << "\t" << fluid_height << std::endl;
                umax_file << std::setprecision(13) << loop_time << "\t" << disp_max << std::endl;
            }

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";

                fillDispMagFromXi(patch_hierarchy, disp_mag_idx, phi_cur_idx, xi0_cur_idx, xi1_cur_idx);

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
        if (!IBTK_MPI::getRank())
        {
            mass_file.close();
            front_file.close();
            height_file.close();
            umax_file.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main

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
    std::snprintf(temp_buf, sizeof(temp_buf), "%05d.samrai.%05d", iteration_num, IBTK_MPI::getRank());
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
