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
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/BrinkmanAdvDiffBcHelper.h>
#include <ibamr/BrinkmanAdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/BrinkmanPenalizationRigidBodyDynamics.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBInterpolantHierarchyIntegrator.h>
#include <ibamr/IBInterpolantMethod.h>
#include <ibamr/IBLevelSetMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include "ibtk/IndexUtilities.h"
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Application specific includes
#include "LevelSetInitialConditionEgg.h"
#include "LevelSetInitialConditionTorus.h"

// Struct to specify the variables required for inhomogeneous Neumann conditions for Brinkman penalization
struct BrinkmanPenalizationCtx
{
    Pointer<BrinkmanAdvDiffSemiImplicitHierarchyIntegrator> adv_diff_hier_integrator;
    Pointer<SideVariable<NDIM, double> > grad_phi_sc_var;
    Pointer<CellVariable<NDIM, double> > grad_phi_cc_var;
    Pointer<CellVariable<NDIM, double> > beta_var;
    Pointer<CellVariable<NDIM, double> > g_cc_var;
    int beta_scratch_idx;
    int grad_phi_sc_idx;
    int grad_phi_cc_idx;
    int g_cc_scratch_idx;
    int interface_cell_patch_idx;
    int num_prop_cells;
    int prop_counter_idx;
};

void
evaluate_brinkman_bc_callback_fcn(int B_idx,
                                  Pointer<CellVariable<NDIM, double> > ls_solid_var,
                                  Pointer<HierarchyMathOps> hier_math_ops,
                                  double time,
                                  void* ctx)
{
    BrinkmanPenalizationCtx* bp_ctx = static_cast<BrinkmanPenalizationCtx*>(ctx);
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int phi_idx =
        var_db->mapVariableAndContextToIndex(ls_solid_var, bp_ctx->adv_diff_hier_integrator->getNewContext());
    const int phi_scratch_idx =
        var_db->mapVariableAndContextToIndex(ls_solid_var, bp_ctx->adv_diff_hier_integrator->getScratchContext());

    // Fill the ghost cells of phi_scratch.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> phi_transaction_comps(1);
    phi_transaction_comps[0] =
        InterpolationTransactionComponent(phi_scratch_idx,
                                          phi_idx,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          false,
                                          "CONSERVATIVE_COARSEN",
                                          "LINEAR",
                                          false,
                                          bp_ctx->adv_diff_hier_integrator->getPhysicalBcCoefs(ls_solid_var));
    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(phi_transaction_comps, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    // Computing the gradient of signed distance at cell center.
    int grad_phi_cc_idx = bp_ctx->grad_phi_cc_idx;
    Pointer<CellVariable<NDIM, double> > grad_phi_cc_var = bp_ctx->grad_phi_cc_var;
    hier_math_ops->grad(grad_phi_cc_idx, grad_phi_cc_var, -1.0, phi_scratch_idx, ls_solid_var, nullptr, time);

    // computing the gradient of signed distance at side center.
    int grad_phi_sc_idx = bp_ctx->grad_phi_sc_idx;
    Pointer<SideVariable<NDIM, double> > grad_phi_sc_var = bp_ctx->grad_phi_sc_var;
    hier_math_ops->grad(grad_phi_sc_idx, grad_phi_sc_var, true, -1.0, phi_scratch_idx, ls_solid_var, nullptr, time);

    int g_cc_scratch_idx = bp_ctx->g_cc_scratch_idx;
    int interface_cell_patch_idx = bp_ctx->interface_cell_patch_idx;
    int prop_counter_idx = bp_ctx->prop_counter_idx;

    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    IntVector<NDIM> ratio = finest_level->getRatio();

    const Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_hierarchy->getGridGeometry();
    const double* const grid_x_lower = grid_geom->getXLower();
    const SAMRAI::hier::Box<NDIM> domain_box =
        SAMRAI::hier::Box<NDIM>::refine(grid_geom->getPhysicalDomain()[0], ratio);
    const hier::Index<NDIM>& grid_lower_idx = domain_box.lower();

    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();

        Pointer<CellData<NDIM, double> > grad_phi_cc_data = patch->getPatchData(grad_phi_cc_idx);
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(phi_scratch_idx);
        Pointer<CellData<NDIM, double> > g_cc_data = patch->getPatchData(g_cc_scratch_idx);
        Pointer<CellData<NDIM, double> > interface_cell_data = patch->getPatchData(interface_cell_patch_idx);
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());
            CellIndex<NDIM> ci_e = ci;
            CellIndex<NDIM> ci_w = ci;
            CellIndex<NDIM> ci_n = ci;
            CellIndex<NDIM> ci_s = ci;
#if (NDIM == 3)
            CellIndex<NDIM> ci_f = ci;
            CellIndex<NDIM> ci_b = ci;
#endif
            ci_e(0) += 1;
            ci_w(0) -= 1;
            ci_n(1) += 1;
            ci_s(1) -= 1;
#if (NDIM == 3)
            ci_f(2) += 1;
            ci_b(2) -= 1;
#endif
            const double ls_c_value = (*ls_data)(ci);
            const double ls_e_value = (*ls_data)(ci_e);
            const double ls_w_value = (*ls_data)(ci_w);
            const double ls_n_value = (*ls_data)(ci_n);
            const double ls_s_value = (*ls_data)(ci_s);
#if (NDIM == 3)
            const double ls_f_value = (*ls_data)(ci_f);
            const double ls_b_value = (*ls_data)(ci_b);
#endif

            // Finding interface cells.
            if ((ls_c_value * ls_e_value <= 0.0) || (ls_c_value * ls_w_value <= 0.0) ||
                (ls_c_value * ls_n_value <= 0.0) || (ls_c_value * ls_s_value <= 0.0)
#if (NDIM == 3)
                || (ls_c_value * ls_f_value <= 0.0) || (ls_c_value * ls_b_value <= 0.0)
#endif
            )
            {
                IBTK::VectorNd coord = IBTK::Vector::Zero();
                IBTK::VectorNd interface_coord = IBTK::Vector::Zero();
                for (int d = 0; d < NDIM; ++d)
                    coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

                // Compute g in the interface cell.
                double g = 0.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    interface_coord[d] = coord[d] + ls_c_value * (*grad_phi_cc_data)(ci, d);
                }
#if (NDIM == 2)
                g = cos(interface_coord[0]) * sin(interface_coord[1]) * (*grad_phi_cc_data)(ci, 0) +
                    sin(interface_coord[0]) * cos(interface_coord[1]) * (*grad_phi_cc_data)(ci, 1);
#elif (NDIM == 3)
                g = cos(interface_coord[1]) * cos(interface_coord[2]) * sin(interface_coord[0]) *
                        (*grad_phi_cc_data)(ci, 0) +
                    cos(interface_coord[0]) * cos(interface_coord[2]) * sin(interface_coord[1]) *
                        (*grad_phi_cc_data)(ci, 1) +
                    cos(interface_coord[0]) * cos(interface_coord[1]) * sin(interface_coord[2]) *
                        (*grad_phi_cc_data)(ci, 2);
#endif
                (*g_cc_data)(ci) = g;
                (*interface_cell_data)(ci) = 1.0;
            }
        }
    }

    // ghost cell filling for n and g. Both will have the ghost cells of n_prop.
    RobinBcCoefStrategy<NDIM>* ls_bc_coef = bp_ctx->adv_diff_hier_integrator->getPhysicalBcCoefs(ls_solid_var).front();
    std::vector<RobinBcCoefStrategy<NDIM>*> ls_bc_coefs(NDIM, ls_bc_coef);
    std::vector<InterpolationTransactionComponent> gn_transaction_comps(3);
    gn_transaction_comps[0] = InterpolationTransactionComponent(
        g_cc_scratch_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, ls_bc_coef);

    gn_transaction_comps[1] = InterpolationTransactionComponent(
        grad_phi_cc_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false, ls_bc_coefs);

    gn_transaction_comps[2] = InterpolationTransactionComponent(interface_cell_patch_idx,
                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                false,
                                                                "CONSERVATIVE_COARSEN",
                                                                "LINEAR",
                                                                false,
                                                                ls_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> g_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    g_hier_bdry_fill->initializeOperatorState(gn_transaction_comps, patch_hierarchy);
    g_hier_bdry_fill->fillData(time);

    const int num_prop_cells = bp_ctx->num_prop_cells;

    // Propagate the interfacial g value to the neighboring cells along the normal direction.
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        Pointer<CellData<NDIM, double> > g_cc_data = patch->getPatchData(g_cc_scratch_idx);
        Pointer<CellData<NDIM, double> > interface_cell_data = patch->getPatchData(interface_cell_patch_idx);
        Pointer<CellData<NDIM, double> > prop_counter_data = patch->getPatchData(prop_counter_idx);
        Pointer<CellData<NDIM, double> > grad_phi_cc_data = patch->getPatchData(grad_phi_cc_idx);
        const Box<NDIM>& ghost_box = g_cc_data->getGhostBox();
        const hier::Index<NDIM>& ghost_lower_idx = ghost_box.lower();
        const hier::Index<NDIM>& ghost_upper_idx = ghost_box.upper();

        for (Box<NDIM>::Iterator it(ghost_box); it; it++)
        {
            CellIndex<NDIM> ci(it());

            if ((*interface_cell_data)(ci) > 0.0)
            {
                std::vector<CellIndex<NDIM> > nbr_cell_normal(num_prop_cells);     // vector stores cell index of
                                                                                   // virtual cells in +n direction
                std::vector<CellIndex<NDIM> > nbr_cell_neg_normal(num_prop_cells); // vector stores cell index of
                                                                                   // virtual cells in -n direction

                // location of the interface cells.
                IBTK::VectorNd coord = IBTK::Vector::Zero();
                IBTK::VectorNd coord_nbr = IBTK::Vector::Zero();
                const double g = (*g_cc_data)(ci);
                for (int d = 0; d < NDIM; ++d)
                    coord[d] = grid_x_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - grid_lower_idx(d)) + 0.5);

                // propagating g in +n direction.
                int counter = 0;
                for (auto& nbr_cell : nbr_cell_normal)
                {
                    counter++;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord_nbr[d] = coord[d] + counter * patch_dx[d] * (*grad_phi_cc_data)(ci, d);
                    }
                    nbr_cell = IndexUtilities::getCellIndex(&coord_nbr[0], grid_geom, ratio);

                    if (nbr_cell[0] >= ghost_lower_idx[0] && nbr_cell[0] <= ghost_upper_idx[0] &&
                        nbr_cell[1] >= ghost_lower_idx[1] && nbr_cell[1] <= ghost_upper_idx[1]
#if (NDIM == 3)
                        && nbr_cell[2] >= ghost_lower_idx[2] && nbr_cell[2] <= ghost_upper_idx[2]
#endif
                    )
                    {
                        (*prop_counter_data)(nbr_cell) = (*prop_counter_data)(nbr_cell) + 1.0;
                        if (std::abs((*g_cc_data)(nbr_cell)) <= std::abs(g))
                        {
                            (*g_cc_data)(nbr_cell) = g;
                        }
                    }
                }

                // propagating g in -n direction.
                counter = 0;
                for (auto& nbr_cell : nbr_cell_neg_normal)
                {
                    counter++;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        coord_nbr[d] = coord[d] - counter * patch_dx[d] * (*grad_phi_cc_data)(ci, d);
                    }
                    nbr_cell = IndexUtilities::getCellIndex(&coord_nbr[0], grid_geom, ratio);

                    if (nbr_cell[0] >= ghost_lower_idx[0] && nbr_cell[0] <= ghost_upper_idx[0] &&
                        nbr_cell[1] >= ghost_lower_idx[1] && nbr_cell[1] <= ghost_upper_idx[1]
#if (NDIM == 3)
                        && nbr_cell[2] >= ghost_lower_idx[2] && nbr_cell[2] <= ghost_upper_idx[2]
#endif
                    )
                    {
                        (*prop_counter_data)(nbr_cell) = (*prop_counter_data)(nbr_cell) + 1.0;
                        if (std::abs((*g_cc_data)(nbr_cell)) <= std::abs(g))
                        {
                            (*g_cc_data)(nbr_cell) = g;
                        }
                    }
                }
            }
        }
    }

    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<SideData<NDIM, double> > grad_phi_sc_data = patch->getPatchData(grad_phi_sc_idx);
        Pointer<SideData<NDIM, double> > B_data = patch->getPatchData(B_idx);
        Pointer<CellData<NDIM, double> > g_cc_data = patch->getPatchData(g_cc_scratch_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);
                // Interpolating g from cell to side.
                const double g_sc = 0.5 * ((*g_cc_data)(si.toCell(0)) + (*g_cc_data)(si.toCell(1)));
                (*B_data)(si) = g_sc * (*grad_phi_sc_data)(si);
            }
        }
    }
    hier_math_ops->interp(bp_ctx->beta_scratch_idx, bp_ctx->beta_var, B_idx, grad_phi_sc_var, nullptr, time, true);
}

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

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

        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        //

        Pointer<AdvDiffHierarchyIntegrator> time_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "SEMI_IMPLICIT");
        time_integrator = new BrinkmanAdvDiffSemiImplicitHierarchyIntegrator(
            "BrinkmanAdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("BrinkmanAdvDiffSemiImplicitHierarchyIntegrator"));

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

        // Set up the advected and diffused quantity.
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        const string& ls_name = "level_set_solid";
        Pointer<CellVariable<NDIM, double> > phi_solid_var = new CellVariable<NDIM, double>(ls_name);
        time_integrator->registerTransportedQuantity(phi_solid_var, true);
        time_integrator->setDiffusionCoefficient(phi_solid_var, 0.0);

        // Origin of the shape.
        IBTK::VectorNd origin;
#if (NDIM == 2)
        origin << M_PI, M_PI;
        // ls initial condition for SimpleEgg
        Pointer<CartGridFunction> phi_solid_init = new LevelSetInitialConditionEgg("ls_init", grid_geometry, origin);
#elif (NDIM == 3)
        origin << M_PI, M_PI, M_PI;

        // ls initial condition for torus.
        IBTK::Vector2d t(1.5, 0.8);
        Pointer<CartGridFunction> phi_solid_init =
            new LevelSetInitialConditionTorus("ls_init", grid_geometry, origin, t);
#endif

        time_integrator->setInitialConditions(phi_solid_var, phi_solid_init);

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
            time_integrator->setPhysicalBcCoef(phi_solid_var, phi_bc_coef);
        }
        Pointer<CellVariable<NDIM, double> > q_var = new CellVariable<NDIM, double>("q");
        time_integrator->registerTransportedQuantity(q_var, true);
        time_integrator->setDiffusionCoefficient(q_var, input_db->getDouble("KAPPA"));

        if (input_db->keyExists("TransportedQuantityInitialConditions"))
        {
            Pointer<CartGridFunction> q_init = new muParserCartGridFunction(
                "q_init", app_initializer->getComponentDatabase("TransportedQuantityInitialConditions"), grid_geometry);
            time_integrator->setInitialConditions(q_var, q_init);
        }

        RobinBcCoefStrategy<NDIM>* q_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("TransportedQuantityBcCoefs"))
        {
            q_bc_coef = new muParserRobinBcCoefs(
                "q_bc_coef", app_initializer->getComponentDatabase("TransportedQuantityBcCoefs"), grid_geometry);
        }
        time_integrator->setPhysicalBcCoef(q_var, q_bc_coef);

        // Variables for computing flux-forcing function.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int num_prop_cells = input_db->getInteger("NUMBER_OF_PROPAGATION_CELLS");
        Pointer<CellVariable<NDIM, double> > g_cc_var = new CellVariable<NDIM, double>("g_cc");
        const int g_cc_scratch_idx =
            var_db->registerVariableAndContext(g_cc_var, var_db->getContext("g_cc_context"), num_prop_cells);
        Pointer<SideVariable<NDIM, double> > grad_phi_sc_var = new SideVariable<NDIM, double>("grad_phi_sc_var");
        Pointer<CellVariable<NDIM, double> > grad_phi_cc_var = new CellVariable<NDIM, double>("grad_phi_cc_var", NDIM);
        const int grad_phi_sc_idx =
            var_db->registerVariableAndContext(grad_phi_sc_var, var_db->getContext("grad_phi_sc_context"), 0);
        const int grad_phi_cc_idx = var_db->registerVariableAndContext(
            grad_phi_cc_var, var_db->getContext("grad_phi_cc_context"), num_prop_cells);
        Pointer<CellVariable<NDIM, double> > beta_var = new CellVariable<NDIM, double>("beta", NDIM);
        const int beta_scratch_idx =
            var_db->registerVariableAndContext(beta_var, var_db->getContext("beta_context"), 0);
        Pointer<CellVariable<NDIM, double> > interface_cell_var = new CellVariable<NDIM, double>("interface_cell_var");
        const int interface_cell_idx = var_db->registerVariableAndContext(
            interface_cell_var, var_db->getContext("interface_cell_context"), num_prop_cells);
        Pointer<CellVariable<NDIM, double> > prop_counter_var = new CellVariable<NDIM, double>("prop_counter");
        const int prop_counter_idx = var_db->registerVariableAndContext(
            prop_counter_var, var_db->getContext("prop_counter_context"), num_prop_cells);

        BrinkmanPenalizationCtx brinkman_ctx;
        brinkman_ctx.adv_diff_hier_integrator = time_integrator;
        brinkman_ctx.grad_phi_sc_var = grad_phi_sc_var;
        brinkman_ctx.grad_phi_sc_idx = grad_phi_sc_idx;
        brinkman_ctx.grad_phi_cc_var = grad_phi_cc_var;
        brinkman_ctx.grad_phi_cc_idx = grad_phi_cc_idx;
        brinkman_ctx.beta_scratch_idx = beta_scratch_idx;
        brinkman_ctx.beta_var = beta_var;
        brinkman_ctx.g_cc_var = g_cc_var;
        brinkman_ctx.g_cc_scratch_idx = g_cc_scratch_idx;
        brinkman_ctx.interface_cell_patch_idx = interface_cell_idx;
        brinkman_ctx.num_prop_cells = num_prop_cells;
        brinkman_ctx.prop_counter_idx = prop_counter_idx;

        const string indicator_function_type = input_db->getString("INDICATOR_FUNCTION_TYPE");
        const double eta = input_db->getDouble("ETA");
        const double num_of_interface_cells = input_db->getDouble("NUMBER_OF_INTERFACE_CELLS");
        Pointer<BrinkmanAdvDiffBcHelper> brinkman_adv_diff =
            new BrinkmanAdvDiffBcHelper("BrinkmanAdvDiffBcHelper", time_integrator);

        // setting inhomogeneous Neumann at the cylinder surface.
        brinkman_adv_diff->registerInhomogeneousBC(q_var,
                                                   phi_solid_var,
                                                   "NEUMANN",
                                                   &evaluate_brinkman_bc_callback_fcn,
                                                   static_cast<void*>(&brinkman_ctx),
                                                   indicator_function_type,
                                                   num_of_interface_cells,
                                                   eta);

        Pointer<BrinkmanAdvDiffSemiImplicitHierarchyIntegrator> bp_adv_diff_hier_integrator = time_integrator;
        bp_adv_diff_hier_integrator->registerBrinkmanAdvDiffBcHelper(brinkman_adv_diff);

        if (input_db->keyExists("TransportedQuantityForcingFunction"))
        {
            Pointer<CellVariable<NDIM, double> > F_var = new CellVariable<NDIM, double>("F", 1);
            Pointer<CartGridFunction> q_forcing_fcn = new muParserCartGridFunction(
                "q_forcing_fcn",
                app_initializer->getComponentDatabase("TransportedQuantityForcingFunction"),
                grid_geometry);
            time_integrator->registerSourceTerm(F_var, true);
            time_integrator->setSourceTermFunction(F_var, q_forcing_fcn);
            time_integrator->setSourceTerm(q_var, F_var);
        }

        Pointer<CartGridFunction> T_ex = new muParserCartGridFunction(
            "q_ex", app_initializer->getComponentDatabase("TransportedQuantityExactSolutions"), grid_geometry);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
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
        // VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int q_idx = var_db->mapVariableAndContextToIndex(q_var, time_integrator->getCurrentContext());
        const int q_ex_cloned_idx = var_db->registerClonedPatchDataIndex(q_var, q_idx);
        const int q_err_cloned_idx = var_db->registerClonedPatchDataIndex(q_var, q_idx);
        const int phi_idx = var_db->mapVariableAndContextToIndex(phi_solid_var, time_integrator->getCurrentContext());
        const int phi_cloned_idx = var_db->registerClonedPatchDataIndex(phi_solid_var, phi_idx);

        // We need a ls var with single ghost cell width.
        const int phi_scratch_idx =
            var_db->registerVariableAndContext(phi_solid_var, var_db->getContext("scratch"), IntVector<NDIM>(1));

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(q_ex_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(q_err_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(phi_scratch_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(phi_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(g_cc_scratch_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(grad_phi_sc_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(grad_phi_cc_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(beta_scratch_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(interface_cell_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(prop_counter_idx, loop_time);
        }
        visit_data_writer->registerPlotQuantity("q_exact", "SCALAR", q_ex_cloned_idx);
        visit_data_writer->registerPlotQuantity("q_error", "SCALAR", q_err_cloned_idx);
        visit_data_writer->registerPlotQuantity("g", "SCALAR", g_cc_scratch_idx);
        visit_data_writer->registerPlotQuantity("prop_counter", "SCALAR", prop_counter_idx);

        for (int i = 0; i < NDIM; i++)
        {
            if (i == 0)
                visit_data_writer->registerPlotQuantity("beta_x", "SCALAR", beta_scratch_idx, i);
            else if (i == 1)
                visit_data_writer->registerPlotQuantity("beta_y", "SCALAR", beta_scratch_idx, i);
        }

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);

        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.setToScalar(g_cc_scratch_idx, 0.0, false /*interior_only*/);
        hier_cc_data_ops.setToScalar(interface_cell_idx, 0.0, false /*interior_only*/);
        hier_cc_data_ops.setToScalar(prop_counter_idx, 0.0, false /*interior_only*/);

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

            // Determine the accuracy of the computed solution.
            pout << "\n"
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                 << "Computing error norms.\n\n";

            // Calculate Heaviside function and mask the error indices.
            for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    const Box<NDIM>& patch_box = patch->getBox();
                    const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
                    Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
                    Pointer<CellData<NDIM, double> > H_data = patch->getPatchData(phi_cloned_idx);

                    for (Box<NDIM>::Iterator it(patch_box); it; it++)
                    {
                        CellIndex<NDIM> ci(it());
                        double phi = (*phi_data)(ci);
                        double Hphi;
                        if (phi > 0.0)
                        {
                            Hphi = 1.0;
                        }
                        else
                        {
                            Hphi = 0.0;
                        }

                        const double chi = 1.0 - Hphi;
                        (*H_data)(ci) = 1.0 - chi;
                    }
                }
            }

            hier_cc_data_ops.multiply(phi_cloned_idx, phi_cloned_idx, wgt_cc_idx);
            T_ex->setDataOnPatchHierarchy(q_ex_cloned_idx, q_var, patch_hierarchy, loop_time);
            hier_cc_data_ops.subtract(q_err_cloned_idx, q_idx, q_ex_cloned_idx);

            pout << "Error in q at time " << loop_time << ":\n"
                 << "  L1-norm:  " << std::setprecision(10) << hier_cc_data_ops.L1Norm(q_err_cloned_idx, phi_cloned_idx)
                 << "\n"
                 << "  L2-norm:  " << hier_cc_data_ops.L2Norm(q_err_cloned_idx, phi_cloned_idx) << "\n"
                 << "  max-norm: " << hier_cc_data_ops.maxNorm(q_err_cloned_idx, phi_cloned_idx) << "\n";

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
            }
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(q_ex_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(q_err_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(phi_cloned_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(phi_scratch_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(g_cc_scratch_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(grad_phi_sc_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(grad_phi_cc_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(beta_scratch_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(interface_cell_idx);
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(prop_counter_idx);
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main
