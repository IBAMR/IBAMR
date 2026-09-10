// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
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
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/MemoryDatabase.h>

#include <LocationIndexRobinBcCoefs.h>
#include <TimeRefinementIntegrator.h>

#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

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
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double>> D_data = patch->getPatchData(D_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
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

namespace
{
enum class RedistributionField
{
    CONSTANT,
    CUTOFF,
    NEAR_CUTOFF,
    STEEP,
    PLANE,
    CURVE
};

void
copy_redistribution_distance(const int idx, Pointer<HierarchyMathOps> ops, double, bool, void* ctx)
{
    const int source_idx = *static_cast<int*>(ctx);
    Pointer<PatchHierarchy<NDIM>> hierarchy = ops->getPatchHierarchy();
    HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, hierarchy->getFinestLevelNumber());
    data_ops.copyData(idx, source_idx);
}

void
check_redistribution(Pointer<PatchHierarchy<NDIM>> hierarchy,
                     Pointer<CellVariable<NDIM, double>> variable,
                     Pointer<VariableContext> context,
                     Pointer<Database> input)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int phi_idx = var_db->mapVariableAndContextToIndex(variable, context);
    int source_idx = var_db->registerClonedPatchDataIndex(variable, phi_idx);
    const int control_idx = var_db->registerClonedPatchDataIndex(variable, phi_idx);
    const int finest_ln = hierarchy->getFinestLevelNumber();
    Pointer<HierarchyMathOps> ops = new HierarchyMathOps("RedistributionMath", hierarchy, 0, finest_ln);
    const int weight_idx = ops->getCellWeightPatchDescriptorIndex();
    HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, finest_ln);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->allocatePatchData(source_idx, 0.0);
        hierarchy->getPatchLevel(ln)->allocatePatchData(control_idx, 0.0);
    }
    const std::string field_name = input->getString("field");
    RedistributionField field = RedistributionField::CONSTANT;
    if (field_name == "cutoff")
    {
        field = RedistributionField::CUTOFF;
    }
    else if (field_name == "near_cutoff")
    {
        field = RedistributionField::NEAR_CUTOFF;
    }
    else if (field_name == "steep")
    {
        field = RedistributionField::STEEP;
    }
    else if (field_name == "plane")
    {
        field = RedistributionField::PLANE;
    }
    else if (field_name == "curve")
    {
        field = RedistributionField::CURVE;
    }
    else if (field_name != "constant")
    {
        TBOX_ERROR("Unknown redistribution field: " << field_name << '\n');
    }
    const std::string scheme = input->getString("scheme");
    const double eps = std::numeric_limits<double>::epsilon();
    Pointer<Database> reference_db;
    if (input->isDatabase("Reference"))
    {
        reference_db = input->getDatabase("Reference");
        // The stored fields use Box iterator order on one serial patch.
        TBOX_ASSERT(IBTK_MPI::getNodes() == 1 && finest_ln == 0);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        TBOX_ASSERT(level->getNumberOfPatches() == 1);
    }
    Pointer<CartesianGridGeometry<NDIM>> geometry = hierarchy->getGridGeometry();
    double length = 0.0, domain_volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        const double extent = geometry->getXUpper()[d] - geometry->getXLower()[d];
        length = std::max(length, extent);
        domain_volume *= extent;
    }
    LocationIndexRobinBcCoefs<NDIM> bc;
    for (int face = 0; face < 2 * NDIM; ++face)
    {
        // The exact-distance plane has compatible normal derivatives; the
        // steep case uses zero-slope boundaries to test absent interface support.
        const double slope = field == RedistributionField::PLANE && face < 2 ? (face == 0 ? -1.0 : 1.0) : 0.0;
        bc.setBoundarySlope(face, slope);
    }
    plog << std::setprecision(16);
    for (const bool mass_constraint : { false, true })
    {
        double initial_volume = 0.0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> source = patch->getPatchData(source_idx);
                Pointer<CellData<NDIM, double>> weight = patch->getPatchData(weight_idx);
                Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
                const Box<NDIM>& box = patch->getBox();
                double dv = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dv *= geom->getDx()[d];
                }
                const double alpha = std::pow(dv, 1.0 / NDIM);
                for (Box<NDIM>::Iterator i(box); i; i++)
                {
                    const double x = geom->getXLower()[0] + (i()(0) - box.lower(0) + 0.5) * geom->getDx()[0];
                    const double y = geom->getXLower()[1] + (i()(1) - box.lower(1) + 0.5) * geom->getDx()[1];
                    double value = -1.0;
                    switch (field)
                    {
                    case RedistributionField::CONSTANT:
                        break;
                    case RedistributionField::CUTOFF:
                        value = -alpha;
                        break;
                    case RedistributionField::NEAR_CUTOFF:
                        value = -alpha * (1.0 - 1.0e-6);
                        break;
                    case RedistributionField::STEEP:
                    case RedistributionField::PLANE:
                        value = (field == RedistributionField::STEEP ? 4.0 : 1.0) * (x - 0.5);
                        break;
                    case RedistributionField::CURVE:
                        value = (1.0 + 0.2 * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y)) *
                                (std::hypot(x - 0.5, y - 0.5) - 0.23);
                        break;
                    default:
                        TBOX_ERROR("Unknown redistribution geometry\n");
                    }
                    (*source)(i()) = value;
                    initial_volume += IBTK::smooth_heaviside(-value, alpha) * (*weight)(i());
                }
            }
        }
        initial_volume = IBTK_MPI::sumReduction(initial_volume);
        for (const bool redistribution : { false, true })
        {
            Pointer<Database> db = new MemoryDatabase("RedistributionInput");
            db->putInteger("max_iterations", 1);
            db->putString("order", "THIRD_ORDER_ENO");
            db->putString("time_stepping_scheme", scheme);
            db->putBool("apply_volume_redistribution", redistribution);
            db->putBool("apply_mass_constraint", false);
            RelaxationLSMethod relaxation("Redistribution", db, false);
            relaxation.registerInterfaceNeighborhoodLocatingFcn(copy_redistribution_distance, &source_idx);
            relaxation.registerPhysicalBoundaryCondition(&bc);
            relaxation.initializeLSData(phi_idx, ops, 0, 0.0, true);
            // A subsequent reinitialization exercises the local mass constraint.
            relaxation.setApplyMassConstraint(mass_constraint);
            relaxation.setReinitializeLSData(true);
            relaxation.initializeLSData(phi_idx, ops, 1, 0.1, false);
            const std::string reference_case = std::to_string(2 * mass_constraint + redistribution);
            std::vector<double> reference_values;
            double reference_volume = 0.0, reference_norm = 0.0, field_error = 0.0;
            std::size_t reference_position = 0;
            if (reference_db)
            {
                const std::string key = "phi_" + reference_case;
                const int size = reference_db->getArraySize(key);
                reference_values.resize(size);
                reference_db->getDoubleArray(key, reference_values.data(), size);
                reference_volume = reference_db->getDouble("volume_" + reference_case);
                if (!std::isfinite(reference_volume))
                {
                    TBOX_ERROR("Redistribution regression: nonfinite H-reference volume\n");
                }
            }
            double volume = 0.0, distance_error = 0.0, difference = 0.0;
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CellData<NDIM, double>> phi = patch->getPatchData(phi_idx);
                    Pointer<CellData<NDIM, double>> control = patch->getPatchData(control_idx);
                    Pointer<CellData<NDIM, double>> weight = patch->getPatchData(weight_idx);
                    Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
                    const Box<NDIM>& box = patch->getBox();
                    double dv = 1.0;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        dv *= geom->getDx()[d];
                    }
                    const double alpha = std::pow(dv, 1.0 / NDIM);
                    for (Box<NDIM>::Iterator i(box); i; i++)
                    {
                        const double value = (*phi)(i());
                        if (!std::isfinite(value))
                        {
                            TBOX_ERROR("Redistribution regression: nonfinite reinitialized field\n");
                        }
                        if (reference_db)
                        {
                            if (reference_position >= reference_values.size())
                            {
                                TBOX_ERROR("Redistribution regression: H-reference field size mismatch\n");
                            }
                            const double expected = reference_values[reference_position++];
                            if (!std::isfinite(expected))
                            {
                                TBOX_ERROR("Redistribution regression: nonfinite H-reference field\n");
                            }
                            reference_norm = std::max(reference_norm, std::abs(expected));
                            field_error = std::max(field_error, std::abs(value - expected));
                        }
                        volume += IBTK::smooth_heaviside(-value, alpha) * (*weight)(i());
                        const double x = geom->getXLower()[0] + (i()(0) - box.lower(0) + 0.5) * geom->getDx()[0];
                        const double y = geom->getXLower()[1] + (i()(1) - box.lower(1) + 0.5) * geom->getDx()[1];
                        const double distance =
                            field == RedistributionField::CURVE ? std::hypot(x - 0.5, y - 0.5) - 0.23 : x - 0.5;
                        if (field == RedistributionField::CURVE || field == RedistributionField::PLANE ||
                            field == RedistributionField::STEEP)
                        {
                            distance_error += std::abs(value - distance) * (*weight)(i());
                        }
                        if (redistribution)
                        {
                            difference = std::max(difference, std::abs(value - (*control)(i())));
                        }
                    }
                }
            }
            volume = IBTK_MPI::sumReduction(volume);
            distance_error = IBTK_MPI::sumReduction(distance_error);
            difference = IBTK_MPI::maxReduction(difference);
            if (reference_db)
            {
                if (reference_position != reference_values.size())
                {
                    TBOX_ERROR("Redistribution regression: H-reference field size mismatch\n");
                }
                if (!(field_error <= 128.0 * eps * std::max(length, reference_norm)))
                {
                    TBOX_ERROR("Redistribution regression: H-reference field mismatch in case " << reference_case
                                                                                                << '\n');
                }
                if (!(std::abs(volume - reference_volume) <= 128.0 * eps * domain_volume))
                {
                    TBOX_ERROR("Redistribution regression: H-reference volume mismatch in case " << reference_case
                                                                                                 << '\n');
                }
            }
            if (field == RedistributionField::CURVE && redistribution && !(difference > 128.0 * eps))
            {
                TBOX_ERROR("Redistribution regression: resolved correction is inactive\n");
            }
            if ((field == RedistributionField::CONSTANT || field == RedistributionField::PLANE) && redistribution &&
                difference > 128.0 * eps)
            {
                TBOX_ERROR("Redistribution regression: absent-source correction changed the field\n");
            }
            if (field == RedistributionField::PLANE && !(distance_error <= 128.0 * eps))
            {
                TBOX_ERROR("Redistribution regression: exact-distance plane changed\n");
            }
            const bool has_distance = field == RedistributionField::CURVE || field == RedistributionField::PLANE ||
                                      field == RedistributionField::STEEP;
            plog << "mass constraint = " << mass_constraint << "; redistribution = " << redistribution
                 << "; volume, drift";
            if (has_distance)
            {
                plog << ", distance L1";
            }
            plog << ", active difference: " << volume << ' ' << volume - initial_volume << ' ';
            if (has_distance)
            {
                plog << distance_error << ' ';
            }
            plog << difference << '\n';
            if (reference_db)
            {
                plog << "normalized H-reference field, volume errors: "
                     << field_error / std::max(length, reference_norm) << ' '
                     << std::abs(volume - reference_volume) / domain_volume << '\n';
            }
            if (!redistribution)
            {
                data_ops.copyData(control_idx, phi_idx);
            }
        }
    }
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        hierarchy->getPatchLevel(ln)->deallocatePatchData(control_idx);
        hierarchy->getPatchLevel(ln)->deallocatePatchData(source_idx);
    }
    var_db->removePatchDataIndex(control_idx);
    var_db->removePatchDataIndex(source_idx);
}
} // namespace

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
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<AdvectorPredictorCorrectorHyperbolicPatchOps> hyp_patch_ops =
            new AdvectorPredictorCorrectorHyperbolicPatchOps(
                "AdvectorPredictorCorrectorHyperbolicPatchOps",
                app_initializer->getComponentDatabase("AdvectorPredictorCorrectorHyperbolicPatchOps"),
                explicit_predictor,
                grid_geometry);
        Pointer<HyperbolicLevelIntegrator<NDIM>> hyp_level_integrator =
            new HyperbolicLevelIntegrator<NDIM>("HyperbolicLevelIntegrator",
                                                app_initializer->getComponentDatabase("HyperbolicLevelIntegrator"),
                                                hyp_patch_ops,
                                                true,
                                                using_refined_timestepping);
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               hyp_level_integrator,
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
        Pointer<TimeRefinementIntegrator<NDIM>> time_integrator =
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
        Pointer<FaceVariable<NDIM, double>> u_var = new FaceVariable<NDIM, double>("u");
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
        Pointer<CellVariable<NDIM, double>> Q_var = new CellVariable<NDIM, double>("Q");
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
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit) hyp_patch_ops->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializeHierarchy();

        if (input_db->isDatabase("RedistributionTests"))
        {
            check_redistribution(patch_hierarchy,
                                 Q_var,
                                 hyp_level_integrator->getCurrentContext(),
                                 input_db->getDatabase("RedistributionTests"));
            return 0;
        }

        // Create inital level set
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        input_db->getDoubleArray("X0", circle.X0.data(), NDIM);

        Pointer<VariableContext> current_ctx = hyp_level_integrator->getCurrentContext();
        Pointer<VariableContext> scratch_ctx = hyp_level_integrator->getScratchContext();
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, current_ctx);
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, scratch_ctx);

        Pointer<CellVariable<NDIM, double>> E_var = new CellVariable<NDIM, double>("E");
        const int E_idx = var_db->registerVariableAndContext(E_var, scratch_ctx);

        // Heaviside
        Pointer<CellVariable<NDIM, double>> H_var = new CellVariable<NDIM, double>("H");
        const int H_idx = var_db->registerVariableAndContext(H_var, scratch_ctx);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
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
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double>> E_data = patch->getPatchData(E_idx);
                Pointer<CellData<NDIM, double>> Q_scratch_data = patch->getPatchData(Q_scratch_idx);
                Pointer<CellData<NDIM, double>> H_data = patch->getPatchData(H_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());

                    // Get physical coordinates
                    IBTK::Vector coord = IBTK::Vector::Zero();
                    Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
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
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double>> D_data = patch->getPatchData(Q_scratch_idx);
                Pointer<CellData<NDIM, double>> E_data = patch->getPatchData(E_idx);
                Pointer<CellData<NDIM, double>> W_data = patch->getPatchData(wgt_cc_idx);
                for (Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    CellIndex<NDIM> ci(it());
                    const double phi = (*D_data)(ci);
                    const double err = (*E_data)(ci);
                    const double dV = (*W_data)(ci);
                    Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
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
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
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
