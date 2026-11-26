// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <FaceData.h>
#include <FaceVariable.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <SAMRAI_config.h>
#include <SideData.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

#include <cmath>
#include <limits>
#include <vector>

#include <ibamr/app_namespaces.h>

class DummyAdvDiffIntegrator : public AdvDiffHierarchyIntegrator
{
public:
    DummyAdvDiffIntegrator(std::string object_name, Pointer<Database> input_db)
        : AdvDiffHierarchyIntegrator(std::move(object_name), input_db, true)
    {
    }

    ~DummyAdvDiffIntegrator() = default;

protected:
    void integrateHierarchySpecialized(double /*current_time*/, double /*new_time*/, int cycle_num = 0) override
    {
        NULL_USE(cycle_num);
        plog << d_object_name << ": integrateHierarchySpecialized()\n";
    }
    void regridHierarchyBeginSpecialized() override
    {
        plog << d_object_name << ": regridHierarchyBeginSpecialized()\n";
    }
    void regridHierarchyEndSpecialized() override
    {
        plog << d_object_name << ": regridHierarchyEndSpecialized()\n";
    }
    double getMinimumTimeStepSizeSpecialized() override
    {
        plog << d_object_name << ": getMinimumTimeStepSizeSpecialized()\n";
        return std::numeric_limits<double>::min();
    }
    double getMaximumTimeStepSizeSpecialized() override
    {
        plog << d_object_name << ": getMaximumTimeStepSizeSpecialized()\n";
        return std::numeric_limits<double>::max();
    }
    void synchronizeHierarchyDataSpecialized(VariableContextType /*ctx_type*/) override
    {
        plog << d_object_name << ": synchronizeHierarchyDataSpecialized()\n";
    }
    void resetTimeDependentHierarchyDataSpecialized(const double /*new_time*/) override
    {
        plog << d_object_name << ": resetTimeDependentHierarchyDataSpecialized()\n";
    }
    void resetIntegratorToPreadvanceStateSpecialized() override
    {
        plog << d_object_name << ": resetIntegratorToPreadvanceStateSpecialized()\n";
    }
    bool atRegridPointSpecialized() const override
    {
        plog << d_object_name << ": atRegridPointSpecialized()\n";
        return false;
    }
    void setupPlotDataSpecialized() override
    {
        plog << d_object_name << ": setupPlotDataSpecialized()\n";
    }
    void initializeCompositeHierarchyDataSpecialized(double /*init_data_time*/, bool /*initial_time*/) override
    {
        plog << d_object_name << ": initializeCompositeHierarchyDataSpecialized()\n";
    }
    void initializeLevelDataSpecialized(Pointer<BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                        int /*level_number*/,
                                        double /*init_data_time*/,
                                        bool /*can_be_refined*/,
                                        bool /*initial_time*/,
                                        Pointer<BasePatchLevel<NDIM>> /*old_level*/,
                                        bool /*allocate_data*/) override
    {
        plog << d_object_name << ": initializeLevelDataSpecialized()\n";
    }
    void resetHierarchyConfigurationSpecialized(Pointer<BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                                int /*coarsest_level*/,
                                                int /*finest_level*/) override
    {
        plog << d_object_name << ": resetHierarchyConfigurationSpecialized()\n";
    }
    void putToDatabaseSpecialized(Pointer<Database> /*db*/) override
    {
        plog << d_object_name << ": putToDatabaseSpecialized()\n";
    }
    void addWorkloadEstimate(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/, const int /*workload_data_idx*/) override
    {
        plog << d_object_name << ": addWorkloadEstimate()\n";
    }
    void executePreprocessIntegrateHierarchyCallbackFcns(double /*current_time*/,
                                                         double /*new_time*/,
                                                         int /*num_cycles*/) override
    {
        plog << d_object_name << ": executePreprocessIntegrateHierarchyCallbackFcns()\n";
    }
};

class CheckedAdvDiffIntegrator : public AdvDiffSemiImplicitHierarchyIntegrator
{
public:
    using CheckFcn = void (*)(CheckedAdvDiffIntegrator*, const std::string&, void*);

    CheckedAdvDiffIntegrator(const std::string& object_name,
                             Pointer<Database> input_db,
                             CheckFcn check_fcn,
                             void* check_ctx,
                             const bool register_for_restart = true)
        : AdvDiffSemiImplicitHierarchyIntegrator(object_name, input_db, register_for_restart),
          d_check_fcn(check_fcn),
          d_check_ctx(check_ctx)
    {
    }

    double getAdvectionTimestepSize()
    {
        return AdvDiffHierarchyIntegrator::getMaximumTimeStepSizeSpecialized();
    }

protected:
    void initializeCompositeHierarchyDataSpecialized(double init_data_time, bool initial_time) override
    {
        AdvDiffSemiImplicitHierarchyIntegrator::initializeCompositeHierarchyDataSpecialized(init_data_time,
                                                                                            initial_time);
        if (d_check_fcn)
        {
            d_check_fcn(this,
                        initial_time ? "initializeCompositeHierarchyData(initial)" :
                                       "initializeCompositeHierarchyData(regrid)",
                        d_check_ctx);
        }
    }

private:
    CheckFcn d_check_fcn = nullptr;
    void* d_check_ctx = nullptr;
};

struct AdvDiffCheckContext
{
    Pointer<INSHierarchyIntegrator> ins_integrator;
    Pointer<hier::Variable<NDIM>> ins_vel_var;
    bool enabled = false;
};

void
check_adv_diff_velocity(CheckedAdvDiffIntegrator* adv_diff, const std::string& stage, void* raw_ctx)
{
    auto* ctx = static_cast<AdvDiffCheckContext*>(raw_ctx);
    if (!ctx || !ctx->enabled)
    {
        return;
    }

    Pointer<PatchHierarchy<NDIM>> hierarchy = adv_diff->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<hier::Variable<NDIM>> ins_vel_var = ctx->ins_vel_var;
    if (!ins_vel_var)
    {
        ins_vel_var = var_db->getVariable(ctx->ins_integrator->getName() + "::u_ADV");
    }
    const int u_adv_diff_current_idx = var_db->mapVariableAndContextToIndex(
        ctx->ins_integrator->getAdvectionVelocityVariable(), adv_diff->getCurrentContext());
    const int ins_vel_current_idx =
        var_db->mapVariableAndContextToIndex(ins_vel_var, ctx->ins_integrator->getCurrentContext());
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(u_adv_diff_current_idx))
        {
            TBOX_ERROR(adv_diff->getName()
                       << " [" << stage << "]: U_adv_diff is not allocated on level " << ln << ".\n");
        }
        if (!level->checkAllocated(ins_vel_current_idx))
        {
            TBOX_ERROR(adv_diff->getName()
                       << " [" << stage << "]: INS velocity is not allocated on level " << ln << ".\n");
        }
    }

    double u_adv_max = 0.0;
    double max_diff = 0.0;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<FaceData<NDIM, double>> u_adv_diff_data = patch->getPatchData(u_adv_diff_current_idx);
            Pointer<SideData<NDIM, double>> u_sc_data = patch->getPatchData(ins_vel_current_idx);
            Pointer<FaceData<NDIM, double>> u_fc_data = patch->getPatchData(ins_vel_current_idx);
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (FaceIterator<NDIM> fi(patch_box, axis); fi; fi++)
                {
                    const double u_adv_diff = (*u_adv_diff_data)(fi());
                    u_adv_max = std::max(u_adv_max, std::abs(u_adv_diff));
                    double u_ins = 0.0;
                    if (u_sc_data)
                    {
                        const SideIndex<NDIM> si(fi().toCell(1), axis, SideIndex<NDIM>::Lower);
                        u_ins = (*u_sc_data)(si);
                    }
                    else if (u_fc_data)
                    {
                        u_ins = (*u_fc_data)(fi());
                    }
                    else
                    {
                        TBOX_ERROR("Unsupported INS velocity centering.\n");
                    }
                    max_diff = std::max(max_diff, std::abs(u_adv_diff - u_ins));
                }
            }
        }
    }

    u_adv_max = IBTK_MPI::maxReduction(u_adv_max);
    max_diff = IBTK_MPI::maxReduction(max_diff);
    const double adv_diff_dt = adv_diff->getAdvectionTimestepSize();
    plog << adv_diff->getName() << " [" << stage << "]: u_adv_max=" << u_adv_max << ", max_diff=" << max_diff
         << ", adv_diff_dt=" << adv_diff_dt << "\n";

    if (!(u_adv_max > 0.0) || !(max_diff <= 1.0e-12) || !(adv_diff_dt < 0.5 * std::numeric_limits<double>::max()))
    {
        TBOX_ERROR(adv_diff->getName() << " [" << stage << "]: unexpected advection velocity (u_adv_max=" << u_adv_max
                                       << ", max_diff=" << max_diff << ", adv_diff_dt=" << adv_diff_dt << ").\n");
    }
}

void
generate_structure(const unsigned int& /*strct_num*/,
                   const int& /*ln*/,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* /*ctx*/)
{
    num_vertices = 1;
    vertex_posn.resize(num_vertices);
    vertex_posn[0] = Point(0.5, 0.5);
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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without silo
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    { // cleanup dynamically allocated objects prior to shutdown
      // prevent a warning about timer initialization
        TimerManager::createManager(nullptr);

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        const std::string solver_type = input_db->getStringWithDefault("SOLVER_TYPE", "STAGGERED");
        const bool use_ib_parent = input_db->getBoolWithDefault("USE_IB_PARENT", true);
        const bool check_adv_diff = input_db->getBoolWithDefault("CHECK_ADV_DIFF_VELOCITY", false);
        const bool test_regrid = input_db->getBoolWithDefault("TEST_REGRID", false);
        const bool use_second_integrator = input_db->getBoolWithDefault("USE_SECOND_INTEGRATOR", false);

        AdvDiffCheckContext check_ctx;
        check_ctx.enabled = check_adv_diff;

        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator_2;
        Pointer<CheckedAdvDiffIntegrator> checked_adv_diff_integrator;
        Pointer<CheckedAdvDiffIntegrator> checked_adv_diff_integrator_2;
        if (check_adv_diff)
        {
            const std::string adv_diff_db_name = input_db->isDatabase("AdvDiffSemiImplicitHierarchyIntegrator") ?
                                                     "AdvDiffSemiImplicitHierarchyIntegrator" :
                                                     "DummyAdvDiff";
            checked_adv_diff_integrator =
                new CheckedAdvDiffIntegrator("CheckedAdvDiff",
                                             app_initializer->getComponentDatabase(adv_diff_db_name),
                                             &check_adv_diff_velocity,
                                             &check_ctx);
            adv_diff_integrator = checked_adv_diff_integrator;
            if (use_second_integrator)
            {
                checked_adv_diff_integrator_2 =
                    new CheckedAdvDiffIntegrator("CheckedAdvDiff2",
                                                 app_initializer->getComponentDatabase(adv_diff_db_name),
                                                 &check_adv_diff_velocity,
                                                 &check_ctx);
                adv_diff_integrator_2 = checked_adv_diff_integrator_2;
            }
        }
        else
        {
            adv_diff_integrator =
                new DummyAdvDiffIntegrator("DummyAdvDiff", app_initializer->getComponentDatabase("DummyAdvDiff"));
            if (use_second_integrator)
            {
                adv_diff_integrator_2 =
                    new DummyAdvDiffIntegrator("DummyAdvDiff2", app_initializer->getComponentDatabase("DummyAdvDiff"));
            }
        }

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        Pointer<INSVCStaggeredHierarchyIntegrator> vc_ins_integrator;
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else if (solver_type == "VC_STAGGERED_NONCONSERVATIVE")
        {
            vc_ins_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
            navier_stokes_integrator = vc_ins_integrator;
        }
        else
        {
            TBOX_ERROR("Unsupported SOLVER_TYPE: " << solver_type << ".\n");
        }

        check_ctx.ins_integrator = navier_stokes_integrator;
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        if (adv_diff_integrator_2)
        {
            navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator_2);
        }

        Pointer<IBMethod> ib_method_ops;
        Pointer<IBHierarchyIntegrator> ib_time_integrator;
        Pointer<HierarchyIntegrator> time_integrator;
        if (use_ib_parent)
        {
            ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
            ib_time_integrator =
                new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                                  app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                  ib_method_ops,
                                                  navier_stokes_integrator);
            time_integrator = ib_time_integrator;
        }
        else
        {
            time_integrator = navier_stokes_integrator;
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

        Pointer<CartGridFunction> u_init;
        if (input_db->getBoolWithDefault("REGISTER_VELOCITY_INITIAL_CONDITIONS", false))
        {
            u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        Pointer<CartGridFunction> q_fcn;
        if (input_db->getBoolWithDefault("REGISTER_FLUID_SOURCE_FUNCTION", false))
        {
            q_fcn = new muParserCartGridFunction(
                "q_fcn", app_initializer->getComponentDatabase("FluidSourceFunction"), grid_geometry);
            navier_stokes_integrator->registerVelocityDivergenceFunction(q_fcn);
        }

        Pointer<CellVariable<NDIM, double>> mu_var;
        Pointer<CartGridFunction> mu_fcn;
        if (vc_ins_integrator && input_db->keyExists("ViscosityFunction"))
        {
            mu_var = new CellVariable<NDIM, double>("mu_var");
            vc_ins_integrator->registerViscosityVariable(mu_var);
            mu_fcn = new muParserCartGridFunction(
                "mu_fcn", app_initializer->getComponentDatabase("ViscosityFunction"), grid_geometry);
            vc_ins_integrator->registerViscosityInitialConditions(mu_fcn);
        }

        // Configure the IB solver when present.
        Pointer<IBRedundantInitializer> ib_initializer;
        if (ib_method_ops)
        {
            ib_initializer = new IBRedundantInitializer(
                "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
            std::vector<std::string> struct_list_vec = { "dummy" };
            ib_initializer->setStructureNamesOnLevel(input_db->getInteger("MAX_LEVELS") - 1, struct_list_vec);
            ib_initializer->registerInitStructureFunction(generate_structure);
            ib_method_ops->registerLInitStrategy(ib_initializer);
        }

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        if (solver_type != "COLLOCATED")
        {
            check_ctx.ins_vel_var = navier_stokes_integrator->getVelocityVariable();
        }
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        if (solver_type == "COLLOCATED")
        {
            check_ctx.ins_vel_var =
                VariableDatabase<NDIM>::getDatabase()->getVariable(navier_stokes_integrator->getName() + "::u_ADV");
        }
        if (checked_adv_diff_integrator)
        {
            check_adv_diff_velocity(
                checked_adv_diff_integrator.getPointer(), "post-initializePatchHierarchy", &check_ctx);
        }
        if (checked_adv_diff_integrator_2)
        {
            check_adv_diff_velocity(
                checked_adv_diff_integrator_2.getPointer(), "post-initializePatchHierarchy", &check_ctx);
        }
        if (test_regrid)
        {
            time_integrator->regridHierarchy();
            if (checked_adv_diff_integrator)
            {
                check_adv_diff_velocity(checked_adv_diff_integrator.getPointer(), "post-regrid", &check_ctx);
            }
            if (checked_adv_diff_integrator_2)
            {
                check_adv_diff_velocity(checked_adv_diff_integrator_2.getPointer(), "post-regrid", &check_ctx);
            }
        }

        // Close the restart manager.
        RestartManager::getManager()->closeRestartFile();

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
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
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
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
        }

        if (checked_adv_diff_integrator)
        {
            check_adv_diff_velocity(checked_adv_diff_integrator.getPointer(), "post-advance", &check_ctx);
        }
        if (checked_adv_diff_integrator_2)
        {
            check_adv_diff_velocity(checked_adv_diff_integrator_2.getPointer(), "post-advance", &check_ctx);
        }

    } // cleanup dynamically allocated objects prior to shutdown

} // main
