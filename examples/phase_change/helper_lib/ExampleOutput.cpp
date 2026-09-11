// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <tbox/MathUtilities.h>
#include <tbox/RestartManager.h>
#include <tbox/TimerManager.h>

#include "ExampleOutput.h"

#include <ibamr/app_namespaces.h>

namespace PhaseChangeExamples
{
ExampleOutput::ExampleOutput(Pointer<AppInitializer> app_initializer)
{
    // Get various standard options set in the input file.
    d_dump_viz_data = app_initializer->dumpVizData();
    d_viz_dump_interval = app_initializer->getVizDumpInterval();
    d_uses_visit = d_dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

    d_dump_restart_data = app_initializer->dumpRestartData();
    d_restart_dump_interval = app_initializer->getRestartDumpInterval();
    d_restart_dump_dirname = app_initializer->getRestartDumpDirectory();

    const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
    const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
    const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
    if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
    {
        Utilities::recursiveMkdir(postproc_data_dump_dirname);
    }

    d_dump_timer_data = app_initializer->dumpTimerData();
    d_timer_dump_interval = app_initializer->getTimerDumpInterval();
}

void
ExampleOutput::registerDataWriter(Pointer<AppInitializer> app_initializer, Pointer<HierarchyIntegrator> time_integrator)
{
    // Set up visualization plot file writers.
    d_visit_data_writer = app_initializer->getVisItDataWriter();
    if (d_uses_visit)
    {
        time_integrator->registerVisItDataWriter(d_visit_data_writer);
    }
}

void
ExampleOutput::writeInitial(Pointer<HierarchyIntegrator> time_integrator,
                            Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                            const int iteration_num,
                            const double loop_time)
{
    if (d_dump_viz_data && d_uses_visit)
    {
        pout << "\n\nWriting visualization files...\n\n";
        time_integrator->setupPlotData();
        d_visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
    }
}

void
ExampleOutput::writeStep(Pointer<HierarchyIntegrator> time_integrator,
                         Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                         const int iteration_num,
                         const double loop_time)
{
    // At specified intervals, write visualization and restart files,
    // print out timer data, and store hierarchy data for post
    // processing.
    const bool last_step = !time_integrator->stepsRemaining();
    if (d_dump_viz_data && d_uses_visit && (iteration_num % d_viz_dump_interval == 0 || last_step))
    {
        pout << "\nWriting visualization files...\n\n";
        time_integrator->setupPlotData();
        d_visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
    }
    if (d_dump_restart_data && (iteration_num % d_restart_dump_interval == 0 || last_step))
    {
        pout << "\nWriting restart files...\n\n";
        RestartManager::getManager()->writeRestartFile(d_restart_dump_dirname, iteration_num);
    }
    if (d_dump_timer_data && (iteration_num % d_timer_dump_interval == 0 || last_step))
    {
        pout << "\nWriting timer data...\n\n";
        TimerManager::getManager()->print(plog);
    }
}

void
run_time_loop(Pointer<HierarchyIntegrator> time_integrator,
              Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
              ExampleOutput& output,
              const std::function<void(double)>& postprocess)
{
    int iteration_num = time_integrator->getIntegratorStep();
    double loop_time = time_integrator->getIntegratorTime();
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

        if (postprocess)
        {
            postprocess(loop_time);
        }

        output.writeStep(time_integrator, patch_hierarchy, iteration_num + 1, loop_time);
    }
}
} // namespace PhaseChangeExamples
