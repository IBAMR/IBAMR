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

#ifndef included_PhaseChangeExamples_ExampleOutput
#define included_PhaseChangeExamples_ExampleOutput

#include <ibamr/config.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyIntegrator.h>

#include <VisItDataWriter.h>

#include <functional>
#include <string>

namespace PhaseChangeExamples
{
class ExampleOutput
{
public:
    explicit ExampleOutput(SAMRAI::tbox::Pointer<IBTK::AppInitializer> app_initializer);
    void registerDataWriter(SAMRAI::tbox::Pointer<IBTK::AppInitializer> app_initializer,
                            SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> time_integrator);
    void writeInitial(SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> time_integrator,
                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                      int iteration_num,
                      double loop_time);
    void writeStep(SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> time_integrator,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                   int iteration_num,
                   double loop_time);

private:
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM>> d_visit_data_writer;
    bool d_dump_viz_data, d_uses_visit, d_dump_restart_data, d_dump_timer_data;
    int d_viz_dump_interval, d_restart_dump_interval, d_timer_dump_interval;
    std::string d_restart_dump_dirname;
};

void run_time_loop(SAMRAI::tbox::Pointer<IBTK::HierarchyIntegrator> time_integrator,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> patch_hierarchy,
                   ExampleOutput& output,
                   const std::function<void(double)>& postprocess = {});
} // namespace PhaseChangeExamples
#endif
