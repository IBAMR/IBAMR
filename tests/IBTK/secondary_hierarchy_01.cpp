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

// Test verifying that we can move workloads around in the expected way:
// 1. Compute a workload with the scratch hierarchy.
// 2. Move that workload back to the primary hierarchy.
// 3. Regrid the primary hierarchy according to some other criterion.
// 4. Reinitialize the scratch hierarchy and load balance it correctly based on
//    the values computed in step 1.

// Headers for major SAMRAI objects
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAIBox.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIterator.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIHierarchyCellDataOpsReal.h>
#include <SAMRAIIndex.h>
#include <SAMRAIIntVector.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariableDatabase.h>
#include <SAMRAIVisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SecondaryHierarchy.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <array>
#include <utility>

#include <ibtk/app_namespaces.h>

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input
        // database. Nearly all SAMRAI applications (at least those in IBAMR)
        // start by setting up the same half-dozen objects.
        Pointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        Pointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        Pointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        // custom workload vars for primary hierarchy (ph) and secondary hierarchy (sh)
        Pointer<SAMRAICellVariable<double>> ph_work_cc_var = new SAMRAICellVariable<double>("ph_work_cc");
        Pointer<SAMRAICellVariable<double>> sh_work_cc_var = new SAMRAICellVariable<double>("sh_work_cc");
        const int ph_work_cc_idx = var_db->registerVariableAndContext(ph_work_cc_var, ctx, SAMRAIIntVector(0));
        const int sh_work_cc_idx = var_db->registerVariableAndContext(sh_work_cc_var, ctx, SAMRAIIntVector(0));

        // setup plotting
        Pointer<SAMRAIVisItDataWriter> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);
        visit_data_writer->registerPlotQuantity(ph_work_cc_var->getName(), "SCALAR", ph_work_cc_idx);
        visit_data_writer->registerPlotQuantity(sh_work_cc_var->getName(), "SCALAR", sh_work_cc_idx);

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

        {
            const int tag_buffer = std::numeric_limits<int>::max();
            int level_number = 0;
            while (gridding_algorithm->levelCanBeRefined(level_number))
            {
                gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
                ++level_number;
            }
        }

        const int finest_level = patch_hierarchy->getFinestLevelNumber();

        // Allocate data for each variable on each level of the patch
        // hierarchy.
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(ph_work_cc_idx, 0.0);
            level->allocatePatchData(sh_work_cc_idx, 0.0);
        }
        SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy, 0, finest_level);
        hier_cc_data_ops.setToScalar(ph_work_cc_idx, 0.0);
        hier_cc_data_ops.setToScalar(sh_work_cc_idx, 0.0);

        // custom workload setup: we just need something nonuniform
        const int n_processes = IBTK_MPI::getNodes();
        const int current_rank = IBTK_MPI::getRank();
        unsigned int index = 0;
        for (int ln = 0; ln <= finest_level; ++ln)
        {
            Pointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                const SAMRAIPatch& patch = *level->getPatch(p());
                const SAMRAIBox& patch_box = patch.getBox();
                Pointer<SAMRAICellData<double>> sh_work_cc_data = patch.getPatchData(sh_work_cc_idx);
                Pointer<SAMRAICellData<double>> ph_work_cc_data = patch.getPatchData(ph_work_cc_idx);
                for (SAMRAICellIterator ic(patch_box); ic; ic++)
                {
                    const SAMRAIIndex& i = ic();
                    (*ph_work_cc_data)(i) = index;
                    (*sh_work_cc_data)(i) = current_rank == n_processes - 1 ? 100 + (index % 50) : 0;
                    ++index;
                }
            }
        }
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

        IBTK::SecondaryHierarchy secondary_hierarchy("secondary",
                                                     app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                     app_initializer->getComponentDatabase("LoadBalancer"));

        secondary_hierarchy.reinit(0, patch_hierarchy->getFinestLevelNumber(), patch_hierarchy, sh_work_cc_idx);

        visit_data_writer->writePlotData(secondary_hierarchy.getSecondaryHierarchy(), 1, 1.0);

        std::array<std::pair<Pointer<SAMRAIPatchHierarchy>, std::string>, 2> data{
            { { patch_hierarchy, "primary" }, { secondary_hierarchy.getSecondaryHierarchy(), "secondary" } }
        };

        for (const auto& pair : data)
        {
            std::vector<double> workload_per_processor(n_processes);
            SAMRAIHierarchyCellDataOpsReal<double> hier_cc_data_ops(pair.first, 0, finest_level);
            workload_per_processor[current_rank] = hier_cc_data_ops.L1Norm(sh_work_cc_idx, IBTK::invalid_index, true);

            const auto right_padding = std::size_t(std::log10(n_processes)) + 1;

            int ierr = MPI_Allreduce(MPI_IN_PLACE,
                                     workload_per_processor.data(),
                                     workload_per_processor.size(),
                                     MPI_DOUBLE,
                                     MPI_SUM,
                                     IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == 0);
            if (current_rank == 0)
            {
                for (int rank = 0; rank < n_processes; ++rank)
                {
                    pout << pair.second << " workload estimate on processor " << std::setw(right_padding) << std::left
                         << rank << " = " << long(workload_per_processor[rank]) << '\n';
                }
            }
        }
    }
} // main()
