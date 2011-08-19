// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_prefix_config.h>
#include <IBTK_prefix_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// C++ namespace delcarations
#include <ibamr/namespaces.h>
using namespace IBAMR;
using namespace IBTK;
using namespace SAMRAI;
using namespace std;
using namespace blitz;

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
main(
    int argc,
    char* argv[])
{
    if (argc == 1)
    {
        pout << "USAGE:  " << argv[0] << " <input filename> <restart dir> <restore number> [options]\n"
             << "OPTIONS: PETSc command line options; use -help for more information\n";
        return -1;
    }

    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Process command line options.
        const string input_filename = argv[1];
        string restart_read_dirname;
        int restore_num = 0;
        bool is_from_restart = false;
        if (argc >= 4)
        {
            // Check whether this appears to be a restarted run.
            FILE* fstream = (SAMRAI_MPI::getRank() == 0 ? fopen(argv[2], "r") : NULL);
            if (SAMRAI_MPI::bcast(fstream != NULL ? 1 : 0, 0) == 1)
            {
                restart_read_dirname = argv[2];
                restore_num = atoi(argv[3]);
                is_from_restart = true;
            }
            if (fstream != NULL)
            {
                fclose(fstream);
            }
        }

        // Create input database and parse all data in input file.
        Pointer<Database> input_db = new InputDatabase("input_db");
        InputManager::getManager()->parseInputFile(input_filename, input_db);

        // Process "Main" section of the input database.
        Pointer<Database> main_db = input_db->getDatabase("Main");

        // Configure logging options.
        const string log_file_name = main_db->getStringWithDefault("log_file_name","IBAMR.log");
        const bool log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",false);
        if (log_all_nodes)
        {
            PIO::logAllNodes(log_file_name);
        }
        else
        {
            PIO::logOnlyNodeZero(log_file_name);
        }

        // Configure visualization options.
        const int viz_dump_interval = main_db->getIntegerWithDefault("viz_dump_interval",0);
        const bool viz_dump_data = viz_dump_interval > 0;
        string viz_dump_dirname;
        bool uses_visit = false;
        int visit_number_procs_per_file = 1;
        if (viz_dump_data)
        {
            Array<string> viz_writer;
            if (main_db->keyExists("viz_writer"))
            {
                viz_writer = main_db->getStringArray("viz_writer");
            }
            for (int i = 0; i < viz_writer.getSize(); i++)
            {
                if (viz_writer[i] == "VisIt") uses_visit = true;
            }

            if (main_db->keyExists("viz_dump_dirname"))
            {
                viz_dump_dirname = main_db->getString("viz_dump_dirname");
                if (viz_dump_dirname.empty())
                {
                    TBOX_ERROR("viz_dump_interval > 0, but `viz_dump_dirname' is empty\n");
                }
            }
            else
            {
                TBOX_ERROR("viz_dump_interval > 0, but key `viz_dump_dirname' not specifed in input file\n");
            }

            if (uses_visit)
            {
                visit_number_procs_per_file = main_db->getIntegerWithDefault("visit_number_procs_per_file",visit_number_procs_per_file);
            }
        }

        // Configure restart options.
        const int restart_interval = main_db->getIntegerWithDefault("restart_interval",0);
        const bool write_restart = restart_interval > 0;
        string restart_write_dirname;
        if (write_restart)
        {
            if (main_db->keyExists("restart_write_dirname"))
            {
                restart_write_dirname = main_db->getString("restart_write_dirname");
                if (restart_write_dirname.empty())
                {
                    TBOX_ERROR("restart_interval > 0, but `restart_write_dirname' is empty\n");
                }
            }
            else
            {
                TBOX_ERROR("restart_interval > 0, but key `restart_write_dirname' not specifed in input file\n");
            }
        }

        // Configure timing options.
        const int timer_dump_interval = main_db->getIntegerWithDefault("timer_dump_interval",0);
        const bool write_timer_data = timer_dump_interval > 0;
        if (write_timer_data)
        {
            TimerManager::createManager(input_db->getDatabase("TimerManager"));
        }

        // Process restart data if this is a restarted run.
        if (is_from_restart)
        {
            RestartManager::getManager()->openRestartFile(
                restart_read_dirname, restore_num, SAMRAI_MPI::getNodes());
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", input_db->getDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<INSHierarchyIntegrator> time_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator", input_db->getDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", time_integrator, input_db->getDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", input_db->getDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", input_db->getDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", input_db->getDatabase("PressureInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);
        time_integrator->registerPressureInitialConditions(p_init);

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        const bool periodic_domain = periodic_shift.min() != 0;
        TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> u_bc_coefs;
        if (!periodic_domain)
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
                    bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry);
            }
            time_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create body force function specification objects (when necessary).
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", input_db->getDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Create source function specification objects (when necessary).
        if (input_db->keyExists("SourceSinkFunction"))
        {
            Pointer<CartGridFunction> q_fcn = new muParserCartGridFunction(
                "q_fcn", input_db->getDatabase("SourceSinkFunction"), grid_geometry);
            time_integrator->registerFluidSourceFunction(q_fcn);
        }

        // Set up visualization plot file writer.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer;
        if (uses_visit)
        {
            visit_data_writer = new VisItDataWriter<NDIM>(
                "VisItDataWriter", viz_dump_dirname, visit_number_procs_per_file);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Close the restart manager.
        RestartManager::getManager()->closeRestartFile();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (viz_dump_data)
        {
            if (uses_visit)
            {
                pout << "Writing visualization files...\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, loop_time, iteration_num);
            }
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout <<                                                    "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";

            dt = time_integrator->getTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout <<                                                    "\n";
            pout << "At end       of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout <<                                                    "\n";

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            if (viz_dump_data &&
                (iteration_num%viz_dump_interval == 0 || !time_integrator->stepsRemaining()))
            {
                if (uses_visit)
                {
                    pout << "\nWriting visualization files...\n\n";
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
            }
            if (write_restart &&
                (iteration_num%restart_interval == 0 || !time_integrator->stepsRemaining()))
            {
                pout << "\nWriting restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_write_dirname, iteration_num);
            }
            if (write_timer_data &&
                (iteration_num%timer_dump_interval == 0 || !time_integrator->stepsRemaining()))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        const Pointer<SideVariable<NDIM,double> > u_var = time_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();

        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);

        const Pointer<CellVariable<NDIM,double> > p_var = time_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();

        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx, loop_time);
        }

        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time-0.5*dt);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
        pout << "Error in u at time " << loop_time << ":\n"
             << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_idx,wgt_sc_idx)  << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx,wgt_sc_idx)  << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx,wgt_sc_idx) << "\n";

        HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        pout << "Error in p at time " << loop_time-0.5*dt << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx,wgt_cc_idx)  << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx,wgt_cc_idx)  << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(p_idx,wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (viz_dump_data)
        {
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num+1, loop_time);
            }
        }

        // Cleanup boundary condition specification objects (when necessary).
        if (!periodic_domain)
        {
            for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        }

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}// main
