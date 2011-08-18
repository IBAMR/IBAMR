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
#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/LSiloDataWriter.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

using namespace IBAMR;
using namespace IBTK;
using namespace SAMRAI;
using namespace std;

// Function prototypes
void
output_data(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy,
    tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
    LDataManager* l_data_manager,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname);

/************************************************************************
 * For each run, the input filename and restart information (if         *
 * needed) must be given on the command line.  For non-restarted case,  *
 * command line is:                                                     *
 *                                                                      *
 *    executable <input file name>                                      *
 *                                                                      *
 * For restarted run, command line is:                                  *
 *                                                                      *
 *    executable <input file name> <restart directory> <restart number> *
 *                                                                      *
 ************************************************************************
 */

int
main(
    int argc,
    char* argv[])
{
    /*
     * Initialize PETSc, MPI, and SAMRAI.
     */
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    tbox::SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    tbox::SAMRAIManager::setMaxNumberPatchDataEntries(1024);
    tbox::SAMRAIManager::startup();

    {// cleanup all smart Pointers prior to shutdown

        /*
         * Process command line and enable logging.
         */
        string input_filename;
        string restart_read_dirname;
        int restore_num = 0;

        bool is_from_restart = false;

        if (argc == 1)
        {
            tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                       << "<restart dir> <restore number> [options]\n"
                       << "  options:\n"
                       << "  PETSc command line options; use -help for more information"
                       << endl;
            tbox::SAMRAI_MPI::abort();
            return -1;
        }
        else
        {
            input_filename = argv[1];
            if (argc >= 4)
            {
                FILE* fstream = NULL;
                if (tbox::SAMRAI_MPI::getRank() == 0)
                {
                    fstream = fopen(argv[2], "r");
                }
                int worked = (fstream ? 1 : 0);
#ifdef HAVE_MPI
                worked = tbox::SAMRAI_MPI::bcast(worked, 0);
#endif
                if (worked)
                {
                    restart_read_dirname = argv[2];
                    restore_num = atoi(argv[3]);
                    is_from_restart = true;
                }
                if (fstream) fclose(fstream);
            }
        }

        tbox::plog << "input_filename = " << input_filename << endl;
        tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
        tbox::plog << "restore_num = " << restore_num << endl;

        /*
         * Create input database and parse all data in input file.
         */
        tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        /*
         * Retrieve "Main" section of the input database.  First, read dump
         * information, which is used for writing plot files.  Second, if proper
         * restart information was given on command line, and the restart
         * interval is non-zero, create a restart database.
         */
        tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

        string log_file_name = "IB.log";
        if (main_db->keyExists("log_file_name"))
        {
            log_file_name = main_db->getString("log_file_name");
        }
        bool log_all_nodes = false;
        if (main_db->keyExists("log_all_nodes"))
        {
            log_all_nodes = main_db->getBool("log_all_nodes");
        }
        if (log_all_nodes)
        {
            tbox::PIO::logAllNodes(log_file_name);
        }
        else
        {
            tbox::PIO::logOnlyNodeZero(log_file_name);
        }

        int viz_dump_interval = 0;
        if (main_db->keyExists("viz_dump_interval"))
        {
            viz_dump_interval = main_db->getInteger("viz_dump_interval");
        }

        tbox::Array<string> viz_writer(1);
        viz_writer[0] = "VisIt";
        string viz_dump_filename;
        string visit_dump_dirname;
        bool uses_visit = false;
        int visit_number_procs_per_file = 1;
        if (viz_dump_interval > 0)
        {
            if (main_db->keyExists("viz_writer"))
            {
                viz_writer = main_db->getStringArray("viz_writer");
            }
            if (main_db->keyExists("viz_dump_filename"))
            {
                viz_dump_filename = main_db->getString("viz_dump_filename");
            }
            string viz_dump_dirname;
            if (main_db->keyExists("viz_dump_dirname"))
            {
                viz_dump_dirname = main_db->getString("viz_dump_dirname");
            }
            for (int i = 0; i < viz_writer.getSize(); i++)
            {
                if (viz_writer[i] == "VisIt") uses_visit = true;
            }
            if (uses_visit)
            {
                visit_dump_dirname = viz_dump_dirname;
            }
            else
            {
                TBOX_ERROR("main(): "
                           << "\nUnrecognized 'viz_writer' entry..."
                           << "\nOptions are 'Vizamrai' and/or 'VisIt'"
                           << endl);
            }
            if (uses_visit)
            {
                if (viz_dump_dirname.empty())
                {
                    TBOX_ERROR("main(): "
                               << "\nviz_dump_dirname is null ... "
                               << "\nThis must be specified for use with VisIt"
                               << endl);
                }
                if (main_db->keyExists("visit_number_procs_per_file"))
                {
                    visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
                }
            }
        }

        const bool viz_dump_data = (viz_dump_interval > 0);

        int restart_interval = 0;
        if (main_db->keyExists("restart_interval"))
        {
            restart_interval = main_db->getInteger("restart_interval");
        }

        string restart_write_dirname;
        if (restart_interval > 0)
        {
            if (main_db->keyExists("restart_write_dirname"))
            {
                restart_write_dirname = main_db->getString("restart_write_dirname");
            }
            else
            {
                TBOX_ERROR("restart_interval > 0, but key `restart_write_dirname'"
                           << " not specified in input file");
            }
        }

        const bool write_restart = restart_interval > 0
            && !restart_write_dirname.empty();

        bool stop_after_writing_restart = false;
        if (write_restart)
        {
            if (main_db->keyExists("stop_after_writing_restart"))
            {
                stop_after_writing_restart = main_db->getBool("stop_after_writing_restart");
            }
        }

        int data_dump_interval = 0;
        if (main_db->keyExists("data_dump_interval"))
        {
            data_dump_interval = main_db->getInteger("data_dump_interval");
        }

        string data_dump_dirname;
        if (data_dump_interval > 0)
        {
            if (main_db->keyExists("data_dump_dirname"))
            {
                data_dump_dirname = main_db->getString("data_dump_dirname");
            }
            else
            {
                TBOX_ERROR("data_dump_interval > 0, but key `data_dump_dirname'"
                           << " not specified in input file");
            }
        }

        const bool write_data = (data_dump_interval > 0) && !(data_dump_dirname.empty());
        if (write_data)
        {
            tbox::Utilities::recursiveMkdir(data_dump_dirname);
        }

        int timer_dump_interval = 0;
        if (main_db->keyExists("timer_dump_interval"))
        {
            timer_dump_interval = main_db->getInteger("timer_dump_interval");
        }

        const bool write_timer_data = (timer_dump_interval > 0);
        if (write_timer_data)
        {
            tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));
        }

        int postprocess_interval = 0;
        if (main_db->keyExists("postprocess_interval"))
        {
            postprocess_interval = main_db->getInteger("postprocess_interval");
        }

        const bool postprocess_data = (postprocess_interval > 0);

        /*
         * Get the restart manager and root restart database.  If run is from
         * restart, open the restart file.
         */
        tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
        if (is_from_restart)
        {
            restart_manager->openRestartFile(restart_read_dirname, restore_num, tbox::SAMRAI_MPI::getNodes());
        }

        /*
         * Create major algorithm and data objects which comprise the
         * application.  Each object will be initialized either from input data
         * or restart files, or a combination of both.  Refer to each class
         * constructor for details.  For more information on the composition of
         * objects for this application, see comments at top of file.
         */
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = new geom::CartesianGridGeometry<NDIM>("CartesianGeometry", input_db->getDatabase("CartesianGeometry"));
        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy = new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator("INSStaggeredHierarchyIntegrator", input_db->getDatabase("INSStaggeredHierarchyIntegrator"));
        tbox::Pointer<IBHierarchyIntegrator> time_integrator = new IBHierarchyIntegrator("IBHierarchyIntegrator", input_db->getDatabase("IBHierarchyIntegrator"), navier_stokes_integrator);
        tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector = new mesh::StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", time_integrator, input_db->getDatabase("StandardTagAndInitialize"));
        tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator = new mesh::BergerRigoutsos<NDIM>();
        tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer = new mesh::LoadBalancer<NDIM>("LoadBalancer", input_db->getDatabase("LoadBalancer"));
        tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm = new mesh::GriddingAlgorithm<NDIM>("GriddingAlgorithm", input_db->getDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        /*
         * Setup IB initializaton.
         */
        tbox::Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer("IBStandardInitializer", input_db->getDatabase("IBStandardInitializer"));
        time_integrator->registerLInitStrategy(ib_initializer);

        /*
         * Setup IB forcing.
         */
        tbox::Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        time_integrator->registerIBLagrangianForceFunction(ib_force_fcn);

        /*
         * Create Eulerian initial condition specification objects.
         */
        tbox::Pointer<CartGridFunction> u_init = new muParserCartGridFunction("u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        /*
         * Create Eulerian boundary condition specification objects (when necessary).
         */
        const hier::IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        const bool periodic_domain = periodic_shift.min() != 0;
        blitz::TinyVector<solv::RobinBcCoefStrategy<NDIM>*,NDIM> u_bc_coefs;
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

                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        /*
         * Create Eulerian body force function specification objects (when
         * necessary).
         */
        if (input_db->keyExists("ForcingFunction"))
        {
            tbox::Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction("f_fcn", input_db->getDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        /*
         * Set up visualization plot file writer.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer;
        tbox::Pointer<LSiloDataWriter> silo_data_writer;
        if (uses_visit)
        {
            visit_data_writer = new appu::VisItDataWriter<NDIM>("VisIt Writer", visit_dump_dirname, visit_number_procs_per_file);
            silo_data_writer = new LSiloDataWriter("LSiloDataWriter", visit_dump_dirname);
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerLSiloDataWriter(silo_data_writer);
        }

        /*
         * Initialize hierarchy configuration and data on all patches.  Then,
         * close restart file and write initial state for visualization.
         */
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        tbox::RestartManager::getManager()->closeRestartFile();

        /*
         * Deallocate the Lagrangian initializer, as it is no longer needed.
         */
        time_integrator->freeLInitStrategy();
        ib_initializer.setNull();

        /*
         * After creating all objects and initializing their state, we print the
         * input database contents to the log file.
         */
        tbox::plog << "\nCheck input data before simulation:" << endl;
        tbox::plog << "Input database..." << endl;
        input_db->printClassData(tbox::plog);

        /*
         * Write initial visualization files.
         */
        double loop_time = time_integrator->getIntegratorTime();
        int iteration_num = time_integrator->getIntegratorStep();
        if (viz_dump_data)
        {
            if (uses_visit)
            {
                tbox::pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, loop_time, iteration_num);
                silo_data_writer->writePlotData(loop_time, iteration_num);
            }
        }

        /*
         * Time step loop.  Note that the step count and integration time are
         * maintained by the time integrator object.
         */
        double loop_time_end = time_integrator->getEndTime();
        double dt_old = 0.0, dt_now = 0.0;
        if (write_data && iteration_num%data_dump_interval == 0)
        {
            output_data(patch_hierarchy,
                        navier_stokes_integrator, time_integrator->getLDataManager(),
                        iteration_num, loop_time, data_dump_dirname);
        }
        while (!tbox::MathUtilities<double>::equalEps(loop_time,loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep() + 1;

            tbox::pout <<                                                        endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            tbox::pout << "At beginning of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                  << endl;

            dt_old = dt_now;
            dt_now = time_integrator->getTimeStepSize();
            time_integrator->advanceHierarchy(dt_now);
            loop_time += dt_now;

            tbox::pout <<                                                        endl;
            tbox::pout << "At end       of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                  << endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            tbox::pout <<                                                        endl;

            /*
             * At specified intervals, write visualization and restart files,
             * print out timer data, write state data, and post process data.
             */
            if (write_timer_data && iteration_num%timer_dump_interval == 0)
            {
                tbox::pout << "\nWriting timer data...\n\n";
                tbox::TimerManager::getManager()->print(tbox::plog);
            }

            if (viz_dump_data && iteration_num%viz_dump_interval == 0)
            {
                if (uses_visit)
                {
                    tbox::pout << "\nWriting visualization files...\n\n";
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                    silo_data_writer->writePlotData(iteration_num, loop_time);
                }
            }

            if (write_restart && iteration_num%restart_interval == 0)
            {
                tbox::pout << "\nWriting restart files...\n\n";
                tbox::RestartManager::getManager()->writeRestartFile(restart_write_dirname, iteration_num);
            }

            if (write_data && iteration_num%data_dump_interval == 0)
            {
                output_data(patch_hierarchy,
                            navier_stokes_integrator, time_integrator->getLDataManager(),
                            iteration_num, loop_time, data_dump_dirname);
            }

            if (postprocess_data && iteration_num%postprocess_interval == 0)
            {
                time_integrator->postProcessData();
            }
        }

        /*
         * Ensure the final timer and visualization data are written out.
         */
        if (write_timer_data && iteration_num%timer_dump_interval != 0)
        {
            tbox::pout << "\nWriting timer data...\n\n";
            tbox::TimerManager::getManager()->print(tbox::plog);
        }

        if (viz_dump_data && iteration_num%viz_dump_interval != 0)
        {
            if (uses_visit)
            {
                tbox::pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
        }

        /*
         * Cleanup Eulerian boundary condition specification objects (when
         * necessary).
         */
        if (!periodic_domain)
        {
            for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        }
    }// cleanup all smart Pointers prior to shutdown

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}// main

void
output_data(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy,
    tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
    LDataManager* l_data_manager,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname)
{
    tbox::plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    tbox::plog << "simulation time is " << loop_time << endl;

    /*
     * Write Cartesian data.
     */
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, tbox::SAMRAI_MPI::getRank());
    file_name += temp_buf;

    tbox::Pointer<tbox::HDFDatabase> hier_db = new tbox::HDFDatabase("hier_db");
    hier_db->create(file_name);

    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    hier::ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(), navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(), navier_stokes_integrator->getCurrentContext()));

    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);

    hier_db->close();

    /*
     * Write Lagrangian data.
     */
    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
    tbox::Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    VecView(X_lag_vec, viewer);
    PetscViewerDestroy(viewer);
    VecDestroy(X_lag_vec);
    return;
}// output_data
