// Config files
#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// Headers for IO routines for file locking
#include <fcntl.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <PatchLevel.h>
#include <VariableDatabase.h>
#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/MathUtilities.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// Headers for major algorithm/data structure objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBImplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibtk/LagSiloDataWriter.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/PETScVecOps.h>
#include <ibtk/muParserDataSetter.h>
#include <ibtk/muParserRobinBcCoefs.h>

using namespace IBAMR;
using namespace IBTK;
using namespace SAMRAI;
using namespace std;

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
                           << " not specifed in input file");
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

        /*
         * Create the lock file.
         */
        string lock_file_name;
        if (main_db->keyExists("lock_file_name"))
        {
            lock_file_name = main_db->getString("lock_file_name");
        }

        const bool enable_lock_file = !lock_file_name.empty();

        if (enable_lock_file && tbox::SAMRAI_MPI::getRank() == 0)
        {
            int fd = open(lock_file_name.c_str(), O_WRONLY | O_CREAT | O_EXCL);
            if (fd < 0)
            {
                tbox::pout << "ERROR: could not create lock file: " << lock_file_name << endl;
                tbox::SAMRAI_MPI::abort();
                return -1;
            }
        }

        /*
         * Get the restart manager and root restart database.  If run is from
         * restart, open the restart file.
         */
        tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
        if (is_from_restart)
        {
            restart_manager->openRestartFile(
                restart_read_dirname, restore_num, tbox::SAMRAI_MPI::getNodes());
        }

        /*
         * Create major algorithm and data objects which comprise the
         * application.  Each object will be initialized either from input data
         * or restart files, or a combination of both.  Refer to each class
         * constructor for details.  For more information on the composition of
         * objects for this application, see comments at top of file.
         */
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new geom::CartesianGridGeometry<NDIM>(
                "CartesianGeometry",
                input_db->getDatabase("CartesianGeometry"));

        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new hier::PatchHierarchy<NDIM>(
                "PatchHierarchy",
                grid_geometry);
#if 0
        tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator =
            new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                input_db->getDatabase("INSStaggeredHierarchyIntegrator"),
                patch_hierarchy);
#endif
        tbox::Pointer<IBSpringForceGen> spring_force_generator = new IBSpringForceGen();
        tbox::Pointer<IBBeamForceGen> beam_force_generator = new IBBeamForceGen();
        tbox::Pointer<IBTargetPointForceGen> target_point_force_generator = new IBTargetPointForceGen();
        tbox::Pointer<IBStandardForceGen> force_generator = new IBStandardForceGen(
            spring_force_generator, beam_force_generator, target_point_force_generator);

        tbox::Pointer<IBImplicitHierarchyIntegrator> time_integrator =
            new IBImplicitHierarchyIntegrator(
                "IBStaggeredHierarchyIntegrator",
                input_db->getDatabase("IBStaggeredHierarchyIntegrator"),
                patch_hierarchy, /*navier_stokes_integrator,*/ force_generator);

        tbox::Pointer<IBStandardInitializer> initializer =
            new IBStandardInitializer(
                "IBStandardInitializer",
                input_db->getDatabase("IBStandardInitializer"));
        time_integrator->registerLNodeInitStrategy(initializer);

        tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
            new mesh::StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize",
                time_integrator,
                input_db->getDatabase("StandardTagAndInitialize"));

        tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator =
            new mesh::BergerRigoutsos<NDIM>();

        tbox::Pointer<mesh::LoadBalancer<NDIM> > load_balancer =
            new mesh::LoadBalancer<NDIM>(
                "LoadBalancer",
                input_db->getDatabase("LoadBalancer"));

        tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
            new mesh::GriddingAlgorithm<NDIM>(
                "GriddingAlgorithm",
                input_db->getDatabase("GriddingAlgorithm"),
                error_detector, box_generator, load_balancer);

        /*
         * Create initial condition specification objects.
         */
        tbox::Pointer<SetDataStrategy> u_init = new muParserDataSetter(
            "u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);

        /*
         * Create boundary condition specification objects (when necessary).
         */
        const hier::IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        const bool periodic_domain = periodic_shift.min() != 0;

        vector<solv::RobinBcCoefStrategy<NDIM>*> u_bc_coefs;
        if (!periodic_domain)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs.push_back(
                    new muParserRobinBcCoefs(
                        bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry));
            }
            time_integrator->registerVelocityPhysicalBcCoefs(u_bc_coefs);
        }

        /*
         * Setup parallel summation options.
         */
        LEInteractor::setFromDatabase(input_db->getDatabase("LEInteractor"));
        LEInteractor::printClassData(tbox::pout);
        PETScVecOps::setFromDatabase(input_db->getDatabase("PETScVecOps"));
        PETScVecOps::printClassData(tbox::pout);

        /*
         * Set up visualization plot file writer.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
            new appu::VisItDataWriter<NDIM>(
                "VisIt Writer",
                visit_dump_dirname, visit_number_procs_per_file);
        tbox::Pointer<LagSiloDataWriter> silo_data_writer =
            new LagSiloDataWriter(
                "LagSiloDataWriter",
                visit_dump_dirname);

        if (uses_visit)
        {
            initializer->registerLagSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerLagSiloDataWriter(silo_data_writer);
        }

        /*
         * Initialize hierarchy configuration and data on all patches.  Then,
         * close restart file and write initial state for visualization.
         */
        time_integrator->initializeHierarchyIntegrator(gridding_algorithm);
        double dt_now = time_integrator->initializeHierarchy();
        tbox::RestartManager::getManager()->closeRestartFile();

        /*
         * Deallocate the Lagrangian initializer, as it is no longer needed.
         */
        time_integrator->freeLNodeInitStrategy();
        initializer.setNull();

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
        if (viz_dump_data)
        {
            if (uses_visit)
            {
                tbox::pout << "\nWriting visualization files...\n\n";
                visit_data_writer->writePlotData(
                    patch_hierarchy,
                    time_integrator->getIntegratorStep(),
                    time_integrator->getIntegratorTime());
                silo_data_writer->writePlotData(
                    time_integrator->getIntegratorStep(),
                    time_integrator->getIntegratorTime());
            }
        }

#if 0
        /*
         * Open files to output the lift and drag coefficients.
         */
        const double radius = input_db->getDouble("R");
        ofstream drag_stream, lift_stream;
        if (tbox::SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.open("C_D.curve", tbox::RestartManager::getManager()->isFromRestart() ? ios::app : ios::out);
            lift_stream.open("C_L.curve", tbox::RestartManager::getManager()->isFromRestart() ? ios::app : ios::out);

            drag_stream.setf(ios_base::scientific);
            drag_stream.setf(ios_base::showpos);
            drag_stream.setf(ios_base::showpoint);
            drag_stream.width(16); drag_stream.precision(15);

            lift_stream.setf(ios_base::scientific);
            lift_stream.setf(ios_base::showpos);
            lift_stream.setf(ios_base::showpoint);
            lift_stream.width(16); lift_stream.precision(15);

            if (!tbox::RestartManager::getManager()->isFromRestart())
            {
                drag_stream << 0.0 << " " << 0.0 << endl;
                lift_stream << 0.0 << " " << 0.0 << endl;
            }
        }
#endif

        /*
         * Time step loop.  Note that the step count and integration time are
         * maintained by the time integrator object.
         */
        double loop_time = time_integrator->getIntegratorTime();
        double loop_time_end = time_integrator->getEndTime();

        int iteration_num = time_integrator->getIntegratorStep();

        while (!tbox::MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep() + 1;

            tbox::pout <<                                                       endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout << "At begining of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                 << endl;

            double dt_new = time_integrator->advanceHierarchy(dt_now);

            loop_time += dt_now;
            dt_now = dt_new;

            tbox::pout <<                                                       endl;
            tbox::pout << "At end      of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                 << endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            tbox::pout <<                                                       endl;
#if 0
            /*
             * Compute the drag and lift coefficients by integrating the
             * components of the Lagrangian force field over the computational
             * domain.
             */
            const int ln = patch_hierarchy->getFinestLevelNumber();
            LDataManager* lag_manager = time_integrator->getLDataManager();
            tbox::Pointer<LNodeLevelData> X_data = lag_manager->getLNodeLevelData("X",ln);
            tbox::Pointer<LNodeLevelData> F_data = lag_manager->createLNodeLevelData("F",ln,NDIM);
            force_generator->computeLagrangianForce(
                F_data, X_data,
                patch_hierarchy, ln, loop_time, lag_manager);

            double F_D = 0.0;
            double F_L = 0.0;
            for (int i = 0; i < F_data->getLocalNodeCount(); ++i)
            {
                F_D -= (*F_data)(i,0);
                F_L -= (*F_data)(i,1);
            }

            F_D = tbox::SAMRAI_MPI::sumReduction(F_D);
            F_L = tbox::SAMRAI_MPI::sumReduction(F_L);

            /*
             * Output the normalized drag and lift coefficients.
             *
             * NOTE: We assume that rho = 1.0, u_oo = 1.0, so:
             *      C_D = F_D/(rho u_oo^2 R) = F_D/R
             *      C_L = F_L/(rho u_oo^2 R) = F_L/R
             */
            if (tbox::SAMRAI_MPI::getRank() == 0)
            {
                drag_stream << loop_time << " " << F_D/radius << endl;
                lift_stream << loop_time << " " << F_L/radius << endl;
            }
#endif
            /*
             * At specified intervals, output timer data and write visualization
             * and restart and files.
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
                    visit_data_writer->writePlotData(
                        patch_hierarchy, iteration_num, loop_time);
                    silo_data_writer->writePlotData(
                        iteration_num, loop_time);
                }
            }

            if (write_restart && iteration_num%restart_interval == 0)
            {
                tbox::pout << "\nWriting restart files...\n\n";
                tbox::RestartManager::getManager()->writeRestartFile(
                    restart_write_dirname, iteration_num);

                if (stop_after_writing_restart) break;
            }
        }

        /*
         * Ensure the final timer data is written out.
         */
        if (write_timer_data && iteration_num%timer_dump_interval != 0)
        {
            tbox::pout << "\nWriting timer data...\n\n";
            tbox::TimerManager::getManager()->print(tbox::plog);
        }

        /*
         * Ensure the last state is written out.
         */
        if (viz_dump_data && iteration_num%viz_dump_interval != 0)
        {
            if (uses_visit)
            {
                tbox::pout << "\nWriting visualization files...\n\n";
                visit_data_writer->writePlotData(
                    patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(
                    iteration_num, loop_time);
            }
        }

#if 0
        /*
         * Close the files used to store the lift and drag coefficients.
         */
        if (tbox::SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.close();
            lift_stream.close();
        }
#endif

        /*
         * Delete the lock file.
         */
        if (enable_lock_file && tbox::SAMRAI_MPI::getRank() == 0)
        {
            if (remove(lock_file_name.c_str()) != 0)
            {
                tbox::pout << "ERROR: could not remove lock file: " << lock_file_name << endl;
                tbox::SAMRAI_MPI::abort();
                return -1;
            }
        }

    }// cleanup all smart Pointers prior to shutdown

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main
