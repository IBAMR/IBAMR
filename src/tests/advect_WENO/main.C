// Config files
#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <PatchLevel.h>
#include <VariableDatabase.h>
#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/MPI.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// Headers for major algorithm/data structure objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <MethodOfLinesIntegrator.h>
#include <PatchHierarchy.h>
#include <StandardTagAndInitialize.h>
#include <TimeRefinementIntegrator.h>
#include <TimeRefinementLevelStrategy.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include "WENOMOLOps.h"

using namespace SAMRAI;
using namespace std;

/************************************************************************
 *                                                                      *
 * This is the main program for an AMR solution of the linear           *
 * advection equation: dQ/dt + div(u*Q) = 0, where "Q" is a             *
 * scalar-valued function and "u" is a possibly time-dependent          *
 * advection velocity.                                                  *
 *                                                                      *
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
    char *argv[])
{
    /*
     * Initialize MPI and SAMRAI, enable logging, and process command line.
     */
    tbox::MPI::init(&argc, &argv);
    tbox::SAMRAIManager::startup();

    string input_filename;
    string restart_read_dirname;
    int restore_num = 0;

    bool is_from_restart = false;

    if ((argc != 2) && (argc != 4))
    {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                   << "<restart dir> <restore number> [options]\n"
                   << "  options:\n"
                   << "  none at this time"
                   << endl;
        tbox::MPI::abort();
        return (-1);
    }
    else
    {
        input_filename = argv[1];
        if (argc == 4)
        {
            TBOX_ERROR("restart not implemented at the present time.\n");
            restart_read_dirname = argv[2];
            restore_num = atoi(argv[3]);

            is_from_restart = true;
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
     * restart information was given on command line, and the restart interval
     * is non-zero, create a restart database.  Third, read in solution
     * algorithm customizations.
     */
    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    string log_file_name = "advect.log";
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
                visit_number_procs_per_file =
                    main_db->getInteger("visit_number_procs_per_file");
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
        tbox::pout << "using subcycled timestepping.\n";
    }
    else
    {
        tbox::pout << "NOT using subcycled timestepping.\n";
    }

    const bool u_is_div_free = main_db->getBoolWithDefault("u_is_div_free", false);
    if (u_is_div_free)
    {
        tbox::pout << "advection velocity u is discretely divergence free.\n";
    }
    else
    {
        tbox::pout << "advection velocity u is NOT discretely divergence free.\n";
    }

    const bool consv_form = main_db->getBoolWithDefault("consv_form", false);
    if (consv_form)
    {
        tbox::pout << "solving the advection equation in CONSERVATION form.\n";
    }
    else
    {
        tbox::pout << "solving the advection equation in NON-CONSERVATION form.\n";
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
     * Get the restart manager and root restart database.  If run is from
     * restart, open the restart file.
     */
    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
    if (is_from_restart)
    {
        restart_manager->openRestartFile(
            restart_read_dirname, restore_num, tbox::MPI::getNodes());
    }

    /*
     * Create major algorithm and data objects which comprise the application.
     * Each object will be initialized either from input data or restart files,
     * or a combination of both.  Refer to each class constructor for details.
     * For more information on the composition of objects for this application,
     * see comments at top of file.
     */
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
        new geom::CartesianGridGeometry<NDIM>(
            "CartesianGeometry",
            input_db->getDatabase("CartesianGeometry"));

    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>(
            "PatchHierarchy",
            grid_geometry);

    WENOMOLOps* weno_mol_ops =
        new WENOMOLOps(
            "WENOMOLOps",
            grid_geometry,
            input_db->getDatabase("WENOMOLOps"));

    tbox::Pointer<algs::MethodOfLinesIntegrator<NDIM> > time_integrator =
        new algs::MethodOfLinesIntegrator<NDIM>(
            "MethodOfLinesIntegrator",
            input_db->getDatabase("MethodOfLinesIntegrator"),
            weno_mol_ops);

    tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
        new mesh::StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize",
            time_integrator,
            input_db->getDatabase("StandardTagAndInitialize"));

    tbox::Pointer<mesh::BergerRigoutsos<NDIM> > box_generator =
        new mesh::BergerRigoutsos<NDIM>();

    tbox::Pointer<mesh::LoadBalancer<NDIM> >  load_balancer =
        new mesh::LoadBalancer<NDIM>(
            "LoadBalancer",
            input_db->getDatabase("LoadBalancer"));

    tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
        new mesh::GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            error_detector, box_generator, load_balancer);

    /*
     * Set up visualization plot file writer.
     */
    tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
        new appu::VisItDataWriter<NDIM>(
            "VisIt Writer", visit_dump_dirname, visit_number_procs_per_file);

    if (uses_visit)
    {
        weno_mol_ops->registerVisItDataWriter(visit_data_writer);
    }

    /*
     * Initialize hierarchy configuration and data on all patches.  Then, close
     * restart file and write initial state for visualization.
     */
    time_integrator->initializeIntegrator(gridding_algorithm);

    tbox::Array<int>* tag_buffer_array = new tbox::Array<int>(gridding_algorithm->getMaxLevels());
    for (int il = 0; il < gridding_algorithm->getMaxLevels(); ++il)
    {
        (*tag_buffer_array)[il] = 2;
        tbox::pout << "il = " << il << " tag_buffer = " << (*tag_buffer_array)[il] << endl;
    }

    double loop_time = 0.0;

    if (tbox::RestartManager::getManager()->isFromRestart())
    {
        patch_hierarchy->getFromRestart(gridding_algorithm->getMaxLevels());

        gridding_algorithm->getTagAndInitializeStrategy()->
            resetHierarchyConfiguration(patch_hierarchy,
                                        0,
                                        patch_hierarchy->getFinestLevelNumber());

    }
    else
    {
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy,loop_time);

        bool done = false;
        bool initial_time = true;
        for (int ln = 0; gridding_algorithm->levelCanBeRefined(ln) && !done; ++ln)
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy,
                                               loop_time,
                                               initial_time,
                                               (*tag_buffer_array)[ln]);
            done = !(patch_hierarchy->finerLevelExists(ln));
        }
    }

    tbox::RestartManager::getManager()->closeRestartFile();

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
                patch_hierarchy, 0, 0.0);
        }
    }

    /*
     * Time step loop.
     */
    double loop_time_end = 1.0;
    int iteration_num = 0;

    while (loop_time < loop_time_end)
    {
        iteration_num = iteration_num + 1;

        tbox::pout <<                                                       endl;
        tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        tbox::pout << "At begining of timestep # " <<  iteration_num - 1 << endl;
        tbox::pout << "Simulation time is " << loop_time                 << endl;

        double dt_now = time_integrator->getTimestep(patch_hierarchy, loop_time);
        if (loop_time+dt_now > loop_time_end) dt_now = loop_time_end - loop_time;
        time_integrator->advanceHierarchy(patch_hierarchy, loop_time, dt_now);
        loop_time += dt_now;

        tbox::pout <<                                                       endl;
        tbox::pout << "At end      of timestep # " <<  iteration_num - 1 << endl;
        tbox::pout << "Simulation time is " << loop_time                 << endl;
        tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        tbox::pout <<                                                       endl;

        /*
         * At specified intervals, output timer data and write visualization and
         * restart and files.
         */
        if (write_timer_data && iteration_num%timer_dump_interval == 0)
        {
            tbox::pout << "\nWriting timer data...\n\n";
            tbox::TimerManager::getManager()->print(tbox::plog);
        }

        if (write_restart && iteration_num%restart_interval == 0)
        {
            tbox::pout << "\nWriting restart files...\n\n";
            tbox::RestartManager::getManager()->writeRestartFile(
                restart_write_dirname, iteration_num);
        }

        if (viz_dump_data && iteration_num%viz_dump_interval == 0)
        {
            if (uses_visit)
            {
                tbox::pout << "\nWriting visualization files...\n\n";
                visit_data_writer->writePlotData(
                    patch_hierarchy, iteration_num, loop_time);
            }
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
     * Ensure the last state is written out for plotting.
     */
    if (viz_dump_data && iteration_num%viz_dump_interval != 0)
    {
        if (uses_visit)
        {
            tbox::pout << "\nWriting visualization files...\n\n";
            visit_data_writer->writePlotData(
                patch_hierarchy, iteration_num, loop_time);
        }
    }

    /*
     * At conclusion of simulation, deallocate objects.
     */
    grid_geometry.setNull();
    patch_hierarchy.setNull();
    box_generator.setNull();
    load_balancer.setNull();
    delete weno_mol_ops;
    time_integrator.setNull();
    error_detector.setNull();
    gridding_algorithm.setNull();
    time_integrator.setNull();
    visit_data_writer.setNull();

    tbox::SAMRAIManager::shutdown();
    tbox::MPI::finalize();

    return 0;
}// main
