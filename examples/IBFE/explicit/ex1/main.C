// Config files
#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic libMesh objects
#include <exodusII_io.h>
#include <mesh.h>
#include <mesh_generation.h>
#include <quadrature.h>
#include <string_to_enum.h>

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
#include <ibamr/IBFEHierarchyIntegrator.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

using namespace IBAMR;
using namespace IBTK;
using namespace SAMRAI;
using namespace std;

namespace
{
// Problem parameters.
static const double R = 0.25;
static const double w = 0.0625;
static const double gamma = 0.1;
static const double mu = 1.0;

// Coordinate mapping function.
Point
coordinate_mapping_function(
    const Point& s,
    void* ctx)
{
    return Point((R+      s(1))*cos(s(0)/R)+0.5,
                 (R+gamma+s(1))*sin(s(0)/R)+0.5);
}// coordinate_mapping_function

// Stress tensor function.
TensorValue<double>
PK1_stress_function(
    const TensorValue<double>& dX_ds,
    const Point& X,
    const Point& s,
    Elem* const elem,
    const double& time,
    void* ctx)
{
    TensorValue<double> P = (mu/w)*dX_ds;
//  P(0,1) = 0.0;
//  P(1,1) = 0.0;
    return P;
}// PK1_stress_function
}

int
main(
    int argc,
    char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    tbox::SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    tbox::SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    tbox::SAMRAIManager::startup();

    {// cleanup all smart Pointers prior to shutdown

        // Process command line and enable logging.
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

        // Create input database and parse all data in input file.
        tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
        tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

        // Retrieve "Main" section of the input database.  First, read dump
        // information, which is used for writing plot files.  Second, if proper
        // restart information was given on command line, and the restart
        // interval is non-zero, create a restart database.
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
        bool uses_exodus = false;
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
                if (viz_writer[i] == "ExodusII") uses_exodus = true;
            }
            if (uses_visit)
            {
                visit_dump_dirname = viz_dump_dirname;
                if (main_db->keyExists("visit_number_procs_per_file"))
                {
                    visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
                }
            }
            else if (uses_exodus)
            {
                // intentionally blank
            }
            else
            {
                TBOX_ERROR("main(): "
                           << "\nUnrecognized 'viz_writer' entry..."
                           << "\nOptions are: 'VisIt', 'ExodusII'"
                           << endl);
            }
            if (uses_visit || uses_exodus)
            {
                if (viz_dump_dirname.empty())
                {
                    TBOX_ERROR("main(): "
                               << "\nviz_dump_dirname is null ... "
                               << "\nThis must be specified for use with VisIt or ExodusII"
                               << endl);
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

        const bool write_restart = restart_interval > 0 && !restart_write_dirname.empty();

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

        const bool use_nonuniform_load_balancer = main_db->getBoolWithDefault("use_nonuniform_load_balancer", false);

        // Create a simple FE mesh.
        Mesh mesh(NDIM);
        const double R = 0.25;
        const double w = 0.0625;
        const int M = input_db->getIntegerWithDefault("M", 4);
        string elem_type = "QUAD4";
        MeshTools::Generation::build_square(mesh,
                                            28*M, M,
                                            0.0, 2.0*M_PI*R,
                                            0.0, w,
                                            Utility::string_to_enum<ElemType>(elem_type));
        ExodusII_IO mesh_writer(mesh);

        PeriodicBoundary pbc;
        pbc.myboundary = 3;
        pbc.pairedboundary = 1;
        VectorValue<double> boundary_translation(2.0*M_PI*R, 0.0, 0.0);
        pbc.translation_vector = boundary_translation;

        // Create the FE data manager used to manage mappings between the FE
        // mesh and the Cartesian grid.
        EquationSystems equation_systems(mesh);
        IBTK::FEDataManager* fe_data_manager = IBTK::FEDataManager::getManager("EP FE Manager");

        const int mesh_level_number = input_db->getInteger("MAX_LEVELS")-1;
        fe_data_manager->setEquationSystems(&equation_systems, mesh_level_number);

        const std::string quad_type = input_db->getStringWithDefault("quad_type", "QGAUSS");
        const std::string quad_order = input_db->getStringWithDefault("quad_order", "SIXTH");
        AutoPtr<QBase> qrule = QBase::build(Utility::string_to_enum<QuadratureType>(quad_type),NDIM,Utility::string_to_enum<Order>(quad_order));
        const std::string weighting_fcn = input_db->getStringWithDefault("weighting_fcn", "IB_4");
        fe_data_manager->setInteractionScheme(weighting_fcn, *qrule);

        const bool use_consistent_mass_matrix = input_db->getBoolWithDefault("use_consistent_mass_matrix", true);
        fe_data_manager->setUseConsistentMassMatrix(use_consistent_mass_matrix);

        // Get the restart manager and root restart database.  If run is from
        // restart, open the restart file.
        tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
        if (is_from_restart)
        {
            restart_manager->openRestartFile(restart_read_dirname, restore_num, tbox::SAMRAI_MPI::getNodes());
        }

        // Create major algorithm and data objects which comprise the
        // application.  Each object will be initialized either from input data
        // or restart files, or a combination of both.  Refer to each class
        // constructor for details.  For more information on the composition of
        // objects for this application, see comments at top of file.
        tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
            new geom::CartesianGridGeometry<NDIM>(
                "CartesianGeometry",
                input_db->getDatabase("CartesianGeometry"));

        tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
            new hier::PatchHierarchy<NDIM>(
                "PatchHierarchy",
                grid_geometry);

        tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator =
            new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                input_db->getDatabase("INSStaggeredHierarchyIntegrator"),
                patch_hierarchy);

        tbox::Pointer<IBFEHierarchyIntegrator> time_integrator =
            new IBFEHierarchyIntegrator(
                "IBFEHierarchyIntegrator",
                input_db->getDatabase("IBFEHierarchyIntegrator"),
                patch_hierarchy, navier_stokes_integrator, fe_data_manager);
        time_integrator->setInitialCoordinateMappingFunction(&coordinate_mapping_function);
        time_integrator->setPK1StressTensorFunction(&PK1_stress_function);

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
        if (use_nonuniform_load_balancer)
        {
            time_integrator->registerLoadBalancer(load_balancer);
        }

        tbox::Pointer<mesh::GriddingAlgorithm<NDIM> > gridding_algorithm =
            new mesh::GriddingAlgorithm<NDIM>(
                "GriddingAlgorithm",
                input_db->getDatabase("GriddingAlgorithm"),
                error_detector, box_generator, load_balancer);

        // Create initial condition specification objects.
        tbox::Pointer<CartGridFunction> u_init = new muParserCartGridFunction("u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        tbox::Pointer<CartGridFunction> p_init = new muParserCartGridFunction("p_init", input_db->getDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Create boundary condition specification objects (when necessary).
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

                u_bc_coefs.push_back(new muParserRobinBcCoefs(bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry));
            }
            time_integrator->registerVelocityPhysicalBcCoefs(u_bc_coefs);
        }

        // Set up visualization plot file writer.
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
            new appu::VisItDataWriter<NDIM>(
                "VisIt Writer",
                visit_dump_dirname, visit_number_procs_per_file);
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        ExodusII_IO exodus_io(mesh);

        // Initialize hierarchy configuration and data on all patches.  Then,
        // close restart file and write initial state for visualization.
        time_integrator->initializeHierarchyIntegrator(gridding_algorithm);
        for (unsigned int k = 0; k < equation_systems.n_systems(); ++k)
        {
            System& system = equation_systems.get_system(k);
            system.get_dof_map().add_periodic_boundary(pbc);
        }
        double dt_now = time_integrator->initializeHierarchy();
        tbox::RestartManager::getManager()->closeRestartFile();

        // After creating all objects and initializing their state, we print the
        // input database contents to the log file.
        tbox::plog << "\nCheck input data before simulation:" << endl;
        tbox::plog << "Input database..." << endl;
        input_db->printClassData(tbox::plog);

        // Write initial visualization files.
        if (viz_dump_data)
        {
            if (uses_visit)
            {
                tbox::pout << "\nWriting visualization files...\n\n";
                visit_data_writer->writePlotData(
                    patch_hierarchy,
                    time_integrator->getIntegratorStep(),
                    time_integrator->getIntegratorTime());
            }
            if (uses_exodus)
            {
                std::ostringstream os;
                os << "output.ex2";
                exodus_io.write_timestep(os.str(), equation_systems, time_integrator->getIntegratorStep()/viz_dump_interval+1, time_integrator->getIntegratorTime());
            }
        }

        // Time step loop.  Note that the step count and integration time are
        // maintained by the time integrator object.
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

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
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
                }
                if (uses_exodus)
                {
                    std::ostringstream os;
                    os << "output.ex2";
                    exodus_io.write_timestep(os.str(), equation_systems, iteration_num/viz_dump_interval+1, loop_time);
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

    }// cleanup all smart Pointers prior to shutdown

    tbox::SAMRAIManager::shutdown();

    // All done.
    return 0;
}// main
