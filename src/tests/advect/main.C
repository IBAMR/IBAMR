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
#include <tbox/Utilities.h>

// Headers for major algorithm/data structure objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HyperbolicLevelIntegrator.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <StandardTagAndInitialize.h>
#include <TimeRefinementIntegrator.h>
#include <TimeRefinementLevelStrategy.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <LocationIndexRobinBcCoefs.h>

#include <ibamr/AdvectHypPatchOps.h>
#include <ibamr/GodunovAdvector.h>

#include <stools/HierarchyMathOps.h>

#include "QInit.h"
#include "USet.h"

using namespace IBAMR;
using namespace SAMRAI;
using namespace std;

/************************************************************************
 *                                                                      *
 * This is the main program for an AMR solution of the linear           *
 * advection equation: dQ/dt + div(u*Q) = 0, where "Q" is a             *
 * scalar-valued function and "u" is a possibly time-dependent          *
 * advection velocity.  This application program is constructed by      *
 * composing several algorithm objects found in the SAMRAI and IBAMR    *
 * libraries with a few that are specific to this application.  A       *
 * brief description of these object follows.                           *
 *                                                                      *
 * There are two main data containment objects.  These are:             *
 *                                                                      *
 *    PatchHierarchy - A container for the AMR patch hierarchy and the  *
 *       data on the grid.                                              *
 *                                                                      *
 *    CartesianGridGeometry - Defines and maintains the Cartesian       *
 *       coordinate system on the grid.  The PatchHierarchy maintains   *
 *       a reference to this object.                                    *
 *                                                                      *
 * A single overarching algorithm object drives the time integration    *
 * and adaptive gridding processes:                                     *
 *                                                                      *
 *    TimeRefinementIntegrator - Coordinates time integration and       *
 *       adaptive gridding procedures for the various levels in the     *
 *       AMR patch hierarchy.  Local time refinement is employed        *
 *       during hierarchy integration; i.e., finer levels are advanced  *
 *       using smaller time increments than coarser level.  Thus, this  *
 *       object also invokes data synchronization procedures which      *
 *       couple the solution on different patch hierarchy levels.       *
 *                                                                      *
 * The time refinement integrator is not specific to the numerical      *
 * methods used and the problem being solved.  It maintains references  *
 * to two other finer grain algorithmic objects that are more specific  *
 * to the problem at hand and with which it is configured when they     *
 * are passed into its constructor.  These finer grain algorithm        *
 * objects are:                                                         *
 *                                                                      *
 *    HyperbolicLevelIntegrator - Defines data management procedures    *
 *       for level integration, data synchronization between levels,    *
 *       and tagging cells for refinement.  These operations are        *
 *       tailored to explicit time integration algorithms used for      *
 *       hyperbolic systems of conservation laws, such as the Euler     *
 *       equations.  This integrator manages data for numerical         *
 *       routines that treat individual patches in the AMR patch        *
 *       hierarchy.  In this particular application, it maintains a     *
 *       pointer to the GodunovAdvector object that defines variables   *
 *       and provides numerical routines for the linear advection       *
 *       problem.                                                       *
 *                                                                      *
 *       GodunovAdvector - Provides the numerical routines necessary    *
 *          to solve the discrete linear advection equation on each     *
 *          patch in the AMR hierarchy.                                 *
 *                                                                      *
 *       AdvectHypPatchOps - Defines variables and numerical routines   *
 *          for the discrete linear advection equation on each patch    *
 *          in the AMR hierarchy and interfaces the GodunovAdvector     *
 *          with the HyperbolicLevelIntegrator.                         *
 *                                                                      *
 *    GriddingAlgorithm - Drives the AMR patch hierarchy generation     *
 *       and regridding procedures.  This object maintains references   *
 *       to three other algorithmic objects with which it is            *
 *       configured when they are passed into its constructor.  They    *
 *       are:                                                           *
 *                                                                      *
 *       BergerRigoutsos - Clusters cells tagged for refinement on a    *
 *          patch level into a collection of logically-rectangular box  *
 *          domains.                                                    *
 *                                                                      *
 *       LoadBalancer - Processes the boxes generated by the            *
 *          BergerRigoutsos algorithm into a configuration from which   *
 *          patches are contructed.  The algorithm used in this class   *
 *          assumes a spatially-uniform workload distribution; thus,    *
 *          it attempts to produce a collection of boxes each of which  *
 *          contains the same number of cells.  The load balancer also  *
 *          assigns patches to processors.                              *
 *                                                                      *
 *       StandardTagAndInitialize - Couples the gridding algorithm      *
 *          to the HyperbolicIntegrator. Selects cells for refinement   *
 *          based on either Gradient detection, Richardson              *
 *          extrapolation, or pre-defined Refine box region.  The       *
 *          object maintains a pointer to the                           *
 *          HyperbolicLevelIntegrator, which is passed into its         *
 *          constructor, for this purpose.                              *
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

int main(int argc, char *argv[])
{
    /*
     * Initialize MPI and SAMRAI, enable logging, and process command
     * line.
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
     * Retrieve "Main" section of the input database.  First, read
     * dump information, which is used for writing plot files.
     * Second, if proper restart information was given on command
     * line, and the restart interval is non-zero, create a restart
     * database.  Third, read in solution algorithm customizations.
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

    bool use_refined_timestepping = false;
    if (main_db->keyExists("timestepping"))
    {
        string timestepping_method = main_db->getString("timestepping");
        if (timestepping_method == "SYNCHRONIZED")
        {
            use_refined_timestepping = false;
        }
        else
        {
            use_refined_timestepping = true;
        }
    }
    if (use_refined_timestepping)
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

    /*
     * Get the restart manager and root restart database.  If run is
     * from restart, open the restart file.
     */
    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
    if (is_from_restart)
    {
        restart_manager->openRestartFile(
            restart_read_dirname, restore_num, tbox::MPI::getNodes());
    }

    /*
     * Create major algorithm and data objects which comprise
     * application.  Each object will be initialized either from input
     * data or restart files, or a combination of both.  Refer to each
     * class constructor for details.  For more information on the
     * composition of objects for this application, see comments at
     * top of file.
     */
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry =
        new geom::CartesianGridGeometry<NDIM>(
            "CartesianGeometry",
            input_db->getDatabase("CartesianGeometry"));

    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>(
            "PatchHierarchy",
            grid_geometry);

    tbox::Pointer<GodunovAdvector> advector =
        new GodunovAdvector(
            "GodunovAdvector",
            input_db->getDatabase("GodunovAdvector"));

    tbox::Pointer<AdvectHypPatchOps> hyp_patch_ops =
        new AdvectHypPatchOps(
            "AdvectHypPatchOps",
            input_db->getDatabase("AdvectHypPatchOps"),
            advector, grid_geometry);

    tbox::Pointer< pdat::FaceVariable<NDIM,double> > u_var =
        new pdat::FaceVariable<NDIM,double>("u");
    USet u_set("USet", grid_geometry, input_db->getDatabase("USet"));
    hyp_patch_ops->registerAdvectionVelocity(
        u_var, u_is_div_free, tbox::Pointer<SetDataStrategy>(&u_set,false));

    tbox::Pointer< pdat::CellVariable<NDIM,double> > Q_var =
        new pdat::CellVariable<NDIM,double>("Q");
    QInit Q_init(
        "QInit", grid_geometry, input_db->getDatabase("QInit"));
    solv::LocationIndexRobinBcCoefs<NDIM> physical_bc_coef(
        "physical_bc_coef",
        input_db->getDatabase("LocationIndexRobinBcCoefs"));
    hyp_patch_ops->registerAdvectedQuantity(
        Q_var, consv_form,
        tbox::Pointer<SetDataStrategy>(&Q_init,false),
        &physical_bc_coef);

    tbox::Pointer<algs::HyperbolicLevelIntegrator<NDIM> > hyp_level_integrator =
        new algs::HyperbolicLevelIntegrator<NDIM>(
            "HyperbolicLevelIntegrator",
            input_db->getDatabase("HyperbolicLevelIntegrator"),
            hyp_patch_ops, true, use_refined_timestepping);

    tbox::Pointer<mesh::StandardTagAndInitialize<NDIM> > error_detector =
        new mesh::StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize",
            hyp_level_integrator,
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

    tbox::Pointer<algs::TimeRefinementIntegrator<NDIM> > time_integrator =
        new algs::TimeRefinementIntegrator<NDIM>(
            "TimeRefinementIntegrator",
            input_db->getDatabase("TimeRefinementIntegrator"),
            patch_hierarchy, hyp_level_integrator, gridding_algorithm);

    /*
     * Set up visualization plot file writer.
     */
    tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
        new appu::VisItDataWriter<NDIM>(
            "VisIt Writer", visit_dump_dirname, visit_number_procs_per_file);

    if (uses_visit)
    {
        hyp_patch_ops->registerVisItDataWriter(visit_data_writer);
    }

    /*
     * Initialize hierarchy configuration and data on all patches.
     * Then, close restart file and write initial state for
     * visualization.
     */
    double dt_now = time_integrator->initializeHierarchy();
    tbox::RestartManager::getManager()->closeRestartFile();

        /*
     * After creating all objects and initializing their state, we
     * print the input database contents to the log file.
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
            visit_data_writer->writePlotData(
                patch_hierarchy,
                time_integrator->getIntegratorStep(),
                time_integrator->getIntegratorTime());
        }
    }

    /*
     * Time step loop.  Note that the step count and integration time
     * are maintained by the algs::TimeRefinementIntegrator object.
     */
    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();

    int iteration_num = time_integrator->getIntegratorStep();

    while (!tbox::Utilities::deq(loop_time,loop_time_end) &&
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

        /*
         * At specified intervals, write restart and visualization files.
         */
        if (write_restart && iteration_num%restart_interval == 0)
        {
            tbox::RestartManager::getManager()->writeRestartFile(
                restart_write_dirname, iteration_num);
        }

        if (viz_dump_data && iteration_num%viz_dump_interval == 0)
        {
            if (uses_visit)
            {
                visit_data_writer->writePlotData(
                    patch_hierarchy, iteration_num, loop_time);
            }
        }
    }

    /*
     * Determine the accuracy of the computed solution.
     */
    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();

    const tbox::Pointer<hier::VariableContext> Q_ctx =
        hyp_level_integrator->getCurrentContext();
    const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, Q_ctx);
    const int Q_cloned_idx = var_db->registerClonedPatchDataIndex(Q_var, Q_idx);

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(Q_cloned_idx, loop_time);
    }

    Q_init.setDataOnPatchHierarchy(
        Q_cloned_idx, Q_var, patch_hierarchy, loop_time);

    STOOLS::HierarchyMathOps hier_math_ops(
        "HierarchyMathOps", patch_hierarchy);
    hier_math_ops.setPatchHierarchy(patch_hierarchy);
    hier_math_ops.resetLevels(coarsest_ln, finest_ln);
    const int wgt_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

    math::HierarchyCellDataOpsReal<NDIM,double> hier_cc_data_ops(
        patch_hierarchy, coarsest_ln, finest_ln);

    hier_cc_data_ops.subtract(Q_cloned_idx, Q_idx, Q_cloned_idx);
    tbox::pout << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
               << "  L1-norm:  " << hier_cc_data_ops.L1Norm(Q_cloned_idx,wgt_idx)  << "\n"
               << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_cloned_idx,wgt_idx)  << "\n"
               << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_cloned_idx,wgt_idx) << "\n";

    /*
     * Ensure the last state is written out for plotting.
     */
    if (viz_dump_data && iteration_num%viz_dump_interval != 0)
    {
        if (uses_visit)
        {
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
    advector.setNull();
    hyp_patch_ops.setNull();
    hyp_level_integrator.setNull();
    error_detector.setNull();
    gridding_algorithm.setNull();
    time_integrator.setNull();
    visit_data_writer.setNull();

    tbox::SAMRAIManager::shutdown();
    tbox::MPI::finalize();

    return 0;
}// main
