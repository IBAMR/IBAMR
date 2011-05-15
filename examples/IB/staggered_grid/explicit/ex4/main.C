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
#include <ibamr/IBAnchorPointSpec.h>
#include <ibamr/IBKirchhoffRodForceGen.h>
#include <ibamr/IBStaggeredHierarchyIntegrator.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LSiloDataWriter.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

using namespace IBAMR;
using namespace IBTK;
using namespace SAMRAI;
using namespace std;

// Function prototypes
void
update_triads(
    const double freq,
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt);

void
output_data(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy,
    LDataManager* l_data_manager,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname);

// Basic 4-point delta function kernel.
inline double
ib4_delta_fcn(
    double r)
{
    r = std::abs(r);
    if (r < 1.0)
    {
        const double t2 = r * r;
        const double t6 = sqrt(-0.4e1 * t2 + 0.4e1 * r + 0.1e1);
        return -r / 0.4e1 + 0.3e1 / 0.8e1 + t6 / 0.8e1;
    }
    else if (r < 2.0)
    {
        const double t2 = r * r;
        const double t6 = sqrt(0.12e2 * r - 0.7e1 - 0.4e1 * t2);
        return -r / 0.4e1 + 0.5e1 / 0.8e1 - t6 / 0.8e1;
    }
    else
    {
        return 0.0;
    }
}// ib4_delta_fcn

// Specified delta-function width.
double W = 4.0;

// Re-scaled 4-point delta function kernel
inline double
scaled_ib4_delta_fcn(
    double r)
{
    return ib4_delta_fcn(r/(W/4.0))/(W/4.0);
}// scaled_ib4_delta_fcn

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
         * Determine the scaling factor to use to re-scale the 4-point delta
         * function, and configure the Lagrangian-Eulerian interaction routines
         * accordingly.
         */
        W = input_db->getDoubleWithDefault("W", W);
        LEInteractor::s_delta_fcn = &scaled_ib4_delta_fcn;
        LEInteractor::s_delta_fcn_stencil_size = std::ceil(W);
        LEInteractor::s_delta_fcn_C = -1.0;  // NOTE: C is not used by the
                                             // explicit IB integrator, but we
                                             // set it to cause an error in case
                                             // it is used elsewhere.

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

        const bool write_data = (data_dump_interval > 0)
            && !(data_dump_dirname.empty());
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

       tbox::Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator =
            new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                input_db->getDatabase("INSStaggeredHierarchyIntegrator"),
                patch_hierarchy);

        tbox::Pointer<IBKirchhoffRodForceGen> force_generator = new IBKirchhoffRodForceGen();

        tbox::Pointer<IBStaggeredHierarchyIntegrator> time_integrator =
            new IBStaggeredHierarchyIntegrator(
                "IBStaggeredHierarchyIntegrator",
                input_db->getDatabase("IBStaggeredHierarchyIntegrator"),
                patch_hierarchy, navier_stokes_integrator, force_generator);

        tbox::Pointer<IBStandardInitializer> initializer =
            new IBStandardInitializer(
                "IBStandardInitializer",
                input_db->getDatabase("IBStandardInitializer"));
        time_integrator->registerLInitStrategy(initializer);

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
        tbox::Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        time_integrator->registerVelocityInitialConditions(u_init);

        /*
         * Read the rotational frequency from the input database.
         */
        if (!input_db->keyExists("FREQ"))
        {
            TBOX_ERROR("value ``FREQ'' must be specified in the input file");
        }
        const double freq = input_db->getDouble("FREQ");

        /*
         * Set up visualization plot file writer.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
            new appu::VisItDataWriter<NDIM>(
                "VisIt Writer",
                visit_dump_dirname, visit_number_procs_per_file);
        tbox::Pointer<LSiloDataWriter> silo_data_writer =
            new LSiloDataWriter(
                "LSiloDataWriter",
                visit_dump_dirname);

        if (uses_visit)
        {
            initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerLSiloDataWriter(silo_data_writer);
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
        time_integrator->freeLInitStrategy();
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

        /*
         * Time step loop.  Note that the step count and integration time are
         * maintained by the time integrator object.
         */
        double loop_time = time_integrator->getIntegratorTime();
        double loop_time_end = time_integrator->getEndTime();
        double dt_old = 0.0;

        int iteration_num = time_integrator->getIntegratorStep();

        /*
         * At specified intervals, write state data for post-processing.
         */
        if (write_data && iteration_num%data_dump_interval == 0)
        {
            output_data(patch_hierarchy, time_integrator->getLDataManager(), iteration_num, loop_time, data_dump_dirname);
        }

        while (!tbox::MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep() + 1;

            tbox::pout <<                                                        endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            tbox::pout << "At beginning of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                  << endl;

            LDataManager* const lag_data_manager = time_integrator->getLDataManager();
            update_triads(freq, patch_hierarchy, lag_data_manager, loop_time, dt_now);

            dt_old = dt_now;
            double dt_new = time_integrator->advanceHierarchy(dt_now);

            loop_time += dt_now;
            dt_now = dt_new;

            tbox::pout <<                                                        endl;
            tbox::pout << "At end       of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                  << endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            tbox::pout <<                                                        endl;

            /*
             * At specified intervals, write visualization and restart files,
             * print out timer data, and write state data.
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
            }

            if (write_data && iteration_num%data_dump_interval == 0)
            {
                output_data(patch_hierarchy, time_integrator->getLDataManager(), iteration_num, loop_time, data_dump_dirname);
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

    }// cleanup all smart Pointers prior to shutdown

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main

void
update_triads(
    const double freq,
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt)
{
    // The angular velocity.
    static const double angw = 2.0*M_PI*freq;

    // The LMesh object provides the set of local Lagrangian nodes.
    const int finest_ln = hierarchy->getFinestLevelNumber();
    const tbox::Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    const std::vector<LNode*>& local_nodes = mesh->getNodes();

    // Update the director triads for all anchored points.
    tbox::Pointer<LData> D_data = l_data_manager->getLData("D", finest_ln);
    blitz::Array<double,2>& D_array = *D_data->getLocalFormVecArray();
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin();
         cit != local_nodes.end(); ++cit)
    {
        const LNode& node_idx = **cit;

        // If spec is NULL, then the IB point is NOT anchored.  Otherwise, the
        // IB point IS anchored.
        tbox::Pointer<IBAnchorPointSpec> spec = node_idx.getNodeData<IBAnchorPointSpec>();
        if (!spec.isNull())
        {
            // `global_idx' is the index of the vertex in the input file.
            const int global_idx = node_idx.getLagrangianIndex();
            (void) global_idx;  // to prevent compiler warnings...

            // `local_petsc_idx' is the index of the vertex in the internal
            // IBAMR data structures.  Generally, global_idx != local_petsc_idx.
            const int local_petsc_idx = node_idx.getLocalPETScIndex();

            double* D1 = &D_array(local_petsc_idx,0);
            D1[0] =  cos(angw*current_time);
            D1[1] =  sin(angw*current_time);
            D1[2] = 0.0;

            double* D2 = &D_array(local_petsc_idx,3);
            D2[0] = -sin(angw*current_time);
            D2[1] =  cos(angw*current_time);
            D2[2] = 0.0;

            double* D3 = &D_array(local_petsc_idx,6);
            D3[0] = 0.0;
            D3[1] = 0.0;
            D3[2] = 1.0;
        }
    }
    D_data->restoreArrays();
    return;
}// update_triads

void
output_data(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy,
    LDataManager* l_data_manager,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname)
{
    tbox::pout << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    tbox::pout << "simulation time is " << loop_time << endl;

    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();

    string file_name;
    char temp_buf[128];

    /*
     * Write Lagrangian data.
     */
    tbox::Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    for (int rank = 0; rank < tbox::SAMRAI_MPI::getNodes(); ++rank)
    {
        if (tbox::SAMRAI_MPI::getRank() == rank)
        {
            ofstream str(file_name.c_str(), rank == 0 ? ios_base::trunc : ios_base::app);

            if (rank == 0)
            {
                str << X_data->getGlobalNodeCount() << "\n";
            }

            int size;
            VecGetLocalSize(X_lag_vec, &size);
            if (size > 0)
            {
                str.precision(12);
                str.setf(ios::scientific);
                str.setf(ios::showpos);
                double* X_lag_arr;
                VecGetArray(X_lag_vec, &X_lag_arr);
                const int depth = NDIM;
                for (int k = 0; k < size/depth; ++k)
                {
                    for (int d = 0; d < depth; ++d)
                    {
                        str << X_lag_arr[depth*k+d];
                        if (d < depth-1)
                        {
                            str << " ";
                        }
                        else
                        {
                            str << "\n";
                        }
                    }
                }
                VecRestoreArray(X_lag_vec, &X_lag_arr);
            }
        }
        tbox::SAMRAI_MPI::barrier();
    }
    VecDestroy(X_lag_vec);

    tbox::Pointer<LData> D_data = l_data_manager->getLData("D", finest_hier_level);
    Vec D_petsc_vec = D_data->getVec();
    Vec D_lag_vec;
    VecDuplicate(D_petsc_vec, &D_lag_vec);
    l_data_manager->scatterPETScToLagrangian(D_petsc_vec, D_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "D.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    for (int rank = 0; rank < tbox::SAMRAI_MPI::getNodes(); ++rank)
    {
        if (tbox::SAMRAI_MPI::getRank() == rank)
        {
            ofstream str(file_name.c_str(), rank == 0 ? ios_base::trunc : ios_base::app);

            if (rank == 0)
            {
                str << D_data->getGlobalNodeCount() << "\n";
            }

            int size;
            VecGetLocalSize(D_lag_vec, &size);
            if (size > 0)
            {
                str.precision(12);
                str.setf(ios::scientific);
                str.setf(ios::showpos);
                double* D_lag_arr;
                VecGetArray(D_lag_vec, &D_lag_arr);
                const int depth = 3;
                for (int k = 0; k < size/depth; ++k)
                {
                    for (int d = 0; d < depth; ++d)
                    {
                        str << D_lag_arr[depth*k+d];
                        if (d < depth-1)
                        {
                            str << " ";
                        }
                        else
                        {
                            str << "\n";
                        }
                    }
                }
                VecRestoreArray(D_lag_vec, &D_lag_arr);
            }
        }
        tbox::SAMRAI_MPI::barrier();
    }
    VecDestroy(D_lag_vec);
    return;
}// output_data
