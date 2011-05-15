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

// C++ standard library includes
#include <limits>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic SAMRAI objects
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
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/muParserRobinBcCoefs.h>

using namespace IBAMR;
using namespace IBTK;
using namespace SAMRAI;
using namespace std;

#define COMPUTE_STREAM_FUNCTION 0

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

        string log_file_name = "navier_stokes.log";
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
                           << " not specified in input file");
            }
        }

        const bool write_restart = restart_interval > 0
            && !restart_write_dirname.empty();

        int hier_dump_interval = 0;
        if (main_db->keyExists("hier_dump_interval"))
        {
            hier_dump_interval = main_db->getInteger("hier_dump_interval");
        }

        string hier_dump_dirname;
        if (hier_dump_interval > 0)
        {
            if (main_db->keyExists("hier_dump_dirname"))
            {
                hier_dump_dirname = main_db->getString("hier_dump_dirname");
            }
            else
            {
                TBOX_ERROR("hier_dump_interval > 0, but key `hier_dump_dirname'"
                           << " not specified in input file");
            }
        }

        const bool write_hier_data = (hier_dump_interval > 0)
            && !(hier_dump_dirname.empty());
        if (write_hier_data)
        {
            tbox::Utilities::recursiveMkdir(hier_dump_dirname);
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

        tbox::Pointer<INSStaggeredHierarchyIntegrator> time_integrator =
            new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                input_db->getDatabase("INSStaggeredHierarchyIntegrator"),
                patch_hierarchy);

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
         * Create boundary condition specification objects.
         */
        muParserRobinBcCoefs u0_bc_coef("u0_bc_coef", input_db->getDatabase("VelocityBcCoefs_0"), grid_geometry);

        solv::LocationIndexRobinBcCoefs<NDIM> u1_bc_coef("u1_bc_coef", tbox::Pointer<tbox::Database>(NULL));
        for (int i = 0; i < 2*NDIM; ++i)
        {
            u1_bc_coef.setBoundaryValue(i, 0.0);
        }
#if (NDIM > 2)
        solv::LocationIndexRobinBcCoefs<NDIM> u2_bc_coef("u2_bc_coef", tbox::Pointer<tbox::Database>(NULL));
        for (int i = 0; i < 2*NDIM; ++i)
        {
            u2_bc_coef.setBoundaryValue(i,0.0);
        }
#endif
        blitz::TinyVector<solv::RobinBcCoefStrategy<NDIM>*,NDIM> U_bc_coefs;
        U_bc_coefs[0] = &u0_bc_coef;
        U_bc_coefs[1] = &u1_bc_coef;
#if (NDIM > 2)
        U_bc_coefs[2] = &u2_bc_coef;
#endif
        time_integrator->registerVelocityPhysicalBcCoefs(U_bc_coefs);

        /*
         * Set up visualization plot file writer.
         */
        tbox::Pointer<appu::VisItDataWriter<NDIM> > visit_data_writer =
            new appu::VisItDataWriter<NDIM>(
                "VisIt Writer",
                visit_dump_dirname, visit_number_procs_per_file);

        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        /*
         * Setup a variable for computing and storing the stream function.
         */
        hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
#if COMPUTE_STREAM_FUNCTION
        tbox::Pointer<hier::VariableContext> main_context = var_db->getContext("main::CONTEXT");
        tbox::Pointer<pdat::NodeVariable<NDIM,double> > psi_var = new pdat::NodeVariable<NDIM,double>("psi");
        const int psi_idx = var_db->registerVariableAndContext(psi_var, main_context, 0);
        if (uses_visit)
        {
            visit_data_writer->registerPlotQuantity(psi_var->getName(), "SCALAR", psi_idx, 0, 1.0);
        }
#endif
        /*
         * Initialize hierarchy configuration and data on all patches.  Then,
         * close restart file and write initial state for visualization.
         */
        time_integrator->initializeHierarchyIntegrator(gridding_algorithm);
        double dt_now = time_integrator->initializeHierarchy();
#if COMPUTE_STREAM_FUNCTION
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(psi_idx))
            {
                level->allocatePatchData(psi_idx, time_integrator->getIntegratorTime());
            }
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                tbox::Pointer<pdat::NodeData<NDIM,double> > psi_data = patch->getPatchData(psi_idx);
                psi_data->fillAll(0.0);
            }
        }
#endif
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
                    patch_hierarchy,
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

        while (!tbox::MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            /*
             * At specified intervals, write state data for post-processing.
             */
            if (write_hier_data && iteration_num%hier_dump_interval == 0)
            {
                hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();

                tbox::plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
                tbox::plog << "simulation time is " << loop_time << endl;

                string file_name = hier_dump_dirname + "/" + "hier_data.";
                char temp_buf[128];
                sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, tbox::SAMRAI_MPI::getRank());
                file_name += temp_buf;

                tbox::Pointer<tbox::HDFDatabase> hier_db = new tbox::HDFDatabase("hier_db");
                hier_db->create(file_name);

                hier::ComponentSelector hier_data;
                hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVar(),
                                                                       time_integrator->getCurrentContext()));
                hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getPressureVar(),
                                                                       time_integrator->getCurrentContext()));

                patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
                hier_db->putDouble("loop_time", loop_time);
                hier_db->putDouble("end_time", loop_time_end);
                hier_db->putDouble("dt", dt_now);
                hier_db->putInteger("iteration_num", iteration_num);
                hier_db->putInteger("hier_dump_interval", hier_dump_interval);

                hier_db->close();
            }

            iteration_num = time_integrator->getIntegratorStep() + 1;

            tbox::pout <<                                                        endl;
            tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++++++"  << endl;
            tbox::pout << "At beginning of timestep # " <<  iteration_num - 1 << endl;
            tbox::pout << "Simulation time is " << loop_time                  << endl;

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
             * and print out timer data.
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
#if COMPUTE_STREAM_FUNCTION
                /*
                 * Compute the stream-function.
                 */
                const tbox::Pointer<pdat::SideVariable<NDIM,double> > u_var = time_integrator->getVelocityVar();
                const tbox::Pointer<hier::VariableContext> u_ctx = time_integrator->getCurrentContext();
                const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);

                TBOX_ASSERT(patch_hierarchy->getFinestLevelNumber() == 0);
                tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(0);
                TBOX_ASSERT(level->getNumberOfPatches() == 1);
                if (!level->checkAllocated(psi_idx))
                {
                    level->allocatePatchData(psi_idx, loop_time);
                }

                for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());

                    const hier::Box<NDIM>& patch_box = patch->getBox();
                    hier::Box<NDIM> lower_box = patch_box;
                    lower_box.upper()(1) = lower_box.lower()(1);
                    hier::Box<NDIM> upper_box = patch_box;
                    upper_box.upper()(1) = upper_box.upper()(1)-1;

                    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                    const double* const dx = pgeom->getDx();
                    tbox::Pointer<pdat::NodeData<NDIM,double> > psi_data = patch->getPatchData(psi_idx);
                    tbox::Pointer<pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_idx);

                    psi_data->fillAll(numeric_limits<double>::quiet_NaN());
                    (*psi_data)(pdat::NodeIndex<NDIM>(hier::Index<NDIM>(0),0)) = 0.0;

                    for (hier::Box<NDIM>::Iterator it(lower_box); it; it++)
                    {
                        const hier::Index<NDIM>& i = it();
                        const pdat::NodeIndex<NDIM> n_i_l(i,0);
                        pdat::NodeIndex<NDIM> n_i_u(n_i_l);
                        n_i_u(0) += 1;
                        const pdat::SideIndex<NDIM> s_i(i,1,0);
                        const double dpsi_dx = -(*u_data)(s_i);
                        (*psi_data)(n_i_u) = (*psi_data)(n_i_l) + dpsi_dx*dx[0];
                    }

                    for (hier::Box<NDIM>::Iterator it(pdat::NodeGeometry<NDIM>::toNodeBox(upper_box)); it; it++)
                    {
                        const hier::Index<NDIM>& i = it();
                        const pdat::NodeIndex<NDIM> n_i_l(i,0);
                        pdat::NodeIndex<NDIM> n_i_u(n_i_l);
                        n_i_u(1) += 1;
                        const pdat::SideIndex<NDIM> s_i(i,0,0);
                        const double dpsi_dy = +(*u_data)(s_i);
                        (*psi_data)(n_i_u) = (*psi_data)(n_i_l) + dpsi_dy*dx[1];
                    }
                }
#endif
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
         * Ensure the last state is written out.
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

#if (NDIM == 2)
        /*
         * Determine the accuracy of the computed solution.
         */
        tbox::plog.precision(12);
        tbox::plog << std::scientific;
        const tbox::Pointer<pdat::SideVariable<NDIM,double> > u_var = time_integrator->getVelocityVar();
        const tbox::Pointer<hier::VariableContext> u_ctx = time_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);

        const tbox::Pointer<pdat::CellVariable<NDIM,double> > p_var = time_integrator->getPressureVar();
        const tbox::Pointer<hier::VariableContext> p_ctx = time_integrator->getCurrentContext();
        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

            tbox::plog << "u on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                tbox::plog << "patch_box = " << patch_box << "\n";

                tbox::Pointer<pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_idx);
                for (hier::Box<NDIM>::Iterator it(pdat::SideGeometry<NDIM>::toSideBox(patch_box,0)); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    const pdat::SideIndex<NDIM> s_i(i,0,pdat::SideIndex<NDIM>::Lower);
                    tbox::plog << (*u_data)(s_i) << "\n";
                }
            }
            tbox::plog << "\n";

            tbox::plog << "v on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                tbox::plog << "patch_box = " << patch_box << "\n";

                tbox::Pointer<pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_idx);
                for (hier::Box<NDIM>::Iterator it(pdat::SideGeometry<NDIM>::toSideBox(patch_box,1)); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    const pdat::SideIndex<NDIM> s_i(i,1,pdat::SideIndex<NDIM>::Lower);
                    tbox::plog << (*u_data)(s_i) << "\n";
                }
            }
            tbox::plog << "\n";

            tbox::plog << "p on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                tbox::plog << "patch_box = " << patch_box << "\n";
                tbox::Pointer<pdat::CellData<NDIM,double> > p_data = patch->getPatchData(p_idx);
                for (hier::Box<NDIM>::Iterator it(patch_box); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    tbox::plog << (*p_data)(i) << "\n";
                }
            }
            tbox::plog << "\n";
        }
#if COMPUTE_STREAM_FUNCTION
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            tbox::Pointer<hier::PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            const hier::Box<NDIM>& domain_box = level->getPhysicalDomain()[0];

            TBOX_ASSERT((domain_box.upper()(0) - domain_box.lower()(0) + 1)%2 == 0);
            const int i_center0 = (domain_box.upper()(0) - domain_box.lower()(0) + 1)/2 + domain_box.lower()(0);
            const hier::Box<NDIM> center_box0(hier::Index<NDIM>(i_center0,domain_box.lower()(1)),hier::Index<NDIM>(i_center0,domain_box.upper()(1)));

            TBOX_ASSERT((domain_box.upper()(1) - domain_box.lower()(1) + 1)%2 == 0);
            const int i_center1 = (domain_box.upper()(1) - domain_box.lower()(1) + 1)/2 + domain_box.lower()(1);
            const hier::Box<NDIM> center_box1(hier::Index<NDIM>(domain_box.lower()(0),i_center1),hier::Index<NDIM>(domain_box.upper()(0),i_center1));

            tbox::plog << "u(y=0.5) on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                tbox::Pointer<pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_idx);
                for (hier::Box<NDIM>::Iterator it(patch_box*center_box0); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    const pdat::SideIndex<NDIM> s_i(i,0,0);
                    tbox::plog << (double(i(1))+0.5)*dx[1] << "\t" << (*u_data)(s_i) << "\n";
                }
            }
            tbox::plog << "\n";

            tbox::plog << "p(y=0.5) on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                tbox::Pointer<pdat::CellData<NDIM,double> > p_data = patch->getPatchData(p_idx);
                for (hier::Box<NDIM>::Iterator it(patch_box*center_box0); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    hier::Index<NDIM> i_ll(i);
                    hier::Index<NDIM> i_l(i);
                    hier::Index<NDIM> i_u(i);
                    hier::Index<NDIM> i_uu(i);
                    i_ll(0) -= 2;
                    i_l (0) -= 1;
                    i_u (0) += 0;
                    i_uu(0) += 1;
                    tbox::plog << (double(i(1))+0.5)*dx[1] << "\t" << 9.0*((*p_data)(i_l) + (*p_data)(i_u))/16.0 - ((*p_data)(i_ll) + (*p_data)(i_uu))/16.0 << "\n";
                }
            }
            tbox::plog << "\n";

            tbox::plog << "v(x=0.5) on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                tbox::Pointer<pdat::SideData<NDIM,double> > u_data = patch->getPatchData(u_idx);
                for (hier::Box<NDIM>::Iterator it(patch_box*center_box1); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    const pdat::SideIndex<NDIM> s_i(i,1,0);
                    tbox::plog << (double(i(0))+0.5)*dx[0] << "\t" << (*u_data)(s_i) << "\n";
                }
            }
            tbox::plog << "\n";

            tbox::plog << "p(x=0.5) on level " << ln << "\n";
            for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
                const hier::Box<NDIM>& patch_box = patch->getBox();
                const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
                const double* const dx = pgeom->getDx();
                tbox::Pointer<pdat::CellData<NDIM,double> > p_data = patch->getPatchData(p_idx);
                for (hier::Box<NDIM>::Iterator it(patch_box*center_box1); it; it++)
                {
                    const hier::Index<NDIM>& i = it();
                    hier::Index<NDIM> i_ll(i);
                    hier::Index<NDIM> i_l(i);
                    hier::Index<NDIM> i_u(i);
                    hier::Index<NDIM> i_uu(i);
                    i_ll(1) -= 2;
                    i_l (1) -= 1;
                    i_u (1) += 0;
                    i_uu(1) += 1;
                    tbox::plog << (double(i(0))+0.5)*dx[0] << "\t" << 9.0*((*p_data)(i_l) + (*p_data)(i_u))/16.0 - ((*p_data)(i_ll) + (*p_data)(i_uu))/16.0 << "\n";
                }
            }
            tbox::plog << "\n";
        }
#endif
#endif

    }// cleanup all smart Pointers prior to shutdown

    tbox::SAMRAIManager::shutdown();
    PetscFinalize();

    return 0;
}// main
