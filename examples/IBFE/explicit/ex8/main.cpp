// Copyright (c) 2002-2013, Boyce Griffith
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
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Elasticity model data.
namespace ModelData
{
static double kappa_s = 1.0e6;
static double kappa_t = 1.0e6;

// Tether (penalty) force function for the solid blocks.
void
block_tether_force_function(
    VectorValue<double>& F,
    const TensorValue<double>& /*FF*/,
    const libMesh::Point& X,
    const libMesh::Point& s,
    Elem* const /*elem*/,
    NumericVector<double>& /*X_vec*/,
    const vector<NumericVector<double>*>& /*system_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    F = kappa_s*(s-X);
    return;
}// block_tether_force_function

// Tether (penalty) force function for the thin beam.
void
beam_tether_force_function(
    VectorValue<double>& F,
    const TensorValue<double>& /*FF*/,
    const libMesh::Point& X,
    const libMesh::Point& s,
    Elem* const /*elem*/,
    NumericVector<double>& /*X_vec*/,
    const vector<NumericVector<double>*>& /*system_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    F.zero();
    static const libMesh::Point s0(0.5,0.5);
    static const libMesh::Point s1(1.5,0.5);
    static const double h = 0.004;
    double d0 = (s-s0).size();
    double d1 = (s-s1).size();
    if (d0 < h)
    {
        F = kappa_t*(1.0-d0/h)*(s-X);
    }
    else if (d1 < h)
    {
        F = kappa_t*(1.0-d1/h)*(s-X);
    }
    return;
}// beam_tether_force_function

// Stress tensor function for the thin beam.
static double mu_s;
void
beam_PK1_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /*X*/,
    const libMesh::Point& /*s*/,
    Elem* const /*elem*/,
    NumericVector<double>& /*X_vec*/,
    const vector<NumericVector<double>*>& /*system_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    const TensorValue<double> CC = FF.transpose()*FF;
    PP = mu_s*(FF-FF_inv_trans);
    return;
}// beam_PK1_stress_function
}
using namespace ModelData;

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
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string block1_exodus_filename = app_initializer->getExodusIIFilename("block1");
        const string block2_exodus_filename = app_initializer->getExodusIIFilename("block2");
        const string   beam_exodus_filename = app_initializer->getExodusIIFilename("beam"  );

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        const double dx = input_db->getDouble("DX");
        const double ds_block = input_db->getDouble("BLOCK_MFAC")*dx;
        const double ds_beam  = input_db->getDouble("BEAM_MFAC" )*dx;

        string block_elem_type = input_db->getString("BLOCK_ELEM_TYPE");
        string  beam_elem_type = input_db->getString("BEAM_ELEM_TYPE" );

        Mesh block1_mesh(NDIM);
        MeshTools::Generation::build_square(block1_mesh,
                                            ceil(0.5/ds_block), ceil(0.5/ds_block),
                                            0.0, 0.5,
                                            0.0, 0.5,
                                            Utility::string_to_enum<ElemType>(block_elem_type));
        Mesh block2_mesh(NDIM);
        MeshTools::Generation::build_square(block2_mesh,
                                            ceil(0.5/ds_block), ceil(0.5/ds_block),
                                            1.5, 2.0,
                                            0.0, 0.5,
                                            Utility::string_to_enum<ElemType>(block_elem_type));

        Mesh beam_mesh(NDIM);
        MeshTools::Generation::build_square(beam_mesh,
                                            ceil(1.0/ds_beam), max(ceil(0.016/ds_beam),4.0),
                                            0.5, 1.5,
                                            0.5-0.008, 0.5+0.008,
                                            Utility::string_to_enum<ElemType>(beam_elem_type));
        block1_mesh.prepare_for_use();
        block2_mesh.prepare_for_use();
        beam_mesh  .prepare_for_use();

        vector<Mesh*> meshes(3);
        meshes[0] = &block1_mesh;
        meshes[1] = &block2_mesh;
        meshes[2] = &  beam_mesh;

        mu_s    = input_db->getDouble("MU_S");
        kappa_s = input_db->getDouble("KAPPA_S");
        kappa_t = input_db->getDouble("KAPPA_T");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBFEMethod> ib_method_ops = new IBFEMethod(
            "IBFEMethod", app_initializer->getComponentDatabase("IBFEMethod"), meshes, app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
            "IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops, navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Configure the IBFE solver.
        ib_method_ops->registerLagBodyForceFunction(&block_tether_force_function, std::vector<unsigned int>(), NULL, 0);
        ib_method_ops->registerLagBodyForceFunction(&block_tether_force_function, std::vector<unsigned int>(), NULL, 1);
        ib_method_ops->registerLagBodyForceFunction(& beam_tether_force_function, std::vector<unsigned int>(), NULL, 2);
        ib_method_ops->registerPK1StressTensorFunction(&beam_PK1_stress_function, std::vector<unsigned int>(), NULL, 2);
        EquationSystems* block1_equation_systems = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* block2_equation_systems = ib_method_ops->getFEDataManager(1)->getEquationSystems();
        EquationSystems*   beam_equation_systems = ib_method_ops->getFEDataManager(2)->getEquationSystems();

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
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
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        AutoPtr<ExodusII_IO> block1_exodus_io(uses_exodus ? new ExodusII_IO(block1_mesh) : NULL);
        AutoPtr<ExodusII_IO> block2_exodus_io(uses_exodus ? new ExodusII_IO(block2_mesh) : NULL);
        AutoPtr<ExodusII_IO>   beam_exodus_io(uses_exodus ? new ExodusII_IO(  beam_mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                block1_exodus_io->write_timestep(block1_exodus_filename, *block1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                block2_exodus_io->write_timestep(block2_exodus_filename, *block2_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                beam_exodus_io  ->write_timestep(  beam_exodus_filename, *  beam_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
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

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout <<                                                    "\n";
            pout << "At end       of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout <<                                                    "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num%viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    block1_exodus_io->write_timestep(block1_exodus_filename, *block1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                    block2_exodus_io->write_timestep(block2_exodus_filename, *block2_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                    beam_exodus_io  ->write_timestep(  beam_exodus_filename, *  beam_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}// main
