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

// Headers for basic libMesh objects
#include <elem.h>
#include <mesh.h>
#include <mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IMPInitializer.h>
#include <ibamr/IMPMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Elasticity model data.
namespace ModelData
{
static const double r0    = 1.250; // radius of thoracic aorta (cm)
static const double w_int = 0.048; // thickness of the intima (cm)
static const double w_med = 0.118; // thickness of the media  (cm)
static const double w_adv = 0.093; // thickness of the adventitia (cm)
static const double w = w_int+w_med+w_adv;

static const int NUM_VARS = 9;
static const int       MU_IDX = 0;
static const int       K1_IDX = 1;
static const int       K2_IDX = 2;
static const int      PHI_IDX = 3;
static const int    KAPPA_IDX = 4;
static const int      R_F_IDX = 5;
static const int      M_F_IDX = 6;
static const int W_F4_MAX_IDX = 7;
static const int W_F6_MAX_IDX = 8;

// Stress tensor function.
void
PK1_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const Point& /*x*/,
    const Point& X,
    subdomain_id_type /*subdomain_id*/,
    std::vector<double>& internal_vars,
    double /*time*/,
    void* /*ctx*/)
{
    VectorValue<double> R = X - Point(0.0,0.0);
    const double r = R.size();
    R = R.unit();
    if (internal_vars.empty())  // ugly hacky hack
    {
        internal_vars.resize(NUM_VARS);
        if (r <= r0+w_int)
        {
            internal_vars[   MU_IDX] = 0.034;  // MPa
            internal_vars[   K1_IDX] = 4.34;   // MPa
            internal_vars[   K2_IDX] = 13.32;  // dimensionless
            internal_vars[  PHI_IDX] = 46.5;   // degrees
            internal_vars[KAPPA_IDX] = 0.20;   // dimensionless
            internal_vars[  R_F_IDX] = 1.0;    // dimensionless
            internal_vars[  M_F_IDX] = 0.014;  // dimensionless
        }
        else if (r <= r0+w_int+w_med)
        {
            internal_vars[   MU_IDX] = 0.028;  // MPa
            internal_vars[   K1_IDX] = 0.14;   // MPa
            internal_vars[   K2_IDX] = 11.90;  // dimensionless
            internal_vars[  PHI_IDX] = 38.4;   // degrees
            internal_vars[KAPPA_IDX] = 0.21;   // dimensionless
            internal_vars[  R_F_IDX] = 1.87;   // dimensionless
            internal_vars[  M_F_IDX] = 0.009;  // dimensionless
        }
        else if (r <= r0+w_int+w_med+w_adv)
        {
            internal_vars[   MU_IDX] = 0.020;  // MPa
            internal_vars[   K1_IDX] = 0.39;   // MPa
            internal_vars[   K2_IDX] = 6.79;   // dimensionless
            internal_vars[  PHI_IDX] = 52.3;   // degrees
            internal_vars[KAPPA_IDX] = 0.23;   // dimensionless
            internal_vars[  R_F_IDX] = 1.15;   // dimensionless
            internal_vars[  M_F_IDX] = 0.022;  // dimensionless
        }
        internal_vars[W_F4_MAX_IDX] = 0.0;
        internal_vars[W_F6_MAX_IDX] = 0.0;
    }

    const double& mu    = internal_vars[   MU_IDX];
    const double& k1    = internal_vars[   K1_IDX];
    const double& k2    = internal_vars[   K2_IDX];
    const double& phi   = internal_vars[  PHI_IDX];
    const double& kappa = internal_vars[KAPPA_IDX];
    const double& r_f   = internal_vars[  R_F_IDX];
    const double& m_f   = internal_vars[  M_F_IDX];
    double& W_f4_max = internal_vars[W_F4_MAX_IDX];
    double& W_f6_max = internal_vars[W_F6_MAX_IDX];

    // setup fiber axes
    static const TensorValue<double> II(1.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0);
    const TensorValue<double> R_cross( 0.0 , -R(2), +R(1),
                                       +R(2),  0.0 , -R(0),
                                       -R(1), +R(0),  0.0 );
    TensorValue<double> RR4 = cos(+phi)*II + sin(+phi)*R_cross + (1.0-cos(+phi))*outer_product(R,R);
    VectorValue<double> f4_0 = RR4*VectorValue<double>(-R(1),R(0),R(2));
    TensorValue<double> RR6 = cos(-phi)*II + sin(-phi)*R_cross + (1.0-cos(-phi))*outer_product(R,R);
    VectorValue<double> f6_0 = RR6*VectorValue<double>(-R(1),R(0),R(2));
    f4_0 = VectorValue<double>(-R(1),R(0),R(2));
    f6_0 = VectorValue<double>(-R(1),R(0),R(2));

    // compute invariants.
    const TensorValue<double> CC = FF.transpose()*FF;
    const double I1 = CC.tr();
    const double I4_star = kappa*I1 + (1.0-3.0*kappa)*f4_0*(CC*f4_0);
    const double I6_star = kappa*I1 + (1.0-3.0*kappa)*f6_0*(CC*f6_0);

    // compute energies and stresses.
    PP = mu*(FF-tensor_inverse_transpose(FF));
    double W_f4 = 0.0;
    if (I4_star-1.0 > 0.0)
    {
        W_f4 = 0.5*k1/k2*(exp(k2*(I4_star-1.0)*(I4_star-1.0))-1.0);
        PP += 2.0*k1*(I4_star-1.0)*exp(k2*(I4_star-1.0)*(I4_star-1.0))*(kappa*FF+(1.0-3.0*kappa)*FF*outer_product(f4_0,f4_0));
        W_f4_max = max(W_f4_max, W_f4);
    }
    double W_f6 = 0.0;
    if (I6_star-1.0 > 0.0)
    {
        W_f6 = 0.5*k1/k2*(exp(k2*(I6_star-1.0)*(I6_star-1.0))-1.0);
        PP += 2.0*k1*(I6_star-1.0)*exp(k2*(I6_star-1.0)*(I6_star-1.0))*(kappa*FF+(1.0-3.0*kappa)*FF*outer_product(f6_0,f6_0));
        W_f6_max = max(W_f6_max, W_f6);
    }
    PP *= 1.0e7;  // convert to CGS units
    return;
}// PK1_stress_function
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
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        Mesh mesh(NDIM);
        const double dx0 = 1.0/64.0;
        const double dx = input_db->getDouble("DX");
        const double ds0 = input_db->getDouble("MFAC")*dx0;
        string elem_type = input_db->getString("ELEM_TYPE");
        MeshTools::Generation::build_square(mesh,
                                            4*static_cast<int>((dx0/dx)*ceil(2.0*M_PI*r0/ds0/4.0)), static_cast<int>((dx0/dx)*ceil(w/ds0)),
                                            0.0, 2.0*M_PI*r0,
                                            0.0, w,
                                            Utility::string_to_enum<ElemType>(elem_type));
        for (MeshBase::node_iterator it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it)
        {
            Node& n = **it;
            Point X((r0+n(1))*cos(n(0)/r0), (r0+n(1))*sin(n(0)/r0),0.0);
            n = X;
        }
        mesh.prepare_for_use();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator("INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IMPMethod> ib_method_ops = new IMPMethod("IMPMethod", app_initializer->getComponentDatabase("IMPMethod"));
        Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops, navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>("CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>("GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Configure the IMP solver.
        Pointer<IMPInitializer> ib_initializer = new IMPInitializer("IMPInitializer", app_initializer->getComponentDatabase("IMPInitializer"), patch_hierarchy, gridding_algorithm);
        ib_initializer->registerMesh(&mesh);
        ib_method_ops->registerPK1StressTensorFunction(PK1_stress_function);
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Create Eulerian boundary condition specification objects.
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
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction("f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }
        if (input_db->keyExists("SourceFunction"))
        {
            Pointer<CartGridFunction> q_fcn = new muParserCartGridFunction("q_fcn", app_initializer->getComponentDatabase("SourceFunction"), grid_geometry);
            navier_stokes_integrator->registerFluidSourceFunction(q_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Open streams to save volume of structure.
        ofstream volume_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            volume_stream.open("volume.curve", ios_base::out | ios_base::trunc);
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
            if (dump_viz_data && uses_visit && (iteration_num%viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
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

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            volume_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}// main
