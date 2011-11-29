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
#include <boundary_info.h>
#include <exodusII_io.h>
#include <mesh.h>
#include <mesh_function.h>
#include <mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Elasticity model data.
namespace ModelData
{
// Stress tensor function for solid block.
void
block_PK1_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& /*FF*/,
    const Point& /*X*/,
    const Point& /*s*/,
    Elem* const /*elem*/,
    NumericVector<double>& /*X_vec*/,
    const vector<NumericVector<double>*>& /*system_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    PP.zero();
    return;
}// block_PK1_stress_function

// Tether (penalty) force function for the solid block.
static double kappa_s = 1.0e6;
void
block_tether_force_function(
    VectorValue<double>& F,
    const TensorValue<double>& /*FF*/,
    const Point& X,
    const Point& s,
    Elem* const /*elem*/,
    NumericVector<double>& /*X_vec*/,
    const vector<NumericVector<double>*>& /*system_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    F = kappa_s*(s-X);
    return;
}// block_tether_force_function

// Stress tensor function for thin beam.
static double mu_s, lambda_s;
void
beam_PK1_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const Point& /*X*/,
    const Point& /*s*/,
    Elem* const /*elem*/,
    NumericVector<double>& /*X_vec*/,
    const vector<NumericVector<double>*>& /*system_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    static const TensorValue<double> II(1.0,0.0,0.0,
                                        0.0,1.0,0.0,
                                        0.0,0.0,1.0);
    const TensorValue<double> CC = FF.transpose()*FF;
    const TensorValue<double> EE = 0.5*(CC - II);
    const TensorValue<double> SS = lambda_s*EE.tr()*II + 2.0*mu_s*EE;
    PP = FF*SS;
    return;
}// beam_PK1_stress_function
}
using namespace ModelData;

// Function prototypes
static ofstream drag_stream, lift_stream, A_x_posn_stream, A_y_posn_stream;
static int stress_diag_idx, stress_odiag_idx;
void
postprocess_data(
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
    Mesh& beam_mesh,
    EquationSystems* beam_equation_systems,
    FEDataManager* beam_fe_data_manager,
    Mesh& block_mesh,
    EquationSystems* block_equation_systems,
    FEDataManager* block_fe_data_manager,
    const int iteration_num,
    const double loop_time,
    const string& data_dump_dirname);

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
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string block_exodus_filename = app_initializer->getExodusIIFilename("block");
        const string  beam_exodus_filename = app_initializer->getExodusIIFilename("beam" );

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC")*dx;

        Mesh block_mesh(NDIM);
        string block_elem_type = input_db->getString("BLOCK_ELEM_TYPE");
        const double R = 0.05;
        if (block_elem_type == "TRI3" || block_elem_type == "TRI6")
        {
            const int num_circum_nodes = ceil(2.0*M_PI*R/ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0*M_PI*static_cast<double>(k)/static_cast<double>(num_circum_nodes);
                block_mesh.add_point(Point(R*cos(theta), R*sin(theta)));
            }
            TriangleInterface triangle(block_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(block_elem_type);
            triangle.desired_area() = sqrt(3.0)/4.0*ds*ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
            block_mesh.prepare_for_use();
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0*M_PI*R/ds;
            const int r = log2(0.25*num_circum_segments);
            MeshTools::Generation::build_sphere(block_mesh, R, r, Utility::string_to_enum<ElemType>(block_elem_type));
        }
        for (MeshBase::node_iterator n_it = block_mesh.nodes_begin(); n_it != block_mesh.nodes_end(); ++n_it)
        {
            Node& n = **n_it;
            n(0) += 0.2;
            n(1) += 0.2;
        }

        Mesh beam_mesh(NDIM);
        string beam_elem_type = input_db->getString("BEAM_ELEM_TYPE");
        MeshTools::Generation::build_square(beam_mesh,
                                            ceil(0.4/ds), ceil(0.02/ds),
                                            0.2, 0.6, 0.19, 0.21,
                                            Utility::string_to_enum<ElemType>(beam_elem_type));
        beam_mesh.prepare_for_use();

        vector<Mesh*> meshes(2);
        meshes[0] = &block_mesh;
        meshes[1] = & beam_mesh;

        mu_s     = input_db->getDouble("MU_S");
        lambda_s = input_db->getDouble("LAMBDA_S");
        kappa_s  = input_db->getDouble("KAPPA_S");
        
        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator", app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n" <<
                       "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops = new IBFEMethod(
            "IBFEMethod", app_initializer->getComponentDatabase("IBFEMethod"), meshes, app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        Pointer<IBHierarchyIntegrator> time_integrator = new IBHierarchyIntegrator(
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
        ib_method_ops->registerPK1StressTensorFunction(&block_PK1_stress_function, std::vector<unsigned int>(), NULL, 0);
        ib_method_ops->registerPK1StressTensorFunction(& beam_PK1_stress_function, std::vector<unsigned int>(), NULL, 1);
        ib_method_ops->registerLagBodyForceFunction(&block_tether_force_function, std::vector<unsigned int>(), NULL, 0);
        FEDataManager* block_fe_data_manager = ib_method_ops->getFEDataManager(0);
        FEDataManager*  beam_fe_data_manager = ib_method_ops->getFEDataManager(1);
        EquationSystems* block_equation_systems = block_fe_data_manager->getEquationSystems();
        EquationSystems*  beam_equation_systems =  beam_fe_data_manager->getEquationSystems();

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
        TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM> u_bc_coefs;
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
        AutoPtr<ExodusII_IO> block_exodus_io(uses_exodus ? new ExodusII_IO(block_mesh) : NULL);
        AutoPtr<ExodusII_IO>  beam_exodus_io(uses_exodus ? new ExodusII_IO( beam_mesh) : NULL);

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
                block_exodus_io->write_timestep(block_exodus_filename, *block_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                beam_exodus_io ->write_timestep( beam_exodus_filename, * beam_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
            }
        }

        // Setup variables to store fluid stress components.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        stress_diag_idx  = var_db->registerVariableAndContext(new CellVariable<NDIM,double>("stress_diag_var" ,2), navier_stokes_integrator->getScratchContext(), IntVector<NDIM>(3));  // specialized for NDIM == 2
        stress_odiag_idx = var_db->registerVariableAndContext(new NodeVariable<NDIM,double>("stress_odiag_var",1), navier_stokes_integrator->getScratchContext(), IntVector<NDIM>(3));
        
        // Open streams to save lift and drag coefficients.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.open("C_D.curve", ios_base::out | ios_base::trunc);
            lift_stream.open("C_L.curve", ios_base::out | ios_base::trunc);
            A_x_posn_stream.open("A_x.curve", ios_base::out | ios_base::trunc);
            A_y_posn_stream.open("A_y.curve", ios_base::out | ios_base::trunc);
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
                    block_exodus_io->write_timestep(block_exodus_filename, *block_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                    beam_exodus_io ->write_timestep( beam_exodus_filename, * beam_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
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
            if (dump_postproc_data && (iteration_num%postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting state data...\n\n";
                postprocess_data(patch_hierarchy, navier_stokes_integrator,
                                 beam_mesh, beam_equation_systems, beam_fe_data_manager,
                                 block_mesh, block_equation_systems, block_fe_data_manager,
                                 iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.close();
            lift_stream.close();
            A_x_posn_stream.close();
            A_y_posn_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
}// main

void
postprocess_data(
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
    Mesh& /*beam_mesh*/,
    EquationSystems* beam_equation_systems,
    FEDataManager* beam_fe_data_manager,
    Mesh& /*block_mesh*/,
    EquationSystems* block_equation_systems,
    FEDataManager* block_fe_data_manager,
    const int /*iteration_num*/,
    const double loop_time,
    const string& /*data_dump_dirname*/)
{
    // Compute the position of the tip of the flexible beam.
    System& X_system = beam_equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    NumericVector<double>* X_vec = X_system.solution.get();
    AutoPtr<NumericVector<Number> > X_serial_vec = NumericVector<Number>::build();
    X_serial_vec->init(X_vec->size(), true, SERIAL);
    X_vec->localize(*X_serial_vec);
    DofMap& X_dof_map = X_system.get_dof_map();
    vector<unsigned int> vars(2);
    vars[0] = 0; vars[1] = 1;
    MeshFunction X_fcn(*beam_equation_systems, *X_serial_vec, X_dof_map, vars);
    X_fcn.init();
    DenseVector<double> X_A(2);
    X_fcn(Point(0.6,0.2,0), 0.0, X_A);
    if (SAMRAI_MPI::getRank() == 0)
    {
        A_x_posn_stream.precision(12);
        A_x_posn_stream.setf(ios::fixed,ios::floatfield);
        A_x_posn_stream << loop_time << " " << X_A(0) << endl;
        A_y_posn_stream.precision(12);
        A_y_posn_stream.setf(ios::fixed,ios::floatfield);
        A_y_posn_stream << loop_time << " " << X_A(1) << endl;
    }

    // Compute the lift and drag coefficients.
    DenseVector<double> F(2);
    
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    
    const Pointer<SideVariable<NDIM,double> > u_var = navier_stokes_integrator->getVelocityVariable();
    const Pointer<CellVariable<NDIM,double> > p_var = navier_stokes_integrator->getPressureVariable();
    const Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
    const Pointer<VariableContext> scratch_ctx = navier_stokes_integrator->getScratchContext();
    
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);
    const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, current_ctx);
    
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(   u_scratch_idx, loop_time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData( stress_diag_idx, loop_time);
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(stress_odiag_idx, loop_time);
    }

    HierarchySideDataOpsReal<NDIM,double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
    hier_sc_data_ops.copyData(u_scratch_idx, u_current_idx);
    
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent u_bc_component(u_scratch_idx);
    HierarchyGhostCellInterpolation u_bc_fill_op;
    u_bc_fill_op.initializeOperatorState(u_bc_component, patch_hierarchy, coarsest_ln, finest_ln);
    u_bc_fill_op.fillData(loop_time);

    const double mu = navier_stokes_integrator->getINSProblemCoefs()->getMu();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            Pointer<SideData<NDIM,double> >            u_data = patch->getPatchData(   u_scratch_idx);
            Pointer<CellData<NDIM,double> >            p_data = patch->getPatchData(   p_current_idx);
            Pointer<CellData<NDIM,double> >  stress_diag_data = patch->getPatchData( stress_diag_idx);
            Pointer<NodeData<NDIM,double> > stress_odiag_data = patch->getPatchData(stress_odiag_idx);
            for (Box<NDIM>::Iterator b(patch_box); b; b++)
            {
                const Index<NDIM>& i = b();
                for (unsigned int d = 0; d < 2; ++d)
                {
                    (*stress_diag_data)(i,d) = -(*p_data)(i) + 2.0*mu*((*u_data)(SideIndex<NDIM>(i,d,SideIndex<NDIM>::Upper)) - (*u_data)(SideIndex<NDIM>(i,d,SideIndex<NDIM>::Lower)))/dx[d];
                }
            }
            for (Box<NDIM>::Iterator b(NodeGeometry<NDIM>::toNodeBox(patch_box)); b; b++)
            {
                const Index<NDIM>& i = b();
                const double du_dy = ((*u_data)(SideIndex<NDIM>(i,0,SideIndex<NDIM>::Lower)) - (*u_data)(SideIndex<NDIM>(i-Index<NDIM>(0,1),0,SideIndex<NDIM>::Lower)))/dx[1];
                const double dv_dx = ((*u_data)(SideIndex<NDIM>(i,1,SideIndex<NDIM>::Lower)) - (*u_data)(SideIndex<NDIM>(i-Index<NDIM>(1,0),1,SideIndex<NDIM>::Lower)))/dx[0];
                const NodeIndex<NDIM> n_i(i,0);
                (*stress_odiag_data)(n_i) = mu*(du_dy+dv_dx);
            }
        }
    }

    std::vector<InterpolationTransactionComponent> stress_bc_components(2);
    stress_bc_components[0] = InterpolationTransactionComponent( stress_diag_idx);
    stress_bc_components[1] = InterpolationTransactionComponent(stress_odiag_idx);
    HierarchyGhostCellInterpolation stress_bc_fill_op;
    stress_bc_fill_op.initializeOperatorState(stress_bc_components, patch_hierarchy, coarsest_ln, finest_ln);
    stress_bc_fill_op.fillData(loop_time);

    const int dim = NDIM;
    AutoPtr<QBase> qrule_face = QBase::build(QGAUSS, dim-1, FIFTH);
    
    // Integrate over the block boundary.
    {
        NumericVector<double>& X_ghost_vec = *block_fe_data_manager->buildGhostedCoordsVector();
        System& system = block_equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
        const DofMap& dof_map = system.get_dof_map();
        blitz::Array<std::vector<unsigned int>,1> side_dof_indices(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) side_dof_indices(d).reserve(9);
        AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
        fe_face->attach_quadrature_rule(qrule_face.get());
        const std::vector<Point>& q_point_face = fe_face->get_xyz();
        const std::vector<double>& JxW_face = fe_face->get_JxW();
        const std::vector<Point>& normal_face = fe_face->get_normals();
        const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
        const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();
        
        const blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_element_map = block_fe_data_manager->getActivePatchElementMap();
        const int level_num = block_fe_data_manager->getLevelNumber();
        VectorValue<double> F_da;
        TensorValue<double> FF, FF_inv_trans, Stress;
        Point X_qp;
        std::vector<Point> elem_X;
        blitz::Array<double,2> X_node, X_node_side;
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(level_num);
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const blitz::Array<Elem*,1>& patch_elems = active_patch_element_map(local_patch_num);
            const int num_active_patch_elems = patch_elems.size();
            if (num_active_patch_elems == 0) continue;
            
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            
            // Setup vectors to store the values of the stress tensor and X at
            // the quadrature points.  We compute a somewhat conservative upper
            // bound on the number of quadrature points to try to avoid
            // reallocations.
            static const unsigned int n_qp_estimate = (NDIM == 2 ? 12 : 12*12);
            std::vector<double> stress_diag_bdry, stress_odiag_bdry;
            stress_diag_bdry .reserve(2*n_qp_estimate*num_active_patch_elems);
            stress_odiag_bdry.reserve(1*n_qp_estimate*num_active_patch_elems);
            std::vector<double> X_bdry;
            X_bdry.reserve(NDIM*n_qp_estimate*num_active_patch_elems);

            // Loop over the elements and compute the values to be spread and
            // the positions of the quadrature points.
            int qp_offset = 0;
            for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems(e_idx);
                                
                // Loop over the element boundaries.
                for (unsigned short int side = 0; side < elem->n_sides(); ++side)
                {
                    // Skip non-physical boundaries.
                    bool at_physical_bdry = elem->neighbor(side) == NULL;
                    if (!at_physical_bdry) continue;
                    
                    AutoPtr<Elem> side_elem = elem->build_side(side);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        dof_map.dof_indices(side_elem.get(), side_dof_indices(d), d);
                    }
                    get_values_for_interpolation(X_node_side, X_ghost_vec, side_dof_indices);
                    
                    fe_face->reinit(elem, side);
                    
                    const unsigned int n_qp = qrule_face->n_points();
                    stress_diag_bdry .resize( stress_diag_bdry.size()+2*n_qp,0.0);
                    stress_odiag_bdry.resize(stress_odiag_bdry.size()+  n_qp,0.0);
                    X_bdry.resize(X_bdry.size()+NDIM*n_qp);
                    for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                    {
                        const int idx = NDIM*qp_offset;
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            X_bdry[idx+i] = X_qp(i);
                        }
                    }
                }
            }
            
            if (qp_offset == 0) continue;
            
            // Interpolate the stresses from the grid.
            const std::string& interp_weighting_fcn = block_fe_data_manager->getInterpWeightingFunction();
            const Box<NDIM> interp_box = patch->getBox();
            Pointer<CellData<NDIM,double> > stress_diag_data = patch->getPatchData(stress_diag_idx);
            LEInteractor::interpolate( stress_diag_bdry, 2, X_bdry, NDIM, stress_diag_data, patch, interp_box, interp_weighting_fcn);
            Pointer<NodeData<NDIM,double> > stress_odiag_data = patch->getPatchData(stress_odiag_idx);
            LEInteractor::interpolate(stress_odiag_bdry, 1, X_bdry, NDIM, stress_odiag_data, patch, interp_box, interp_weighting_fcn);

            // Accumulate values.
            qp_offset = 0;
            for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems(e_idx);
                                
                // Loop over the element boundaries.
                for (unsigned short int side = 0; side < elem->n_sides(); ++side)
                {
                    // Skip non-physical boundaries.
                    bool at_physical_bdry = elem->neighbor(side) == NULL;
                    if (!at_physical_bdry) continue;
                    
                    AutoPtr<Elem> side_elem = elem->build_side(side);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        dof_map.dof_indices(side_elem.get(), side_dof_indices(d), d);
                    }
                    get_values_for_interpolation(X_node_side, X_ghost_vec, side_dof_indices);
                    
                    fe_face->reinit(elem, side);
                    
                    const unsigned int n_qp = qrule_face->n_points();
                    for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                    {
                        const Point& s_qp = q_point_face[qp];
                        interpolate(X_qp,qp,X_node,phi_face);
                        jacobian(FF,qp,X_node,dphi_face);
                        const double J = std::abs(FF.det());
                        tensor_inverse_transpose(FF_inv_trans,FF,NDIM);

                        Stress.zero();
                        Stress(0,0) =  stress_diag_bdry[2*qp_offset  ];
                        Stress(1,1) =  stress_diag_bdry[2*qp_offset+1];
                        Stress(0,1) = stress_odiag_bdry[  qp_offset  ];
                        Stress(1,0) = stress_odiag_bdry[  qp_offset  ];

                        F_da = Stress * J * JxW_face[qp] * FF_inv_trans * normal_face[qp];

                        if (s_qp(0) < 0.2 || !(s_qp(1) <= 0.21 && s_qp(1) >= 0.19)) for (unsigned int d = 0; d < 2; ++d) F(d) += F_da(d);
                    }
                }
            }
        }
    }
    
    // Integrate over the beam boundary.
    {
        NumericVector<double>& X_ghost_vec = *beam_fe_data_manager->buildGhostedCoordsVector();
        System& system = beam_equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
        const DofMap& dof_map = system.get_dof_map();
        blitz::Array<std::vector<unsigned int>,1> side_dof_indices(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) side_dof_indices(d).reserve(9);
        AutoPtr<FEBase> fe_face(FEBase::build(dim, dof_map.variable_type(0)));
        fe_face->attach_quadrature_rule(qrule_face.get());
        const std::vector<Point>& q_point_face = fe_face->get_xyz();
        const std::vector<double>& JxW_face = fe_face->get_JxW();
        const std::vector<Point>& normal_face = fe_face->get_normals();
        const std::vector<std::vector<double> >& phi_face = fe_face->get_phi();
        const std::vector<std::vector<VectorValue<double> > >& dphi_face = fe_face->get_dphi();
        
        const blitz::Array<blitz::Array<Elem*,1>,1>& active_patch_element_map = beam_fe_data_manager->getActivePatchElementMap();
        const int level_num = beam_fe_data_manager->getLevelNumber();
        VectorValue<double> F_da;
        TensorValue<double> FF, FF_inv_trans, Stress;
        Point X_qp;
        std::vector<Point> elem_X;
        blitz::Array<double,2> X_node, X_node_side;
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(level_num);
        int local_patch_num = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
        {
            // The relevant collection of elements.
            const blitz::Array<Elem*,1>& patch_elems = active_patch_element_map(local_patch_num);
            const int num_active_patch_elems = patch_elems.size();
            if (num_active_patch_elems == 0) continue;
            
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            
            // Setup vectors to store the values of the stress tensor and X at
            // the quadrature points.  We compute a somewhat conservative upper
            // bound on the number of quadrature points to try to avoid
            // reallocations.
            static const unsigned int n_qp_estimate = (NDIM == 2 ? 12 : 12*12);
            std::vector<double> stress_diag_bdry, stress_odiag_bdry;
            stress_diag_bdry .reserve(2*n_qp_estimate*num_active_patch_elems);
            stress_odiag_bdry.reserve(1*n_qp_estimate*num_active_patch_elems);
            std::vector<double> X_bdry;
            X_bdry.reserve(NDIM*n_qp_estimate*num_active_patch_elems);

            // Loop over the elements and compute the values to be spread and
            // the positions of the quadrature points.
            int qp_offset = 0;
            for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems(e_idx);
                                
                // Loop over the element boundaries.
                for (unsigned short int side = 0; side < elem->n_sides(); ++side)
                {
                    // Skip non-physical boundaries.
                    bool at_physical_bdry = elem->neighbor(side) == NULL;
                    if (!at_physical_bdry) continue;
                    
                    AutoPtr<Elem> side_elem = elem->build_side(side);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        dof_map.dof_indices(side_elem.get(), side_dof_indices(d), d);
                    }
                    get_values_for_interpolation(X_node_side, X_ghost_vec, side_dof_indices);
                    
                    fe_face->reinit(elem, side);
                    
                    const unsigned int n_qp = qrule_face->n_points();
                    stress_diag_bdry .resize( stress_diag_bdry.size()+2*n_qp,0.0);
                    stress_odiag_bdry.resize(stress_odiag_bdry.size()+  n_qp,0.0);
                    X_bdry.resize(X_bdry.size()+NDIM*n_qp);
                    for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                    {
                        const int idx = NDIM*qp_offset;
                        for (unsigned int i = 0; i < NDIM; ++i)
                        {
                            X_bdry[idx+i] = X_qp(i);
                        }
                    }
                }
            }
            
            if (qp_offset == 0) continue;
            
            // Interpolate the stresses from the grid.
            const std::string& interp_weighting_fcn = beam_fe_data_manager->getInterpWeightingFunction();
            const Box<NDIM> interp_box = patch->getBox();
            Pointer<CellData<NDIM,double> > stress_diag_data = patch->getPatchData(stress_diag_idx);
            LEInteractor::interpolate( stress_diag_bdry, 2, X_bdry, NDIM, stress_diag_data, patch, interp_box, interp_weighting_fcn);
            Pointer<NodeData<NDIM,double> > stress_odiag_data = patch->getPatchData(stress_odiag_idx);
            LEInteractor::interpolate(stress_odiag_bdry, 1, X_bdry, NDIM, stress_odiag_data, patch, interp_box, interp_weighting_fcn);

            // Accumulate values.
            qp_offset = 0;
            for (int e_idx = 0; e_idx < num_active_patch_elems; ++e_idx)
            {
                Elem* const elem = patch_elems(e_idx);
                                
                // Loop over the element boundaries.
                for (unsigned short int side = 0; side < elem->n_sides(); ++side)
                {
                    // Skip non-physical boundaries.
                    bool at_physical_bdry = elem->neighbor(side) == NULL;
                    if (!at_physical_bdry) continue;
                    
                    AutoPtr<Elem> side_elem = elem->build_side(side);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        dof_map.dof_indices(side_elem.get(), side_dof_indices(d), d);
                    }
                    get_values_for_interpolation(X_node_side, X_ghost_vec, side_dof_indices);
                    
                    fe_face->reinit(elem, side);
                    
                    const unsigned int n_qp = qrule_face->n_points();
                    for (unsigned int qp = 0; qp < n_qp; ++qp, ++qp_offset)
                    {
                        const Point& s_qp = q_point_face[qp];
                        interpolate(X_qp,qp,X_node,phi_face);
                        jacobian(FF,qp,X_node,dphi_face);
                        const double J = std::abs(FF.det());
                        tensor_inverse_transpose(FF_inv_trans,FF,NDIM);

                        Stress.zero();
                        Stress(0,0) =  stress_diag_bdry[2*qp_offset  ];
                        Stress(1,1) =  stress_diag_bdry[2*qp_offset+1];
                        Stress(0,1) = stress_odiag_bdry[  qp_offset  ];
                        Stress(1,0) = stress_odiag_bdry[  qp_offset  ];

                        F_da = Stress * J * JxW_face[qp] * FF_inv_trans * normal_face[qp];

                        Point r(0.2,0.2);
                        r = r - s_qp;
                        if (r.size() > 0.05) for (unsigned int d = 0; d < 2; ++d) F(d) += F_da(d);
                    }
                }
            }
        }
    }

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(   u_scratch_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData( stress_diag_idx);
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(stress_odiag_idx);
    }

    if (SAMRAI_MPI::getRank() == 0)
    {
        drag_stream.precision(12);
        drag_stream.setf(ios::fixed,ios::floatfield);
        drag_stream << loop_time << " " << F(0) << endl;
        lift_stream.precision(12);
        lift_stream.setf(ios::fixed,ios::floatfield);
        lift_stream << loop_time << " " << F(1) << endl;
    }    
    return;
}// postprocess_data 
