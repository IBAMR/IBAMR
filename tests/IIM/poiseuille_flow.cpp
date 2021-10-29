//  2016-2021, created by Ebrahim (Amin) Kolahdouz
// Fully developed flow in a 2D (slanted) channel
// See Sec. 4.1.1 of "An immersed interface method for discrete surfaces"
//  by Ebrahim M. Kolahdouz et al., Journal of Computational Physics 400 (2020) 108854
// This code was run with a tight solver/preconditioner tolerance
// through the command line flag: -stokes_ksp_rtol 1.0e-10 -ksp_rtol 1.0e-10

// This code was tested against the commit 57fb379454ea3f8f50c476f10226cb6b520a11a0
// which is currently on branch iim-1 at https://github.com/drwells/IBAMR

// Config files

#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

//#include <FeedbackForcer.h>


// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include "libmesh/edge_edge2.h"
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>


// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/SpongeLayerForceFunction.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibamr/StaggeredStokesOpenBoundaryStabilizer.h>

static const int FEDIM = 1;

// Elasticity model data for thin body.
namespace ModelData
{
    // Tether (penalty) force function for the thin body.
    static double kappa_s_thin = 0.0;
    static double X1_lower = 0.0;
    static double X1_upper = 0.0;
    static double eta_s_thin = 0.0;
    static double theta = 0.0;
    static double FFAC = 0.0;
    static double fac = 0.0;
    static double y_position_low = 0.0;
    static double y_position_up = 0.0;
    static double L = 0.0;
    static double D = 0.0;
    static double p_e =  0.0;
    static double MU = 0.0;
    static double Re = 0.0;
    
    
void tether_body_force_function_thin(VectorValue<double>& F,
                      const VectorValue<double>& n,
                      const VectorValue<double>& /*N*/,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x,
                      const libMesh::Point& X,
                      Elem* const /*elem*/,
                      const unsigned short /*side*/,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double /*time*/,
                      void* /*ctx*/)
{
    const std::vector<double>& U = *var_data[0];
 
    double u_bndry_n = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
		u_bndry_n += n(d) * U[d];
	}

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
         F(d) =  kappa_s_thin * (X(d) - x(d)) + eta_s_thin * (0.0 - u_bndry_n) * n(d);
    }
    return; 
}
    static double aa = 0.4;
    static double bb = 0.6;
    static double x_loc;
    static double x_loc_min = 0.0;
    static double x_loc_max = 0.0;
    static double y_loc_min = 0.0;
    static double y_loc_max = 0.0;
}
using namespace ModelData;

static ofstream  p_norm_stream;

void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);



void compute_velocity_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int u_idx,
                              const double data_time,
                              const string& data_dump_dirname);
                              
void compute_pressure_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int p_idx,
                              const double data_time,
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
main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    
    pout << std::setprecision(10);
    plog << std::setprecision(10);
    
    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-10");
    PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-10");

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus =
        dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string lower_exodus_filename = app_initializer->getExodusIIFilename("lower");
        const string upper_exodus_filename = app_initializer->getExodusIIFilename("upper");

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval =
        app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname =
        app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) &&
            !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();
        
        // Create a hybrid FE mesh.
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        MU = input_db->getDouble("MU");
        fac = input_db->getDouble("FAC");
        Re = input_db->getDouble("Re");
        
        D = input_db->getDouble("D"); // channel parameter (cm)
        L = input_db->getDouble("L"); // channel length (cm)
        p_e = input_db->getDouble("P_E"); // channel length (cm)
        const double H = input_db->getDouble("H"); // Height (cm)
        FFAC = input_db->getDouble("FFAC"); // Height (cm)
        theta = input_db->getDouble("THETA");
        // Thin structure part of the mesh.
        string elem_type = input_db->getString("ELEM_TYPE");
        X1_lower = input_db->getDouble("X1_lower");
        X1_upper = input_db->getDouble("X1_upper");
        x_loc_max = aa*L;
        x_loc_min = bb*L;
        const unsigned int  nn = ceil(L/(cos(theta)*ds));
        
        
        Mesh lower_mesh(init.comm(), FEDIM);
       
        y_position_low = 0.5*H - 0.5*D/cos(theta) -0.5*L*tan(theta);
        int node_id = 0;
     
        lower_mesh.reserve_nodes (nn + 1);
        lower_mesh.reserve_elem (nn);
		for (unsigned int i = 0; i <= nn; i++)
        {     			
			lower_mesh.add_point (libMesh::Point(L*static_cast<Real>(i)/static_cast<Real>(nn), y_position_low + tan(theta)*L*static_cast<Real>(i)/static_cast<Real>(nn), 0.0), node_id++);

		}
        
        BoundaryInfo& boundary_info_lower = lower_mesh.get_boundary_info(); 
        
        for (unsigned int i=0; i<nn; i++)
        {
                     Elem * elem = lower_mesh.add_elem (new Edge2);
                     elem->set_node(0) = lower_mesh.node_ptr(i);
                    elem->set_node(1) = lower_mesh.node_ptr(i+1);
                    if (i == 0)
                     boundary_info_lower.add_side(elem, 0, 0); 
					if (i == (nn-1))
                      boundary_info_lower.add_side(elem, 1, 1);
        }
        lower_mesh.prepare_for_use();

        Mesh upper_mesh(init.comm(), FEDIM);

        y_position_up = y_position_low + D/cos(theta);
        node_id = 0;
        upper_mesh.reserve_nodes (nn + 1);
        upper_mesh.reserve_elem (nn);

		for (unsigned int i = 0; i <= nn; i++)
        {     
			upper_mesh.add_point (libMesh::Point(L*static_cast<Real>(i)/static_cast<Real>(nn), y_position_up + tan(theta)*L*static_cast<Real>(i)/static_cast<Real>(nn) , 0.0), node_id++);
		}

        
        BoundaryInfo& boundary_info_upper = upper_mesh.get_boundary_info(); 
        for (unsigned int i=0; i<nn; i++)
        {
                     Elem * elem = upper_mesh.add_elem (new Edge2);
                     elem->set_node(0) = upper_mesh.node_ptr(i);
                    elem->set_node(1) = upper_mesh.node_ptr(i+1);
                    if (i == 0)
                     boundary_info_upper.add_side(elem, 0, 0); 
                    if (i == (nn-1))
                      boundary_info_upper.add_side(elem, 1, 1);
        }
        
        
        
        upper_mesh.prepare_for_use();
        const int LOWER_MESH_ID = 0;
        const int UPPER_MESH_ID = 1;      
        vector<MeshBase*> meshes(2);
        meshes[LOWER_MESH_ID] = &lower_mesh;
        meshes[UPPER_MESH_ID] = &upper_mesh;


        // Setup data for imposing constraints.
        kappa_s_thin = input_db->getDoubleWithDefault("KAPPA_S_THIN", kappa_s_thin);
        eta_s_thin = input_db->getDoubleWithDefault("ETA_S_THIN", eta_s_thin);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator =
        new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        
        Pointer<IBAMR::IIMethod> ib_method_ops =
        new IIMethod("IIMethod", app_initializer->getComponentDatabase("IIMethod"),
                       meshes, app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
                       
        Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
            "IBHierarchyIntegrator",
            app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops,
            navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
        new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
        new StandardTagAndInitialize<NDIM>(
                "StandardTagAndInitialize", time_integrator,
                app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"),
            error_detector, box_generator, load_balancer);

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1);
        sys_data[0] = SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars);
        IIMethod::LagSurfaceForceFcnData  tether_body_force_thin_data(tether_body_force_function_thin, sys_data);
 
        
        pair<int,int> thin_id_range = make_pair(0,2);
        for (int id = thin_id_range.first; id<thin_id_range.second; id++)
        {
            ib_method_ops->registerLagSurfaceForceFunction(tether_body_force_thin_data, id);

        }
        
        // Create Eulerian initial condition specification objects.

		Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
																		"u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"),
																		grid_geometry);
		navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        

		Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
			"p_init", app_initializer->getComponentDatabase("PressureInitialConditions"),
			grid_geometry);
		navier_stokes_integrator->registerPressureInitialConditions(p_init);

        
        // Create Eulerian boundary condition specification objects.
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(
            NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));
        const bool periodic_domain = grid_geometry->getPeriodicShift().min() > 0;
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
                u_bc_coefs[d] = new muParserRobinBcCoefs(
                                                         bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name),
                                                         grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
            
            
            if (input_db->keyExists("BoundaryStabilization"))
            {
                time_integrator->registerBodyForceFunction(new StaggeredStokesOpenBoundaryStabilizer(
                    "BoundaryStabilization",
                    app_initializer->getComponentDatabase("BoundaryStabilization"),
                    navier_stokes_integrator,
                    grid_geometry));
            }
        }
        // Create Eulerian body force function specification objects.


        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer =
        app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        UniquePtr<ExodusII_IO>  lower_exodus_io(uses_exodus ? new ExodusII_IO(lower_mesh) : NULL);
        UniquePtr<ExodusII_IO>  upper_exodus_io(uses_exodus ? new ExodusII_IO(upper_mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();
        
        
        // Set up locations to get velocity data.
        x_loc = input_db->getDoubleWithDefault("X_LOC", L / 2);
        y_loc_min = input_db->getDoubleWithDefault("Y_LOC_MIN", 0 );
        y_loc_max = input_db->getDoubleWithDefault("Y_LOC_MAX", H );
        
        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);
        
        // Write out initial visualization data.
        EquationSystems* lower_equation_systems =
        ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* upper_equation_systems =
        ib_method_ops->getFEDataManager(1)->getEquationSystems();
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
                lower_exodus_io->write_timestep(lower_exodus_filename, *lower_equation_systems,
                                                iteration_num / viz_dump_interval + 1,
                                                loop_time);
                upper_exodus_io->write_timestep(upper_exodus_filename, *upper_equation_systems,
                                                iteration_num / viz_dump_interval + 1,
                                                loop_time);
            }
        }
        
    

        // Open streams to save lift and drag coefficients and the norms of the
        // velocity.
        if (SAMRAI_MPI::getRank() == 0)
        {
			p_norm_stream.open("p_norm.curve");
        }
        
        // Main time step loop.

        double loop_time_end = time_integrator->getEndTime();
        
        
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";
            
            
            
			//****** Check to see if the solution is reached steady state ***********//           
			
            Pointer<hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
            Pointer<VariableContext> current_ctx = time_integrator->getCurrentContext();
            
            
            const Pointer<hier::Variable<NDIM> > p_var = time_integrator->getPressureVariable();
			const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();

            // ************************************************************************//
            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num,
                                                     loop_time);
                }
                if (uses_exodus)
                {
                    lower_exodus_io->write_timestep(
                                                    lower_exodus_filename, *lower_equation_systems,
                                                    iteration_num / viz_dump_interval + 1, loop_time);
                    upper_exodus_io->write_timestep(
                                                    upper_exodus_filename, *upper_equation_systems,
                                                    iteration_num / viz_dump_interval + 1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname,
                                                               iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            pout << "\nWriting state data...\n\n";
                   
                   
			postprocess_data(patch_hierarchy,
							 navier_stokes_integrator,
							 lower_mesh,
							 lower_equation_systems,
							  iteration_num,
							 loop_time,
							  postproc_data_dump_dirname);
     
            if (dump_postproc_data &&
            (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {

                                         
            }
        }

     // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.resetLevels(finest_ln, finest_ln);
		 Pointer<hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
        
        

        const Pointer<hier::Variable<NDIM> > p_var = time_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();

        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
        
        
        
        const int coarsest_ln = 0;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx, loop_time);
        }
        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);
        
        
       compute_pressure_profile(patch_hierarchy, p_idx, loop_time,
                                         postproc_data_dump_dirname);
                                         
       compute_velocity_profile(patch_hierarchy, u_idx, loop_time,
                                         postproc_data_dump_dirname);
        
        
        
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
            HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
            pout << std::setprecision(16) << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
                  << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
             
             
             
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        pout << "Error in p at time " << loop_time - 0.5 * dt << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                 
                 
             pout<< " MU = "<< MU <<"\n"
                  << "  dx:  " << dx << "\n"
                  << " theta = "<<theta<<"\n"
					<< "  dt: " << dt << "\n";

             if  (input_db->getBool("USE_VELOCITY_JUMP_CONDITIONS"))
					pout<<" Using the jump condition"<<"\n";
			else
				    pout<<" Using regular IB"<<"\n";
				    

       
        if (dump_viz_data && uses_visit)
        {
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }
        
        if (SAMRAI_MPI::getRank() == 0)
        {
            p_norm_stream.close();

        }
        
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        
    } // cleanup dynamically allocated objects prior to shutdown
    
    SAMRAIManager::shutdown();
    return 0;
} // main

void compute_velocity_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int u_idx,
                              const double data_time,
                              const string& data_dump_dirname)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    
    HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
    hier_math_ops.resetLevels(finest_ln, finest_ln);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const double X_min[2] = { 0.2*L , 0.0 };
    const double X_max[2] = { 0.8*L  , L};
    
    double u_Eulerian_L2_norm = 0.0;
    double u_Eulerian_max_norm = 0.0;
    int N_max = 0;
    vector<double> pos_values;
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();

            const double* const patch_dx = patch_geom->getDx();
            

            // Entire box containing the required data.
            Box<NDIM> box(IndexUtilities::getCellIndex(&X_min[0], patch_x_lower, patch_x_upper,
                                                       patch_dx, patch_lower, patch_upper),
                          IndexUtilities::getCellIndex(&X_max[0], patch_x_lower, patch_x_upper,
                                                       patch_dx, patch_lower, patch_upper));
            // Part of the box on this patch
            Box<NDIM> trim_box = patch_box * box;
            BoxList<NDIM> iterate_box_list = trim_box;
            
            // Trim the box covered by the finer region
            BoxList<NDIM> covered_boxes;
            if (ln < finest_ln)
            {
                BoxArray<NDIM> refined_region_boxes;
                Pointer<PatchLevel<NDIM> > next_finer_level =
                patch_hierarchy->getPatchLevel(ln + 1);
                refined_region_boxes = next_finer_level->getBoxes();
                refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> covered_box = trim_box * refined_box;
                    covered_boxes.unionBoxes(covered_box);
                }
            }
            iterate_box_list.removeIntersections(covered_boxes);
            
            // Loop over the boxes and store the location and interpolated value.
            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);
            const Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);

            for (BoxList<NDIM>::Iterator lit(iterate_box_list); lit; lit++)
            {
                const Box<NDIM>& iterate_box = *lit;
                for (Box<NDIM>::Iterator bit(iterate_box); bit; bit++)
                {
                    const CellIndex<NDIM>& lower_idx = *bit;
                   // CellIndex<NDIM> upper_idx = lower_idx;
                   // upper_idx(0) += 1;
                    const double yu =
                    patch_x_lower[1] + patch_dx[1] * (lower_idx(1) - patch_lower(1) + 0.5);
                    const double xu =
                    patch_x_lower[0] + patch_dx[0] * (lower_idx(0) - patch_lower(0));
                   // upper_idx(0) += 1;
                    const double yv =
                    patch_x_lower[1] + patch_dx[1] * (lower_idx(1) - patch_lower(1) );
                    const double xv =
                    patch_x_lower[0] + patch_dx[0] * (lower_idx(0) - patch_lower(0) + 0.5);
                    
                    
					double u_ex, v_ex;
                    if (yu < (X1_lower + xu*tan(theta)) || yu > (X1_upper + xu*tan(theta) ))
                    {
                        u_ex = 0.0;
                  //      v_ex = 0.0;
                    }    
                    else
                    { 
                        u_ex = cos(theta)*(p_e*D/(MU*(L/cos(theta) + D*tan(theta))))*(-xu*sin(theta) + (yu-X1_lower)*cos(theta))*(1 - (-xu*sin(theta) + (yu-X1_lower)*cos(theta))/D);
                     }
                     
                     if (yv < (X1_lower + xv*tan(theta)) || yv > (X1_upper + xv*tan(theta) ))
                    {
						v_ex = 0.0;
                    }    
                    else
                    { 
						v_ex = sin(theta)*(p_e*D/(MU*(L/cos(theta) + D*tan(theta))))*(-xv*sin(theta) + (yv-X1_lower)*cos(theta))*(1 - (-xv*sin(theta) + (yv-X1_lower)*cos(theta))/D);
                     }
                     
                    //~ const double x1 = x0 + patch_dx[0];
                    const double u0 =
                    (*u_data)(SideIndex<NDIM>(lower_idx, 0, SideIndex<NDIM>::Lower));
                     const double v0 =
                    (*u_data)(SideIndex<NDIM>(lower_idx, 1, SideIndex<NDIM>::Lower));

                    //~ const double u1 =
                    //~ (*u_data)(SideIndex<NDIM>(upper_idx, 0, SideIndex<NDIM>::Lower));
                    //~ pos_values.push_back(y);
                    //~ pos_values.push_back(u0 + (u1 - u0) * (x_loc - x0) / (x1 - x0));

                   if (xu>0.4*L && xv>0.4*L && xu<0.6*L && xv<0.6*L)
                  {
					N_max +=1;                   
                    u_Eulerian_L2_norm += std::abs(u0 - u_ex) * std::abs(u0 - u_ex)*(*wgt_cc_data)(lower_idx);
                    u_Eulerian_L2_norm += std::abs(v0 - v_ex) * std::abs(v0 - v_ex)*(*wgt_cc_data)(lower_idx);

                    u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(u0 - u_ex));
                    u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(v0 - v_ex));
                  }  
                }
            }
        }
    }  
    
         SAMRAI_MPI::sumReduction(&N_max, 1);
         SAMRAI_MPI::sumReduction(&u_Eulerian_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&u_Eulerian_max_norm, 1);

        u_Eulerian_L2_norm = sqrt(u_Eulerian_L2_norm/static_cast<Real>(N_max));
		//u_Eulerian_L2_norm = sqrt(u_Eulerian_L2_norm);
         
         pout << " u_Eulerian_L2_norm = " << u_Eulerian_L2_norm <<"\n\n";
         pout << " u_Eulerian_max_norm = " << u_Eulerian_max_norm <<"\n\n";

    
    return;
} // compute_velocity_profile




void compute_pressure_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int p_idx,
                              const double data_time,
                              const string& data_dump_dirname)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    
            HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.resetLevels(finest_ln, finest_ln);
            const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

    const double X_min[2] = { 0.1*L , 0.1*L };
    const double X_max[2] = { 0.9*L  ,0.9*L };
    double p_Eulerian_L2_norm = 0.0;
    double p_Eulerian_max_norm = 0.0;
    int N_max = 0;
    for (int ln = finest_ln; ln >= coarsest_ln; --ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const CellIndex<NDIM>& patch_lower = patch_box.lower();
            const CellIndex<NDIM>& patch_upper = patch_box.upper();

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom =
            patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();
            const double* const patch_dx = patch_geom->getDx();
            

            
            // Entire box containing the required data.
            Box<NDIM> box(IndexUtilities::getCellIndex(&X_min[0], patch_x_lower, patch_x_upper,
                                                       patch_dx, patch_lower, patch_upper),
                          IndexUtilities::getCellIndex(&X_max[0], patch_x_lower, patch_x_upper,
                                                       patch_dx, patch_lower, patch_upper));
            // Part of the box on this patch
            Box<NDIM> trim_box = patch_box * box;
            BoxList<NDIM> iterate_box_list = trim_box;
            
            // Trim the box covered by the finer region
            BoxList<NDIM> covered_boxes;
            if (ln < finest_ln)
            {
                BoxArray<NDIM> refined_region_boxes;
                Pointer<PatchLevel<NDIM> > next_finer_level =
                patch_hierarchy->getPatchLevel(ln + 1);
                refined_region_boxes = next_finer_level->getBoxes();
                refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
                for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                {
                    const Box<NDIM> refined_box = refined_region_boxes[i];
                    const Box<NDIM> covered_box = trim_box * refined_box;
                    covered_boxes.unionBoxes(covered_box);
                }
            }
            iterate_box_list.removeIntersections(covered_boxes);
            
            // Loop over the boxes and store the location and interpolated value.
            //~ Pointer<Data<NDIM, double> > p_data = patch->getPatchData(p_idx);
            const Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_idx);
            const Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);
            for (BoxList<NDIM>::Iterator lit(iterate_box_list); lit; lit++)
            {
                const Box<NDIM>& iterate_box = *lit;
                for (Box<NDIM>::Iterator bit(iterate_box); bit; bit++)
                {
                    const CellIndex<NDIM>& cell_idx = *bit;
                    
                    const double y =
                    patch_x_lower[1] + patch_dx[1] * (cell_idx(1) - patch_lower(1) + 0.5);
                    const double x =
                    patch_x_lower[0] + patch_dx[0] * (cell_idx(0) - patch_lower(0) + 0.5);
                    //~ double p_ex_qp = -2.*p_e*x/L + p_e;
                    
                    const double p1 = (*p_data)(cell_idx);
					double p_ex_qp;
                    if (y < (X1_lower + x*tan(theta) - fac*patch_dx[1]) || y > (X1_upper + x*tan(theta) + fac*patch_dx[1] ))
                    { 
                        p_ex_qp = 0.0;
					} 
                    else if (std::abs(y - (X1_lower + x*tan(theta)))<= fac*patch_dx[1] || std::abs(y - (X1_upper + x*tan(theta)))<= fac*patch_dx[1])
                    {
						p_ex_qp = p1;
					}
                    else
                    {
                        p_ex_qp = (-2.*p_e*(x/cos(theta) + (y-X1_lower-x*tan(theta))*sin(theta))/(L/cos(theta) + D*tan(theta)) + p_e);
                        
                    }

                    
                     if (x>0.4*L && x<0.6*L)
                     {
						 N_max +=1;
						p_Eulerian_L2_norm += std::abs(p1 - p_ex_qp) * std::abs(p1 - p_ex_qp) * (*wgt_cc_data)(cell_idx);
						p_Eulerian_max_norm = std::max(p_Eulerian_max_norm, std::abs(p1 - p_ex_qp));
					}
                }
            }
        }
    }
		SAMRAI_MPI::sumReduction(&N_max, 1);
        SAMRAI_MPI::sumReduction(&p_Eulerian_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&p_Eulerian_max_norm, 1);

        p_Eulerian_L2_norm = sqrt(p_Eulerian_L2_norm/static_cast<Real>(N_max));
		//p_Eulerian_L2_norm = sqrt(p_Eulerian_L2_norm);
         
         pout << " p_Eulerian_L2_norm = " << p_Eulerian_L2_norm <<"\n\n";
         pout << " p_Eulerian_max_norm = " << p_Eulerian_max_norm <<"\n\n";


    return;
} // compute_pressure_profile






void
postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
                 const int /*iteration_num*/,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{

    const unsigned int dim = mesh.mesh_dimension();
    double F_integral[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;

    System& x_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system(IIMethod::VELOCITY_SYSTEM_NAME);
    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>& X0_vec = x_system.get_vector("INITIAL_COORDINATES");
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    NumericVector<double>* U_vec = U_system.solution.get();
    NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
    U_vec->localize(*U_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    UniquePtr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    UniquePtr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<libMesh::Point>& q_point = fe->get_xyz();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();
    const std::vector<std::vector<double> >& dphi_dxi = fe->get_dphidxi();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF_qp;
    boost::multi_array<double, 2> x_node, X_node, U_node, P_o_node, P_i_node, P_j_node;
    VectorValue<double> F_qp, U_qp, x_qp, X_qp, N,n;

    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            jacobian(FF_qp, qp, x_node, dphi);
            interpolate(U_qp, qp, U_node, phi);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U_qp(d);
            }
            tether_body_force_function_thin(F_qp, n, N, FF_qp, x_qp, q_point[qp], elem, 0, var_data, grad_var_data, loop_time, force_fcn_ctx);

                
                
            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F_qp(d) * JxW[qp];
            }
        }
  
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    static const double D = 1.0;

    

    {
 
        double WSS_L2_norm = 0.0, WSS_max_norm = 0.0;
        double U_L2_norm = 0.0, U_max_norm = 0.0;
        double P_L2_norm = 0.0, P_max_norm = 0.0;
        double disp_L2_norm = 0.0, disp_max_norm = 0.0;
        System& U_system = equation_systems->get_system<System>(IIMethod::VELOCITY_SYSTEM_NAME);
		System& WSS_system = equation_systems->get_system<System>(IIMethod::WSS_IN_SYSTEM_NAME);
		System& P_o_system = equation_systems->get_system<System>(IIMethod::PRESSURE_OUT_SYSTEM_NAME);
		System& P_i_system = equation_systems->get_system<System>(IIMethod::PRESSURE_IN_SYSTEM_NAME);
		System& P_j_system = equation_systems->get_system<System>(IIMethod::PRESSURE_JUMP_SYSTEM_NAME);
		
        NumericVector<double>* U_vec = U_system.solution.get();
        NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
        U_vec->localize(*U_ghost_vec);
        DofMap& U_dof_map = U_system.get_dof_map();
        std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);
        

        NumericVector<double>* WSS_vec = WSS_system.solution.get();
        NumericVector<double>* WSS_ghost_vec = WSS_system.current_local_solution.get();
        WSS_vec->localize(*WSS_ghost_vec);
        DofMap& WSS_dof_map = WSS_system.get_dof_map();
        std::vector<std::vector<unsigned int> > WSS_dof_indices(NDIM);
        UniquePtr<FEBase> fe(FEBase::build(dim, WSS_dof_map.variable_type(0)));
        
        
        
        NumericVector<double>* P_o_vec = P_o_system.solution.get();
        NumericVector<double>* P_o_ghost_vec = P_o_system.current_local_solution.get();
        P_o_vec->localize(*P_o_ghost_vec);
        DofMap& P_o_dof_map = P_o_system.get_dof_map();
         std::vector<unsigned int> P_o_dof_indices;
         
        NumericVector<double>* P_i_vec = P_i_system.solution.get();
        NumericVector<double>* P_i_ghost_vec = P_i_system.current_local_solution.get();
        P_i_vec->localize(*P_i_ghost_vec);
        DofMap& P_i_dof_map = P_i_system.get_dof_map();
         std::vector<unsigned int> P_i_dof_indices;
		
        NumericVector<double>* P_j_vec = P_j_system.solution.get();
        NumericVector<double>* P_j_ghost_vec = P_j_system.current_local_solution.get();
        P_j_vec->localize(*P_j_ghost_vec);
        DofMap& P_j_dof_map = P_j_system.get_dof_map();
        std::vector<unsigned int> P_j_dof_indices;
       
		
		
		
		VectorValue<double> U_qp, WSS_qp;
        double  P_o_qp,P_j_qp,P_i_qp;
        VectorValue<double> tau1, tau2;
        int qp_tot = 0;
        boost::multi_array<double, 2> U_node, WSS_node;
         boost::multi_array<double, 1> P_o_node, P_i_node, P_j_node;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
				dof_map.dof_indices(elem, dof_indices[d], d);
                U_dof_map.dof_indices(elem, U_dof_indices[d], d);
                WSS_dof_map.dof_indices(elem, WSS_dof_indices[d], d);
                
            }
            P_j_dof_map.dof_indices(elem, P_j_dof_indices);
            P_o_dof_map.dof_indices(elem, P_o_dof_indices);
            P_i_dof_map.dof_indices(elem, P_i_dof_indices);
            const int n_qp = qrule->n_points();
            get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
            get_values_for_interpolation(WSS_node, *WSS_ghost_vec, WSS_dof_indices);
            get_values_for_interpolation(P_j_node, *P_j_ghost_vec, P_j_dof_indices);
            get_values_for_interpolation(P_i_node, *P_i_ghost_vec, P_i_dof_indices);
            get_values_for_interpolation(P_o_node, *P_o_ghost_vec, P_o_dof_indices);
			get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
			get_values_for_interpolation(X_node, X0_vec, dof_indices);

			
            for (int qp = 0; qp < n_qp; ++qp)
            {

				interpolate(x_qp, qp, x_node, phi);
				interpolate(X_qp, qp, X_node, phi);
                interpolate(U_qp, qp, U_node, phi);
                interpolate(WSS_qp, qp, WSS_node, phi);
                interpolate(P_o_qp, qp, P_o_node, phi);
                interpolate(P_j_qp, qp, P_j_node, phi);
                interpolate(P_i_qp, qp, P_i_node, phi);
                P_i_qp = -(P_j_qp - P_o_qp);
                
                double p_ex_qp =(-2.0*p_e*(X_qp(0)/cos(theta) + (X_qp(1)-X1_lower-X_qp(0)*tan(theta))*sin(theta))/(L/cos(theta) + D*tan(theta)) + p_e);

                interpolate(&tau1(0), qp, X_node, dphi_dxi);
   
                tau2 = VectorValue<double>(0.0, 0.0, 1.0);

                N = tau1.cross(tau2);
                N = N.unit();
                
                double ex_wss[NDIM];
                double ex_U[NDIM];
				ex_U[1] = 0.0; 
				ex_U[0] = 0.0;
				ex_wss[0] =  -p_e*D* cos(theta)/((L/cos(theta) + D*tan(theta)));
				ex_wss[1] = -p_e*D* sin(theta)/(L/cos(theta) + D*tan(theta)); 
				

				if (x_qp(0)>0.2*L && x_qp(0)<0.8*L)
				{
					qp_tot +=1;
					for (unsigned int d = 0; d < NDIM; ++d)
					{ 
						
						WSS_L2_norm += (WSS_qp(d) - ex_wss[d])* (WSS_qp(d) - ex_wss[d]) * JxW[qp];
						WSS_max_norm = std::max(WSS_max_norm, std::abs(WSS_qp(d) - ex_wss[d]));
					   U_L2_norm  +=  (U_qp(d)- ex_U[d])* (U_qp(d)- ex_U[d]) * JxW[qp];
					   U_max_norm = std::max(U_max_norm, std::abs(U_qp(d)- ex_U[d]));
					   disp_L2_norm += std::abs(X_qp(d) - x_qp(d)) * std::abs(X_qp(d) - x_qp(d)) * JxW[qp];
					   disp_max_norm = std::max(disp_max_norm, std::abs(X_qp(d) - x_qp(d)));
					
					}
					P_L2_norm += std::abs(P_i_qp - p_ex_qp) * std::abs(P_i_qp - p_ex_qp) * JxW[qp];
					P_max_norm = std::max(P_max_norm, std::abs(P_i_qp - p_ex_qp));
                }
				
			}
		}

		SAMRAI_MPI::sumReduction(&qp_tot, 1);
		
        SAMRAI_MPI::sumReduction(&WSS_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&WSS_max_norm, 1);
        
        SAMRAI_MPI::sumReduction(&U_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&U_max_norm, 1);
        
        SAMRAI_MPI::sumReduction(&P_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&P_max_norm, 1);
        
        SAMRAI_MPI::sumReduction(&disp_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&disp_max_norm, 1);
        
        
        // U_L2_norm = sqrt(U_L2_norm/static_cast<Real>(qp_tot));
        // WSS_L2_norm = sqrt(WSS_L2_norm/static_cast<Real>(qp_tot));
        // P_L2_norm = sqrt(P_L2_norm/static_cast<Real>(qp_tot));
        // disp_L2_norm = sqrt(disp_L2_norm/static_cast<Real>(qp_tot));
        
        U_L2_norm = sqrt(U_L2_norm);
        WSS_L2_norm = sqrt(WSS_L2_norm);
        disp_L2_norm = sqrt(disp_L2_norm);
        P_L2_norm = sqrt(P_L2_norm);
        
         pout << " WSS_L2_norm = " << WSS_L2_norm <<"\n\n";
         pout << " WSS_max_norm = " << WSS_max_norm <<"\n\n";
         
         
         pout << " U_L2_norm = " << U_L2_norm <<"\n\n";
         pout << " U_max_norm = " << U_max_norm <<"\n\n";
         
         
         pout << " P_L2_norm = " << P_L2_norm <<"\n\n";
         pout << " P_max_norm = " << P_max_norm <<"\n\n";
         
         pout << " disp_L2_norm = " << disp_L2_norm <<"\n\n";
         pout << " disp_max_norm = " << disp_max_norm <<"\n\n";
         
         
		if (SAMRAI_MPI::getRank() == 0)
		{
			p_norm_stream.precision(12);
			p_norm_stream.setf(ios::fixed, ios::floatfield);
			p_norm_stream << loop_time << "\t" << P_L2_norm << "\t" << P_max_norm << endl;			
		}
         
    }


    return;
} // postprocess_data

