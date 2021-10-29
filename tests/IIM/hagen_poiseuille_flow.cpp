//  2018-2021, created by  Amin Kolahdouz
// 3D Hagen–Poiseuille flow in a pipe
// See Sec. 4.1.2 of "An immersed interface method for discrete surfaces" 
//  by Ebrahim M. Kolahdouz et al., Journal of Computational Physics 400 (2020) 108854
// To get the right convergence results the code was run with a tight solver/preconditioner tolerance
// with the command line flag: -stokes_ksp_rtol 1.0e-10 -ksp_rtol 1.0e-10

// This code was tested against the commit 57fb379454ea3f8f50c476f10226cb6b520a11a0
// which is currently on branch iim-1 at https://github.com/drwells/IBAMR

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Application-specific includes.

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>
#include <hagen_poiseuille_FeedbackForcer.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/face_quad4.h>
#include <libmesh/face_quad.h>
// Headers for application-specific algorithm/data structure objects

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/SpongeLayerForceFunction.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibamr/StaggeredStokesOpenBoundaryStabilizer.h>


namespace ModelData
{
// Elasticity model data for thin body.
// Tether (penalty) force function for thin body.
static double kappa_s_thin = 1.0;
static double eta_s_thin = 1.0;

static double L = 0.0;
static double D = 0.0;
static double H = 0.0;
static double p_e =  0.0;
static double MU = 0.0;
static double theta = 0.0;
static double zi = 0.0;
static double yi = 0.0;
static double yo = 0.0;
static double U_MAX = 0.0;
static double zo = 0.0;
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
    double x_kappa_n = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
		x_kappa_n += n(d) * kappa_s_thin * (X(d) - x(d));
		u_bndry_n += n(d) * U[d];
	}

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
       F(d) = kappa_s_thin * (X(d) - x(d)) - eta_s_thin * (U[d] - u_bndry_n * n(d));
    }
    return; 
}

inline
unsigned int idx(const unsigned int nr,
                 const unsigned int i,
                 const unsigned int j)
{

    return i + j*nr;

  return libMesh::invalid_uint;
}

static double x_loc_min, x_loc_max, z_loc_min, z_loc_max, y_loc_min, y_loc_max;
}
using namespace ModelData;

// Function prototypes
void compute_velocity_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int u_idx,
                              const double data_time,
                              const string& data_dump_dirname);
                              
    void compute_pressure_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                   const int p_idx,
                   const double data_time,
					const string& data_dump_dirname);
                              
                              
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);

void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<IBHierarchyIntegrator> ins_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);
void compute_velocity_profile(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                              const int u_idx,
                              const double data_time,
                              const string& data_dump_dirname);
                              
void compute_flow_rate(const double dt,
                       const Pointer<PatchHierarchy<NDIM> > hierarchy,
                       const int U_idx,
                       const double loop_time,
                       ostream& flow_rate_stream,
                       const int wgt_sc_idx);
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
int main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    
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
        const string thin_exodus_filename = app_initializer->getExodusIIFilename("cylinder_thin");

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

        // Create FE mesh.
        D = input_db->getDouble("D");         // channel height (cm)
        L = input_db->getDouble("L");         // channel length (cm)
        H = input_db->getDouble("H");         // channel length (cm)
        const double H = input_db->getDouble("H");         // wall thickness (cm)
        

        zo = input_db->getDouble("Zo");
        yo = input_db->getDouble("Yo");
        yi = input_db->getDouble("Yi"); 
        zi = input_db->getDouble("Zi"); 
        U_MAX = input_db->getDouble("U_MAX"); 
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        MU = input_db->getDouble("MU");
        theta = input_db->getDouble("THETA");
        p_e = input_db->getDouble("P_E");


		Mesh cylinder_mesh_thin(init.comm(), NDIM-1);
          
        cylinder_mesh_thin.boundary_info->clear_boundary_node_ids();
        const unsigned int  NXi_elem = ceil(L/ds);
        const unsigned int NRi_elem = ceil(M_PI*D/ds);
        int node_id = 0;
        cylinder_mesh_thin.reserve_nodes (NRi_elem*(NXi_elem + 1));
        cylinder_mesh_thin.reserve_elem (NRi_elem*NXi_elem);
    
		for (unsigned int j = 0; j <= NXi_elem; j++)
		{              
			for (unsigned int i = 0; i <= NRi_elem -1; i++)
			{ 
	   
				   const double theta = 2.0 * M_PI * static_cast<Real>(i) / static_cast<Real>(NRi_elem);
				   cylinder_mesh_thin.add_point(libMesh::Point(L*static_cast<Real>(j)/static_cast<Real>(NXi_elem),0.5*H + 0.5*D*cos(theta), 0.5*H + 0.5*D*sin(theta)), node_id++);
			}
		}
        
      
        
        for (unsigned int j = 0; j <= NXi_elem-1; j++)
			{  
        for (unsigned int i = 0; i <= NRi_elem -2 ; i++)
        { 
   
                     Elem * elem = cylinder_mesh_thin.add_elem (new Quad4);
					 elem->set_node(0) = cylinder_mesh_thin.node_ptr(idx(NRi_elem,i,j));
                     elem->set_node(1) = cylinder_mesh_thin.node_ptr(idx(NRi_elem,i+1,j));
                     elem->set_node(2) = cylinder_mesh_thin.node_ptr(idx(NRi_elem,i+1,j+1));
                     elem->set_node(3) = cylinder_mesh_thin.node_ptr(idx(NRi_elem,i,j+1));
           }
	    }
	      
	    
	    for (unsigned int j = 0; j <= NXi_elem-1; j++)
		{ 	
			  Elem * elem = cylinder_mesh_thin.add_elem (new Quad4);
		      elem->set_node(0) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, NRi_elem-1,j));
              elem->set_node(1) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, 0, j));
              elem->set_node(2) = cylinder_mesh_thin.node_ptr(idx(NRi_elem,0, j+1));
              elem->set_node(3) = cylinder_mesh_thin.node_ptr(idx(NRi_elem, NRi_elem-1,j+1));		
		}
		

        MeshBase::const_element_iterator el_end = cylinder_mesh_thin.elements_end();
        for (MeshBase::const_element_iterator el = cylinder_mesh_thin.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    BoundaryInfo* boundary_info_cylinder = cylinder_mesh_thin.boundary_info.get();
                    boundary_info_cylinder->add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID);
                }
            }
        }
        cylinder_mesh_thin.prepare_for_use();
          

        
        vector<MeshBase*> meshes(1);
        meshes[0] = &cylinder_mesh_thin;

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
        
        Pointer<IIMethod> ib_method_ops =
        new IIMethod("IIMethod", app_initializer->getComponentDatabase("IIMethod"),
                       meshes, app_initializer->getComponentDatabase("GriddingAlgorithm")
                       ->getInteger("max_levels"));
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
                                                               
                                                              
        ib_method_ops->registerLagSurfaceForceFunction(tether_body_force_thin_data, 0);

       

        ib_method_ops->initializeFEEquationSystems();
        // Add target points
     

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

        }
        // Create Eulerian body force function specification objects.

        time_integrator->registerBodyForceFunction(
            new hagen_poiseuille_FeedbackForcer(H, D, navier_stokes_integrator, patch_hierarchy));

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer =
        app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        UniquePtr<ExodusII_IO> exodus_io_thin(uses_exodus ? new ExodusII_IO(cylinder_mesh_thin) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Set up locations to get velocity data.
        x_loc_min = input_db->getDoubleWithDefault("X_LOC_MIN", 0.1*L);
        x_loc_max = input_db->getDoubleWithDefault("X_LOC_MAX", 0.9*L);
        z_loc_min = input_db->getDoubleWithDefault("Z_LOC_MIN", 0.0);
        y_loc_min = input_db->getDoubleWithDefault("Y_LOC_MIN", 0.0 );
        y_loc_max = input_db->getDoubleWithDefault("Y_LOC_MAX", L );
        z_loc_max = input_db->getDoubleWithDefault("Z_LOC_MAX", L );

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);
        
        // Write out initial visualization data.
        EquationSystems* equation_systems_thin =
            ib_method_ops->getFEDataManager(0)->getEquationSystems();
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                exodus_io_thin->write_timestep(thin_exodus_filename, *equation_systems_thin,
                                          iteration_num / viz_dump_interval + 1, loop_time);
            }
        }
        
		ofstream flow_rate_stream;
        if (SAMRAI_MPI::getRank() == 0)
        {
            flow_rate_stream.open("flow_rate.curve");
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
            Pointer<hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
            Pointer<VariableContext> current_ctx = time_integrator->getCurrentContext();

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
                                                                         exodus_io_thin->write_timestep(thin_exodus_filename, *equation_systems_thin,
                                              iteration_num / viz_dump_interval + 1,
                                              loop_time);
                                              
					exodus_io_thin->write_timestep(thin_exodus_filename, *equation_systems_thin,
                                              iteration_num / viz_dump_interval + 1,
                                              loop_time);
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
            
            
            {
                 VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
				Pointer<hier::Variable<NDIM> > u_var = time_integrator->getVelocityVariable();
				const Pointer<VariableContext> u_ctx = time_integrator->getCurrentContext();
				const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
				const int coarsest_ln = 0;
				const int finest_ln = patch_hierarchy->getFinestLevelNumber();
				HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
				hier_math_ops.setPatchHierarchy(patch_hierarchy);
				hier_math_ops.resetLevels(coarsest_ln, finest_ln);
				const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
				 
				compute_flow_rate(dt, patch_hierarchy, u_idx, loop_time, flow_rate_stream, wgt_sc_idx);  
			}

			pout << "\nWriting state data...\n\n";
			 //  output_data(patch_hierarchy, navier_stokes_integrator, iteration_num, loop_time, postproc_data_dump_dirname);

			 postprocess_data(patch_hierarchy,
						 navier_stokes_integrator,
						  cylinder_mesh_thin,
						  equation_systems_thin,
						  iteration_num,
						 loop_time,
						  postproc_data_dump_dirname);

        }

     // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
		 Pointer<hier::Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const Pointer<VariableContext> u_ctx = navier_stokes_integrator->getCurrentContext();
        const int u_idx = var_db->mapVariableAndContextToIndex(u_var, u_ctx);
        const int u_cloned_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();

      const Pointer<hier::Variable<NDIM> > p_var = time_integrator->getPressureVariable();
        const Pointer<VariableContext> p_ctx = time_integrator->getCurrentContext();

        const int p_idx = var_db->mapVariableAndContextToIndex(p_var, p_ctx);
        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
        
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
        
        
         HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
		const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
            pout << std::setprecision(16) << "Error in u at time " << loop_time << ":\n"
                 << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
                 << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n"
                  << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
             
             
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(p_idx, p_idx, p_cloned_idx);
        pout << "Error in p at time " << loop_time - 0.5 * dt << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(p_idx, wgt_cc_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(p_idx, wgt_cc_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(p_idx, wgt_cc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
                 
                 
             pout<< " MU = "<< MU <<"\n"
                  << "  dx:  " << dx << "\n"
					<< "  dt: " << dt << "\n";

             if  (input_db->getBool("USE_VELOCITY_JUMP_CONDITIONS"))
					pout<<" Using the jump condition"<<"\n";
			else
				    pout<<" Using regular IB"<<"\n";
        if (dump_viz_data && uses_visit)
        {
            navier_stokes_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
            
             exodus_io_thin->write_timestep(thin_exodus_filename, *equation_systems_thin,
                                                iteration_num / viz_dump_interval + 1,
                                                loop_time);
        }


       if (SAMRAI_MPI::getRank() == 0)
        {
            flow_rate_stream.close();
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
    const double X_min[3] = { 0.1*L , y_loc_min, z_loc_min };
    const double X_max[3] = { 0.9*L  , y_loc_max, y_loc_max };
    int qp_tot=0;
    double u_Eulerian_L2_norm = 0.0;
    double u_Eulerian_max_norm = 0.0;
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

                    const double yu =
                    patch_x_lower[1] + patch_dx[1] * (lower_idx(1) - patch_lower(1) + 0.5);
                    const double zu =
                    patch_x_lower[2] + patch_dx[2] * (lower_idx(2) - patch_lower(2) + 0.5);
                    const double xu =
                    patch_x_lower[0] + patch_dx[0] * (lower_idx(0) - patch_lower(0));


					double u_ex, v_ex, w_ex;
					
                    if (sqrt((yu - 0.5*H)*(yu - 0.5*H) + (zu - 0.5*H)*(zu - 0.5*H)) >= D/2)
                    {
                        u_ex = 0.0;
                        v_ex = 0.0;
                        w_ex = 0.0;
                    }    
                    else
                    { 
                        u_ex =U_MAX*(1.0 - 4.0*((yu - 0.5*H)*(yu - 0.5*H) + (zu - 0.5*H)*(zu - 0.5*H))/(D*D));
                        v_ex = 0.0;
                        w_ex = 0.0;
                     }
                     
                    const double u0 =
                    (*u_data)(SideIndex<NDIM>(lower_idx, 0, SideIndex<NDIM>::Lower));
                     const double v0 =
                    (*u_data)(SideIndex<NDIM>(lower_idx, 1, SideIndex<NDIM>::Lower));
                     const double w0 =
                    (*u_data)(SideIndex<NDIM>(lower_idx, 2, SideIndex<NDIM>::Lower));

                    if (xu > 0.4*L && xu < 0.6*L) 
                    {
						qp_tot +=1;
						u_Eulerian_L2_norm += std::abs(u0 - u_ex) * std::abs(u0 - u_ex)*(*wgt_cc_data)(lower_idx);
						u_Eulerian_L2_norm += std::abs(v0 - v_ex) * std::abs(v0 - v_ex)*(*wgt_cc_data)(lower_idx);
						u_Eulerian_L2_norm += std::abs(w0 - w_ex) * std::abs(w0 - w_ex)*(*wgt_cc_data)(lower_idx);

						u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(u0 - u_ex));
						u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(v0 - v_ex));
						u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(w0 - w_ex));
					}
                    
                }
            }
        }
    }  
    
		SAMRAI_MPI::sumReduction(&qp_tot, 1);
        SAMRAI_MPI::sumReduction(&u_Eulerian_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&u_Eulerian_max_norm, 1);

		u_Eulerian_L2_norm = sqrt(u_Eulerian_L2_norm);
         
         pout << " u_Eulerian_L2_norm = " << u_Eulerian_L2_norm <<"\n\n";
         pout << " u_Eulerian_max_norm = " << u_Eulerian_max_norm <<"\n\n";

    
    return;
} // compute_velocity_profile

void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(
        var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                             navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(
        var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                             navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data



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
            boost::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
        dphi_dxi[0] = &fe->get_dphidxi();
        if (NDIM > 2) dphi_dxi[1] = &fe->get_dphideta();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF_qp;
    boost::multi_array<double, 2> x_node, X_node, U_node;
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
	     int qp_tot = 0;
         double WSS_L2_norm = 0.0, WSS_max_norm = 0.0;
        double U_L2_norm = 0.0, U_max_norm = 0.0;
        double disp_L2_norm = 0.0, disp_max_norm = 0.0;
        double P_L2_norm = 0.0, P_max_norm = 0.0;
        System& U_system = equation_systems->get_system<System>(IIMethod::VELOCITY_SYSTEM_NAME);
		System& WSS_system = equation_systems->get_system<System>(IIMethod::WSS_IN_SYSTEM_NAME);
		System& P_o_system = equation_systems->get_system<System>(IIMethod::PRESSURE_OUT_SYSTEM_NAME);
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
        boost::array<VectorValue<double>, 2> dX_dxi, dx_dxi;
        
        
        NumericVector<double>* P_o_vec = P_o_system.solution.get();
        NumericVector<double>* P_o_ghost_vec = P_o_system.current_local_solution.get();
        P_o_vec->localize(*P_o_ghost_vec);
        DofMap& P_o_dof_map = P_o_system.get_dof_map();
         std::vector<unsigned int> P_o_dof_indices;
		
        NumericVector<double>* P_j_vec = P_j_system.solution.get();
        NumericVector<double>* P_j_ghost_vec = P_j_system.current_local_solution.get();
        P_j_vec->localize(*P_j_ghost_vec);
        DofMap& P_j_dof_map = P_j_system.get_dof_map();
        std::vector<unsigned int> P_j_dof_indices;

		
		
		
		
        VectorValue<double> U_qp, WSS_qp;
        double  P_o_qp,P_j_qp,P_i_qp;
        boost::multi_array<double, 2> U_node, WSS_node;
        boost::multi_array<double, 1> P_o_node, P_j_node;
        const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            // fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
				dof_map.dof_indices(elem, dof_indices[d], d);
                U_dof_map.dof_indices(elem, U_dof_indices[d], d);
                WSS_dof_map.dof_indices(elem, WSS_dof_indices[d], d);
                
            }
            P_j_dof_map.dof_indices(elem, P_j_dof_indices);
            P_o_dof_map.dof_indices(elem, P_o_dof_indices);
            const int n_qp = qrule->n_points();
            get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
            get_values_for_interpolation(WSS_node, *WSS_ghost_vec, WSS_dof_indices);
			get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
			get_values_for_interpolation(X_node, X0_vec, dof_indices);
			get_values_for_interpolation(P_j_node, *P_j_ghost_vec, P_j_dof_indices);
            get_values_for_interpolation(P_o_node, *P_o_ghost_vec, P_o_dof_indices);

			
            for (int qp = 0; qp < n_qp; ++qp)
            {

				interpolate(x_qp, qp, x_node, phi);
				interpolate(X_qp, qp, X_node, phi);
                interpolate(U_qp, qp, U_node, phi);
                interpolate(WSS_qp, qp, WSS_node, phi);
                interpolate(P_o_qp, qp, P_o_node, phi);
                interpolate(P_j_qp, qp, P_j_node, phi);
                P_i_qp = -(P_j_qp - P_o_qp);
                
                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi[k]);
                }
                if (NDIM == 2)
                {
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }
                n = (dx_dxi[0].cross(dx_dxi[1])).unit();
                
                double p_ex_qp = -2.*p_e*x_qp(0)/(L) + p_e;
                

                double ex_wss[NDIM];
                double ex_U[NDIM];
				ex_wss[0] = -0.5*p_e*D/L; //P_E*D/(MU*L) ; 
				ex_wss[1] = 0.0;
				ex_wss[2] = 0.0;
				ex_U[0] = 0.0;
				ex_U[1] = 0.0;
				ex_U[2] = 0.0;		
				if (x_qp(0)>0.4*L && x_qp(0)<0.6*L)
				{
					qp_tot += 1;
					for (unsigned int d = 0; d < NDIM; ++d)
					{ 
						U_L2_norm += (U_qp(d) - ex_U[d])*(U_qp(d) - ex_U[d])* JxW[qp];
						U_max_norm = std::max(U_max_norm, std::abs(U_qp(d) - ex_U[d]));
						WSS_L2_norm += (WSS_qp(d) - ex_wss[d])* (WSS_qp(d) - ex_wss[d]) * JxW[qp];
						WSS_max_norm = std::max(WSS_max_norm, std::abs(WSS_qp(d) - ex_wss[d]));
						disp_L2_norm += (X_qp(d) - x_qp(d))* (X_qp(d) - x_qp(d)) * JxW[qp];
						disp_max_norm = std::max(disp_max_norm, std::abs(X_qp(d) - x_qp(d)));
					}
					P_L2_norm += std::abs(P_i_qp - p_ex_qp) * std::abs(P_i_qp -p_ex_qp) * JxW[qp];
					P_max_norm = std::max(P_max_norm, std::abs(P_i_qp - p_ex_qp));
                }
				
			}
		}

		SAMRAI_MPI::sumReduction(&qp_tot, 1);
		
        SAMRAI_MPI::sumReduction(&WSS_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&WSS_max_norm, 1);
        SAMRAI_MPI::sumReduction(&U_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&U_max_norm, 1);
        SAMRAI_MPI::sumReduction(&disp_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&disp_max_norm, 1);
        SAMRAI_MPI::sumReduction(&P_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&P_max_norm, 1);

        U_L2_norm = sqrt(U_L2_norm);
        WSS_L2_norm = sqrt(WSS_L2_norm);
        disp_L2_norm = sqrt(disp_L2_norm);
        P_L2_norm = sqrt(P_L2_norm);
       
         pout << " WSS_L2_norm = " << WSS_L2_norm <<"\n\n";
         pout << " WSS_max_norm = " << WSS_max_norm <<"\n\n";
         
         
         pout << " U_L2_norm = " << U_L2_norm <<"\n\n";
         pout << " U_max_norm = " << U_max_norm <<"\n\n";
         
         
         pout << " disp_L2_norm = " << disp_L2_norm <<"\n\n";
         pout << " disp_max_norm = " << disp_max_norm <<"\n\n";
         
         pout << " P_L2_norm = " << P_L2_norm <<"\n\n";
         pout << " P_max_norm = " << P_max_norm <<"\n\n";
    }

    return;
} // postprocess_data



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

    const double X_min[3] = { x_loc_min , y_loc_min, z_loc_min };
    const double X_max[3] = { x_loc_max , y_loc_max, z_loc_max };
    //vector<double> pos_values;
    double p_Eulerian_L2_norm = 0.0;
    double p_Eulerian_max_norm = 0.0;

    int qp_tot = 0;
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
            const Pointer<CellData<NDIM, double> > p_data = patch->getPatchData(p_idx);
            const Pointer<CellData<NDIM, double> > wgt_cc_data = patch->getPatchData(wgt_cc_idx);
        
            for (BoxList<NDIM>::Iterator lit(iterate_box_list); lit; lit++)
            {
                const Box<NDIM>& iterate_box = *lit;
                for (Box<NDIM>::Iterator bit(iterate_box); bit; bit++)
                {
                    const CellIndex<NDIM>& cell_idx = *bit;
                    const double p1 = (*p_data)(cell_idx);
                    const double y =
                    patch_x_lower[1] + patch_dx[1] * (cell_idx(1) - patch_lower(1) + 0.5);
                    const double z =
                    patch_x_lower[2] + patch_dx[2] * (cell_idx(2) - patch_lower(2) + 0.5);
                    const double x =
                    patch_x_lower[0] + patch_dx[0] * (cell_idx(0) - patch_lower(0) + 0.5);
                    
                    double p_ex_qp;
                    if (sqrt((y - 0.5*H)*(y - 0.5*H) + (z - 0.5*H)*(z - 0.5*H)) > D/2 + patch_dx[1])
                    {
						p_ex_qp = 0.0;
                    }    
                    else if (sqrt((y - 0.5*H)*(y - 0.5*H) + (z - 0.5*H)*(z - 0.5*H)) < D/2 - patch_dx[1])
                    { 
						p_ex_qp = -2.0*p_e*x /L + p_e;
						
					}
					else
					{
						p_ex_qp = p1;
					}

                    
                   
                    if (x > 0.4*L && x < 0.6*L) 
					{
						qp_tot +=1;
						p_Eulerian_L2_norm += std::abs(p1 - p_ex_qp) * std::abs(p1 - p_ex_qp) * (*wgt_cc_data)(cell_idx);
						p_Eulerian_max_norm = std::max(p_Eulerian_max_norm, std::abs(p1 - p_ex_qp));
					}
                }
            }
                    
        }
    }
		SAMRAI_MPI::sumReduction(&qp_tot, 1);
        SAMRAI_MPI::sumReduction(&p_Eulerian_L2_norm, 1);
        SAMRAI_MPI::maxReduction(&p_Eulerian_max_norm, 1);
 
        p_Eulerian_L2_norm = sqrt(p_Eulerian_L2_norm);
        
         pout << " p_Eulerian_L2_norm = " << p_Eulerian_L2_norm <<"\n\n";
         pout << " p_Eulerian_max_norm = " << p_Eulerian_max_norm <<"\n\n";


    return;
} // compute_pressure_profile


void compute_flow_rate(const double dt,
                       const Pointer<PatchHierarchy<NDIM> > hierarchy,
                       const int U_idx,
                       const double loop_time,
                       ostream& flow_rate_stream,
                       const int wgt_sc_idx)
{
	vector<double> qsrc;
	qsrc.resize(2);
	std::fill(qsrc.begin(), qsrc.end(), 0.0);
    const double posni[3] = { 0.0 , yi , zi };
    const double posno[3] = { L , yo , zo };
    const double rsrc[2] = {0.5*D, 0.5*D};
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            if (pgeom->getTouchesRegularBoundary())
            {
                Pointer<SideData<NDIM, double> > U_data = patch->getPatchData(U_idx);
                Pointer<SideData<NDIM, double> > wgt_sc_data = patch->getPatchData(wgt_sc_idx);
                const Box<NDIM>& patch_box = patch->getBox();
                const double* const x_lower = pgeom->getXLower();
                const double* const dx = pgeom->getDx();
                double dV = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    dV *= dx[d];
                }
                double X[NDIM];
                static const int axis = 0;
                for (int side = 0; side <= 1; ++side)
                {
                    const bool is_lower = side == 0;
                    if (pgeom->getTouchesRegularBoundary(axis, side))
                    {
                       // const double rsrc = d_rsrc[side];
						//
                        Vector n;
                        for (int d = 0; d < NDIM; ++d)
                        {
                            n[d] = axis == d ? (is_lower ? -1.0 : +1.0) : 0.0;
                        }
                        Box<NDIM> side_box = patch_box;
                        if (is_lower)
                        {
                            side_box.lower(axis) = patch_box.lower(axis);
                            side_box.upper(axis) = patch_box.lower(axis);
                        }
                        else
                        {
                            side_box.lower(axis) = patch_box.upper(axis) + 1;
                            side_box.upper(axis) = patch_box.upper(axis) + 1;
                        }
                        for (Box<NDIM>::Iterator b(side_box); b; b++)
                        {
                            const hier::Index<NDIM>& i = b();
                            double r_sq = 0.0;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                X[d] =
                                    x_lower[d] + dx[d] * (double(i(d) - patch_box.lower(d)) + (d == axis ? 0.0 : 0.5));
                                if (d != axis && side ==0) r_sq += pow(X[d] - posni[d], 2.0);
                                if (d != axis && side ==1) r_sq += pow(X[d] - posno[d], 2.0);
                            }
                            const double r = sqrt(r_sq);
                            if (r <= rsrc[side])
                            {
                                const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                if ((*wgt_sc_data)(i_s) > std::numeric_limits<double>::epsilon())
                                {
                                    double dA = n[axis] * dV / dx[axis];
                                    qsrc[side] += (*U_data)(i_s)*dA;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(&qsrc[0], 2);

    if (SAMRAI_MPI::getRank() == 0)
    {
        flow_rate_stream.precision(12);
        flow_rate_stream.setf(ios::fixed, ios::floatfield);
        flow_rate_stream << loop_time << " " << qsrc[0] <<" " << qsrc[1] << endl;
   
    }

    
}


