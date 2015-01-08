// Filename: main.cpp
// Created on 6 Jan 2015 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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
#include <libmesh/centroid_partitioner.h>
#include <libmesh/dense_matrix.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/fe.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/parallel.h>
#include <libmesh/quadrature.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBStrategySet.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibtk/muParserCartGridFunction.h>


namespace
{
	
static double R;
static double gamma = 0.0;
// Coordinate mapping function.
void
coordinate_mapping_function(libMesh::Point& X,
							const libMesh::Point& s,
							void* /*ctx*/)
{
	X(0) = (R + s(1)) * cos(s(0) / R) + 0.0;
	X(1) = (R + gamma + s(1)) * sin(s(0) / R) + 0.0;
	return;
} // coordinate_mapping_function
	
struct RigidBodyCtx
{
	VectorValue<double> V, W;     // rigid translational and rotational velocity
	double dt;                    // present time step size
};

// Set the rigid body velocity to zero.
void
update_rigid_body_velocity(libMesh::NumericVector<double>& U_b,
						   libMesh::NumericVector<double>& /*U*/,
						   libMesh::NumericVector<double>& /*X*/,
						   const Eigen::Vector3d& /*X_com*/,
						   Eigen::Vector3d& U_com,
						   Eigen::Vector3d& W_com,
						   libMesh::EquationSystems* /*equation_systems*/,
						   double /*data_time*/,
						   void* ctx)
{
	// Get housing data.
	RigidBodyCtx* rb_ctx = static_cast<RigidBodyCtx*>(ctx);
	VectorValue<double>& W = rb_ctx->W;
	VectorValue<double>& V = rb_ctx->V;
	
	// Keep body stationary.
	W.zero();
	V.zero();
	U_com.Zero();
	W_com.Zero();
	U_b.zero();
	
	// Communicate ghost values.
	U_b.close();
	
	return;
}// update_rigid_body_velocity
}

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
		const string exodus_filename = app_initializer->getExodusIIFilename();

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

        // Create the FE mesh.
        pout << "Creating the mesh...\n";
		Mesh mesh(NDIM);
		
		const double dx    = input_db->getDouble("DX");
		const double MFAC  = input_db->getDoubleWithDefault("MFAC", 1.0);
		const double ds = MFAC*dx;
		R = dx/2.0;
		const double w = 0.5 - dx;
		string elem_type = input_db->getString("ELEM_TYPE");
		const int n_x = ceil(2.0 * M_PI * R / ds);
		const int n_y = ceil(w / (ds / 2.0));
		MeshTools::Generation::build_square(
			mesh, n_x, n_y, 0.0, 2.0 * M_PI * R, 0.0, w, Utility::string_to_enum<ElemType>(elem_type));

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator("INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBFEMethod> ibfe_method_ops = new IBFEMethod("IBFEMethod", app_initializer->getComponentDatabase("IBFEMethod"), &mesh, app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        vector<Pointer<IBStrategy> > ib_ops_vec;
        ib_ops_vec.push_back(ibfe_method_ops);
        Pointer<IBStrategySet> ib_ops_set = new IBStrategySet(ib_ops_vec.begin(), ib_ops_vec.end());
        Pointer<IBExplicitHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_ops_set, navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>("CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>("GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Register the constraint parts.
        pout << "\nRegistering constraint part...\n";
		RigidBodyCtx rb_ctx;
		rb_ctx.dt = numeric_limits<double>::quiet_NaN();
		ibfe_method_ops->registerConstrainedPart(0);
		ibfe_method_ops->registerConstrainedVelocityFunction(
			(IBFEMethod::ConstrainedVelocityFcnPtr)&update_rigid_body_velocity,(void*)&rb_ctx,0);
		ibfe_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);

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
        AutoPtr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize FE data.
		EquationSystems* rb_systems = ibfe_method_ops->getFEDataManager(0)->getEquationSystems();
        ibfe_method_ops->initializeFEData();

        // Initialize hierarchy configuration and data on all patches.
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
                exodus_io ->write_timestep(exodus_filename, *rb_systems,
									   iteration_num/viz_dump_interval+1, loop_time);
            }
        }

		return 0;
		
        // Main time step loop.
        pout << "Entering main time step loop...\n";
        const double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();

            pout <<                                                       endl;
            pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            pout << "At beginning of timestep # " << iteration_num     << endl;
            pout << "Simulation time is " << loop_time                 << endl;

            dt = time_integrator->getMaximumTimeStepSize();
			rb_ctx.dt  = dt;
            time_integrator->advanceHierarchy(dt);

			loop_time += dt;

            pout <<                                                       endl;
            pout << "At end       of timestep # " << iteration_num     << endl;
            pout << "Simulation time is " << loop_time                 << endl;
            pout << "++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            pout <<                                                       endl;

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
                    exodus_io ->write_timestep(exodus_filename, *rb_systems,
											   iteration_num/viz_dump_interval+1, loop_time);
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

        for (int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
    }

    // Shutdown SAMRAI.
    SAMRAIManager::shutdown();
    return 0;
}// main
