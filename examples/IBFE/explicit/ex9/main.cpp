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
#include <libmesh/system.h>
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
	
// Coordinate mapping function.
void
coordinate_mapping_function(libMesh::Point& X,
							const libMesh::Point& s,
							void* /*ctx*/)
{
	X(0) = s(0) + 0.0;
	X(1) = s(1) + 0.0;
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
						   Eigen::Vector3d& u_com,
						   Eigen::Vector3d& w_com,
						   libMesh::EquationSystems* /*equation_systems*/,
						   double /*data_time*/,
						   void* ctx)
{
	// Get rigid body data.
	RigidBodyCtx* rb_ctx = static_cast<RigidBodyCtx*>(ctx);
	VectorValue<double>& W = rb_ctx->W;
	VectorValue<double>& V = rb_ctx->V;
	
	// Keep body stationary.
	W.zero();
	V.zero();
	u_com.fill(0.0);
	w_com.fill(0.0);
	U_b.zero();
	
	// Communicate ghost values.
	U_b.close();
	
	return;
}// update_rigid_body_velocity
	
}

// Function prototypes
static ofstream drag_stream, lift_stream, U_L1_norm_stream, U_L2_norm_stream, U_max_norm_stream;

void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
					  Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
					  Mesh& mesh,
					  EquationSystems* equation_systems,
					  const int iteration_num,
					  const double dt,
					  const double loop_time,
					  const string& data_dump_dirname);
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

		// Create a simple FE mesh.
		pout << "Creating the mesh...\n";
		Mesh mesh(NDIM);
		const double dx = input_db->getDouble("DX");
		const double ds = input_db->getDouble("MFAC") * dx;
		string elem_type = input_db->getString("ELEM_TYPE");
		const double R = 0.5;
		const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
		for (int k = 0; k < num_circum_nodes; ++k)
		{
			const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
			mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
		}
		TriangleInterface triangle(mesh);
		triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
		triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
		triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
		triangle.insert_extra_points() = true;
		triangle.smooth_after_generating() = true;
		triangle.triangulate();
		mesh.prepare_for_use();
		
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
		
		// Open streams to save lift and drag coefficients and the norms of the
		// velocity.
		if (SAMRAI_MPI::getRank() == 0)
		{
			drag_stream.open((postproc_data_dump_dirname + "/C_D.curve").c_str(),
							 ios_base::out | ios_base::trunc);
			lift_stream.open((postproc_data_dump_dirname + "/C_L.curve").c_str(),
							 ios_base::out | ios_base::trunc);
			U_L1_norm_stream.open((postproc_data_dump_dirname + "/U_L1.curve").c_str(),
								   ios_base::out | ios_base::trunc);
			U_L2_norm_stream.open((postproc_data_dump_dirname  + "/U_L2.curve").c_str(),
								  ios_base::out | ios_base::trunc);
			U_max_norm_stream.open((postproc_data_dump_dirname + "/U_max.curve").c_str(),
								   ios_base::out | ios_base::trunc);
			
			drag_stream.precision(10);
			lift_stream.precision(10);
			U_L1_norm_stream.precision(10);
			U_L2_norm_stream.precision(10);
			U_max_norm_stream.precision(10);
		}

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
			
			if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
			{
				postprocess_data(patch_hierarchy,
								 navier_stokes_integrator,
								 mesh,
								 rb_systems,
								 iteration_num,
								 dt,
								 loop_time,
								 postproc_data_dump_dirname);
			}
			
        }

		// Close the logging streams.
		if (SAMRAI_MPI::getRank() == 0)
		{
			drag_stream.close();
			lift_stream.close();
			U_L1_norm_stream.close();
			U_L2_norm_stream.close();
			U_max_norm_stream.close();
		}
		
        for (int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
    }

    // Shutdown SAMRAI.
    SAMRAIManager::shutdown();
    return 0;
}// main

void postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
					  Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
					  Mesh& mesh,
					  EquationSystems* equation_systems,
					  const int /*iteration_num*/,
					  const double dt,
					  const double loop_time,
					  const string& /*data_dump_dirname*/)
{
	const unsigned int dim = mesh.mesh_dimension();
	{
		double F_integral[NDIM];
		for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
		System& F_system = equation_systems->get_system<System>(IBFEMethod::FORCE_SYSTEM_NAME);
		NumericVector<double>* F_vec = F_system.solution.get();
		NumericVector<double>* F_ghost_vec = F_system.current_local_solution.get();
		F_vec->localize(*F_ghost_vec);
		DofMap& F_dof_map = F_system.get_dof_map();
		std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);
		AutoPtr<FEBase> fe(FEBase::build(dim, F_dof_map.variable_type(0)));
		AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
		fe->attach_quadrature_rule(qrule.get());
		const std::vector<std::vector<double> >& phi = fe->get_phi();
		const std::vector<double>& JxW = fe->get_JxW();
		boost::multi_array<double, 2> F_node;
		const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
		for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
		{
			Elem* const elem = *el_it;
			fe->reinit(elem);
			for (unsigned int d = 0; d < NDIM; ++d)
			{
				F_dof_map.dof_indices(elem, F_dof_indices[d], d);
			}
			const int n_qp = qrule->n_points();
			const int n_basis = (int)F_dof_indices[0].size();
			get_values_for_interpolation(F_node, *F_ghost_vec, F_dof_indices);
			for (int qp = 0; qp < n_qp; ++qp)
			{
				for (int k = 0; k < n_basis; ++k)
				{
					for (int d = 0; d < NDIM; ++d)
					{
						F_integral[d] += F_node[k][d] * phi[k][qp] * JxW[qp];
					}
				}
			}
		}
		SAMRAI_MPI::sumReduction(F_integral, NDIM);
		static const double U_max = 1.0;
		static const double D = 1.0;
		if (SAMRAI_MPI::getRank() == 0)
		{
			drag_stream << loop_time << " " << -F_integral[0] / (0.5 * dt * U_max * U_max * D) << endl;
			lift_stream << loop_time << " " << -F_integral[1] / (0.5 * dt * U_max * U_max * D) << endl;
		}
	}
	
	{
		double U_L1_norm = 0.0, U_L2_norm = 0.0, U_max_norm = 0.0;
		System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
		NumericVector<double>* U_vec = U_system.solution.get();
		NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
		U_vec->localize(*U_ghost_vec);
		DofMap& U_dof_map = U_system.get_dof_map();
		std::vector<std::vector<unsigned int> > U_dof_indices(NDIM);
		AutoPtr<FEBase> fe(FEBase::build(dim, U_dof_map.variable_type(0)));
		AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, FIFTH);
		fe->attach_quadrature_rule(qrule.get());
		const std::vector<std::vector<double> >& phi = fe->get_phi();
		const std::vector<double>& JxW = fe->get_JxW();
		VectorValue<double> U_qp;
		boost::multi_array<double, 2> U_node;
		const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
		const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
		for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
		{
			Elem* const elem = *el_it;
			fe->reinit(elem);
			for (unsigned int d = 0; d < NDIM; ++d)
			{
				U_dof_map.dof_indices(elem, U_dof_indices[d], d);
			}
			const int n_qp = qrule->n_points();
			get_values_for_interpolation(U_node, *U_ghost_vec, U_dof_indices);
			for (int qp = 0; qp < n_qp; ++qp)
			{
				interpolate(U_qp, qp, U_node, phi);
				for (unsigned int d = 0; d < NDIM; ++d)
				{
					U_L1_norm += std::abs(U_qp(d)) * JxW[qp];
					U_L2_norm += U_qp(d) * U_qp(d) * JxW[qp];
					U_max_norm = std::max(U_max_norm, std::abs(U_qp(d)));
				}
			}
		}
		SAMRAI_MPI::sumReduction(&U_L1_norm, 1);
		SAMRAI_MPI::sumReduction(&U_L2_norm, 1);
		SAMRAI_MPI::maxReduction(&U_max_norm, 1);
		U_L2_norm = sqrt(U_L2_norm);
		if (SAMRAI_MPI::getRank() == 0)
		{
			U_L1_norm_stream << loop_time << " " << U_L1_norm << endl;
			U_L2_norm_stream << loop_time << " " << U_L2_norm << endl;
			U_max_norm_stream << loop_time << " " << U_max_norm << endl;
		}
	}
	return;
} // postprocess_data
