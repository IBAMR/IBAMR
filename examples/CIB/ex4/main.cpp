// Active squimers suspension 

// Filename main.cpp
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

#include <Eigen/Geometry>
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
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/CIBStandardInitializer.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/KrylovMobilitySolver.h>
#include <ibamr/app_namespaces.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibamr/RNG.h>

//////////////////////////////////////////////////////////////////////////////
static double B1,B2;
void SlipVelocity(const unsigned part, Vec U_k,  const RigidDOFVector& U, Vec X, const Eigen::Vector3d& Xin_com, double data_time, void* ctx, CIBMethod* cib_method)
{
    Vec Vslip;
    VecDuplicate(U_k, &Vslip);
    VecSet(Vslip, 0.0);

    const int struct_ln = cib_method->getStructuresLevelNumber();
    // Get the position data
    double* V_array=NULL; 
    const unsigned num_of_nodes = cib_method->getNumberOfNodes(part);
    unsigned size = num_of_nodes*NDIM;

    static std::vector<Eigen::Vector3d> center_of_mass_initial;
    static bool initializeCOM=true;
    if (initializeCOM)
    {
	const int num_rigid_parts = cib_method->getNumberOfRigidStructures();
	center_of_mass_initial.resize(num_rigid_parts, Eigen::Vector3d::Zero());
	cib_method->computeInitialCOMOfStructures(center_of_mass_initial);
	initializeCOM=false;
    }

    const Eigen::Vector3d& X_com  =  center_of_mass_initial[part];
    Eigen::Quaterniond* Q = cib_method->getBodyQuaternion(part,true);
    Eigen::Matrix3d  body_rot_matrix=Q->toRotationMatrix();

    
    if (!SAMRAI_MPI::getRank())
    {
	V_array = new double[size];
	for (unsigned i=0; i<num_of_nodes; ++i)
	{
	    const IBTK::Point& X = cib_method->getStandardInitializer()->getVertexPosn(struct_ln,part,i);
	    
	    Eigen::Vector3d coord;
	    for (unsigned d=0; d<NDIM; ++d) coord[d]=X[d]- X_com[d];
	    
	    double r=coord.norm();
	    double r_proj=sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
	    double cos_theta=coord[0]/r;
	    double sin_theta=r_proj/r;
	    double v_theta=(B1*sin_theta+B2*cos_theta*sin_theta);
	    
	    Eigen::Vector3d temp_vec = Eigen::Vector3d::Zero();
	    temp_vec[0] = -v_theta*sin_theta;
	    if (r_proj>1e-9)
	    {
		temp_vec[1] = v_theta*cos_theta*coord[1]/r_proj;
		temp_vec[2] = v_theta*cos_theta*coord[2]/r_proj;
	    }
	    Eigen::Vector3d rot_vec=body_rot_matrix*temp_vec;
	    
	    for (unsigned d=0; d<NDIM; ++d) V_array[i*NDIM+d] = rot_vec[d];
	}
    }
    cib_method->copyArrayToVec(Vslip, V_array, (std::vector<unsigned int>(1, part)), NDIM, 0);
    

    // Wrap the PETSc V into LData
    std::vector<int> nonlocal_indices;
    LData V_data("V", U_k, nonlocal_indices, false);

    boost::multi_array_ref<double, 2>& V_data_array = *V_data.getLocalFormVecArray();
    
    const Pointer<LMesh> mesh = cib_method->getLDataManager()->getLMesh(struct_ln);
    const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
    const std::pair<int, int>& part_idx_range = cib_method->getLDataManager()->getLagrangianStructureIndexRange(part, struct_ln);
    
    for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
    {
	const LNode* const node_idx = *cit;
	const int lag_idx = node_idx->getLagrangianIndex();
	Eigen::Vector3d dr;
	Eigen::Vector3d R_dr;
	
	
	if (part_idx_range.first <= lag_idx && lag_idx < part_idx_range.second)
	{
	    const int local_idx = node_idx->getLocalPETScIndex();
	    double* const V_node = &V_data_array[local_idx][0];
	    
	    const IBTK::Point& X =  cib_method->getStandardInitializer()->getVertexPosn(struct_ln,part,lag_idx-part_idx_range.first);
	    for (unsigned int d = 0; d < NDIM; ++d)  dr[d] = X[d] - X_com[d];
	    
	    R_dr = body_rot_matrix*dr; 
	    
	    V_node[0] = U[0] + U[4] * R_dr[2] - U[5] * R_dr[1];
	    V_node[1] = U[1] + U[5] * R_dr[0] - U[3] * R_dr[2];
	    V_node[2] = U[2] + U[3] * R_dr[1] - U[4] * R_dr[0];
	}
    }
    
    // Restore underlying arrays
    V_data.restoreArrays();
    
    VecAXPY(U_k, 1.0, Vslip);

    VecDestroy(&Vslip);
    if (!SAMRAI_MPI::getRank()) delete[] V_array;
	
    return;
}

void NetExternalForceTorque(double /*data_time*/, Eigen::Vector3d& F_ext, Eigen::Vector3d& T_ext)
{
    double F[NDIM], T[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)   RNG::genrand(F+d);
    for (unsigned int d = 0; d < NDIM; ++d)   RNG::genrand(T+d);

    F_ext << B1*(F[0]-0.5), B1*(F[1]-0.5), B1*(F[2]-0.5);
    T_ext << B2*(T[0]-0.5), B2*(T[1]-0.5), B2*(T[2]-0.5);

    return;
} // NetExternalForceTorque

void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 LDataManager* l_data_manager,
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
int main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    SAMRAIManager::setMaxNumberPatchDataEntries(2054);
    RNG::parallel_seed(1);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

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

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.

        // INS integrator
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

	//**********************************************
	//parameters for structures configuration   
	unsigned int num_structures=0;
	//First get number of structure using modified IBStandardInitializer
	const unsigned num_structs_types = app_initializer->getComponentDatabase("CIBStandardInitializer")->getArraySize("structure_names");

	std::vector<int> structs_clones_num(num_structs_types);
	app_initializer->getComponentDatabase("CIBStandardInitializer")->getIntegerArray("structs_clones_num", &structs_clones_num[0], num_structs_types);

	std::vector<std::string> structure_type_names(num_structs_types);
	app_initializer->getComponentDatabase("CIBStandardInitializer")->getStringArray("structure_names", &structure_type_names[0], num_structs_types);

	for(unsigned itype=0;itype<num_structs_types;itype++)  num_structures +=structs_clones_num[itype];

        // CIB method
//        const unsigned int num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        Pointer<CIBMethod> ib_method_ops =
            new CIBMethod("CIBMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);

        // Krylov solver for INS integrator that solves for [u,p,U,L]
        Pointer<CIBStaggeredStokesSolver> CIBSolver =
            new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver", input_db->getDatabase("CIBStaggeredStokesSolver"),
                                         navier_stokes_integrator, ib_method_ops, "SP_");

        // Register the Krylov solver with INS integrator
        navier_stokes_integrator->setStokesSolver(CIBSolver);

        Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
            "IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops,
            navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector, box_generator, load_balancer);
        // Configure the IB solver.
        Pointer<CIBStandardInitializer> ib_initializer = new CIBStandardInitializer(
            "CIBStandardInitializer", app_initializer->getComponentDatabase("CIBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
	ib_method_ops->registerStandardInitializer(ib_initializer);

        std::string example_to_run = input_db->getString("run_example");
	if (example_to_run=="SQUIMER")
	{
	    //Read  constants B1 and B2 for active slip
	    B1=input_db->getDouble("B1_squimer");
	    B2=input_db->getDouble("B2_squimer");
	}
	else if (example_to_run=="RANDOM_FORCE")
	{
	    //Read  constants B1 and B2 for active slip
	    B1=input_db->getDouble("DRAG_SCALE");
	    B2=input_db->getDouble("TORQUE_SCALE");
	} else
	{
	    TBOX_ERROR("Only two possible examples SQUIMER or RANDOM_FORCE" << std::endl);
	}

        // Specify kinematics os a suspension as all free moving structures
	for(unsigned i=0;i<num_structures;i++)
	{
	    FreeRigidDOFVector free_dofs;
	    free_dofs << 1, 1, 1, 1, 1, 1;
	    ib_method_ops->setSolveRigidBodyVelocity(i, free_dofs);
	    //Register slip velocity function for saddle point solver RHS 
	    //register slip velocity function
	    if (example_to_run=="SQUIMER")
	    {
		ib_method_ops->registerConstrainedVelocityFunction(&SlipVelocity, NULL, NULL, i);
	    }
	    else if (example_to_run=="RANDOM_FORCE")
	    {
		ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorque, NULL, i);
	    }
	}


        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create boundary condition specification objects (when necessary).
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

                Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

	//First initialize orientation vectors
	//ib_method_ops->updateBodyRotationMatrix();
	// Register mobility matrices (if needed)
        std::string mobility_solver_type = input_db->getString("MOBILITY_SOLVER_TYPE");
        if (mobility_solver_type == "DIRECT")
        {
            DirectMobilitySolver* direct_solvers = NULL;
            KrylovMobilitySolver* krylov_solvers = NULL;
            CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(&krylov_solvers,
                                                                                          &direct_solvers);
	    const unsigned num_nodes = SAMRAI_MPI::getNodes();
	    unsigned prototype_ID=0;
	    unsigned counter =0;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *
 *  Here are  examples of construction mobility matrix in parallel
 *
 *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1*/

#if 0
	    //Each clone has its own copy of mobility matrix
	    //This is example that all structures have their own mobility matrix 
	    //distributed parallel
	    for(unsigned itype=0;itype<num_structs_types;itype++)
	    {
		//use all availabel processors for efficiency
		const unsigned num_clones = structs_clones_num[itype];
		
		for(unsigned istruct=0; istruct<num_clones; istruct++) 
		{
		    //different mat names over all processors
		    std::stringstream convert;
		    convert<<itype;
		    std::string mat_name = "struct-"+convert.str();
		    convert<<istruct;
		    mat_name += convert.str();
		    
		    direct_solvers->registerMobilityMat(mat_name, prototype_ID+istruct, EMPIRICAL, LAPACK_SVD, counter%num_nodes);
		    
		    std::vector<std::vector<unsigned> > struct_ids;
		    struct_ids.push_back(std::vector<unsigned int>(1, prototype_ID+istruct));
		    direct_solvers->registerStructIDsWithMobilityMat(mat_name, struct_ids);
		    counter++;
		}
		prototype_ID +=num_clones;
	    }
	
#else
	    //More memory efficient way to manage mobility matrices.
	    //Uses only one copy of mobility matrix for each prototype on each processor
	    //coresponding to prototype structure (which is the first in subset) distributed parallel
	    
	    for(unsigned itype=0;itype<num_structs_types;itype++)
	    {
		//use all availabel processors for efficiency
		const unsigned num_clones = structs_clones_num[itype];
		
		for(unsigned inode=0; inode < std::min(num_nodes, num_clones); inode++) 
		{
		    //different mat names over all processors
		    std::stringstream convert;
		    convert<<itype;
		    std::string mat_name = "struct-"+convert.str();
		    convert<<inode;
		    mat_name += convert.str();
		    
		    direct_solvers->registerMobilityMat(mat_name, prototype_ID, EMPIRICAL, LAPACK_SVD, counter%num_nodes);
		    counter++;
		    
		    std::vector<std::vector<unsigned> > struct_ids;
		    std::vector<unsigned> local_struct_ids;
		    for(unsigned istruct=inode; istruct<num_clones; istruct += num_nodes) 		    
		    {
			
			struct_ids.push_back(std::vector<unsigned int>(1, prototype_ID+istruct));
		    }
		    
		    direct_solvers->registerStructIDsWithMobilityMat(mat_name, struct_ids);
		}
		prototype_ID +=num_clones;
	    }
#endif
	}

	// Deallocate initialization objects.
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
        if (dump_postproc_data)
        {
            output_data(patch_hierarchy, ib_method_ops->getLDataManager(), iteration_num, loop_time,
                        postproc_data_dump_dirname);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();

            pout << "Advancing hierarchy by timestep size dt = " << dt << "\n";
	    
            if (time_integrator->atRegridPoint()) navier_stokes_integrator->setStokesSolverNeedsInit();
		    
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy, ib_method_ops->getLDataManager(), iteration_num, loop_time,
                            postproc_data_dump_dirname);
            }
        }

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main

void output_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 LDataManager* /*l_data_manager*/,
                 const int iteration_num,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    return;
} // output_data
