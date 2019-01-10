// Filename main.cpp
// Created on 26 Jul 2016 by Amneet Bhalla
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
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/CIBStochasticMethod.h>
#include <ibamr/CIBStochasticStokesSolver.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredStochasticForcing.h>
#include <ibamr/KrylovMobilitySolver.h>
#include <ibamr/RNG.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

//////////////////////////////////////////////////////////////////////////////

struct StructureCtx
{
    std::string name;
    int part;
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry;
    Pointer<CIBStochasticMethod> ib_method_ops;
    void (*ExternalLagForce)(Eigen::Vector3d& F_ext,
                             Eigen::Vector3d& T_ext,
                             const double* const X,
                             Eigen::Vector3d& domain_corner,
                             double dx);
    void (*LagLagForce)(Eigen::Vector3d& F_ext, Eigen::Vector3d& R, double dx);
}; // StructureCtx

struct StructureFT
{
    IBTK::Vector F;
    IBTK::Vector T;

}; // StructureBdry

// Center of mass velocity
void
ConstrainedCOMVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    U_com[1] = 1.0;

    return;
} // ConstrainedCOMOuterVel

void
NetExternalForceTorque(double data_time, Eigen::Vector3d& F_ext, Eigen::Vector3d& T_ext, void* ctx)
{
    static StructureFT* struct_FT = static_cast<StructureFT*>(ctx);
    F_ext = struct_FT->F;
    T_ext = struct_FT->T;
    return;
} // NetExternalForceTorque

void
ComputeLagrangianForceTorque(double data_time, Eigen::Vector3d& F_ext, Eigen::Vector3d& T_ext, StructureCtx* struct_ctx)
{
    F_ext.setZero();
    T_ext.setZero();

    // static StructureCtx* struct_ctx = static_cast<StructureCtx*>(ctx);
    const double* const domain_x_upper = struct_ctx->grid_geometry->getXUpper();
    const IntVector<NDIM>& periodic_shift = struct_ctx->grid_geometry->getPeriodicShift();
    const double* const DX = struct_ctx->grid_geometry->getDx();
    double dx = DX[0];

    int my_part = struct_ctx->part;

    Eigen::Vector3d domain_corner;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        domain_corner(d) = domain_x_upper[d];
        if (periodic_shift[d])
        {
            domain_corner(d) = -1.0;
        }
    }

    ///////////////////////////////////
    //                               //
    //         Lag--Wall Force       //
    //                               //
    ///////////////////////////////////

    const int coarsest_ln = 0;
    const int finest_ln = 0;

    // Eigen::Vector3d F_sum, T_sum;

    int s_max_free_dofs = 6;

    LDataManager* d_l_data_manager = struct_ctx->ib_method_ops->getLDataManager();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_l_data_manager->levelContainsLagrangianData(ln)) continue;
        const boost::multi_array_ref<double, 2>& X_array =
            *(d_l_data_manager->getLData("X", ln)->getLocalFormVecArray());
        const Pointer<LMesh> mesh = d_l_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        // Get position info.
        const boost::multi_array_ref<double, 2>& X0_array =
            *(d_l_data_manager->getLData("X0_unshifted", ln)->getLocalFormVecArray());
        Eigen::Matrix3d rotation_mat;
        Eigen::Vector3d dr = Eigen::Vector3d::Zero();
        Eigen::Vector3d R_dr = Eigen::Vector3d::Zero();

        // Get structures on this level.
        const std::vector<int> structIDs = d_l_data_manager->getLagrangianStructureIDs(ln);
        const unsigned structs_on_this_ln = static_cast<unsigned>(structIDs.size());

        std::vector<int> idx(s_max_free_dofs);
        std::vector<PetscScalar> F_body(s_max_free_dofs);

        Eigen::Vector3d F_blob, T_blob, COM;
        Eigen::Quaterniond Q;

        struct_ctx->ib_method_ops->getConfiguration(my_part, COM, Q);
        rotation_mat = Q.toRotationMatrix();
        struct_ctx->ib_method_ops->getInitialCOM(my_part, COM);

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int lag_idx = node_idx->getLagrangianIndex();
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_array[local_idx][0];
            const double* const X0 = &X0_array[local_idx][0];

            int struct_handle = 0;
            if (structs_on_this_ln > 1) struct_handle = struct_ctx->ib_method_ops->getStructureHandle(lag_idx);
            if (my_part != struct_handle) continue;

            if (struct_ctx->ExternalLagForce)
            {
                struct_ctx->ExternalLagForce(F_blob, T_blob, X, domain_corner, dx);

                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    dr[d] = X0[d] - COM(d);
                }
                R_dr = rotation_mat * dr;
                T_blob = -1.0 * F_blob.cross(R_dr);
            }
            else
            {
                F_blob.setZero();
                T_blob.setZero();
            }

            F_ext += F_blob;
            T_ext += T_blob;
        }
        d_l_data_manager->getLData("X", ln)->restoreArrays();
        ;
        d_l_data_manager->getLData("X0_unshifted", ln)->restoreArrays();
    }

    SAMRAI_MPI::sumReduction(&F_ext[0], 3);
    SAMRAI_MPI::sumReduction(&T_ext[0], 3);

    return;
}

void
ExternalLagForce(Eigen::Vector3d& F_ext,
                 Eigen::Vector3d& T_ext,
                 const double* const X,
                 Eigen::Vector3d& domain_corner,
                 double dx)
{
    F_ext.setZero();
    T_ext.setZero();

    double g = 0.0001539384;
    double blob_radius = 1.47 * dx;
    double repulsion_strength = 0.0165677856;
    double debye_length = 0.1 * blob_radius;

    double top;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        top = domain_corner(d);
        if (top <= 0)
        {
            continue;
        }
        if (X[d] > blob_radius)
        {
            F_ext(d) += (repulsion_strength / debye_length) * exp(-(X[d] - blob_radius) / debye_length);
        }

        if (X[d] < blob_radius)
        {
            F_ext(d) += (repulsion_strength / debye_length);
        }

        if (X[d] < (top - blob_radius))
        {
            F_ext(d) -= (repulsion_strength / debye_length) * exp(-((top - X[d]) - blob_radius) / debye_length);
        }

        if (X[d] > (top - blob_radius))
        {
            F_ext(d) -= (repulsion_strength / debye_length);
        }
    }

    F_ext(2) -= g;

    return;
}

void
LagLagForce(Eigen::Vector3d& F_ext, Eigen::Vector3d& R, double dx)
{
    double blob_radius = 1.47 * dx;
    double repulsion_strength = 0.0165677856;
    double debye_length = 0.1 * blob_radius;
    double r_norm = R.norm();
    double r_fact = 0;

    if (MathUtilities<double>::equalEps(0.0, r_norm))
    {
        r_norm = 1e-8;
    }

    if (r_norm > 2 * blob_radius)
    {
        r_fact = -((repulsion_strength / debye_length) * exp(-(r_norm - 2 * blob_radius) / debye_length) / r_norm);
    }
    else
    {
        r_fact = -((repulsion_strength / debye_length) / r_norm);
    }
    // r_fact = 0;
    F_ext[0] = r_fact * R[0];
    F_ext[1] = r_fact * R[1];
    F_ext[2] = r_fact * R[2];

    return;
}

void
ConstrainedNodalVel(Vec /*U_k*/, const RigidDOFVector& /*U*/, const Eigen::Vector3d& /*X_com*/, void* /*ctx*/)
{
    // intentionally left blank
    return;
} // ConstrainedNodalVel

ofstream U_stream;
ofstream Cfg_stream;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

// void
// CQout(Eigen::Vector3d& COM, Eigen::Quaterniond& Quat);

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
bool
run_example(int argc, char* argv[], double& end_time, double& end_u)
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    SAMRAIManager::setMaxNumberPatchDataEntries(2054);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "CIB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Read default Petsc options
        if (input_db->keyExists("petsc_options_file"))
        {
            std::string petsc_options_file = input_db->getString("petsc_options_file");
            PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, petsc_options_file.c_str(), PETSC_TRUE);
        }

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

        // CIB method
        const unsigned int num_structures = input_db->getIntegerWithDefault("num_structures", 1);

        Pointer<CIBStochasticMethod> ib_method_ops = new CIBStochasticMethod(
            "CIBStochasticMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);

        /*	Pointer<CIBMethod> ib_method_ops =
                    new CIBMethod("CIBMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);
        */
        double kT = input_db->getDouble("kT");
        double L_scale = input_db->getDouble("L_SCALE");
        ib_method_ops->setkT(kT);
        ib_method_ops->setLScale(L_scale);

        // Krylov solver for INS integrator that solves for [u,p,U,L]
        Pointer<CIBStochasticStokesSolver> CIBSolver =
            new CIBStochasticStokesSolver("CIBStochasticStokesSolver",
                                          input_db->getDatabase("CIBStaggeredStokesSolver"),
                                          navier_stokes_integrator,
                                          ib_method_ops,
                                          "SP_");

        /*        Pointer<CIBStaggeredStokesSolver> CIBSolver =
                    new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver",
                                                 input_db->getDatabase("CIBStaggeredStokesSolver"),
                                                 navier_stokes_integrator,
                                                 ib_method_ops,
                                                 "SP_");
        */
        // Register the Krylov solver with INS integrator
        navier_stokes_integrator->setStokesSolver(CIBSolver);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Specify structure kinematics
        FreeRigidDOFVector struct_0_free_dofs;
        struct_0_free_dofs << 1, 1, 1, 1, 1, 1;
        ib_method_ops->setSolveRigidBodyVelocity(0, struct_0_free_dofs);

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

        // Create stochastic forcing function specification object.
        Pointer<INSStaggeredStochasticForcing> f_fcn =
            new INSStaggeredStochasticForcing("INSStaggeredStochasticForcing",
                                              app_initializer->getComponentDatabase("INSStaggeredStochasticForcing"),
                                              navier_stokes_integrator);

        ///////////////////////////////////////////////////////////////////////
        //	                                                             //
        time_integrator->registerBodyForceFunction(f_fcn); //
                                                           //								     //
        ///////////////////////////////////////////////////////////////////////

        // Seed the random number generator.
        int seed = 0;
        if (input_db->keyExists("SEED"))
        {
            seed = input_db->getInteger("SEED");
        }
        else
        {
            TBOX_ERROR("Key data `seed' not found in input.");
        }
        RNG::parallel_seed(seed);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Set physical boundary operator used in spreading.
        ib_method_ops->setVelocityPhysBdryOp(time_integrator->getVelocityPhysBdryOp());

        // Register mobility matrices (if needed)
        std::string mobility_solver_type = input_db->getString("MOBILITY_SOLVER_TYPE");
        if (mobility_solver_type == "DIRECT")
        {
            std::string mat_name1 = "struct-1";
            std::vector<std::vector<unsigned> > struct_ids1;
            std::vector<unsigned> prototype_structs1;

            // Dense matrix
            prototype_structs1.push_back(0);

            struct_ids1.push_back(prototype_structs1);

            DirectMobilitySolver* direct_solvers = NULL;
            CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(NULL, &direct_solvers, NULL);

            direct_solvers->registerMobilityMat(
                mat_name1, prototype_structs1, EMPIRICAL, std::make_pair(LAPACK_SVD, LAPACK_SVD), 0);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name1, struct_ids1);
        }

        navier_stokes_integrator->setStokesSolverNeedsInit();

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
            output_data(patch_hierarchy,
                        ib_method_ops->getLDataManager(),
                        iteration_num,
                        loop_time,
                        postproc_data_dump_dirname);
        }

        if (SAMRAI_MPI::getRank() == 0)
        {
            char U_file[128];
            char C_file[128];
            strcpy(U_file, "./Lambda/U_");
            strcpy(C_file, "./Lambda/Config_");
            stringstream strs;
            strs << seed;
            string temp_str = strs.str();
            char const* file_no = temp_str.c_str();
            strcat(U_file, file_no);
            strcat(U_file, ".txt");
            strcat(C_file, file_no);
            strcat(C_file, ".txt");
            pout << U_file << "\n";
            pout << C_file << "\n";
            U_stream.open(U_file, std::ios_base::out | ios_base::trunc);
            U_stream.precision(16);
            Cfg_stream.open(C_file, std::ios_base::out | ios_base::trunc);
            Cfg_stream.precision(16);
        }

        bool reject = false;

        int p_0 = 0;

        StructureCtx struct0;
        struct0.grid_geometry = grid_geometry;
        struct0.ib_method_ops = ib_method_ops;
        struct0.ExternalLagForce = ExternalLagForce;
        struct0.part = p_0;
        struct0.name = "boom_0";

        StructureFT structFT;
        ComputeLagrangianForceTorque(0, structFT.F, structFT.T, &struct0);
        ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorque, &structFT, 0);

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        double time_reject = 0.0;
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
            if (ib_method_ops->flagRegrid())
            {
                time_integrator->regridHierarchy();
                navier_stokes_integrator->setStokesSolverNeedsInit();
            }

            ComputeLagrangianForceTorque(loop_time, structFT.F, structFT.T, &struct0);
            ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorque, &structFT, 0);

            pout << "reject value is: " << reject << "\n";
            RDV U0;
            Eigen::Vector3d COM;
            Eigen::Quaterniond Quat;
            if (!reject)
            {
                U_stream << time_reject << std::endl;
                Cfg_stream << time_reject << std::endl;
                for (unsigned itype = 0; itype < num_structures; itype++)
                {
                    ib_method_ops->getCurrentRigidBodyVelocity(itype, U0);
                    U_stream << loop_time << "\t" << U0(0) << "\t" << U0(1) << "\t" << U0(2) << "\t" << U0(3) << "\t"
                             << U0(4) << "\t" << U0(5) << std::endl;

                    ib_method_ops->getConfiguration(itype, COM, Quat);
                    Cfg_stream << loop_time << "\t" << COM(0) << "\t" << COM(1) << "\t" << COM(2) << "\t" << Quat.w()
                               << "\t" << Quat.x() << "\t" << Quat.y() << "\t" << Quat.z() << std::endl;
                }
            }

            // store velocity and loop time in variables passed to run_example as reference
            // used to test accuracy relative to past preformance in test_main.cpp
            end_u = std::abs(U0(2));
            end_time = loop_time;
            pout << "z-Velocity of rigid body is " << std::setprecision(10) << end_u << "\n";

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
                output_data(patch_hierarchy,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }

        if (SAMRAI_MPI::getRank() == 0)
        {
            U_stream.close();
            Cfg_stream.close();
        }

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return true;
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
            LDataManager* /*l_data_manager*/,
            const int iteration_num,
            const double loop_time,
            const string& /*data_dump_dirname*/)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    return;
} // output_data
