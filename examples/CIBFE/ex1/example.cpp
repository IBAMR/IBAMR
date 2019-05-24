// Filename main.cpp
// Created on 17 Dec 2014 by Amneet Bhalla
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
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/petsc_matrix.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/CIBFEMethod.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Various model parameters and functions.
namespace ModelData
{
// Coordinate mapping function.
// This determines the initial position of the structure
// "s" is the Lagrangian coordinate system of the solid body.
// The code does not really care if s corresponds to body coordinates, or if they are general curvilinear coordinates.
// Here for a cube they are body coordinates (x = X(s,t) = R(t)*s + D(t))

static double SHIFT = 0.0;

// Nodal velocity function
void
ConstrainedNodalVel(libMesh::NumericVector<double>& U_k,
                    const RigidDOFVector& U,
                    libMesh::NumericVector<double>& X,
                    const Eigen::Vector3d& X_com,
                    libMesh::EquationSystems* equation_systems,
                    double /*data_time*/,
                    void* /*ctx*/)
{
    MeshBase& mesh = equation_systems->get_mesh();
    const unsigned int total_local_nodes = mesh.n_nodes_on_proc(SAMRAI_MPI::getRank());
    System& X_system = equation_systems->get_system<System>(CIBFEMethod::COORDS_SYSTEM_NAME);
    System& U_system = equation_systems->get_system<System>(CIBFEMethod::CONSTRAINT_VELOCITY_SYSTEM_NAME);
    const unsigned int X_sys_num = X_system.number();
    const unsigned int U_sys_num = U_system.number();

    std::vector<std::vector<unsigned int> > nodal_X_indices(NDIM), nodal_U_indices(NDIM);
    std::vector<std::vector<double> > nodal_X_values(NDIM), nodal_U_values(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        nodal_X_indices[d].reserve(total_local_nodes);
        nodal_U_indices[d].reserve(total_local_nodes);
        nodal_X_values[d].reserve(total_local_nodes);
        nodal_U_values[d].reserve(total_local_nodes);
    }

    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
    {
        const Node* const n = *it;
        if (n->n_vars(X_sys_num) && n->n_vars(U_sys_num))
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM && n->n_vars(U_sys_num) == NDIM);
#endif
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                nodal_X_indices[d].push_back(n->dof_number(X_sys_num, d, 0));
                nodal_U_indices[d].push_back(n->dof_number(U_sys_num, d, 0));
            }
        }
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        X.get(nodal_X_indices[d], nodal_X_values[d]);
    }

    // Set the cross-product matrix
    Eigen::Matrix3d W(Eigen::Matrix3d::Zero());
#if (NDIM == 2)
    W(0, 1) = -U[2];
    W(1, 0) = U[2];
#elif (NDIM == 3)
    W(0, 1) = -U[5];
    W(1, 0) = U[5];
    W(0, 2) = U[4];
    W(2, 0) = -U[4];
    W(1, 2) = -U[3];
    W(2, 1) = U[3];
#endif

    Eigen::Vector3d R(Eigen::Vector3d::Zero());
    Eigen::Vector3d WxR(Eigen::Vector3d::Zero());
    for (unsigned int k = 0; k < total_local_nodes; ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            R[d] = nodal_X_values[d][k] - X_com[d];
        }

        WxR = W * R;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            U_k.set(nodal_U_indices[d][k], U[d] + WxR[d]);
        }
    }

    U_k.close();
    return;
} // ConstrainedNodalVel

// Center of mass velocity
void
ConstrainedCOMVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com)
{
    U_com.setZero();
    W_com.setZero();

    return;
} // ConstrainedCOMVel
} // namespace ModelData
using namespace ModelData;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 Mesh& mesh,
                 EquationSystems* equation_systems,
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
bool
run_example(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::setMaxNumberPatchDataEntries(2048);
    SAMRAIManager::startup();

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
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string exodus_filename = app_initializer->getExodusIIFilename();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        // Create mesh based upon input file
        const std::string mesh_type = input_db->getString("MESH_TYPE");
        const bool use_boundary_mesh = input_db->getBoolWithDefault("USE_BOUNDARY_MESH", false);
        const double DX = input_db->getDouble("DX");
        SHIFT = input_db->getDouble("SHIFT");
        const double MFAC = input_db->getDouble("MFAC");
        const double DS = DX * MFAC;
        string elem_type = input_db->getString("ELEM_TYPE");

        // Create solid mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        if (mesh_type == "CUBIC")
        {
            const tbox::Array<double> cube_extents = input_db->getDoubleArray("CUBE_EXTENTS");
            std::vector<int> num_elems(3, 0);
            for (unsigned d = 0; d < NDIM; ++d)
            {
                num_elems[d] = round(cube_extents[d] / DS);
                pout << "Constructing cubic FEM mesh with " << num_elems[d] << " cells along dim " << d << "\n";
            }
            MeshTools::Generation::build_cube(solid_mesh,
                                              num_elems[0],
                                              num_elems[1],
                                              num_elems[2],
                                              -cube_extents[0] / 2.0,
                                              cube_extents[0] / 2.0,
                                              -cube_extents[1] / 2.0,
                                              cube_extents[1] / 2.0,
                                              -cube_extents[2] / 2.0,
                                              cube_extents[2] / 2.0,
                                              Utility::string_to_enum<ElemType>(elem_type));
        }
        else if (mesh_type == "SPHERICAL")
        {
            const double R = input_db->getDouble("RADIUS");
            if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
            {
                const int num_circum_nodes = ceil(2.0 * M_PI * R / DS);
                for (int k = 0; k < num_circum_nodes; ++k)
                {
                    const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                    solid_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
                }
                TriangleInterface triangle(solid_mesh);
                triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
                triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
                triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * DS * DS;
                triangle.insert_extra_points() = true;
                triangle.smooth_after_generating() = true;
                triangle.triangulate();
            }
            else
            {
                // NOTE: number of segments along boundary is 4*2^r.
                const double num_circum_segments = 2.0 * M_PI * R / DS;
                const int r = log2(0.25 * num_circum_segments);
                MeshTools::Generation::build_sphere(solid_mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
            }
        }
        else
        {
            TBOX_ERROR("Only CUBIC and SPHERICAL mesh types are supported...\n");
        }
        solid_mesh.prepare_for_use();

        // Create boundary mesh (if necessary).
        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        solid_mesh.boundary_info->sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        Mesh& mesh = use_boundary_mesh ? boundary_mesh : solid_mesh;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<CIBFEMethod> ib_method_ops =
            new CIBFEMethod("CIBFEMethod",
                            app_initializer->getComponentDatabase("CIBFEMethod"),
                            &mesh,
                            app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        // Krylov solver for extended INS system that solves for [u,p,L,U]
        Pointer<CIBStaggeredStokesSolver> c_stokes_solver =
            new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver",
                                         input_db->getDatabase("CIBStaggeredStokesSolver"),
                                         navier_stokes_integrator,
                                         ib_method_ops,
                                         "SP_");
        navier_stokes_integrator->setStokesSolver(c_stokes_solver);

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

        // Configure the IBFE solver.
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
        for (MeshBase::node_iterator it = mesh.nodes_begin(); it != mesh.nodes_end(); ++it)
        {
            Node* const n = *it;
            libMesh::Point& X = *n;

            for (int d = 0; d < NDIM; ++d)
            {
                X(d) = X(d) + SHIFT; // No need to convert here
            }
        }

        FreeRigidDOFVector solve_dofs;
        solve_dofs.setZero();
        ib_method_ops->setSolveRigidBodyVelocity(0, solve_dofs);
        ib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedCOMVel);

        // Create Eulerian initial condition specification objects.  These
        // objects also are used to specify exact solution values for error
        // analysis.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

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
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

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

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Get equation system for FE data.
        EquationSystems* equation_systems = fe_data_manager->getEquationSystems();

        // Setup variable data used to determine mobility matrix and visualization.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("wide_ctx");
        Pointer<SideVariable<NDIM, double> > u_var = new SideVariable<NDIM, double>("u_var", 1);
        Pointer<SideVariable<NDIM, double> > f_var = new SideVariable<NDIM, double>("f_var", 1);
        Pointer<CellVariable<NDIM, double> > p_var = new CellVariable<NDIM, double>("p_var", 1);
        Pointer<CellVariable<NDIM, double> > h_var = new CellVariable<NDIM, double>("h_var", 1);

        Pointer<CellVariable<NDIM, double> > u_plot_var = new CellVariable<NDIM, double>("u_plot_var", NDIM);
        Pointer<CellVariable<NDIM, double> > f_plot_var = new CellVariable<NDIM, double>("f_plot_var", NDIM);

        const IntVector<NDIM> ib_width = ib_method_ops->getMinimumGhostCellWidth();
        const IntVector<NDIM> ghost_width = 1;
        const IntVector<NDIM> no_width = 0;

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, ib_width);
        const int f_idx = var_db->registerVariableAndContext(f_var, ctx, ib_width);
        const int p_idx = var_db->registerVariableAndContext(p_var, ctx, ghost_width);
        const int h_idx = var_db->registerVariableAndContext(h_var, ctx, ghost_width);

        const int u_plot_idx = var_db->registerVariableAndContext(u_plot_var, ctx, no_width);
        const int f_plot_idx = var_db->registerVariableAndContext(f_plot_var, ctx, no_width);

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_idx, 0.0);
            level->allocatePatchData(f_idx, 0.0);
            level->allocatePatchData(p_idx, 0.0);
            level->allocatePatchData(h_idx, 0.0);
            level->allocatePatchData(u_plot_idx, 0.0);
            level->allocatePatchData(f_plot_idx, 0.0);
        }

        // Setup Eulerian vectors.
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const double vol_domain = hier_math_ops.getVolumeOfPhysicalDomain();
        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        SAMRAIVectorReal<NDIM, double> x_wide("wide_x", patch_hierarchy, coarsest_ln, finest_ln);
        SAMRAIVectorReal<NDIM, double> b_wide("wide_b", patch_hierarchy, coarsest_ln, finest_ln);
        x_wide.addComponent(u_var, u_idx, wgt_sc_idx);
        x_wide.addComponent(p_var, p_idx, wgt_cc_idx);
        b_wide.addComponent(f_var, f_idx, wgt_sc_idx);
        b_wide.addComponent(h_var, h_idx, wgt_cc_idx);
        x_wide.setToScalar(0.0);
        b_wide.setToScalar(0.0);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_writer = app_initializer->getVisItDataWriter();
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);
        if (dump_viz_data && uses_visit)
        {
            visit_writer->registerPlotQuantity("U", "VECTOR", u_plot_idx, 0);
            visit_writer->registerPlotQuantity("S_F", "VECTOR", f_plot_idx, 0);
            for (unsigned d = 0; d < NDIM; ++d)
            {
                if (d == 0)
                {
                    visit_writer->registerPlotQuantity("U_x", "SCALAR", u_plot_idx, d);
                    visit_writer->registerPlotQuantity("S_F_x", "SCALAR", f_plot_idx, d);
                }
                if (d == 1)
                {
                    visit_writer->registerPlotQuantity("U_y", "SCALAR", u_plot_idx, d);
                    visit_writer->registerPlotQuantity("S_F_y", "SCALAR", f_plot_idx, d);
                }
                if (d == 2)
                {
                    visit_writer->registerPlotQuantity("U_z", "SCALAR", u_plot_idx, d);
                    visit_writer->registerPlotQuantity("S_F_z", "SCALAR", f_plot_idx, d);
                }
            }
        }

        // Stokes solver and ghost fill schedule.
        const bool normalize_force = c_stokes_solver->getSaddlePointSolver()->getNormalizeSpreadForce();
        Pointer<StaggeredStokesSolver> LInv = c_stokes_solver->getSaddlePointSolver()->getStokesSolver();
        RefineAlgorithm<NDIM> ghost_fill_alg;
        ghost_fill_alg.registerRefine(u_idx, u_idx, u_idx, NULL);
        Pointer<RefineSchedule<NDIM> > ghost_fill_schd =
            ghost_fill_alg.createSchedule(patch_hierarchy->getPatchLevel(0));

        // Initialize time integrator.
        time_integrator->initializeHierarchyIntegrator(patch_hierarchy, gridding_algorithm);
        double dt = time_integrator->getMaximumTimeStepSize();
        double current_time = 0.0, new_time = dt, half_time = 0.5 * (current_time + new_time);
        time_integrator->preprocessIntegrateHierarchy(current_time, new_time);

        // Setup constraint force vector.
        Vec L;
        Vec* vL;
        ib_method_ops->getConstraintForce(&L, 0.0);
        VecNestGetSubVecs(L, NULL, &vL);
        PetscInt global_size;
        VecGetSize(vL[0], &global_size);
        const int total_nodes = mesh.n_nodes();
        TBOX_ASSERT(global_size == total_nodes * NDIM);
        pout << "\n\n Total number of nodes in the mesh   = " << total_nodes << "\n";
        pout << " Total number of elems in the mesh       = " << mesh.n_elem() << "\n";
        pout << " Total DOFs in the force/velocity vector = " << global_size << "\n";

        // Interpolated velocity vector
        Vec V;
        Vec* vV;
        VecDuplicate(L, &V);
        VecNestGetSubVecs(V, NULL, &vV);

        // Index information.
        System& V_sys = equation_systems->get_system<System>(CIBFEMethod::CONSTRAINT_VELOCITY_SYSTEM_NAME);
        const int V_sys_num = V_sys.number();
        PetscVector<double>& V_vec = *dynamic_cast<PetscVector<double>*>(V_sys.solution.get());
        PetscInt local_size;
        VecGetLocalSize(vV[0], &local_size);
        std::vector<PetscInt> mm_indices;
        mm_indices.reserve(local_size);
        std::vector<PetscScalar> mm_values(local_size);
        for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
        {
            const Node* const n = *it;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                mm_indices.push_back(n->dof_number(V_sys_num, d, 0));
            }
        }
        TBOX_ASSERT((PetscInt)mm_indices.size() == local_size);

        // PETSc dense mobility matrix.
        Mat MM;
        MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, global_size, global_size, NULL, &MM);
        MatMPIDenseSetPreallocation(MM, NULL);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Fill the dense matrix.
        pout << "\n\nCreating mobility matrix...\n\n";
        for (int col = 0; col < global_size; ++col)
        {
            pout << "\n\nColumn " << col << "...\n\n";
            VecSet(vL[0], 0.0);
            VecSetValue(vL[0], col, 1.0, INSERT_VALUES);
            VecAssemblyBegin(vL[0]);
            VecAssemblyEnd(vL[0]);
            b_wide.setToScalar(0.0);

            ib_method_ops->setConstraintForce(L, half_time);
            ib_method_ops->spreadForce(f_idx, NULL, std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
            if (normalize_force)
            {
                pout << "Subtracting mean force ...\n";
                ib_method_ops->subtractMeanConstraintForce(L, f_idx);
                const double mean_force = (1 / vol_domain) * hier_sc_data_ops.integral(f_idx, wgt_sc_idx);
                pout << "Mean of the spread force  = " << mean_force << "\n\n";
            }
            LInv->solveSystem(x_wide, b_wide);
            ghost_fill_schd->fillData(half_time);
            ib_method_ops->setInterpolatedVelocityVector(V, half_time);
            ib_method_ops->interpolateVelocity(u_idx,
                                               std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                               std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                               half_time);
            ib_method_ops->getInterpolatedVelocity(V, half_time);

            PetscScalar* v_array;
            VecGetArray(vV[0], &v_array);
            for (int j = 0; j < local_size; ++j)
            {
                mm_values[j] = v_array[V_vec.map_global_to_local_index(mm_indices[j])];
            }
            MatSetValues(MM, local_size, &mm_indices[0], 1, &col, &mm_values[0], INSERT_VALUES);
            VecRestoreArray(vV[0], &v_array);

            if (dump_viz_data && (col % viz_dump_interval == 0))
            {
                hier_math_ops.interp(u_plot_idx, u_plot_var, u_idx, u_var, NULL, 0.0, false);
                hier_math_ops.interp(f_plot_idx, f_plot_var, f_idx, f_var, NULL, 0.0, false);

                pout << "\n\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    visit_writer->writePlotData(patch_hierarchy,
                                                /*iteration_num*/ col,
                                                /*time*/ col);
                }
                if (uses_exodus)
                {
                    System& L_system = equation_systems->get_system<System>(CIBFEMethod::FORCE_SYSTEM_NAME);
                    PetscVector<double>* L_petsc_vec = dynamic_cast<PetscVector<double>*>(L_system.solution.get());
                    VecCopy(vL[0], L_petsc_vec->vec());
                    exodus_io->write_timestep(exodus_filename,
                                              *equation_systems,
                                              /*iteration_num*/ col / viz_dump_interval + 1,
                                              /*time*/ col);
                }
            }
        }
        MatAssemblyBegin(MM, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(MM, MAT_FINAL_ASSEMBLY);

        // Output mobility matrix for MATLAB viewing.
        pout << "\n\nWriting mobility matrix file in MATLAB format...\n\n";
        PetscViewer matlab_viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mobility_mat.dat", FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerPushFormat(matlab_viewer, PETSC_VIEWER_NATIVE);
        PetscObjectSetName((PetscObject)MM, "MobilityMatrix");
        MatView(MM, matlab_viewer);
        PetscViewerDestroy(&matlab_viewer);

        // Get the mass matrix
        std::pair<LinearSolver<double>*, SparseMatrix<double>*> filter =
            fe_data_manager->buildL2ProjectionSolver(CIBFEMethod::VELOCITY_SYSTEM_NAME);
        Mat mm = static_cast<PetscMatrix<double>*>(filter.second)->mat();

        // Output mass matrix for MATLAB viewing.
        pout << "\n\nWriting mass matrix file in MATLAB format...\n\n";
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "mass_mat.dat", FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerPushFormat(matlab_viewer, PETSC_VIEWER_NATIVE);
        PetscObjectSetName((PetscObject)mm, "MassMatrix");
        MatView(mm, matlab_viewer);
        PetscViewerDestroy(&matlab_viewer);

        // Cleanup memory for Lagrangian data.
        VecDestroy(&V);
        MatDestroy(&MM);
        ib_method_ops->postprocessIntegrateData(current_time, new_time, 1);

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return true;
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            Mesh& mesh,
            EquationSystems* equation_systems,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    file_name = data_dump_dirname + "/" + "fe_mesh.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    file_name += ".xda";
    mesh.write(file_name);
    file_name = data_dump_dirname + "/" + "fe_equation_systems.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
    return;
} // output_data
