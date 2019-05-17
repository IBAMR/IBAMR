// Filename: example.cpp
// Created by Amneet Bhalla.

// Copyright (c) 2018, Amneet Bhalla
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
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEDirectForcingKinematics.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

inline double
kernel(double x)
{
    x += 4.;
    const double x2 = x * x;
    const double x3 = x * x2;
    const double x4 = x * x3;
    const double x5 = x * x4;
    const double x6 = x * x5;
    const double x7 = x * x6;
    if (x <= 0.)
        return 0.;
    else if (x <= 1.)
        return .1984126984126984e-3 * x7;
    else if (x <= 2.)
        return .1111111111111111e-1 * x6 - .1388888888888889e-2 * x7 - .3333333333333333e-1 * x5 +
               .5555555555555556e-1 * x4 - .5555555555555556e-1 * x3 + .3333333333333333e-1 * x2 -
               .1111111111111111e-1 * x + .1587301587301587e-2;
    else if (x <= 3.)
        return .4333333333333333 * x5 - .6666666666666667e-1 * x6 + .4166666666666667e-2 * x7 - 1.500000000000000 * x4 +
               3.055555555555556 * x3 - 3.700000000000000 * x2 + 2.477777777777778 * x - .7095238095238095;
    else if (x <= 4.)
        return 9. * x4 - 1.666666666666667 * x5 + .1666666666666667 * x6 - .6944444444444444e-2 * x7 -
               28.44444444444444 * x3 + 53. * x2 - 54.22222222222222 * x + 23.59047619047619;
    else if (x <= 5.)
        return 96. * x3 - 22.11111111111111 * x4 + 3. * x5 - .2222222222222222 * x6 + .6944444444444444e-2 * x7 -
               245.6666666666667 * x2 + 344. * x - 203.9650793650794;
    else if (x <= 6.)
        return 483.5000000000000 * x2 - 147.0555555555556 * x3 + 26.50000000000000 * x4 - 2.833333333333333 * x5 +
               .1666666666666667 * x6 - .4166666666666667e-2 * x7 - 871.2777777777778 * x + 664.0904761904762;
    else if (x <= 7.)
        return 943.1222222222222 * x - 423.7000000000000 * x2 + 104.9444444444444 * x3 - 15.50000000000000 * x4 +
               1.366666666666667 * x5 - .6666666666666667e-1 * x6 + .1388888888888889e-2 * x7 - 891.1095238095238;
    else if (x <= 8.)
        return 416.1015873015873 - 364.0888888888889 * x + 136.5333333333333 * x2 - 28.44444444444444 * x3 +
               3.555555555555556 * x4 - .2666666666666667 * x5 + .1111111111111111e-1 * x6 - .1984126984126984e-3 * x7;
    else
        return 0.;
} // kernel

void
cylinder_kinematics(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();

    return;
} // cylinder_kinematics

// Function prototypes
static ofstream drag_stream, lift_stream;
static double R;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
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
run_example(int argc, char** argv)
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Setup user-defined kernel function.
        LEInteractor::s_kernel_fcn = &kernel;
        LEInteractor::s_kernel_fcn_stencil_size = 8;

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
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        R = input_db->getDouble("R");
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        // Ensure nodes on the surface are on the analytic boundary.
        MeshBase::element_iterator el_end = solid_mesh.elements_end();
        for (MeshBase::element_iterator el = solid_mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (!at_mesh_bdry) continue;
                for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                {
                    if (!elem->is_node_on_side(k, side)) continue;
                    Node& n = elem->node_ref(k);
                    n = R * n.unit();
                }
            }
        }
        solid_mesh.prepare_for_use();
        Mesh& mesh = solid_mesh;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
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

        // Create IBFE direct forcing kinematics object.
        Pointer<IBFEDirectForcingKinematics> df_kinematics_ops = new IBFEDirectForcingKinematics(
            "cylinder_dfk",
            app_initializer->getComponentDatabase("CylinderIBFEDirectForcingKinematics"),
            ib_method_ops,
            /*part*/ 0,
            /*register_for_restart*/ true);
        ib_method_ops->registerDirectForcingKinematics(df_kinematics_ops, /*part*/ 0);

        // Specify structure kinematics
        FreeRigidDOFVector solve_dofs;
        solve_dofs.setZero();
        df_kinematics_ops->setSolveRigidBodyVelocity(solve_dofs);
        df_kinematics_ops->registerKinematicsFunction(&cylinder_kinematics, NULL);

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

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

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

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
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save lift and drag coefficients.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.open("C_D.curve", ios_base::out | ios_base::trunc);
            lift_stream.open("C_L.curve", ios_base::out | ios_base::trunc);

            drag_stream.precision(10);
            lift_stream.precision(10);
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
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_io->write_timestep(
                        exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
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
                postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 mesh,
                                 equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.close();
            lift_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return true;
} // run_example

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

    System& F_system = equation_systems->get_system(IBFEMethod::FORCE_SYSTEM_NAME);
    NumericVector<double>* F_vec = F_system.solution.get();
    NumericVector<double>* F_ghost_vec = F_system.current_local_solution.get();
    F_vec->localize(*F_ghost_vec);
    const DofMap& dof_map = F_system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    std::unique_ptr<QBase> qrule = fe_type.default_quadrature_rule(dim);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double> >& phi = fe->get_phi();

    std::vector<std::vector<unsigned int> > dof_indices(NDIM);
    boost::multi_array<double, 2> F_node;
    VectorValue<double> F;

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
        get_values_for_interpolation(F_node, *F_ghost_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(F, qp, F_node, phi);
            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F(d) * JxW[qp];
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    static const double rho = 1.0;
    static const double U_max = 1.0;
    static const double D = 2.0 * R;
    if (SAMRAI_MPI::getRank() == 0)
    {
        drag_stream << loop_time << " " << -F_integral[0] / (0.5 * rho * U_max * U_max * D) << endl;
        lift_stream << loop_time << " " << -F_integral[1] / (0.5 * rho * U_max * U_max * D) << endl;
    }
    return;
} // postprocess_data
