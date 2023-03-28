// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Taylor Couette flow
// See Sec. 4.2 of "An immersed interface method for discrete surfaces"
// by Ebrahim M. Kolahdouz et al., Journal of Computational Physics 400 (2020) 108854
// To get the right convergence results this code was run with tight Stokes solver tolerance
// with the command line flag: -stokes_ksp_rtol 1.0e-10 -ksp_rtol 1.0e-10

// This code was tested against the commit 57fb379454ea3f8f50c476f10226cb6b520a11a0
// which is currently on branch iim-1 at https://github.com/drwells/IBAMR

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/face_quad.h>
#include <libmesh/face_quad4.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
// Tether (penalty) force functions.

static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double AA = 0.0;
static double BB = 0.0;
static double MU = 0.0;
static double R1 = 0.0;
static double R2 = 0.0;
static double OMEGA1 = 0.0;
static double OMEGA2 = 0.0;
static double shift = 0.0;
static double fac = 0.0;
static double L = 0.0;

void
tether_force_function_inner(VectorValue<double>& F,
                            const VectorValue<double>& n,
                            const VectorValue<double>& /*N*/,
                            const TensorValue<double>& /*FF*/,
                            const libMesh::Point& x,
                            const libMesh::Point& X,
                            Elem* const /*elem*/,
                            const unsigned short /*side*/,
                            const vector<const vector<double>*>& var_data,
                            const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                            double time,
                            void* /*ctx*/)
{
    const std::vector<double>& U = *var_data[0];

    double u_bndry_n = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        u_bndry_n += n(d) * U[d];
    }
    F(0) = kappa_s * (X(0) * cos(OMEGA1 * time) - X(1) * sin(OMEGA1 * time) - x(0));
    F(1) = kappa_s * (X(0) * sin(OMEGA1 * time) + X(1) * cos(OMEGA1 * time) - x(1));
    F(2) = 0.0;

    return;
} // tether_force_function

inline unsigned int
idx(const unsigned int nr, const unsigned int i, const unsigned int j)
{
    return i + j * nr;

    return libMesh::invalid_uint;
}

} // namespace ModelData
using namespace ModelData;

void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);

// Function prototypes
void pressure_convergence(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
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

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-16");
    PetscOptionsSetValue(nullptr, "-stokes_ksp_atol", "1e-10");

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Setup user-defined kernel function.

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string inner_exodus_filename = app_initializer->getExodusIIFilename("inner");

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
        OMEGA1 = input_db->getDouble("OMEGA1"); // radius of the inner circle
        OMEGA2 = input_db->getDouble("OMEGA2"); // radius of the outer circle

        AA = input_db->getDouble("AA"); // radius of the inner circle
        BB = input_db->getDouble("BB"); // radius of the outer circle
        shift = input_db->getDouble("SHIFT");
        fac = input_db->getDouble("FAC");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        R1 = input_db->getDouble("R1"); // radius of the inner circle
        R2 = input_db->getDouble("R2"); // radius of the outer circle
        L = input_db->getDouble("L");   // length in the z direction
        MU = input_db->getDouble("MU"); // length in the z direction

        Mesh inner_mesh(init.comm(), NDIM - 1);

        BoundaryInfo& boundary_info = inner_mesh.get_boundary_info();
        boundary_info.clear_boundary_node_ids();
        const unsigned int NXi_elem = ceil(L / ds);
        const unsigned int NRi_elem = ceil(2.0 * M_PI * R1 / ds);
        int node_id = 0;
        inner_mesh.reserve_nodes(NRi_elem * (NXi_elem + 1));
        inner_mesh.reserve_elem(NRi_elem * NXi_elem);

        for (unsigned int j = 0; j <= NXi_elem; j++)
        {
            for (unsigned int i = 0; i <= NRi_elem - 1; i++)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(NRi_elem);
                inner_mesh.add_point(libMesh::Point(R1 * cos(theta),
                                                    R1 * sin(theta),
                                                    -0.5 * L + L * static_cast<Real>(j) / static_cast<Real>(NXi_elem)),
                                     node_id++);
            }
        }

        for (unsigned int j = 0; j <= NXi_elem - 1; j++)
        {
            for (unsigned int i = 0; i <= NRi_elem - 2; i++)
            {
                Elem* elem = inner_mesh.add_elem(new Quad4);
                elem->set_node(0) = inner_mesh.node_ptr(idx(NRi_elem, i, j));
                elem->set_node(1) = inner_mesh.node_ptr(idx(NRi_elem, i + 1, j));
                elem->set_node(2) = inner_mesh.node_ptr(idx(NRi_elem, i + 1, j + 1));
                elem->set_node(3) = inner_mesh.node_ptr(idx(NRi_elem, i, j + 1));
            }
        }

        for (unsigned int j = 0; j <= NXi_elem - 1; j++)
        {
            Elem* elem = inner_mesh.add_elem(new Quad4);
            elem->set_node(0) = inner_mesh.node_ptr(idx(NRi_elem, NRi_elem - 1, j));
            elem->set_node(1) = inner_mesh.node_ptr(idx(NRi_elem, 0, j));
            elem->set_node(2) = inner_mesh.node_ptr(idx(NRi_elem, 0, j + 1));
            elem->set_node(3) = inner_mesh.node_ptr(idx(NRi_elem, NRi_elem - 1, j + 1));
        }

        const auto el_end = inner_mesh.elements_end();
        for (auto el = inner_mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    if (boundary_info.has_boundary_id(elem, side, 1) || boundary_info.has_boundary_id(elem, side, 3))
                    {
                        boundary_info.add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID);
                    }
                }
            }
        }
        inner_mesh.prepare_for_use();

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<IIMethod> ib_method_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         &inner_mesh,
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

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));

        IIMethod::LagSurfaceForceFcnData tether_force_inner_data(tether_force_function_inner, sys_data);

        ib_method_ops->registerLagSurfaceForceFunction(tether_force_inner_data);

        EquationSystems* inner_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);

        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));
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
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }
        // Create Eulerian body force function specification objects.
        // if (input_db->keyExists("ForcingFunction"))
        // {
        //    Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
        //         "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
        //    time_integrator->registerBodyForceFunction(f_fcn);
        // }

        // Set up visualization plot file writers.

        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> inner_exodus_io(uses_exodus ? new ExodusII_IO(inner_mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

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
                inner_exodus_io->write_timestep(
                    inner_exodus_filename, *inner_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
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
                    inner_exodus_io->write_timestep(inner_exodus_filename,
                                                    *inner_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                }
            }

            postprocess_data(patch_hierarchy,
                             navier_stokes_integrator,
                             inner_mesh,
                             inner_equation_systems,
                             iteration_num,
                             loop_time,
                             postproc_data_dump_dirname);
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
        }

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

        pressure_convergence(patch_hierarchy, p_idx, loop_time, postproc_data_dump_dirname);

        const int p_cloned_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
        const int coarsest_ln = 0;
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_cloned_idx, loop_time);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_cloned_idx, loop_time);
        }
        u_init->setDataOnPatchHierarchy(u_cloned_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_cloned_idx, p_var, patch_hierarchy, loop_time - 0.5 * dt);

        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_sc_data_ops.subtract(u_idx, u_idx, u_cloned_idx);
        pout << std::setprecision(16) << "Error in u at time " << loop_time << ":\n"
             << "  L1-norm:  " << hier_sc_data_ops.L1Norm(u_idx, wgt_sc_idx) << "\n"
             << "  L2-norm:  " << hier_sc_data_ops.L2Norm(u_idx, wgt_sc_idx) << "\n"
             << "  max-norm: " << hier_sc_data_ops.maxNorm(u_idx, wgt_sc_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        pout << " MU = " << MU << "\n"
             << "  dx:  " << dx << "\n"
             << "  dt: " << dt << "\n";

        if (input_db->getBool("USE_VELOCITY_JUMP_CONDITIONS"))
            pout << " Using the jump condition"
                 << "\n";
        else
            pout << " Using regular IB"
                 << "\n";
        if (dump_viz_data && uses_visit)
        {
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main

void
velocity_convergence(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                     const int u_idx,
                     const double /*data_time*/,
                     const string& /*data_dump_dirname*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
    hier_math_ops.resetLevels(finest_ln, finest_ln);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
    const double X_min[3] = { -0.45 * L, -0.45 * L, -0.45 * L };
    const double X_max[3] = { 0.45 * L, 0.45 * L, 0.45 * L };

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

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();

            const double* const patch_dx = patch_geom->getDx();

            // Entire box containing the required data.
            Box<NDIM> box(IndexUtilities::getCellIndex(
                              &X_min[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper),
                          IndexUtilities::getCellIndex(
                              &X_max[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper));
            // Part of the box on this patch
            Box<NDIM> trim_box = patch_box * box;
            BoxList<NDIM> iterate_box_list = trim_box;

            // Trim the box covered by the finer region
            BoxList<NDIM> covered_boxes;
            if (ln < finest_ln)
            {
                BoxArray<NDIM> refined_region_boxes;
                Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
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

                    const double y = patch_x_lower[1] + patch_dx[1] * (lower_idx(1) - patch_lower(1) + 0.5);
                    const double x = patch_x_lower[0] + patch_dx[0] * (lower_idx(0) - patch_lower(0));

                    double u_ex, v_ex;
                    if (sqrt(x * x + y * y) < R1)
                    {
                        u_ex = -OMEGA1 * y;
                        v_ex = OMEGA1 * x;
                    }
                    else if (sqrt(x * x + y * y) > R2)
                    {
                        u_ex = 0.0;
                        v_ex = 0.0;
                    }
                    else
                    {
                        u_ex = -y * (AA + BB / (x * x + y * y));
                        v_ex = x * (AA + BB / (x * x + y * y));
                    }

                    const double u0 = (*u_data)(SideIndex<NDIM>(lower_idx, 0, SideIndex<NDIM>::Lower));
                    const double v0 = (*u_data)(SideIndex<NDIM>(lower_idx, 1, SideIndex<NDIM>::Lower));

                    N_max += 1;
                    u_Eulerian_L2_norm += std::abs(u0 - u_ex) * std::abs(u0 - u_ex) * (*wgt_cc_data)(lower_idx);
                    u_Eulerian_L2_norm += std::abs(v0 - v_ex) * std::abs(v0 - v_ex) * (*wgt_cc_data)(lower_idx);

                    u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(u0 - u_ex));
                    u_Eulerian_max_norm = std::max(u_Eulerian_max_norm, std::abs(v0 - v_ex));
                }
            }
        }
    }

    SAMRAI_MPI::sumReduction(&N_max, 1);
    SAMRAI_MPI::sumReduction(&u_Eulerian_L2_norm, 1);
    SAMRAI_MPI::maxReduction(&u_Eulerian_max_norm, 1);

    u_Eulerian_L2_norm = sqrt(u_Eulerian_L2_norm);

    pout << " u_Eulerian_L2_norm = " << u_Eulerian_L2_norm << "\n\n";
    pout << " u_Eulerian_max_norm = " << u_Eulerian_max_norm << "\n\n";

    return;
} // velocity_convergence

void
pressure_convergence(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                     const int p_idx,
                     const double /*data_time*/,
                     const string& /*data_dump_dirname*/)
{
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
    hier_math_ops.resetLevels(finest_ln, finest_ln);
    const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

    const double X_min[3] = { -0.45 * L, -0.45 * L, -0.45 * L };
    const double X_max[3] = { 0.45 * L, 0.45 * L, 0.45 * L };
    // vector<double> pos_values;
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

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_x_lower = patch_geom->getXLower();
            const double* const patch_x_upper = patch_geom->getXUpper();

            const double* const patch_dx = patch_geom->getDx();

            // Entire box containing the required data.
            Box<NDIM> box(IndexUtilities::getCellIndex(
                              &X_min[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper),
                          IndexUtilities::getCellIndex(
                              &X_max[0], patch_x_lower, patch_x_upper, patch_dx, patch_lower, patch_upper));
            // Part of the box on this patch
            Box<NDIM> trim_box = patch_box * box;
            BoxList<NDIM> iterate_box_list = trim_box;

            // Trim the box covered by the finer region
            BoxList<NDIM> covered_boxes;
            if (ln < finest_ln)
            {
                BoxArray<NDIM> refined_region_boxes;
                Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
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

                    const double y = patch_x_lower[1] + patch_dx[1] * (cell_idx(1) - patch_lower(1) + 0.5);
                    const double x = patch_x_lower[0] + patch_dx[0] * (cell_idx(0) - patch_lower(0) + 0.5);
                    //~ double p_ex_qp = -2.*p_e*x/L + p_e;

                    const double p1 = (*p_data)(cell_idx);
                    double p_ex_qp;

                    N_max += 1;

                    if (sqrt(x * x + y * y) <= R1 - fac * patch_dx[0])
                    {
                        p_ex_qp = 0.5 * OMEGA1 * OMEGA1 * (x * x + y * y) + shift; // p1;
                    }
                    else if (sqrt(x * x + y * y) > (R1 - fac * patch_dx[0]) &&
                             sqrt(x * x + y * y) < (R1 + fac * patch_dx[0]))
                    {
                        p_ex_qp = p1;
                    }
                    else
                    {
                        p_ex_qp = 0.5 * AA * AA * (x * x + y * y) - 0.5 * BB * BB / (x * x + y * y) +
                                  AA * BB * log(x * x + y * y) + 0.5 * OMEGA1 * OMEGA1 * R1 * R1 -
                                  (0.5 * AA * AA * (R1 * R1) - 0.5 * BB * BB / (R1 * R1) + AA * BB * log(R1 * R1)) +
                                  shift;
                    }

                    p_Eulerian_L2_norm += std::abs(p1 - p_ex_qp) * std::abs(p1 - p_ex_qp) * (*wgt_cc_data)(cell_idx);
                    p_Eulerian_max_norm = std::max(p_Eulerian_max_norm, std::abs(p1 - p_ex_qp));
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(&N_max, 1);
    SAMRAI_MPI::sumReduction(&p_Eulerian_L2_norm, 1);
    SAMRAI_MPI::maxReduction(&p_Eulerian_max_norm, 1);

    p_Eulerian_L2_norm = sqrt(p_Eulerian_L2_norm);

    pout << " p_Eulerian_L2_norm = " << p_Eulerian_L2_norm << "\n\n";
    pout << " p_Eulerian_max_norm = " << p_Eulerian_max_norm << "\n\n";

    return;
} // pressure_convergence

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

    std::unique_ptr<FEBase> fe(FEBase::build(dim, dof_map.variable_type(0)));
    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<libMesh::Point>& q_point = fe->get_xyz();
    const vector<vector<double> >& phi = fe->get_phi();
    const vector<vector<VectorValue<double> > >& dphi = fe->get_dphi();
    // const std::vector<std::vector<double> >& dphi_dxi = fe->get_dphidxi();
    //  const std::vector<std::vector<double> >& dphi_deta = fe->get_dphideta();
    boost::array<VectorValue<double>, 2> dx_dxi;

    boost::array<const std::vector<std::vector<double> >*, NDIM - 1> dphi_dxi;
    dphi_dxi[0] = &fe->get_dphidxi();
    if (NDIM > 2) dphi_dxi[1] = &fe->get_dphideta();

    std::vector<double> U_qp_vec(NDIM);
    std::vector<const std::vector<double>*> var_data(1);
    var_data[0] = &U_qp_vec;
    std::vector<const std::vector<libMesh::VectorValue<double> >*> grad_var_data;
    void* force_fcn_ctx = NULL;

    TensorValue<double> FF_qp;
    boost::multi_array<double, 2> x_node, X_node, U_node, P_o_node, P_j_node;
    VectorValue<double> F_qp, U_qp, x_qp, X_qp, N, n;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        auto elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(X_node, X0_vec, dof_indices);
        get_values_for_interpolation(U_node, *U_ghost_vec, dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            interpolate(X_qp, qp, X_node, phi);
            jacobian(FF_qp, qp, x_node, dphi);
            interpolate(U_qp, qp, U_node, phi);

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U_qp_vec[d] = U_qp(d);
            }
            tether_force_function_inner(
                F_qp, n, N, FF_qp, x_qp, q_point[qp], elem, 0, var_data, grad_var_data, loop_time, force_fcn_ctx);

            for (int d = 0; d < NDIM; ++d)
            {
                F_integral[d] += F_qp(d) * JxW[qp];
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);

    {
        double WSS_L2_norm = 0.0, WSS_max_norm = 0.0;
        double U_L2_norm = 0.0, U_max_norm = 0.0;
        double P_L2_norm = 0.0, P_max_norm = 0.0;
        double disp_L2_norm = 0.0, disp_max_norm = 0.0;
        System& U_system = equation_systems->get_system<System>(IIMethod::VELOCITY_SYSTEM_NAME);
        System& WSS_system = equation_systems->get_system<System>(IIMethod::WSS_OUT_SYSTEM_NAME);
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
        std::unique_ptr<FEBase> fe(FEBase::build(dim, WSS_dof_map.variable_type(0)));

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
        double P_o_qp, P_j_qp;
        VectorValue<double> tau1, tau2;
        int qp_tot = 0;
        boost::multi_array<double, 2> U_node, WSS_node;
        boost::multi_array<double, 1> P_o_node, P_j_node;
        const auto el_begin = mesh.active_local_elements_begin();
        const auto el_end = mesh.active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            const auto elem = *el_it;
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
            get_values_for_interpolation(P_j_node, *P_j_ghost_vec, P_j_dof_indices);
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

                for (unsigned int k = 0; k < NDIM - 1; ++k)
                {
                    interpolate(dx_dxi[k], qp, x_node, *dphi_dxi[k]);
                }
                if (NDIM == 2)
                {
                    dx_dxi[1] = VectorValue<double>(0.0, 0.0, 1.0);
                }
                n = (dx_dxi[0].cross(dx_dxi[1])).unit();

                double ex_wss[NDIM];
                double ex_U[NDIM];

                ex_U[0] = (-x_qp(1) / sqrt(x_qp(0) * x_qp(0) + x_qp(1) * x_qp(1))) * R1 * OMEGA1;
                ex_U[1] = (x_qp(0) / sqrt(x_qp(0) * x_qp(0) + x_qp(1) * x_qp(1))) * R1 * OMEGA1;
                ex_U[2] = 0.0;
                ex_wss[0] = (-x_qp(1) / sqrt(x_qp(0) * x_qp(0) + x_qp(1) * x_qp(1))) * MU * (AA - BB / (R1 * R1));
                ex_wss[1] = (x_qp(0) / sqrt(x_qp(0) * x_qp(0) + x_qp(1) * x_qp(1))) * MU * (AA - BB / (R1 * R1));
                ex_wss[2] = 0.0;
                libMesh::Point X = q_point[qp];
                double p_ex_qp = 0.5 * OMEGA1 * OMEGA1 * R1 * R1 + shift;
                qp_tot += 1;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    U_L2_norm += (U_qp(d) - ex_U[d]) * (U_qp(d) - ex_U[d]) * JxW[qp];
                    U_max_norm = std::max(U_max_norm, std::abs(U_qp(d) - ex_U[d]));

                    WSS_L2_norm += (WSS_qp(d) - ex_wss[d]) * (WSS_qp(d) - ex_wss[d]) * JxW[qp];
                    WSS_max_norm = std::max(WSS_max_norm, std::abs(WSS_qp(d) - ex_wss[d]));
                }
                P_L2_norm += std::abs(P_o_qp - p_ex_qp) * std::abs(P_o_qp - p_ex_qp) * JxW[qp];
                P_max_norm = std::max(P_max_norm, std::abs(P_o_qp - p_ex_qp));

                disp_L2_norm += (X_qp(0) * cos(OMEGA1 * loop_time) - X_qp(1) * sin(OMEGA1 * loop_time) - x_qp(0)) *
                                (X_qp(0) * cos(OMEGA1 * loop_time) - X_qp(1) * sin(OMEGA1 * loop_time) - x_qp(0)) *
                                JxW[qp];
                disp_L2_norm += (X_qp(0) * sin(OMEGA1 * loop_time) + X_qp(1) * cos(OMEGA1 * loop_time) - x_qp(1)) *
                                (X_qp(0) * sin(OMEGA1 * loop_time) + X_qp(1) * cos(OMEGA1 * loop_time) - x_qp(1)) *
                                JxW[qp];
                disp_L2_norm += (X_qp(2) - x_qp(2)) * (X_qp(2) - x_qp(2)) * JxW[qp];
                disp_max_norm =
                    std::max(disp_max_norm,
                             std::abs(X_qp(0) * cos(OMEGA1 * loop_time) - X_qp(1) * sin(OMEGA1 * loop_time) - x_qp(0)));
                disp_max_norm =
                    std::max(disp_max_norm,
                             std::abs(X_qp(0) * sin(OMEGA1 * loop_time) + X_qp(1) * cos(OMEGA1 * loop_time) - x_qp(1)));
                disp_max_norm = std::max(disp_max_norm, std::abs(X_qp(2) - x_qp(2)));
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

        pout << " Lagrangian WSS_L2_norm = " << WSS_L2_norm << "\n\n";
        pout << " Lagrangian WSS_max_norm = " << WSS_max_norm << "\n\n";

        pout << " Lagrangian U_L2_norm = " << U_L2_norm << "\n\n";
        pout << " Lagrangian U_max_norm = " << U_max_norm << "\n\n";

        pout << " Lagrangian disp_L2_norm = " << disp_L2_norm << "\n\n";
        pout << " Lagrangian disp_max_norm = " << disp_max_norm << "\n\n";

        pout << "Lagrangian P_L2_norm = " << P_L2_norm << "\n\n";
        pout << "Lagrangian P_max_norm = " << P_max_norm << "\n\n";
    }

    return;
} // postprocess_data
