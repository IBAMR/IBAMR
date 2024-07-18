// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/matlab_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/BrinkmanPenalizationRigidBodyDynamics.h>
#include <ibamr/FESurfaceDistanceEvaluator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBInterpolantHierarchyIntegrator.h>
#include <ibamr/IBInterpolantMethod.h>
#include <ibamr/IBLevelSetMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application specific includes.
#include "LSLocateGasInterface.h"

int coarsest_ln, max_finest_ln;
double dx, ds;

// Struct to maintain the properties of the circular interface
struct CircularInterface
{
    Eigen::Vector3d X0;
    double R;
    double rho_solid;
    double g_y;
};
CircularInterface circle;

// Struct to reset solid level set
struct SolidLevelSetResetter
{
    SAMRAIPointer<IBInterpolantMethod> ib_interp_ops;
    SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
    SAMRAIPointer<CellVariableNd<double> > ls_solid_var;
    SAMRAIPointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd;
};

void
reset_solid_level_set_callback_fcn(double current_time, double new_time, int /*cycle_num*/, void* ctx)
{
    SolidLevelSetResetter* resetter = static_cast<SolidLevelSetResetter*>(ctx);
    resetter->ib_interp_ops->copyEulerianDataToIntegrator(new_time);

    // Get the new centroid of the body
    const double dt = new_time - current_time;
    Eigen::Vector3d XCOM_current = resetter->bp_rbd->getCurrentCOMPosn();
    Eigen::Vector3d XCOM_new = XCOM_current + dt * (resetter->bp_rbd->getNewCOMTransVelocity());

    // Set a large value away from the solid body.
    SAMRAIPointer<PatchHierarchyNd> patch_hier = resetter->adv_diff_integrator->getPatchHierarchy();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
    const int ls_solid_idx =
        var_db->mapVariableAndContextToIndex(resetter->ls_solid_var, resetter->adv_diff_integrator->getNewContext());

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        SAMRAIPointer<PatchLevelNd> patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(patch_level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = patch_level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            const SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::IndexNd& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            SAMRAIPointer<CellDataNd<double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                const hier::IndexNd& ci = it();
                Eigen::Vector3d coord = Eigen::Vector3d::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                const double distance =
                    std::sqrt(std::pow((coord[0] - XCOM_new(0)), 2.0) + std::pow((coord[1] - XCOM_new(1)), 2.0)
#if (NDIM == 3)
                              + std::pow((coord[2] - XCOM_new(2)), 2.0)
#endif
                                  ) -
                    circle.R;

                (*ls_solid_data)(ci) = distance;
                /* if (distance >= 7 * patch_dx[0])
                 {
                     (*ls_solid_data)(ci) = 123.0;
                 }*/
            }
        }
    }

    return;
}

void
calculate_distance_analytically(SAMRAIPointer<PatchHierarchyNd> patch_hierarchy, int E_idx)
{
    int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
    {
        SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::IndexNd& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            SAMRAIPointer<CellDataNd<double> > E_data = patch->getPatchData(E_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                CellIndexNd ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                const double distance =
                    std::sqrt(std::pow((coord[0] - circle.X0(0)), 2.0) + std::pow((coord[1] - circle.X0(1)), 2.0)
#if (NDIM == 3)
                              + std::pow((coord[2] - circle.X0(2)), 2.0)
#endif
                    );

                (*E_data)(ci) = distance - circle.R;
            }
        }
    }

    return;
} // calculate_distance_analytically

void
calculate_error_near_band(SAMRAIPointer<PatchHierarchyNd> patch_hierarchy,
                          int E_idx,
                          int d_idx,
                          int W_idx,
                          double& E_interface,
                          int& num_interface_pts,
                          double& volume_near_interface)

{
    E_interface = 0.0;
    num_interface_pts = 0;
    volume_near_interface = 0.0;
    int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
    {
        SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            SAMRAIPointer<CellDataNd<double> > D_data = patch->getPatchData(d_idx);
            SAMRAIPointer<CellDataNd<double> > E_data = patch->getPatchData(E_idx);
            SAMRAIPointer<CellDataNd<double> > W_data = patch->getPatchData(W_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                CellIndexNd ci(it());
                const double phi = (*D_data)(ci);
                const double err = (*E_data)(ci);
                const double dV = (*W_data)(ci);
                SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();

                if (std::abs(err) <= 4.0 * patch_dx[0])
                {
                    if (std::abs(phi) <= 1.2 * patch_dx[0])
                    {
                        E_interface += std::abs(err - phi) * dV;
                        volume_near_interface += dV;
                        num_interface_pts++;
                    }
                    (*E_data)(ci) = std::abs(err - phi);
                }
                else
                {
                    (*E_data)(ci) = 0.0;
                }
            }
        }
    }
    num_interface_pts = IBTK_MPI::sumReduction(num_interface_pts);
    E_interface = IBTK_MPI::sumReduction(E_interface);
    volume_near_interface = IBTK_MPI::sumReduction(volume_near_interface);

    return;
} // calculate_error_near_band

void
generate_interp_mesh(const unsigned int& /*strct_num*/,
                     const int& ln,
                     int& num_vertices,
                     std::vector<IBTK::Point>& vertex_posn,
                     void* /*ctx*/)
{
    if (ln != max_finest_ln)
    {
        num_vertices = 0;
        vertex_posn.resize(num_vertices);
    }

    // Create a block mesh around the FE mesh.
    double X_lower[2] = { circle.X0[0] - 2.5 * circle.R, circle.X0[1] - 2.5 * circle.R };
    double X_upper[2] = { circle.X0[0] + 2.5 * circle.R, circle.X0[1] + 2.5 * circle.R };

    int Nx = static_cast<int>((X_upper[0] - X_lower[0]) / dx);
    int Ny = static_cast<int>((X_upper[1] - X_lower[1]) / dx);
    num_vertices = Nx * Ny;
    vertex_posn.resize(num_vertices);

    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            IBTK::Point& X = vertex_posn[j * Nx + i];
            X(0) = X_lower[0] + i * dx;
            X(1) = X_lower[1] + j * dx;
        }
    }

    return;
} // generate_interp_mesh

void
imposed_kinematics(double /*data_time*/,
                   int /*cycle_num*/,
                   Eigen::Vector3d& U_com,
                   Eigen::Vector3d& W_com,
                   void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // imposed_kinematics

void
external_force_torque(double /*data_time*/, int /*cycle_num*/, Eigen::Vector3d& F, Eigen::Vector3d& T, void* /*ctx*/)
{
    F.setZero();
    F[1] = circle.rho_solid * M_PI * std::pow(circle.R, 2) * circle.g_y;
    T.setZero();
    return;
} // imposed_kinematics

static double shift_x, shift_y;
#if (NDIM == 3)
static double shift_z;
#endif

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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IBLevelSet.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();
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
        if (dump_restart_data && (restart_dump_interval > 0) && !restart_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(restart_dump_dirname);
        }

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Setup solid information
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);

        // Create mesh based upon input file
        dx = input_db->getDouble("DX");
        ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const bool use_vol_extracted_bdry_mesh = input_db->getBool("USE_VOLUME_EXTRACTED_BOUNDARY_MESH");
        shift_x = circle.X0[0];
        shift_y = circle.X0[1];
#if (NDIM == 3)
        shift_z = circle.X0[2];
#endif
        if (input_db->keyExists("XDA_FILENAME"))
        {
            TBOX_ASSERT(elem_type == "TRI3" || elem_type == "TRI6");

            std::string filename = input_db->getString("XDA_FILENAME");
            MatlabIO distmesh(solid_mesh);
            distmesh.read(filename);

            if (elem_type == "TRI6") solid_mesh.all_second_order();
        }
        else if (input_db->keyExists("GMSH_FILENAME"))
        {
            TBOX_ASSERT(elem_type == "TRI3" || elem_type == "TRI6");

            std::string filename = input_db->getString("GMSH_FILENAME");
            solid_mesh.read(filename);

            if (elem_type == "TRI6") solid_mesh.all_second_order();
        }
        else if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * circle.R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(circle.R * cos(theta), circle.R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();

            if (elem_type == "TRI6") solid_mesh.all_second_order();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * circle.R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, circle.R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        // Translate the mesh to a new location.
        // This can also be achieved by libMesh's canned routine MeshTools::Modification::translate().
        for (MeshBase::node_iterator it = solid_mesh.nodes_begin(); it != solid_mesh.nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& x = *n;
            x(0) += shift_x;
            x(1) += shift_y;
#if (NDIM == 3)
            x(2) += shift_z;
#endif
        }
        solid_mesh.prepare_for_use();

        // Create boundary mesh
        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        Mesh& mesh = solid_mesh;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        SAMRAIPointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        // Set up the advection diffusion hierarchy integrator
        SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "PREDICTOR_CORRECTOR");
        if (adv_diff_solver_type == "PREDICTOR_CORRECTOR")
        {
            SAMRAIPointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps(
                "AdvectorExplicitPredictorPatchOps",
                app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
            adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator(
                "AdvDiffPredictorCorrectorHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"),
                predictor);
        }
        else if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: PREDICTOR_CORRECTOR, SEMI_IMPLICIT");
        }
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        SAMRAIPointer<IBFEMethod> ibfe_method_ops = nullptr;
        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (!from_restart)
        {
            ibfe_method_ops =
                new IBFEMethod("IBFEMethod",
                               app_initializer->getComponentDatabase("IBFEMethod"),
                               &mesh,
                               app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                               /*register_for_restart*/ false);
        }
        SAMRAIPointer<IBInterpolantMethod> ib_interpolant_method_ops = new IBInterpolantMethod(
            "IBInterpolantMethod", app_initializer->getComponentDatabase("IBInterpolantMethod"));
        SAMRAIPointer<IBLevelSetMethod> ib_level_set_method_ops =
            new IBLevelSetMethod(ib_interpolant_method_ops, ibfe_method_ops);

        SAMRAIPointer<IBHierarchyIntegrator> time_integrator = new IBInterpolantHierarchyIntegrator(
            "IBInterpolantHierarchyIntegrator",
            app_initializer->getComponentDatabase("IBInterpolantHierarchyIntegrator"),
            ib_level_set_method_ops,
            navier_stokes_integrator);

        SAMRAIPointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);

        SAMRAIPointer<StandardTagAndInitializeNd> error_detector =
            new StandardTagAndInitializeNd("StandardTagAndInitialize",
                                           time_integrator,
                                           app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        SAMRAIPointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Create level sets for solid interface.
        const string& ls_name_solid = "level_set_solid";
        SAMRAIPointer<CellVariableNd<double> > phi_var_solid = new CellVariableNd<double>(ls_name_solid);

        // Create level sets for gas/liquid interface.
        const double fluid_height = input_db->getDouble("GAS_LS_INIT");
        const string& ls_name_gas = "level_set_gas";
        SAMRAIPointer<CellVariableNd<double> > phi_var_gas = new CellVariableNd<double>(ls_name_gas);
        SAMRAIPointer<RelaxationLSMethod> level_set_gas_ops =
            new RelaxationLSMethod(ls_name_gas, app_initializer->getComponentDatabase("LevelSet_Gas"));
        LSLocateGasInterface setLSLocateGasInterface(
            "LSLocateGasInterface", adv_diff_integrator, phi_var_gas, fluid_height);
        level_set_gas_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateGasInterfaceCallbackFunction,
                                                                    static_cast<void*>(&setLSLocateGasInterface));

        // Register the level sets with advection diffusion integrator.
        adv_diff_integrator->registerTransportedQuantity(phi_var_solid);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_solid, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_solid,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        adv_diff_integrator->registerTransportedQuantity(phi_var_gas);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_gas, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_gas,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Register the reinitialization functions for the level set variables
        IBAMR::LevelSetUtilities::SetLSProperties setSetLSProperties("SetLSProperties", level_set_gas_ops);
        adv_diff_integrator->registerResetFunction(
            phi_var_gas, &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy, static_cast<void*>(&setSetLSProperties));
        SolidLevelSetResetter solid_level_set_resetter;
        solid_level_set_resetter.ib_interp_ops = ib_interpolant_method_ops;
        solid_level_set_resetter.adv_diff_integrator = adv_diff_integrator;
        solid_level_set_resetter.ls_solid_var = phi_var_solid;
        adv_diff_integrator->registerIntegrateHierarchyCallback(&reset_solid_level_set_callback_fcn,
                                                                static_cast<void*>(&solid_level_set_resetter));

        // Setup the advected and diffused fluid quantities.
        SAMRAIPointer<CellVariableNd<double> > mu_var = new CellVariableNd<double>("mu");
        SAMRAIPointer<hier::VariableNd> rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariableNd<double>("rho");
        }
        else
        {
            rho_var = new CellVariableNd<double>("rho");
        }
        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_fluid = input_db->getDouble("RHO_F");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const double mu_fluid = input_db->getDouble("MU_F");
        const double mu_gas = input_db->getDouble("MU_G");
        const double mu_solid = input_db->getDoubleWithDefault("MU_S", std::numeric_limits<double>::quiet_NaN());
        const bool set_mu_solid = input_db->getBool("SET_MU_S");
        const int num_solid_interface_cells = input_db->getDouble("NUM_SOLID_INTERFACE_CELLS");
        const int num_gas_interface_cells = input_db->getDouble("NUM_GAS_INTERFACE_CELLS");
        circle.rho_solid = rho_solid;
        IBAMR::VCINSUtilities::SetFluidProperties setSetFluidProperties("SetFluidProperties",
                                                                        adv_diff_integrator,
                                                                        phi_var_gas,
                                                                        phi_var_solid,
                                                                        rho_fluid,
                                                                        rho_gas,
                                                                        rho_solid,
                                                                        mu_fluid,
                                                                        mu_gas,
                                                                        mu_solid,
                                                                        num_gas_interface_cells,
                                                                        num_solid_interface_cells,
                                                                        set_mu_solid);
        navier_stokes_integrator->registerResetFluidDensityFcn(&IBAMR::VCINSUtilities::callSetDensityCallbackFunction,
                                                               static_cast<void*>(&setSetFluidProperties));
        navier_stokes_integrator->registerResetFluidViscosityFcn(
            &IBAMR::VCINSUtilities::callSetViscosityCallbackFunction, static_cast<void*>(&setSetFluidProperties));

        // Register callback function for tagging refined cells for level set data
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        const double tag_min_value = -tag_thresh;
        const double tag_max_value = tag_thresh;
        IBAMR::LevelSetUtilities::TagLSRefinementCells ls_tagger(
            adv_diff_integrator, phi_var_gas, tag_min_value, tag_max_value);
        time_integrator->registerApplyGradientDetectorCallback(&IBAMR::LevelSetUtilities::tagLSCells,
                                                               static_cast<void*>(&ls_tagger));

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVectorNd& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategyNd*> u_bc_coefs(NDIM);
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

        RobinBcCoefStrategyNd* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategyNd* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategyNd* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_gas, phi_bc_coef);
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef);

        // LS reinit boundary conditions, which is set to be the same as the BCs
        // for advection
        RobinBcCoefStrategyNd* ls_reinit_bcs = phi_bc_coef;
        level_set_gas_ops->registerPhysicalBoundaryCondition(ls_reinit_bcs);

        // Body forces.
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        circle.g_y = grav_const[1];
        SAMRAIPointer<CartGridFunction> grav_force =
            new IBAMR::VCINSUtilities::GravityForcing("GravityForcing",
                                                      adv_diff_integrator,
                                                      phi_var_gas,
                                                      app_initializer->getComponentDatabase("FlowGravityForcing"),
                                                      grav_const);

        SAMRAIPointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var_gas);

        SAMRAIPointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(grav_force);
        eul_forces->addFunction(surface_tension_force);
        time_integrator->registerBodyForceFunction(eul_forces);

        // Configure the IBFE solver.
        EquationSystems* equation_systems = nullptr;
        if (ibfe_method_ops)
        {
            ibfe_method_ops->initializeFEEquationSystems();
            equation_systems = ibfe_method_ops->getFEDataManager(/*part*/ 0)->getEquationSystems();
        }

        // Configure the IBInterpolant solver.
        SAMRAIPointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        std::vector<std::string> struct_list_vec(1, "InterpolationMesh");
        coarsest_ln = 0;
        max_finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
        ib_initializer->setStructureNamesOnLevel(max_finest_ln, struct_list_vec);
        ib_initializer->registerInitStructureFunction(generate_interp_mesh);
        ib_interpolant_method_ops->registerLInitStrategy(ib_initializer);
        ib_interpolant_method_ops->registerVariableAndHierarchyIntegrator(
            ls_name_solid, /*depth*/ 1, phi_var_solid, adv_diff_integrator);

        // Configure the Brinkman penalization object to do the rigid body dynamics.
        SAMRAIPointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd =
            new BrinkmanPenalizationRigidBodyDynamics("Brinkman Body",
                                                      phi_var_solid,
                                                      adv_diff_integrator,
                                                      navier_stokes_integrator,
                                                      app_initializer->getComponentDatabase("BrinkmanPenalization"),
                                                      /*register_for_restart*/ true);
        FreeRigidDOFVector free_dofs;
        free_dofs << 0, 1, 0;
        Eigen::Vector3d U_i = Eigen::Vector3d::Zero();
        const double mass = rho_solid * M_PI * std::pow(circle.R, 2);
        bp_rbd->setSolveRigidBodyVelocity(free_dofs);
        bp_rbd->registerKinematicsFunction(&imposed_kinematics);
        bp_rbd->registerExternalForceTorqueFunction(&external_force_torque);
        bp_rbd->setInitialConditions(circle.X0, U_i, U_i, mass);
        navier_stokes_integrator->registerBrinkmanPenalizationStrategy(bp_rbd);
        solid_level_set_resetter.bp_rbd = bp_rbd;

        // Set up visualization plot file writers.
        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        SAMRAIPointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_interpolant_method_ops->registerLSiloDataWriter(silo_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        if (ibfe_method_ops)
        {
            ibfe_method_ops->initializeFEData();
        }
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Compute distance function from FE mesh.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<CellVariableNd<double> > n_var = new CellVariableNd<double>("num_elements", 1);
        SAMRAIPointer<CellVariableNd<double> > d_var = new CellVariableNd<double>("distance", 1);
        const IntVectorNd no_width = 0;
        SAMRAIPointer<VariableContext> main_ctx = var_db->getContext("Main");
        const int n_idx = var_db->registerVariableAndContext(n_var, main_ctx, no_width);
        const int d_idx = var_db->registerVariableAndContext(d_var, main_ctx, no_width);
        int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(n_idx, 0.0);
            level->allocatePatchData(d_idx, 0.0);
        }
        HierarchyCellDataOpsRealNd<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, hier_finest_ln);
        hier_cc_data_ops.setToScalar(n_idx, 0.0);
        hier_cc_data_ops.setToScalar(d_idx, 5 * dx);
        visit_data_writer->registerPlotQuantity("num_elements", "SCALAR", n_idx);
        visit_data_writer->registerPlotQuantity("distance", "SCALAR", d_idx);

        const int gcw = input_db->getIntegerWithDefault("GCW", 1);
        FESurfaceDistanceEvaluator surface_distance_eval("FESurfaceDistanceEvaluator",
                                                         patch_hierarchy,
                                                         mesh,
                                                         boundary_mesh,
                                                         /*gcw*/ gcw,
                                                         use_vol_extracted_bdry_mesh);
        pout << "Started mapping intersections" << std::endl;
        surface_distance_eval.mapIntersections();
        pout << "Finished mapping intersections" << std::endl;
        pout << "Computing face normal" << std::endl;
        surface_distance_eval.calculateSurfaceNormals();
        pout << "Finished calculation of face normal" << std::endl;
        pout << "Computing distances" << std::endl;
        surface_distance_eval.computeSignedDistance(n_idx, d_idx);
        pout << "Finished computing distances" << std::endl;
        pout << "Updating sign" << std::endl;
        surface_distance_eval.updateSignAwayFromInterface(d_idx, patch_hierarchy, 5 * dx);
        pout << "Finished updating sign" << std::endl;

        // Compute the error.
        SAMRAIPointer<CellVariableNd<double> > E_var = new CellVariableNd<double>("E");
        const int E_idx = var_db->registerVariableAndContext(E_var, main_ctx);
        for (int ln = 0; ln <= hier_finest_ln; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(E_idx)) level->allocatePatchData(E_idx, time_integrator->getIntegratorTime());
        }

        SAMRAIPointer<HierarchyMathOps> hier_math_ops =
            new HierarchyMathOps("HierarchyMathOps", patch_hierarchy, coarsest_ln, hier_finest_ln);
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

        double E_interface = 0.0;
        int num_interface_pts = 0;
        double volume_near_interface = 0.0;
        calculate_distance_analytically(patch_hierarchy, E_idx);
        calculate_error_near_band(
            patch_hierarchy, E_idx, d_idx, wgt_cc_idx, E_interface, num_interface_pts, volume_near_interface);
        pout << "Error in D near interface after level set initialization:" << std::endl
             << "L1-norm:  " << std::setprecision(10) << E_interface / volume_near_interface << std::endl;
        pout << "Number of points within the interface (used to compute interface error):" << std::endl
             << num_interface_pts << std::endl;

        // Register for plotting
        visit_data_writer->registerPlotQuantity("Error", "SCALAR", E_idx);

        // Deallocate initialization objects.
        ib_interpolant_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        // Copy the distance function into level set of solid.
        int phi_solid_idx =
            var_db->mapVariableAndContextToIndex(phi_var_solid, adv_diff_integrator->getCurrentContext());
        hier_cc_data_ops.copyData(phi_solid_idx, d_idx);
        if (!RestartManager::getManager()->isFromRestart())
        {
            navier_stokes_integrator->initializeCompositeHierarchyData(loop_time, /*initial_time*/ true);
        }

        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (!from_restart && uses_exodus)
            {
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Deactivate IBFEMethod as we don't need it anymore.
        if (ibfe_method_ops)
        {
            ib_level_set_method_ops->deactivateIBFEMethod();
            ibfe_method_ops.setNull();
        }

        // Open streams to save position and velocity of the structure.
        ofstream rbd_stream;
        if (IBTK_MPI::getRank() == 0)
        {
            rbd_stream.open("rbd.curve", ios_base::out | ios_base::app);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            pout << "Advancing hierarchy with timestep size dt = " << dt << "\n";
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "Writing visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "Writing restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "Writing timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            if (IBTK_MPI::getRank() == 0)
            {
                const Eigen::Vector3d& rbd_posn = bp_rbd->getCurrentCOMPosn();
                const Eigen::Vector3d& rbd_trans_vel = bp_rbd->getCurrentCOMTransVelocity();

                rbd_stream.precision(12);
                rbd_stream.setf(ios::fixed, ios::floatfield);
                rbd_stream << loop_time << "\t" << rbd_posn[1] << "\t" << rbd_trans_vel[1] << std::endl;
            }
        }

        // Close the logging streams.
        if (IBTK_MPI::getRank() == 0)
        {
            rbd_stream.close();
        }

        // Delete dumb pointers.
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
