// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/BrinkmanPenalizationRigidBodyDynamics.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application specific includes.
#include "LevelSetInitialCondition.h"
#include "SetFluidSolidDensity.h"
#include "SetFluidSolidViscosity.h"
#include "TagLSRefinementCells.h"

// Struct to maintain the properties of the circular interface
struct CircularInterface
{
    Eigen::Vector2d X0;
    double R;
};

CircularInterface circle;
RectangularInterface rectangle;

// Struct to reset solid level set
struct SolidLevelSetResetter
{
    SolidLevelSetResetter(Pointer<AdvDiffHierarchyIntegrator> integrator,
                          Pointer<CellVariable<NDIM, double> > var,
                          Pointer<BrinkmanPenalizationRigidBodyDynamics> bp)
        : adv_diff_integrator(integrator), ls_solid_var(var), bp_rbd(bp)
    {
        return;
    }

    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
    Pointer<CellVariable<NDIM, double> > ls_solid_var;
    Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd;
};

void
reset_solid_level_set_callback_fcn(double current_time, double new_time, int /*cycle_num*/, void* ctx)
{
    SolidLevelSetResetter* resetter = static_cast<SolidLevelSetResetter*>(ctx);

    // Get the new centroid of the body
    const double dt = new_time - current_time;
    Eigen::Vector3d XCOM_current = resetter->bp_rbd->getCurrentCOMPosn();
    Eigen::Vector3d XCOM_new = XCOM_current + dt * (resetter->bp_rbd->getNewCOMTransVelocity());
    IBTK::VectorNd origin(XCOM_new[0], XCOM_new[1]);

    // Set a large value away from the solid body.
    Pointer<PatchHierarchy<NDIM> > patch_hier = resetter->adv_diff_integrator->getPatchHierarchy();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_solid_idx =
        var_db->mapVariableAndContextToIndex(resetter->ls_solid_var, resetter->adv_diff_integrator->getNewContext());

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > ls_solid_data = patch->getPatchData(ls_solid_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const hier::Index<NDIM>& ci = it();
                IBTK::VectorNd coord = IBTK::VectorNd::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                const IBTK::VectorNd p = coord - origin;
                const IBTK::VectorNd& b = rectangle.S;
                IBTK::VectorNd d = p.cwiseAbs() - b;
                const double mm = std::max(d[0], d[1]);
                d[0] = std::max(d[0], 0.0);
                d[1] = std::max(d[1], 0.0);
                const double distance = (d.norm() + std::min(0.0, mm));
                (*ls_solid_data)(ci) = distance;
            }
        }
    }

    return;
}

void
generate_interp_mesh(const unsigned int& /*strct_num*/,
                     const int& /*ln*/,
                     int& /*num_vertices*/,
                     std::vector<IBTK::Point>& /*vertex_posn*/)
{
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
    U_com[1] = 1.0;
    W_com.setZero();
    return;
} // imposed_kinematics

void
external_force_torque(double /*data_time*/, int /*cycle_num*/, Eigen::Vector3d& F, Eigen::Vector3d& T, void* /*ctx*/)
{
    F.setZero();
    T.setZero();
    return;
} // imposed_kinematics

// IBFE Elasticity model data.
namespace ModelData
{
// Problem parameters.
static const double mu = 10.0; // Elastic modulus

// Stress tensor functions.
void
PK1_dev_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const std::vector<const std::vector<double>*>& /*var_data*/,
                        const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = mu * FF;
    return;
} // PK1_dev_stress_function

void
PK1_dil_stress_function(TensorValue<double>& PP,
                        const TensorValue<double>& FF,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* const /*elem*/,
                        const std::vector<const std::vector<double>*>& /*var_data*/,
                        const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                        double /*time*/,
                        void* /*ctx*/)
{
    PP = -mu * tensor_inverse_transpose(FF, NDIM);
    return;
} // PK1_dil_stress_function
} // namespace ModelData
using namespace ModelData;

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IBLevelSet.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

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

        // Setup disk information
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM_DISK");
        circle.X0[1] = input_db->getDouble("YCOM_DISK");

        // Set up rectangle information
        rectangle.S[0] = 0.5 * input_db->getDouble("WIDTH_RECTANGLE");
        rectangle.S[1] = 0.5 * input_db->getDouble("HEIGHT_RECTANGLE");
        rectangle.X0[0] = input_db->getDouble("XCOM_RECTANGLE");
        rectangle.X0[1] = input_db->getDouble("YCOM_RECTANGLE");

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);

        // Create mesh based upon input file
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");

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
            x(0) += circle.X0[0];
            x(1) += circle.X0[1];
        }
        solid_mesh.prepare_for_use();
        Mesh& mesh = solid_mesh;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator =
            new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));

        // Set up the advection diffusion hierarchy integrator
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator",
            app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        Pointer<IBFEMethod> ibfe_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           /*register_for_restart*/ true);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ibfe_method_ops,
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

        // Create level sets for solid interface.
        const string& ls_name_solid = "level_set_solid";
        Pointer<CellVariable<NDIM, double> > phi_var_solid = new CellVariable<NDIM, double>(ls_name_solid);

        // Register the level sets with advection diffusion integrator.
        adv_diff_integrator->registerTransportedQuantity(phi_var_solid);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_solid, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_solid,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Solid level set initial condition
        Pointer<CartGridFunction> phi_solid_init = new LevelSetInitialCondition("phi_solid_init", rectangle);
        adv_diff_integrator->setInitialConditions(phi_var_solid, phi_solid_init);

        SolidLevelSetResetter solid_level_set_resetter(adv_diff_integrator, phi_var_solid, /*bp_rbd*/ nullptr);
        adv_diff_integrator->registerIntegrateHierarchyCallback(&reset_solid_level_set_callback_fcn,
                                                                static_cast<void*>(&solid_level_set_resetter));

        // Setup the advected and diffused fluid quantities.
        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        Pointer<hier::Variable<NDIM> > rho_var = new SideVariable<NDIM, double>("rho");
        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_fluid = input_db->getDouble("RHO_F");
        SetFluidSolidDensity* ptr_setFluidSolidDensity = new SetFluidSolidDensity("SetFluidSolidDensity", rho_fluid);
        navier_stokes_integrator->registerResetFluidDensityFcn(&callSetFluidSolidDensityCallbackFunction,
                                                               static_cast<void*>(ptr_setFluidSolidDensity));

        const double mu_fluid = input_db->getDouble("MU_F");
        SetFluidSolidViscosity* ptr_setFluidSolidViscosity =
            new SetFluidSolidViscosity("SetFluidSolidViscosity", mu_fluid);
        navier_stokes_integrator->registerResetFluidViscosityFcn(&callSetFluidSolidViscosityCallbackFunction,
                                                                 static_cast<void*>(ptr_setFluidSolidViscosity));

        // Register callback function for tagging refined cells for level set data
        const double tag_value = input_db->getDouble("LS_TAG_VALUE");
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        TagLSRefinementCells ls_tagger(adv_diff_integrator, phi_var_solid, tag_value, tag_thresh);
        time_integrator->registerApplyGradientDetectorCallback(&callTagSolidLSRefinementCellsCallbackFunction,
                                                               static_cast<void*>(&ls_tagger));

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

        RobinBcCoefStrategy<NDIM>* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategy<NDIM>* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef);

        // Configure the IBFE solver.
        ibfe_method_ops->initializeFEEquationSystems();
        EquationSystems* equation_systems = ibfe_method_ops->getFEDataManager(/*part*/ 0)->getEquationSystems();
        IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
        IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
        PK1_dev_stress_data.quad_order =
            Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        PK1_dil_stress_data.quad_order =
            Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "FIRST"));
        ibfe_method_ops->registerPK1StressFunction(PK1_dev_stress_data);
        ibfe_method_ops->registerPK1StressFunction(PK1_dil_stress_data);

        // Configure the Brinkman penalization object to do the rigid body dynamics of the rectangle.
        Pointer<BrinkmanPenalizationRigidBodyDynamics> bp_rbd =
            new BrinkmanPenalizationRigidBodyDynamics("Brinkman Body",
                                                      phi_var_solid,
                                                      adv_diff_integrator,
                                                      navier_stokes_integrator,
                                                      app_initializer->getComponentDatabase("BrinkmanPenalization"),
                                                      /*register_for_restart*/ true);
        FreeRigidDOFVector free_dofs;
        free_dofs << 0, 0, 0;
        Eigen::Vector3d U_i = Eigen::Vector3d::Zero();
        const double mass = 4 * rho_fluid * rectangle.S[0] * rectangle.S[1];
        bp_rbd->setSolveRigidBodyVelocity(free_dofs);
        bp_rbd->registerKinematicsFunction(&imposed_kinematics);
        bp_rbd->registerExternalForceTorqueFunction(&external_force_torque);
        bp_rbd->setInitialConditions(Eigen::Vector3d(rectangle.X0[0], rectangle.X0[1], 0.0), U_i, U_i, mass);
        navier_stokes_integrator->registerBrinkmanPenalizationStrategy(bp_rbd);
        solid_level_set_resetter.bp_rbd = bp_rbd;

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ibfe_method_ops->initializeFEData();
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

        // Open streams to save position and velocity of the structure.
        ofstream rbd_stream;
        if (IBTK_MPI::getRank() == 0)
        {
            rbd_stream.open("rbd.curve", ios_base::out | ios_base::app);
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

                if (uses_exodus)
                {
                    exodus_io->write_timestep(
                        exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
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
        delete ptr_setFluidSolidDensity;
        delete ptr_setFluidSolidViscosity;
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
