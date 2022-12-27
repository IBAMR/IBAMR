// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/FEMechanicsExplicitIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesOpenBoundaryStabilizer.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/MathUtilities.h>
#include <tbox/Utilities.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/dense_vector.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/enum_solver_package.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/fe.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/sparse_matrix.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Application includes
#include "FeedbackForcer.h"
#include "VelocityBcCoefs.h"

// Elasticity model data.
namespace ModelData
{
static const unsigned int HOUSING_PART = 0;
static const unsigned int BEAM_PART = 1;

static const double g = -980.665; // cm/s^2
static double t_load = 0.5;
static double tg_load = 0.5;
static double kappa_FSI_housing = 0.0;
static double eta_FSI_housing = 0.0;
static double kappa_FSI_beam = 0.0;
static double eta_damping = 0.0;
static double eta_FSI_beam = 0.0;
static double rho_s = 0.0;
static double rho_f = 0.0;
static double lambda_s = 0.0;
static double mu_s = 0.0;
static double bulk_mod = 0.0;
static double shear_mod = 0.0;
static double dx = 0.0;

System* x_new_beam_system;
System* u_new_beam_system;
System* x_new_beam_surface_system;
System* Tau_new_beam_surface_system;

EquationSystems* bndry_beam_G_systems;

static const double TOL = sqrt(std::numeric_limits<double>::epsilon());
static BoundaryInfo* vol_beam_bndry_info;

// function prototype
double time_ramp(double t);

// housing functions
// Tether (penalty) stress function.

void
tether_FSI_force_function_housing(VectorValue<double>& F,
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
    // tether_force_function() is called on elements of the boundary mesh.  Here
    // we look up the element in the solid mesh that the current boundary
    // element was extracted from.

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.

    const std::vector<double>& U = *var_data[0];

    double u_bndry_n = 0.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        u_bndry_n += n(d) * U[d];
    }

    // The tether force is proportional to the mismatch between the positions
    // and velocities.

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_FSI_housing * (X(d) - x(d)) + eta_FSI_housing * (0.0 - u_bndry_n) * n(d);
    }

    return;
} // FSI_tether_force_function

// beam functions
// Tether (penalty) stress function.
void
PK1_dev_stress_function_beam(TensorValue<double>& PP,
                             const TensorValue<double>& FF,
                             const libMesh::Point& /*x*/,
                             const libMesh::Point& /*X*/,
                             Elem* const /*elem*/,
                             const vector<const vector<double>*>& /*var_data*/,
                             const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    // NH material model with modified invariants and volumetric term
    double J = FF.det();
    double I1 = (FF.transpose() * FF).tr();
    TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> EE = 0.5 * (CC - II);
    //  Unmodified St. Venant-Kirchhoff Model
#if 0
	PP = FF * (2 * mu_s * EE + lambda_s * EE.tr() * II);
#endif

// Neo-Hookean Model
#if 0
    PP = mu_s * (FF - FF_inv_trans);
#endif

//  Modified Neo-Hookean Model
#if 1
    PP = shear_mod * pow(J, -2.0 / 3.0) * (FF - (I1 / 3.0) * FF_inv_trans);
#endif

    return;
} // PK1_dev_stress_function_beam

void
PK1_dil_stress_function_beam(TensorValue<double>& PP,
                             const TensorValue<double>& FF,
                             const libMesh::Point& /*X*/,
                             const libMesh::Point& /*s*/,
                             Elem* const /*elem*/,
                             const std::vector<const std::vector<double>*>& /*var_data*/,
                             const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
                             double /*time*/,
                             void* /*ctx*/)
{
    double J = FF.det();
    TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);
    PP = bulk_mod * J * log(J) * FF_inv_trans;
}

// Tether (penalty) force functions.
void
tether_force_function_beam(VectorValue<double>& F,
                           const VectorValue<double>& /*n*/,
                           const VectorValue<double>& /*N*/,
                           const TensorValue<double>& /*FF*/,
                           const libMesh::Point& /*x*/,
                           const libMesh::Point& X,
                           Elem* const elem,
                           const unsigned short side,
                           const vector<const vector<double>*>& /*var_data*/,
                           const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                           double /*time*/,
                           void* /*ctx*/)
{
    //~ const BeamData* const beam_data = reinterpret_cast<BeamData*>(ctx);

    MeshBase& mesh_bndry = bndry_beam_G_systems->get_mesh();

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const auto el_begin = mesh_bndry.active_local_elements_begin();
        const auto el_end = mesh_bndry.active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* elem_bndry = *el_it;

            if ((elem_bndry->contains_point(X)) && !((vol_beam_bndry_info->has_boundary_id(elem, side, 0))))
            {
                F(d) = Tau_new_beam_surface_system->point_value(d, X, elem_bndry); //&side_elem);
            }
        }
    }

    return;
} // tether_force_function_beam

void
tether_FSI_force_function_beam(VectorValue<double>& F,
                               const VectorValue<double>& /*n*/,
                               const VectorValue<double>& /*N*/,
                               const TensorValue<double>& /*FF*/,
                               const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the boundary
                                                              // mesh
                               const libMesh::Point& X_bndry, // X_bndry gives reference coordinates on the boundary
                                                              // mesh
                               Elem* const elem,
                               const unsigned short /*side*/,
                               const vector<const vector<double>*>& var_data,
                               const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                               double /*time*/,
                               void* /*ctx*/)
{
    // tether_force_function() is called on elements of the boundary mesh.  Here
    // we look up the element in the solid mesh that the current boundary
    // element was extracted from.
    const Elem* const interior_parent = elem->interior_parent();
    const libMesh::Point cp_elem = elem->centroid();

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.
    std::vector<double> x_solid(NDIM, 0.0);
    std::vector<double> u_solid(NDIM, 0.0);

    double disp = 0.0;

    const std::vector<double>& U = *var_data[0];

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_beam_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_beam_system->point_value(d, X_bndry, interior_parent);
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        disp += (x_solid[d] - x_bndry(d)) * (x_solid[d] - x_bndry(d));
    }
    disp = sqrt(disp);
    TBOX_ASSERT(disp < 4.0 * dx);

    // The tether force is proportional to the mismatch between the positions
    // and velocities.

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_FSI_beam * (x_solid[d] - x_bndry(d)) + eta_FSI_beam * (u_solid[d] - U[d]);
    }

    return;
} // FSI_tether_force_function

// Tether (penalty) force functions.
void
body_force_function_beam(VectorValue<double>& F,
                         const TensorValue<double>& /*FF*/,
                         const libMesh::Point& /*x*/,
                         const libMesh::Point& /*X*/,
                         Elem* const /*elem*/,
                         const vector<const vector<double>*>& var_data,
                         const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                         double time,
                         void* /*ctx*/)
{
    const std::vector<double>& U = *var_data[0];
    F.zero();
    //~ const BeamData* const beam_data = reinterpret_cast<BeamData*>(ctx);
    for (int d = 0; d < NDIM; ++d) F(d) = -eta_damping * U[d];
    F(1) += time_ramp(time) * g * (rho_s - rho_f);
    return;
} //

double
time_ramp(double t)
{
    return (t < tg_load) ? -16.0 * t * t * t + 12.0 * t * t : 1.0;
}

} // namespace ModelData
using namespace ModelData;

// Function prototypes
static ofstream x_stream, y_stream, z_stream;
void record_position(const FEMechanicsExplicitIntegrator* fem_solver,
                     EquationSystems* beam_systems,
                     vector<libMesh::Point> evaluation_point,
                     const double loop_time);

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

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        tbox::Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        tbox::Pointer<tbox::Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();

        const string viz_dump_dirname = app_initializer->getVizDumpDirectory();
        const string beam_filename = viz_dump_dirname + "/beam.ex2";
        const string bndry_beam_filename = viz_dump_dirname + "/bndry_beam.ex2";
        const string bndry_housing_filename = viz_dump_dirname + "/bndry_housing.ex2";

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

        string housing_elem_type = input_db->getStringWithDefault("HOUSING_ELEM_TYPE", "HEX8");
        string beam_elem_type = input_db->getString("BEAM_ELEM_TYPE");
        double mfac_b = input_db->getDouble("MFAC");
        const double n_cycles = input_db->getDouble("NCYCLE");
        if (beam_elem_type == "TET10" || beam_elem_type == "HEX27") mfac_b *= 2.0;
        dx = input_db->getDouble("DX");
        const double ds_b = mfac_b * dx;

        lambda_s = input_db->getDouble("LAMBDA_S");

        mu_s = input_db->getDouble("MU_S");
        bulk_mod = input_db->getDouble("BULK_MOD");
        shear_mod = input_db->getDouble("SHEAR_MOD");

        rho_s = input_db->getDouble("RHO_S");

        rho_f = input_db->getDouble("RHO_F");

        kappa_FSI_housing = input_db->getDouble("KAPPA_FSI_HOUSING");

        kappa_FSI_beam = input_db->getDouble("KAPPA_FSI_BEAM");

        eta_FSI_beam = input_db->getDouble("ETA_FSI_BEAM");

        eta_damping = input_db->getDouble("ETA_DAMPING");

        eta_FSI_housing = input_db->getDouble("ETA_FSI_HOUSING");

        double width = 1.1;
        double height = 0.2;
        double length = 6.5;
        double offset = 0.10245;

        Mesh beam_mesh(init.comm(), NDIM);

        MeshTools::Generation::build_cube(beam_mesh,
                                          static_cast<int>(ceil(width / ds_b)),
                                          static_cast<int>(ceil(height / ds_b)),
                                          static_cast<int>(ceil(length / ds_b)),
                                          -0.5 * width,
                                          0.5 * width,
                                          -0.5 * height,
                                          0.5 * height,
                                          offset,
                                          length + offset,
                                          Utility::string_to_enum<ElemType>(beam_elem_type));

        const auto end_el = beam_mesh.elements_end();
        for (auto el = beam_mesh.elements_begin(); el != end_el; ++el)
        {
            const Elem* elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    BoundaryInfo& boundary_info = beam_mesh.get_boundary_info();
                    if (boundary_info.has_boundary_id(elem, side, 0))
                    {
                        boundary_info.add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID);
                    }
                }
            }
        }
        beam_mesh.prepare_for_use();
        vol_beam_bndry_info = &beam_mesh.get_boundary_info();

        Mesh housing_bndry_mesh(beam_mesh.comm(), NDIM - 1);
        housing_bndry_mesh.read(input_db->getString("HOUSING_MESH_FILENAME"));

        MeshRefinement mesh_refinement_housing(housing_bndry_mesh);

        // to scale mesh from units in mm to cm
        using MeshTools::Modification::scale;
        scale(housing_bndry_mesh, 0.10245, 0.10245, 0.10245);

        BoundaryMesh beam_bndry_mesh(beam_mesh.comm(), beam_mesh.mesh_dimension() - 1);

        beam_mesh.get_boundary_info().sync(beam_bndry_mesh);
        beam_bndry_mesh.prepare_for_use();
        // to scale mesh from units in mm to cm

        vector<MeshBase*> bndry_meshes(2);
        bndry_meshes[HOUSING_PART] = &housing_bndry_mesh;
        bndry_meshes[BEAM_PART] = &beam_bndry_mesh;

        // points to calculate position of beam
        double incr = 0.05;
        int num_pts = static_cast<int>(length / incr) + 1;
        vector<libMesh::Point> points(num_pts);
        for (int i = 0; i < num_pts; i++) points[i] = libMesh::Point(0.0, 0.0, offset + i * incr);

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        tbox::Pointer<IIMethod> ib_method_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         bndry_meshes,
                         app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        tbox::Pointer<FEMechanicsExplicitIntegrator> fem_solver = new FEMechanicsExplicitIntegrator(
            "FEMechanicsExplicitIntegrator",
            app_initializer->getComponentDatabase("FEMechanicsExplicitIntegrator"),
            &beam_mesh,
            app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        tbox::Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        tbox::Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        tbox::Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        tbox::Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        tbox::Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        tbox::Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IBFE solver.

        BcData bc_data(input_db->getDatabase("BcCoefs"));
        t_load = bc_data.t_load;
        tg_load = bc_data.tg_load;
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));

        vector<SystemData> velocity_data(1);
        velocity_data[0] = SystemData(fem_solver->getVelocitySystemName(), vars);

        const bool USE_DISCON_ELEMS = input_db->getBool("USE_DISCON_ELEMS");
        const bool USE_NORMALIZED_PRESSURE_JUMP = input_db->getBool("USE_NORMALIZED_PRESSURE_JUMP");

        if (USE_DISCON_ELEMS) ib_method_ops->registerDisconElemFamilyForJumps(BEAM_PART);
        if (USE_NORMALIZED_PRESSURE_JUMP) ib_method_ops->registerPressureJumpNormalization(BEAM_PART);

        ib_method_ops->initializeFEEquationSystems();

        IIMethod::LagSurfaceForceFcnData surface_FSI_fcn_data_beam(tether_FSI_force_function_beam, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_FSI_fcn_data_beam, BEAM_PART);

        IIMethod::LagSurfaceForceFcnData surface_FSI_fcn_data_housing(tether_FSI_force_function_housing, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_FSI_fcn_data_housing, HOUSING_PART);

        EquationSystems* bndry_housing_systems = ib_method_ops->getFEDataManager(HOUSING_PART)->getEquationSystems();
        EquationSystems* bndry_beam_systems = ib_method_ops->getFEDataManager(BEAM_PART)->getEquationSystems();

        // beam stress and surface/body forces

        FEMechanicsBase::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function_beam, velocity_data);
        PK1_dev_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "THIRD"));
        fem_solver->registerPK1StressFunction(PK1_dev_stress_data);

        FEMechanicsBase::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function_beam, velocity_data);
        PK1_dil_stress_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER", "CONSTANT"));
        fem_solver->registerPK1StressFunction(PK1_dil_stress_data);

        FEMechanicsBase::LagSurfaceForceFcnData surface_tether_force_data(tether_force_function_beam, velocity_data);
        fem_solver->registerLagSurfaceForceFunction(surface_tether_force_data);

        fem_solver->initializeFEEquationSystems();
        EquationSystems* beam_systems = fem_solver->getEquationSystems();

        // buoyancy force due to Boussinesq-like approximation + body damping
        FEMechanicsBase::LagBodyForceFcnData body_fcn_data_beam(body_force_function_beam, velocity_data);
        //~ body_fcn_data_beam.ctx = beam_data_ptr;
        fem_solver->registerLagBodyForceFunction(body_fcn_data_beam);

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            tbox::Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            tbox::Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            tbox::Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // register bc's
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (int d = 0; d < NDIM; ++d)
            u_bc_coefs[d] = new VelocityBcCoefs(navier_stokes_integrator, bc_data, d); // VelocityBcCoefs extends
                                                                                       // RobinBcCoefStrategy
        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        tbox::Pointer<FeedbackForcer> feedback_forcer =
            new FeedbackForcer(navier_stokes_integrator, patch_hierarchy, bc_data);
        time_integrator->registerBodyForceFunction(feedback_forcer);

        if (input_db->keyExists("BoundaryStabilization"))
        {
            time_integrator->registerBodyForceFunction(new StaggeredStokesOpenBoundaryStabilizer(
                "BoundaryStabilization",
                app_initializer->getComponentDatabase("BoundaryStabilization"),
                navier_stokes_integrator,
                grid_geometry));
        }

        // Set up visualization plot file writers.
        tbox::Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        std::unique_ptr<ExodusII_IO> beam_io(uses_exodus ? new ExodusII_IO(beam_mesh) : NULL);

        std::unique_ptr<ExodusII_IO> bndry_housing_io(uses_exodus ? new ExodusII_IO(*bndry_meshes[HOUSING_PART]) :
                                                                    NULL);
        std::unique_ptr<ExodusII_IO> bndry_beam_io(uses_exodus ? new ExodusII_IO(*bndry_meshes[BEAM_PART]) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        fem_solver->initializeFEData();

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
                beam_io->write_timestep(beam_filename, *beam_systems, iteration_num / viz_dump_interval + 1, loop_time);

                bndry_housing_io->write_timestep(
                    bndry_housing_filename, *bndry_housing_systems, iteration_num / viz_dump_interval + 1, loop_time);

                bndry_beam_io->write_timestep(
                    bndry_beam_filename, *bndry_beam_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // Open streams to save lift and drag coefficients and position of point A on beam
        if (SAMRAI_MPI::getRank() == 0)
        {
            x_stream.open(postproc_data_dump_dirname + "/x_posn.txt", ios_base::out | ios_base::trunc);
            y_stream.open(postproc_data_dump_dirname + "/y_posn.txt", ios_base::out | ios_base::trunc);
            z_stream.open(postproc_data_dump_dirname + "/z_posn.txt", ios_base::out | ios_base::trunc);
            x_stream.precision(10);
            y_stream.precision(10);
            z_stream.precision(10);
        }

        std::cout << "About to start time loop... \n";

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

            bndry_beam_G_systems = bndry_beam_systems;

            Tau_new_beam_surface_system = &bndry_beam_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_beam_surface_system = &bndry_beam_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            x_new_beam_system = &beam_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());

            u_new_beam_system = &beam_systems->get_system<System>(fem_solver->getVelocitySystemName());

            for (int ii = 0; ii < static_cast<int>(n_cycles); ii++)
            {
                fem_solver->preprocessIntegrateData(loop_time + (0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                    loop_time + (0.5 * static_cast<double>(ii + 1)) * dt / n_cycles,
                                                    /*num_cycles*/ 1);
                fem_solver->modifiedTrapezoidalStep(loop_time + (0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                    loop_time + (0.5 * static_cast<double>(ii + 1)) * dt / n_cycles);
                fem_solver->postprocessIntegrateData(loop_time + (0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                     loop_time + (0.5 * static_cast<double>(ii + 1)) * dt / n_cycles,
                                                     /*num_cycles*/ 1);
            }

            time_integrator->advanceHierarchy(dt);

            bndry_beam_G_systems = bndry_beam_systems;

            Tau_new_beam_surface_system = &bndry_beam_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_beam_surface_system = &bndry_beam_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            for (int ii = 0; ii < static_cast<int>(n_cycles); ii++)
            {
                fem_solver->preprocessIntegrateData(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                    loop_time +
                                                        (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles,
                                                    /*num_cycles*/ 1);
                fem_solver->modifiedTrapezoidalStep(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                    loop_time +
                                                        (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles);
                fem_solver->postprocessIntegrateData(loop_time + (0.5 + 0.5 * static_cast<double>(ii)) * dt / n_cycles,
                                                     loop_time +
                                                         (0.5 + 0.5 * static_cast<double>(ii + 1)) * dt / n_cycles,
                                                     /*num_cycles*/ 1);
            }

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
                    beam_io->write_timestep(
                        beam_filename, *beam_systems, iteration_num / viz_dump_interval + 1, loop_time);

                    bndry_housing_io->write_timestep(bndry_housing_filename,
                                                     *bndry_housing_systems,
                                                     iteration_num / viz_dump_interval + 1,
                                                     loop_time);

                    bndry_beam_io->write_timestep(
                        bndry_beam_filename, *bndry_beam_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }

                record_position(fem_solver, beam_systems, points, loop_time);
            }

            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

        } // end of time loop

        x_new_beam_system->clear();
        u_new_beam_system->clear();
        x_new_beam_surface_system->clear();
        Tau_new_beam_surface_system->clear();

        bndry_beam_systems->clear();
        beam_systems->clear();
        bndry_housing_systems->clear();

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            x_stream.close();
            y_stream.close();
            z_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
} // main

void
record_position(const FEMechanicsExplicitIntegrator* const fem_solver,
                EquationSystems* beam_systems,
                vector<libMesh::Point> evaluation_points,
                const double loop_time)
{
    System& x_system = beam_systems->get_system(fem_solver->getCurrentCoordinatesSystemName());
    DofMap& x_dof_map = x_system.get_dof_map();
    std::vector<unsigned int> vars;
    x_system.get_all_variable_numbers(vars);

    NumericVector<double>* x_vec = x_system.solution.get();
    std::unique_ptr<NumericVector<Number> > x_serial_vec = NumericVector<Number>::build(x_vec->comm());
    x_serial_vec->init(x_vec->size(), true, SERIAL);
    x_vec->localize(*x_serial_vec);

    MeshFunction mesh_fcn(*beam_systems, *x_serial_vec, x_dof_map, vars, 0);
    mesh_fcn.init();
    vector<DenseVector<Number> > x(evaluation_points.size());
    for (unsigned int i = 0; i < evaluation_points.size(); i++) mesh_fcn(evaluation_points[i], loop_time, x[i]);

    if (SAMRAI_MPI::getRank() == 0)
    {
        x_stream << loop_time;
        y_stream << loop_time;
        z_stream << loop_time;
        for (unsigned int i = 0; i < evaluation_points.size(); i++)
        {
            x_stream << " " << x[i](0);
            y_stream << " " << x[i](1);
            z_stream << " " << x[i](2);
        }
        x_stream << "\n";
        y_stream << "\n";
        z_stream << "\n";
    }

    return;
}
