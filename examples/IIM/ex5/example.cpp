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

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/MathUtilities.h>
#include <tbox/Utilities.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

#include <boost/multi_array.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double grav_const[3] = { 0.0, -981, 0.0 };
static const double TOL = sqrt(std::numeric_limits<double>::epsilon());
static VectorValue<double> COM;

System *x_new_solid_system, *u_new_solid_system;
System *x_half_solid_system, *u_half_solid_system;
void
tether_force_function(VectorValue<double>& F,
                      const VectorValue<double>& /*n*/,
                      const VectorValue<double>& /*N*/,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the boundary mesh
                      const libMesh::Point& X_bndry, // X_bndry gives reference coordinates on the boundary mesh
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

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.
    std::vector<double> x_solid(NDIM, 0.0), u_solid(NDIM, 0.0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_solid_system->point_value(d, X_bndry, interior_parent);
    }
    // Look up the velocity of the boundary mesh.
    const std::vector<double>& u_bndry = *var_data[0];

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s * (x_solid[d] - x_bndry(d)) + eta_s * (u_solid[d] - u_bndry[d]);
    }
    return;
} // tether_force_function
} // namespace ModelData
using namespace ModelData;

static ofstream w_new_stream, v_new_stream, x_com_stream, c1_traction_stream, c2_traction_stream, c3_traction_stream;
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

void
calculateGeomQuantitiesOfStructure(const double* params,
                                   double& M,                  // mass of the body
                                   TensorValue<double>& I_w,   // moment of inertia tensor
                                   VectorValue<double>& x_com, // new center of the mass
                                   EquationSystems* solid_equation_systems)
{
    // Get the structure mesh for codim-0 solid.
    // For now the eqs are setup only for one part but this will be extended
    // to multiple parts

    // Extract the FE system and DOF map, and setup the FE object.

    System& X_system = solid_equation_systems->get_system("position_new");
    X_system.solution->localize(*X_system.current_local_solution);
    MeshBase& mesh = solid_equation_systems->get_mesh();

    const double rho_s = params[0];
    const unsigned int dim = mesh.mesh_dimension();

    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);

    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType fe_type = X_dof_map.variable_type(0);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    PetscVector<double>& X_petsc = dynamic_cast<PetscVector<double>&>(*X_system.current_local_solution.get());
    X_petsc.close();
    Vec X_global_vec = X_petsc.vec();
    Vec X_local_ghost_vec;
    VecGhostGetLocalForm(X_global_vec, &X_local_ghost_vec);
    double* X_local_ghost_soln;
    VecGetArray(X_local_ghost_vec, &X_local_ghost_soln);

    // 3D identity tensor.
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    M = 0.0;
    I_w.zero();
    x_com.zero();
    double vol = 0.0; // volume of the body
    boost::multi_array<double, 2> X_node;

    VectorValue<double> X_qp, R_qp;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi);

            x_com += X_qp * JxW[qp];

            vol += JxW[qp];
        }
    }
    SAMRAI_MPI::sumReduction(&x_com(0), NDIM);
    SAMRAI_MPI::sumReduction(&vol);

    x_com /= vol;

    M = rho_s * vol;

    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }
        get_values_for_interpolation(X_node, X_petsc, X_local_ghost_soln, X_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp, qp, X_node, phi);

            R_qp = X_qp - x_com;

            // Accumulate the inertia tensor:
            I_w += rho_s * ((R_qp * R_qp) * II - outer_product(R_qp, R_qp)) * JxW[qp];
        }
    }

    SAMRAI_MPI::sumReduction(&I_w(0, 0), NDIM * NDIM);

    VecRestoreArray(X_local_ghost_vec, &X_local_ghost_soln);
    VecGhostRestoreLocalForm(X_global_vec, &X_local_ghost_vec);

    X_system.solution->close();

    return;

} // calculateGeomQuantitiesOfStructure

void
calculateFluidForceAndTorque(VectorValue<double>& F, // net force  acting on the body
                             VectorValue<double>& T, // net torque acting on the body
                             VectorValue<double> x_com_current,
                             const double /*loop_time*/,
                             EquationSystems* equation_systems)
{
    MeshBase& mesh = equation_systems->get_mesh();

    const unsigned int dim = mesh.mesh_dimension();
    F.zero();
    T.zero();

    System& x_system = equation_systems->get_system(IIMethod::COORDS_SYSTEM_NAME);
    System& TAU_system = equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);

    NumericVector<double>* TAU_vec = TAU_system.solution.get();
    NumericVector<double>* TAU_ghost_vec = TAU_system.current_local_solution.get();
    TAU_vec->localize(*TAU_ghost_vec);
    DofMap& TAU_dof_map = TAU_system.get_dof_map();
    std::vector<std::vector<unsigned int> > TAU_dof_indices(NDIM);
    std::unique_ptr<FEBase> fe(FEBase::build(dim, TAU_dof_map.variable_type(0)));

    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);

    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double> >& phi = fe->get_phi();

    boost::multi_array<double, 2> x_node, TAU_node;
    VectorValue<double> F_qp, x_qp, W_qp, TAU_qp, R_qp;

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
            TAU_dof_map.dof_indices(elem, TAU_dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(TAU_node, *TAU_ghost_vec, TAU_dof_indices);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            interpolate(TAU_qp, qp, TAU_node, phi);

            R_qp = x_qp - x_com_current;

            F += TAU_qp * JxW[qp];
            T += R_qp.cross(TAU_qp) * JxW[qp];
        }
    }
    SAMRAI_MPI::sumReduction(&F(0), NDIM);
    SAMRAI_MPI::sumReduction(&T(0), NDIM);

    x_ghost_vec->close();
    TAU_ghost_vec->close();

    return;

} // calculateFluidForceAndTorque

void
calculateGravitationalForce(const double* params, VectorValue<double>& F_g, EquationSystems* solid_equation_systems)
{
    MeshBase& mesh = solid_equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    std::unique_ptr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);

    const double rho_s = params[0];
    const double rho_f = params[1];
    // Extract the FE system and DOF map, and setup the FE object.
    System& Y_system = solid_equation_systems->get_system("position_new");

    Y_system.solution->localize(*Y_system.current_local_solution);
    DofMap& Y_dof_map = Y_system.get_dof_map();
    std::vector<std::vector<unsigned int> > Y_dof_indices(NDIM);

    FEType fe_type = Y_dof_map.variable_type(0);

    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

    // Extract the FE system and DOF map, and setup the FE object.
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();

    // Zero out the F_g force.
    F_g.zero();

    // Loop over the local elements to compute the local integrals.

    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            for (int d = 0; d < NDIM; ++d) F_g(d) += (rho_s - rho_f) * grav_const[d] * JxW[qp];
        }
    }
    SAMRAI_MPI::sumReduction(&F_g(0), NDIM);

    Y_system.solution->close();

    return;
} // calculateGravitationalForce

void
getSkewSymmetricAngVelTensor(TensorValue<double>& Omega, VectorValue<double> W)
{
    TBOX_ASSERT(NDIM == 3); // The code is currently setup only for 3D cases //

    Omega.zero();

    Omega(0, 1) = -W(2);
    Omega(0, 2) = W(1);
    Omega(1, 0) = W(2);
    Omega(1, 2) = -W(0);
    Omega(2, 0) = -W(1);
    Omega(2, 1) = W(0);

    return;
}

void
Solve6DOFSystemofEquations(const double* params,
                           double dt,
                           VectorValue<double>& V_new,     // linear velocity of the body
                           VectorValue<double>& W_new,     // angular velocity of the body
                           VectorValue<double>& x_com_new, // New position of the body
                           TensorValue<double>& Q_new,     // Rotation Matrix
                           double M,                       // Calculated Mass
                           TensorValue<double>& I_w_new,   // moment of inertia tensor
                           TensorValue<double> I_w_0,      // initial moment of inertia tensor
                           VectorValue<double> F_b,        // total external body force
                           VectorValue<double> F_s,        // total external surface force
                           VectorValue<double> T)          // Torque applied on the surface
{
    // This time-stepping scheme is implemented from the paper by Akkerman et al., J of Applied Mechanics,2012

    TensorValue<double> Q_current, Q_new_iter, I_w_current, I_w_new_iter, Omega_current, Omega_new;
    VectorValue<double> W_current, V_current, x_com_current;
    Q_new_iter.zero();
    Omega_current.zero();
    Omega_new.zero();
    I_w_new_iter.zero();

    const double rho_f = params[1];

    V_current = V_new;
    W_current = W_new;
    Q_current = Q_new;
    I_w_current = I_w_new;
    x_com_current = x_com_new;

    V_new = dt * (F_b / M + rho_f * F_s / M) + V_current;
    x_com_new = 0.5 * dt * (V_new + V_current) + x_com_current;

    int iter = 0;

    while ((Q_new_iter - Q_new).norm() > TOL || (I_w_new_iter - I_w_new).norm() > TOL)
    {
        Q_new_iter = Q_new;
        I_w_new_iter = I_w_new;
        getSkewSymmetricAngVelTensor(Omega_current, W_current);
        I_w_new = Q_new * I_w_0 * Q_new.transpose();
        W_new = I_w_new.inverse() * (dt * T + I_w_current * W_current);
        getSkewSymmetricAngVelTensor(Omega_new, W_new);
        Q_new = Q_current + 0.25 * dt * (Omega_new + Omega_current) * (Q_new + Q_current);
        ++iter;
    }

    pout << " Number of 6DOF iterations = " << iter << "\n\n";

    return;
} // Solve6DOFSystemofEquations

void
updateVelocityAndPositionOfSolidPoints(const double* params,
                                       VectorValue<double>& x_com,
                                       double M,
                                       VectorValue<double>& V, // linear velocity of the body
                                       VectorValue<double>& W, // angular velocity of the body
                                       TensorValue<double>& Q,
                                       TensorValue<double>& I_w,
                                       TensorValue<double> I_w_0,
                                       VectorValue<double> XCOM,
                                       VectorValue<double> F_b,
                                       EquationSystems* bndry_equation_systems,
                                       EquationSystems& solid_equation_systems,
                                       const double current_time,
                                       const double dts)
{
    VectorValue<double> F_s, Torque;
    F_s.zero();
    Torque.zero();

    calculateFluidForceAndTorque(F_s, Torque, x_com, current_time, bndry_equation_systems);

    Solve6DOFSystemofEquations(params, dts, V, W, x_com, Q, M, I_w, I_w_0, F_b, F_s, Torque);

    VectorValue<double> SS, RR, WxR, X_new;

    RR.zero();
    SS.zero();
    WxR.zero();
    X_new.zero();

    //	setRotationMatrix(W, q_new, q_half, rotation_mat, dt);

    MeshBase& mesh = solid_equation_systems.get_mesh();
    System& X_system = solid_equation_systems.get_system("position_new");
    const unsigned int X_sys_num = X_system.number();

    NumericVector<double>& X_coords = *X_system.solution;
    System& U_system = solid_equation_systems.get_system("velocity_new");
    const unsigned int U_sys_num = U_system.number();
    NumericVector<double>& U_coords = *U_system.solution;

    System& X_half_system = solid_equation_systems.get_system("position_half");
    NumericVector<double>& X_half_coords = *X_half_system.solution;

    const auto node_end = mesh.local_nodes_end();
    for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != node_end; ++it)
    {
        Node* n = *it;
        if (n->n_vars(X_sys_num))
        {
            TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
            const libMesh::Point& X = *n;
            SS = Q * (X - XCOM);
            WxR = W.cross(SS);
            SS += x_com;

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const int u_dof_index = n->dof_number(U_sys_num, d, 0);
                const int x_dof_index = n->dof_number(X_sys_num, d, 0);
                X_coords.set(x_dof_index, SS(d));
                U_coords.set(u_dof_index, V(d) + WxR(d));

                X_half_coords.set(x_dof_index, 0.5 * (SS(d) + X_new(d)));
            }
        }
    }
    X_coords.close();
    X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
    X_system.solution->localize(*X_system.current_local_solution);

    X_half_coords.close();
    X_half_system.get_dof_map().enforce_constraints_exactly(X_half_system, &X_coords);
    X_half_system.solution->localize(*X_half_system.current_local_solution);

    U_coords.close();
    U_system.get_dof_map().enforce_constraints_exactly(U_system, &U_coords);
    U_system.solution->localize(*U_system.current_local_solution);

    return;

} // updateVelocityAndPositionOfSolidPoints

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        const double rho_s = input_db->getDouble("RHO_S");
        const double dx = input_db->getDouble("DX");
        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string exodus_solid_filename = "solid_BC_output_dx_" + std::to_string(double(dx)) + "_rho_" +
                                             std::to_string(double(rho_s)) + "_kappa_" +
                                             std::to_string(double(kappa_s)) + ".ex2";
        const string exodus_bndry_filename = app_initializer->getExodusIIFilename();

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
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double rho_f = input_db->getDouble("RHO");

        input_db->getDoubleArray("R_COM", &COM(0), LIBMESH_DIM);
        double params[2] = { rho_s, rho_f };

        //~ const double grav_const =input_db->getDouble("RHO");

        const double R = input_db->getDouble("R");
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

        const auto node_end = solid_mesh.nodes_end();
        for (auto it = solid_mesh.nodes_begin(); it != node_end; ++it)
        {
            Node* n = *it;
            libMesh::Point& X = *n;
            X(0) += COM(0);
            X(1) += COM(1);
            X(2) += COM(2);
        }

        solid_mesh.prepare_for_use();

        BoundaryMesh bndry_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        solid_mesh.get_boundary_info().sync(bndry_mesh);
        bndry_mesh.prepare_for_use();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
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
        Pointer<IIMethod> ib_method_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         &bndry_mesh,
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
        //~ IBFEMethod::LagForceFcnData body_fcn_data(tether_force_function, sys_data);
        //~ ib_method_ops->registerLagForceFunction(body_fcn_data);

        IIMethod::LagSurfaceForceFcnData surface_fcn_data(tether_force_function, sys_data);
        ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data);
        EquationSystems* bndry_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        // Setup solid systems.
        libMesh::EquationSystems* solid_equation_systems(new EquationSystems(solid_mesh));
        x_new_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("position_new");
        u_new_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("velocity_new");
        x_half_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("position_half");
        u_half_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("velocity_half");

        Order order = FIRST;
        FEFamily family = LAGRANGE;
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_new_" << d;
            x_new_solid_system->add_variable(os.str(), order, family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_new_" << d;
            u_new_solid_system->add_variable(os.str(), order, family);
        }

        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_half_" << d;
            x_half_solid_system->add_variable(os.str(), order, family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_half_" << d;
            u_half_solid_system->add_variable(os.str(), order, family);
        }

        solid_equation_systems->init();

        // This is a horrible hack to set up the position vector.
        {
            MeshBase& mesh = solid_equation_systems->get_mesh();
            System& X_new_system = solid_equation_systems->get_system("position_new");
            const unsigned int X_new_sys_num = X_new_system.number();
            NumericVector<double>& X_new_coords = *X_new_system.solution;
            const auto node_end = mesh.local_nodes_end();
            for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != node_end; ++it)
            {
                Node* n = *it;
                if (n->n_vars(X_new_sys_num))
                {
                    TBOX_ASSERT(n->n_vars(X_new_sys_num) == NDIM);
                    const libMesh::Point& X = *n;
                    libMesh::Point x = X;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int dof_index = n->dof_number(X_new_sys_num, d, 0);
                        X_new_coords.set(dof_index, x(d));
                    }
                }
            }
            X_new_coords.close();
            X_new_system.get_dof_map().enforce_constraints_exactly(X_new_system, &X_new_coords);
            X_new_system.solution->localize(*X_new_system.current_local_solution);

            System& X_half_system = solid_equation_systems->get_system("position_half");
            const unsigned int X_half_sys_num = X_half_system.number();
            NumericVector<double>& X_half_coords = *X_half_system.solution;

            for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != node_end; ++it)
            {
                Node* n = *it;
                if (n->n_vars(X_half_sys_num))
                {
                    TBOX_ASSERT(n->n_vars(X_half_sys_num) == NDIM);
                    const libMesh::Point& X = *n;
                    libMesh::Point x = X;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int dof_index = n->dof_number(X_half_sys_num, d, 0);
                        X_half_coords.set(dof_index, x(d));
                    }
                }
            }
            X_half_coords.close();
            X_half_system.get_dof_map().enforce_constraints_exactly(X_half_system, &X_half_coords);
            X_half_system.solution->localize(*X_half_system.current_local_solution);
        }

        x_new_solid_system->assemble_before_solve = false;
        x_new_solid_system->assemble();

        u_new_solid_system->assemble_before_solve = false;
        u_new_solid_system->assemble();

        x_half_solid_system->assemble_before_solve = false;
        x_half_solid_system->assemble();

        u_half_solid_system->assemble_before_solve = false;
        u_half_solid_system->assemble();

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

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_solid_io(uses_exodus ? new ExodusII_IO(solid_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_io(uses_exodus ? new ExodusII_IO(bndry_mesh) : NULL);

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

        TensorValue<double> I_w_0, I_w;

        VectorValue<double> V, W, F_b, F_s, Torque, W_hat, V_hat;
        VectorValue<double> x_com, x_com_hat;
        double M;

        TensorValue<double> Q(std::cos(TOL), -std::sin(TOL), 0.0, std::sin(TOL), std::cos(TOL), 0.0, 0.0, 0.0, 1.0);
        I_w.zero();
        x_com.zero();
        x_com_hat.zero();
        W.zero();
        V.zero();
        W_hat.zero();
        V_hat.zero();
        Torque.zero();
        F_s.zero();
        F_b.zero();
        //****************************** Initialize RBD parameters **************************************//
        calculateGeomQuantitiesOfStructure(params, M, I_w, x_com, solid_equation_systems);
        I_w_0 = I_w;
        calculateGravitationalForce(params, F_b, solid_equation_systems);
        calculateFluidForceAndTorque(F_s, Torque, x_com, loop_time, bndry_equation_systems);

        //******************************************************************************//
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
                exodus_solid_io->write_timestep(
                    exodus_solid_filename, *solid_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_bndry_io->write_timestep(
                    exodus_bndry_filename, *bndry_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        if (SAMRAI_MPI::getRank() == 0)
        {
            x_com_stream.open("xcom_BC_3D_FSI_dx_" + std::to_string(double(dx)) + "_rho_" +
                                  std::to_string(double(rho_s)) + "_kappa_" + std::to_string(double(kappa_s)) +
                                  "_eta_" + std::to_string(double(eta_s)) + "_glass.curve",
                              ios_base::out | ios_base::trunc);
            v_new_stream.open("v_BC_3D_FSI_dx_" + std::to_string(double(dx)) + "_rho_" + std::to_string(double(rho_s)) +
                                  "_kappa_" + std::to_string(double(kappa_s)) + "_eta_" +
                                  std::to_string(double(eta_s)) + "_glass.curve",
                              ios_base::out | ios_base::trunc);
            w_new_stream.open("w_BC_3D_FSI_dx_" + std::to_string(double(dx)) + "_rho_" + std::to_string(double(rho_s)) +
                                  "_kappa_" + std::to_string(double(kappa_s)) + "_eta_" +
                                  std::to_string(double(eta_s)) + "_glass.curve",
                              ios_base::out | ios_base::trunc);
            x_com_stream.precision(10);
            v_new_stream.precision(10);
            w_new_stream.precision(10);
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

            //****************************** RBD code **************************************//
            calculateGravitationalForce(params, F_b, solid_equation_systems);

            updateVelocityAndPositionOfSolidPoints(params,
                                                   x_com,
                                                   M,
                                                   V,
                                                   W,
                                                   Q,
                                                   I_w,
                                                   I_w_0,
                                                   COM,
                                                   F_b,
                                                   bndry_equation_systems,
                                                   *solid_equation_systems,
                                                   loop_time,
                                                   0.5 * dt);

            time_integrator->advanceHierarchy(dt);

            calculateGravitationalForce(params, F_b, solid_equation_systems);

            updateVelocityAndPositionOfSolidPoints(params,
                                                   x_com,
                                                   M,
                                                   V,
                                                   W,
                                                   Q,
                                                   I_w,
                                                   I_w_0,
                                                   COM,
                                                   F_b,
                                                   bndry_equation_systems,
                                                   *solid_equation_systems,
                                                   loop_time + 0.5 * dt,
                                                   0.5 * dt);

            pout << "\n\n"
                 << " Translationl Velocity = " << V << "\n";

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
                    exodus_solid_io->write_timestep(exodus_solid_filename,
                                                    *solid_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
                    exodus_bndry_io->write_timestep(exodus_bndry_filename,
                                                    *bndry_equation_systems,
                                                    iteration_num / viz_dump_interval + 1,
                                                    loop_time);
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

            if (SAMRAI_MPI::getRank() == 0)
            {
                x_com_stream << loop_time << "\t " << x_com(0) << "\t " << x_com(1) << "\t " << x_com(2) << endl;
                v_new_stream << loop_time << "\t " << V(0) << "\t " << V(1) << "\t " << V(2) << endl;
                w_new_stream << loop_time << "\t " << W(0) << "\t " << W(1) << "\t " << W(2) << endl;
                //~ w_new_stream << loop_time << " " << W_new(0) << " " << W_new(1) << " " << W_new(2) << endl;
            }
        }
        x_new_solid_system->clear();
        u_new_solid_system->clear();
        x_half_solid_system->clear();
        u_half_solid_system->clear();
        solid_equation_systems->clear();
        bndry_equation_systems->clear();

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            x_com_stream.close();
            v_new_stream.close();
            w_new_stream.close();
        }
        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main
