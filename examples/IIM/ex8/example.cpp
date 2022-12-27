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

#include "ibamr/FEMechanicsBase.h"
#include <ibamr/FEMechanicsExplicitIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IIMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
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
#include <libmesh/edge_edge2.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/explicit_system.h>
#include <libmesh/fe.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/solver_configuration.h>

#include <boost/multi_array.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Application includes
#include "FeedbackForcer.h"

// Elasticity model data.
namespace ModelData
{
static double kappa_s_line = 1.0e6;
static double eta_s_line = 1.0e6;
static double kappa_s_FSI_lower = 1.0e6;
static double eta_FSI_lower = 1.0e6;
static double kappa_s_FSI_upper = 1.0e6;
static double eta_FSI_upper = 1.0e6;
static double DX = 0.01;
static double D = 1.0;
static double L = 6.0;
static double H0 = 2.0;
static double mu_s_lower, lambda_s_lower, mu_s_upper, lambda_s_upper;
static std::string stress_funtion;

System* x_new_tube_lower_solid_system;
System* u_new_tube_lower_solid_system;
System* x_new_tube_lower_surface_system;
System* Tau_new_tube_lower_surface_system;

System* x_new_tube_upper_solid_system;
System* u_new_tube_upper_solid_system;
System* x_new_tube_upper_surface_system;
System* Tau_new_tube_upper_surface_system;

EquationSystems* boundary_tube_lower_systems;
EquationSystems* boundary_tube_upper_systems;

static BoundaryInfo* tube_lower_copy_info;
static BoundaryInfo* tube_upper_copy_info;

void
FSI_tether_line_force_function(VectorValue<double>& F,
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
        F(d) = kappa_s_line * (X(d) - x(d)) + eta_s_line * (0.0 - u_bndry_n) * n(d);
    }
    return;
}

void
FSI_tether_tube_lower_force_function(VectorValue<double>& F,
                                     const VectorValue<double>& /*n*/,
                                     const VectorValue<double>& /*N*/,
                                     const TensorValue<double>& /*FF*/,
                                     const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the
                                                                    // boundary mesh
                                     const libMesh::Point& X_bndry, // X_bndry gives reference coordinates on the
                                                                    // boundary mesh
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

    const std::vector<double>& U = *var_data[0];

    double disp = 0.0;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_tube_lower_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_tube_lower_solid_system->point_value(d, X_bndry, interior_parent);
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        disp += (x_solid[d] - x_bndry(d)) * (x_solid[d] - x_bndry(d));
    }
    disp = sqrt(disp);
    TBOX_ASSERT(disp < DX);

    // The tether force is proportional to the mismatch between the positions
    // and velocities.

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s_FSI_lower * (x_solid[d] - x_bndry(d)) + eta_FSI_lower * (u_solid[d] - U[d]);
    }

    return;
} // FSI_tether_force_function

void
FSI_tether_tube_upper_force_function(VectorValue<double>& F,
                                     const VectorValue<double>& /*n*/,
                                     const VectorValue<double>& /*N*/,
                                     const TensorValue<double>& /*FF*/,
                                     const libMesh::Point& x_bndry, // x_bndry gives current   coordinates on the
                                                                    // boundary mesh
                                     const libMesh::Point& X_bndry, // X_bndry gives reference coordinates on the
                                                                    // boundary mesh
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

    const std::vector<double>& U = *var_data[0];

    double disp = 0.0;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_tube_upper_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_tube_upper_solid_system->point_value(d, X_bndry, interior_parent);
    }

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        disp += (x_solid[d] - x_bndry(d)) * (x_solid[d] - x_bndry(d));
    }
    disp = sqrt(disp);

    // The tether force is proportional to the mismatch between the positions
    // and velocities.

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s_FSI_upper * (x_solid[d] - x_bndry(d)) + eta_FSI_upper * (u_solid[d] - U[d]);
    }

    return;
} // FSI_tether_tube_upper_force_function

void
solid_surface_force_tube_upper_function(VectorValue<double>& F,
                                        const VectorValue<double>& /*n*/,
                                        const VectorValue<double>& /*N*/,
                                        const TensorValue<double>& /*FF*/,
                                        const libMesh::Point& /*x*/,
                                        const libMesh::Point& X,
                                        Elem* const elem,
                                        const unsigned short int side,
                                        const vector<const vector<double>*>& /*var_data*/,
                                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                        double /*time*/,
                                        void* /*ctx*/)
{
    MeshBase& mesh_bndry = boundary_tube_upper_systems->get_mesh();
    std::unique_ptr<Elem> side_elem = elem->side_ptr(side);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const auto el_begin = mesh_bndry.active_local_elements_begin();
        const auto el_end = mesh_bndry.active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* elem_bndry = *el_it;

            if ((elem_bndry->contains_point(X)) && (tube_upper_copy_info->has_boundary_id(elem, side, 6)))
            {
                F(d) = Tau_new_tube_upper_surface_system->point_value(d, X, elem_bndry); //&side_elem);
            }
        }
    }

    return;
} // solid_surface_force_tube_upper_function

void
solid_surface_force_tube_lower_function(VectorValue<double>& F,
                                        const VectorValue<double>& /*n*/,
                                        const VectorValue<double>& /*N*/,
                                        const TensorValue<double>& /*FF*/,
                                        const libMesh::Point& /*x*/,
                                        const libMesh::Point& X,
                                        Elem* const elem,
                                        const unsigned short int side,
                                        const vector<const vector<double>*>& /*var_data*/,
                                        const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                        double /*time*/,
                                        void* /*ctx*/)
{
    MeshBase& mesh_bndry = boundary_tube_lower_systems->get_mesh();
    std::unique_ptr<Elem> side_elem = elem->side_ptr(side);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const auto el_begin = mesh_bndry.active_local_elements_begin();
        const auto el_end = mesh_bndry.active_local_elements_end();
        for (auto el_it = el_begin; el_it != el_end; ++el_it)
        {
            const Elem* elem_bndry = *el_it;

            if ((elem_bndry->contains_point(X)) && (tube_lower_copy_info->has_boundary_id(elem, side, 5)))
            {
                F(d) = Tau_new_tube_lower_surface_system->point_value(d, X, elem_bndry); //&side_elem);
            }
        }
    }

    return;
} // solid_surface_force_tube_lower_function

void
PK1_dev_stress_tube_upper_function(TensorValue<double>& PP,
                                   const TensorValue<double>& FF,
                                   const libMesh::Point& /*x*/,
                                   const libMesh::Point& /*X*/,
                                   Elem* const /*elem*/,
                                   const vector<const vector<double>*>& /*var_data*/,
                                   const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                   double /*time*/,
                                   void* /*ctx*/)
{
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    static const TensorValue<double> IO(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);

    std::vector<double> x_surface(NDIM, 0.0);

    //  Unmodified St. Venant-Kirchhoff Model

    const TensorValue<double> EE = 0.5 * (CC - II);
    PP = FF * (2 * mu_s_upper * EE + lambda_s_upper * EE.tr() * II);

    return;
} // PK1_dev_stress_function

void
PK1_dev_stress_tube_lower_function(TensorValue<double>& PP,
                                   const TensorValue<double>& FF,
                                   const libMesh::Point& /*x*/,
                                   const libMesh::Point& /*X*/,
                                   Elem* const /*elem*/,
                                   const vector<const vector<double>*>& /*var_data*/,
                                   const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                                   double /*time*/,
                                   void* /*ctx*/)
{
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    static const TensorValue<double> IO(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF, NDIM);

    std::vector<double> x_surface(NDIM, 0.0);

    //  Unmodified St. Venant-Kirchhoff Model

    const TensorValue<double> EE = 0.5 * (CC - II);
    PP = FF * (2 * mu_s_lower * EE + lambda_s_lower * EE.tr() * II);

    return;
} // PK1_dev_stress_function

} // namespace ModelData
using namespace ModelData;

static ofstream dx_posn_stream;

void postprocess_data(tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      tbox::Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      const FEMechanicsExplicitIntegrator* fem_solver,
                      ReplicatedMesh& tube_mesh,
                      EquationSystems* tube_equation_systems,
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
int
main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        tbox::Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        tbox::Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string viz_dump_dirname = app_initializer->getVizDumpDirectory();
        const string exodus_tube_lower_filename = viz_dump_dirname + "/tube_lower.ex2";
        const string exodus_tube_upper_filename = viz_dump_dirname + "/tube_upper.ex2";
        const string exodus_bndry_tube_lower_filename = viz_dump_dirname + "/bndry_tube_lower.ex2";
        const string exodus_bndry_tube_upper_filename = viz_dump_dirname + "/bndry_tube_upper.ex2";
        const string exodus_bndry_line1_filename = viz_dump_dirname + "/bndry_line1.ex2";
        const string exodus_bndry_line2_filename = viz_dump_dirname + "/bndry_line2.ex2";
        const string exodus_bndry_line3_filename = viz_dump_dirname + "/bndry_line3.ex2";
        const string exodus_bndry_line4_filename = viz_dump_dirname + "/bndry_line4.ex2";

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
        DX = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * DX;
        D = input_db->getDouble("D");
        L = input_db->getDouble("L");
        H0 = input_db->getDouble("H0");
        const double n_cycles = input_db->getDouble("NCYCLE");

        libMesh::RealVectorValue extrusion_vec = { 0, 0, L };

        ReplicatedMesh tube_lower_mesh(init.comm(), NDIM);

        string tube_elem_type = input_db->getString("TUBE_ELEM_TYPE");
        tube_lower_mesh.read(input_db->getString("TUBE_LOWER_MESH_FILENAME"), NULL);

        using MeshTools::Modification::translate;

        translate(tube_lower_mesh, 0, H0 - 1, 0);

        // Imposing Dirichlet BC at the two ends
        const auto end_lower_el = tube_lower_mesh.elements_end();
        for (auto el = tube_lower_mesh.elements_begin(); el != end_lower_el; ++el)
        {
            const Elem* elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    BoundaryInfo& boundary_info = tube_lower_mesh.get_boundary_info();
                    if ((boundary_info.has_boundary_id(elem, side, 3)) ||
                        (boundary_info.has_boundary_id(elem, side, 4)))
                    {
                        boundary_info.add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID);
                    }
                }
            }
        }
        tube_lower_mesh.prepare_for_use();

        // Side 1 appears to be the inner part of the tube
        // Side 2 appears to be the outer part of the tube (and maybe the end)
        // Side 3 appears to be the enterance!
        // Side 4 appears to be the outlet!
        tube_lower_copy_info = &tube_lower_mesh.get_boundary_info();
        BoundaryMesh bndry_tube_lower_mesh(init.comm(), NDIM - 1);
        tube_lower_mesh.get_boundary_info().sync({ 5 }, bndry_tube_lower_mesh);
        bndry_tube_lower_mesh.prepare_for_use();

        // ************************************************************************//
        ReplicatedMesh tube_upper_mesh(init.comm(), NDIM);

        tube_upper_mesh.read(input_db->getString("TUBE_UPPER_MESH_FILENAME"), NULL);

        translate(tube_upper_mesh, 0, H0 - 1, 0);

        // Imposing Dirichlet BC at the two ends
        const auto end_upper_el = tube_upper_mesh.elements_end();
        for (auto el = tube_upper_mesh.elements_begin(); el != end_upper_el; ++el)
        {
            const Elem* elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (at_mesh_bdry)
                {
                    BoundaryInfo& boundary_info = tube_upper_mesh.get_boundary_info();
                    if ((boundary_info.has_boundary_id(elem, side, 1)) ||
                        (boundary_info.has_boundary_id(elem, side, 2)))
                    {
                        boundary_info.add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XY_BDRY_ID);
                    }
                }
            }
        }
        tube_upper_mesh.prepare_for_use();

        tube_upper_copy_info = &tube_upper_mesh.get_boundary_info();
        BoundaryMesh bndry_tube_upper_mesh(init.comm(), NDIM - 1);
        tube_upper_mesh.get_boundary_info().sync({ 6 }, bndry_tube_upper_mesh);
        bndry_tube_upper_mesh.prepare_for_use();

        // ************************************************************************//
        const unsigned int nn = ceil(H0 / ds);

        Mesh line1_mesh(init.comm(), NDIM - 1);

        double xc1_position = input_db->getDouble("XC1");
        int node_id = 0;

        line1_mesh.reserve_nodes(nn + 1);
        line1_mesh.reserve_elem(nn);
        for (unsigned int i = 0; i <= nn; i++)
        {
            line1_mesh.add_point(
                libMesh::Point(xc1_position - 0.5 * D, H0 * static_cast<Real>(i) / static_cast<Real>(nn)), node_id++);
        }

        BoundaryInfo& boundary_info_line1 = line1_mesh.get_boundary_info();

        for (unsigned int i = 0; i < nn; i++)
        {
            Elem* elem = line1_mesh.add_elem(new Edge2);
            elem->set_node(0) = line1_mesh.node_ptr(i);
            elem->set_node(1) = line1_mesh.node_ptr(i + 1);
            if (i == 0) boundary_info_line1.add_side(elem, 0, 0);
            if (i == (nn - 1)) boundary_info_line1.add_side(elem, 1, 1);
        }
        line1_mesh.prepare_for_use();
        // ************************************************************************//
        Mesh line2_mesh(init.comm(), NDIM - 1);

        node_id = 0;

        line2_mesh.reserve_nodes(nn + 1);
        line2_mesh.reserve_elem(nn);
        for (unsigned int i = 0; i <= nn; i++)
        {
            line2_mesh.add_point(
                libMesh::Point(xc1_position + 0.5 * D, H0 * static_cast<Real>(i) / static_cast<Real>(nn)), node_id++);
        }

        BoundaryInfo& boundary_info_line2 = line2_mesh.get_boundary_info();

        for (unsigned int i = 0; i < nn; i++)
        {
            Elem* elem = line2_mesh.add_elem(new Edge2);
            elem->set_node(0) = line2_mesh.node_ptr(i);
            elem->set_node(1) = line2_mesh.node_ptr(i + 1);
            if (i == 0) boundary_info_line2.add_side(elem, 0, 0);
            if (i == (nn - 1)) boundary_info_line2.add_side(elem, 1, 1);
        }
        line2_mesh.prepare_for_use();
        // ************************************************************************//
        Mesh line3_mesh(init.comm(), NDIM - 1);

        double xc2_position = input_db->getDouble("XC2");
        node_id = 0;

        line3_mesh.reserve_nodes(nn + 1);
        line3_mesh.reserve_elem(nn);
        for (unsigned int i = 0; i <= nn; i++)
        {
            line3_mesh.add_point(
                libMesh::Point(xc2_position - 0.5 * D, H0 * static_cast<Real>(i) / static_cast<Real>(nn)), node_id++);
        }

        BoundaryInfo& boundary_info_line3 = line3_mesh.get_boundary_info();

        for (unsigned int i = 0; i < nn; i++)
        {
            Elem* elem = line3_mesh.add_elem(new Edge2);
            elem->set_node(0) = line3_mesh.node_ptr(i);
            elem->set_node(1) = line3_mesh.node_ptr(i + 1);
            if (i == 0) boundary_info_line3.add_side(elem, 0, 0);
            if (i == (nn - 1)) boundary_info_line3.add_side(elem, 1, 1);
        }
        line3_mesh.prepare_for_use();
        // ************************************************************************//
        Mesh line4_mesh(init.comm(), NDIM - 1);
        node_id = 0;

        line4_mesh.reserve_nodes(nn + 1);
        line4_mesh.reserve_elem(nn);
        for (unsigned int i = 0; i <= nn; i++)
        {
            line4_mesh.add_point(
                libMesh::Point(xc2_position + 0.5 * D, H0 * static_cast<Real>(i) / static_cast<Real>(nn)), node_id++);
        }

        BoundaryInfo& boundary_info_line4 = line4_mesh.get_boundary_info();

        for (unsigned int i = 0; i < nn; i++)
        {
            Elem* elem = line4_mesh.add_elem(new Edge2);
            elem->set_node(0) = line4_mesh.node_ptr(i);
            elem->set_node(1) = line4_mesh.node_ptr(i + 1);
            if (i == 0) boundary_info_line4.add_side(elem, 0, 0);
            if (i == (nn - 1)) boundary_info_line4.add_side(elem, 1, 1);
        }
        line4_mesh.prepare_for_use();
        // ************************************************************************//
        vector<MeshBase*> bndry_meshes(6);
        bndry_meshes[0] = &bndry_tube_lower_mesh;
        bndry_meshes[1] = &bndry_tube_upper_mesh;
        bndry_meshes[2] = &line1_mesh;
        bndry_meshes[3] = &line2_mesh;
        bndry_meshes[4] = &line3_mesh;
        bndry_meshes[5] = &line4_mesh;

        vector<MeshBase*> meshes(2);
        meshes[0] = &tube_lower_mesh;
        meshes[1] = &tube_upper_mesh;

        mu_s_lower = input_db->getDouble("MU_S_LOWER");
        lambda_s_lower = input_db->getDouble("LAMBDA_S_LOWER");
        mu_s_upper = input_db->getDouble("MU_S_UPPER");
        lambda_s_upper = input_db->getDouble("LAMBDA_S_UPPER");

        // Setup the model parameters.
        kappa_s_line = input_db->getDouble("KAPPA_S_LINE");
        eta_s_line = input_db->getDouble("ETA_S_LINE");
        kappa_s_FSI_lower = input_db->getDouble("KAPPA_S_FSI_LOWER");
        eta_FSI_lower = input_db->getDouble("ETA_FSI_LOWER");
        kappa_s_FSI_upper = input_db->getDouble("KAPPA_S_FSI_UPPER");
        eta_FSI_upper = input_db->getDouble("ETA_FSI_UPPER");

        // Setup the time stepping parameters.
        const double loop_time_end = input_db->getDouble("END_TIME");
        double dt = input_db->getDouble("DT");

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

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.

        tbox::Pointer<IIMethod> ibfe_bndry_ops =
            new IIMethod("IIMethod",
                         app_initializer->getComponentDatabase("IIMethod"),
                         bndry_meshes,
                         app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        tbox::Pointer<FEMechanicsExplicitIntegrator> fem_solver = new FEMechanicsExplicitIntegrator(
            "FEMechanicsExplicitIntegrator",
            app_initializer->getComponentDatabase("FEMechanicsExplicitIntegrator"),
            meshes,
            app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));

        tbox::Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ibfe_bndry_ops,
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

        // attach velocity

        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> velocity_data(1);
        velocity_data[0] = SystemData(fem_solver->getVelocitySystemName(), vars);

        ibfe_bndry_ops->initializeFEEquationSystems();

        vector<SystemData> sys_data(1, SystemData(IIMethod::VELOCITY_SYSTEM_NAME, vars));

        // Configure the FE solver.

        FEMechanicsBase::LagSurfaceForceFcnData solid_surface_force_tube_lower_data(
            solid_surface_force_tube_lower_function, velocity_data);
        fem_solver->registerLagSurfaceForceFunction(solid_surface_force_tube_lower_data, 0);

        FEMechanicsBase::LagSurfaceForceFcnData solid_surface_force_tube_upper_data(
            solid_surface_force_tube_upper_function, velocity_data);
        fem_solver->registerLagSurfaceForceFunction(solid_surface_force_tube_upper_data, 1);

        FEMechanicsBase::PK1StressFcnData PK1_dev_stress_tube_lower_data(PK1_dev_stress_tube_lower_function,
                                                                         velocity_data);
        PK1_dev_stress_tube_lower_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
        fem_solver->registerPK1StressFunction(PK1_dev_stress_tube_lower_data, 0);

        FEMechanicsBase::PK1StressFcnData PK1_dev_stress_tube_upper_data(PK1_dev_stress_tube_upper_function,
                                                                         velocity_data);
        PK1_dev_stress_tube_upper_data.quad_order =
            Utility::string_to_enum<libMesh::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER", "FIFTH"));
        fem_solver->registerPK1StressFunction(PK1_dev_stress_tube_upper_data, 1);

        IIMethod::LagSurfaceForceFcnData surface_FSI_tube_lower_fcn_data(FSI_tether_tube_lower_force_function,
                                                                         sys_data);
        IIMethod::LagSurfaceForceFcnData surface_FSI_tube_upper_fcn_data(FSI_tether_tube_upper_force_function,
                                                                         sys_data);
        IIMethod::LagSurfaceForceFcnData surface_FSI_line_fcn_data(FSI_tether_line_force_function, sys_data);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_tube_lower_fcn_data, 0);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_tube_upper_fcn_data, 1);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_line_fcn_data, 2);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_line_fcn_data, 3);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_line_fcn_data, 4);
        ibfe_bndry_ops->registerLagSurfaceForceFunction(surface_FSI_line_fcn_data, 5);

        EquationSystems* bndry_tube_lower_equation_systems = ibfe_bndry_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* bndry_tube_upper_equation_systems = ibfe_bndry_ops->getFEDataManager(1)->getEquationSystems();
        EquationSystems* bndry_line1_equation_systems = ibfe_bndry_ops->getFEDataManager(2)->getEquationSystems();
        EquationSystems* bndry_line2_equation_systems = ibfe_bndry_ops->getFEDataManager(3)->getEquationSystems();
        EquationSystems* bndry_line3_equation_systems = ibfe_bndry_ops->getFEDataManager(4)->getEquationSystems();
        EquationSystems* bndry_line4_equation_systems = ibfe_bndry_ops->getFEDataManager(5)->getEquationSystems();

        fem_solver->initializeFEEquationSystems();
        EquationSystems* tube_lower_equation_systems = fem_solver->getFEData(0)->getEquationSystems();
        EquationSystems* tube_upper_equation_systems = fem_solver->getFEData(1)->getEquationSystems();

        // **************Create Eulerian initial condition specification objects.**************** //
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

        if (input_db->keyExists("ForcingFunction"))
        {
            tbox::Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
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

            tbox::Pointer<FeedbackForcer> feedback_forcer =
                new FeedbackForcer(xc1_position, xc2_position, D, navier_stokes_integrator, patch_hierarchy);
            time_integrator->registerBodyForceFunction(feedback_forcer);
        }
        // Set up visualization plot file writers.

        tbox::Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        std::unique_ptr<ExodusII_IO> exodus_tube_lower_io(uses_exodus ? new ExodusII_IO(tube_lower_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_tube_upper_io(uses_exodus ? new ExodusII_IO(tube_upper_mesh) : NULL);

        std::unique_ptr<ExodusII_IO> exodus_bndry_tube_lower_io(uses_exodus ? new ExodusII_IO(bndry_tube_lower_mesh) :
                                                                              NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_tube_upper_io(uses_exodus ? new ExodusII_IO(bndry_tube_upper_mesh) :
                                                                              NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_line1_io(uses_exodus ? new ExodusII_IO(line1_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_line2_io(uses_exodus ? new ExodusII_IO(line2_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_line3_io(uses_exodus ? new ExodusII_IO(line3_mesh) : NULL);
        std::unique_ptr<ExodusII_IO> exodus_bndry_line4_io(uses_exodus ? new ExodusII_IO(line4_mesh) : NULL);

        ibfe_bndry_ops->initializeFEData();

        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Initialize solver data.
        fem_solver->initializeFEData();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = 0;
        double loop_time = 0.0;
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
                exodus_tube_lower_io->write_timestep(exodus_tube_lower_filename,
                                                     *tube_lower_equation_systems,
                                                     iteration_num / viz_dump_interval + 1,
                                                     loop_time);

                exodus_tube_upper_io->write_timestep(exodus_tube_upper_filename,
                                                     *tube_upper_equation_systems,
                                                     iteration_num / viz_dump_interval + 1,
                                                     loop_time);

                exodus_bndry_tube_lower_io->write_timestep(exodus_bndry_tube_lower_filename,
                                                           *bndry_tube_lower_equation_systems,
                                                           iteration_num / viz_dump_interval + 1,
                                                           loop_time);

                exodus_bndry_tube_upper_io->write_timestep(exodus_bndry_tube_upper_filename,
                                                           *bndry_tube_upper_equation_systems,
                                                           iteration_num / viz_dump_interval + 1,
                                                           loop_time);

                exodus_bndry_line1_io->write_timestep(exodus_bndry_line1_filename,
                                                      *bndry_line1_equation_systems,
                                                      iteration_num / viz_dump_interval + 1,
                                                      loop_time);

                exodus_bndry_line2_io->write_timestep(exodus_bndry_line2_filename,
                                                      *bndry_line2_equation_systems,
                                                      iteration_num / viz_dump_interval + 1,
                                                      loop_time);

                exodus_bndry_line3_io->write_timestep(exodus_bndry_line3_filename,
                                                      *bndry_line3_equation_systems,
                                                      iteration_num / viz_dump_interval + 1,
                                                      loop_time);

                exodus_bndry_line4_io->write_timestep(exodus_bndry_line4_filename,
                                                      *bndry_line4_equation_systems,
                                                      iteration_num / viz_dump_interval + 1,
                                                      loop_time);
            }
        }

        // Open streams to save lift and drag coefficients.

        if (SAMRAI_MPI::getRank() == 0)
        {
            dx_posn_stream.open("dx_MFAC2_QUAD9_DS2_128_dt4em5_H02.curve", ios_base::out | ios_base::trunc);
        }

        // Open streams to save volume of structure.

        // Main time step loop.
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end))
        {
            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            boundary_tube_lower_systems = bndry_tube_lower_equation_systems;
            System& X_tube_lower_system =
                tube_lower_equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
            x_new_tube_lower_solid_system = &X_tube_lower_system;

            System& U_tube_lower_system =
                tube_lower_equation_systems->get_system<System>(fem_solver->getVelocitySystemName());
            u_new_tube_lower_solid_system = &U_tube_lower_system;

            boundary_tube_upper_systems = bndry_tube_upper_equation_systems;
            System& X_tube_upper_system =
                tube_upper_equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
            x_new_tube_upper_solid_system = &X_tube_upper_system;

            System& U_tube_upper_system =
                tube_upper_equation_systems->get_system<System>(fem_solver->getVelocitySystemName());
            u_new_tube_upper_solid_system = &U_tube_upper_system;

            Tau_new_tube_lower_surface_system =
                &bndry_tube_lower_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_tube_lower_surface_system =
                &bndry_tube_lower_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            Tau_new_tube_upper_surface_system =
                &bndry_tube_upper_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_tube_upper_surface_system =
                &bndry_tube_upper_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            dt = time_integrator->getMaximumTimeStepSize();

            for (int ii = 0; ii < static_cast<int>(n_cycles); ii++)
            {
                fem_solver->preprocessIntegrateData(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                    loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles,
                                                    /*num_cycles*/ 1);
                fem_solver->modifiedTrapezoidalStep(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                    loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles);
                fem_solver->postprocessIntegrateData(loop_time + 0.5 * static_cast<double>(ii) * dt / n_cycles,
                                                     loop_time + 0.5 * static_cast<double>(ii + 1) * dt / n_cycles,
                                                     /*num_cycles*/ 1);
            }

            x_new_tube_lower_solid_system =
                &tube_lower_equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());

            u_new_tube_lower_solid_system =
                &tube_lower_equation_systems->get_system<System>(fem_solver->getVelocitySystemName());

            x_new_tube_upper_solid_system =
                &tube_upper_equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());

            u_new_tube_upper_solid_system =
                &tube_upper_equation_systems->get_system<System>(fem_solver->getVelocitySystemName());

            time_integrator->advanceHierarchy(dt); // FSI solution (IIMethod)

            boundary_tube_lower_systems = bndry_tube_lower_equation_systems;
            boundary_tube_upper_systems = bndry_tube_upper_equation_systems;

            Tau_new_tube_lower_surface_system =
                &bndry_tube_lower_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_tube_lower_surface_system =
                &bndry_tube_lower_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

            Tau_new_tube_upper_surface_system =
                &bndry_tube_upper_equation_systems->get_system<System>(IIMethod::TAU_OUT_SYSTEM_NAME);
            x_new_tube_upper_surface_system =
                &bndry_tube_upper_equation_systems->get_system<System>(IIMethod::COORDS_SYSTEM_NAME);

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
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0))
            {
                pout << "\nWriting visualization files...\n\n";

                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }

                if (uses_exodus)
                {
                    exodus_tube_lower_io->write_timestep(exodus_tube_lower_filename,
                                                         *tube_lower_equation_systems,
                                                         iteration_num / viz_dump_interval + 1,
                                                         loop_time);

                    exodus_tube_upper_io->write_timestep(exodus_tube_upper_filename,
                                                         *tube_upper_equation_systems,
                                                         iteration_num / viz_dump_interval + 1,
                                                         loop_time);

                    exodus_bndry_tube_lower_io->write_timestep(exodus_bndry_tube_lower_filename,
                                                               *bndry_tube_lower_equation_systems,
                                                               iteration_num / viz_dump_interval + 1,
                                                               loop_time);

                    exodus_bndry_tube_upper_io->write_timestep(exodus_bndry_tube_upper_filename,
                                                               *bndry_tube_upper_equation_systems,
                                                               iteration_num / viz_dump_interval + 1,
                                                               loop_time);

                    exodus_bndry_line1_io->write_timestep(exodus_bndry_line1_filename,
                                                          *bndry_line1_equation_systems,
                                                          iteration_num / viz_dump_interval + 1,
                                                          loop_time);

                    exodus_bndry_line2_io->write_timestep(exodus_bndry_line2_filename,
                                                          *bndry_line2_equation_systems,
                                                          iteration_num / viz_dump_interval + 1,
                                                          loop_time);

                    exodus_bndry_line3_io->write_timestep(exodus_bndry_line3_filename,
                                                          *bndry_line3_equation_systems,
                                                          iteration_num / viz_dump_interval + 1,
                                                          loop_time);

                    exodus_bndry_line4_io->write_timestep(exodus_bndry_line4_filename,
                                                          *bndry_line4_equation_systems,
                                                          iteration_num / viz_dump_interval + 1,
                                                          loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }

            pout << "\nWriting state data...\n\n";
            postprocess_data(patch_hierarchy,
                             navier_stokes_integrator,
                             fem_solver,
                             tube_lower_mesh,
                             tube_lower_equation_systems,
                             iteration_num,
                             loop_time,
                             postproc_data_dump_dirname);

        } // end of time loop

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            dx_posn_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main

void
postprocess_data(tbox::Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                 tbox::Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                 const FEMechanicsExplicitIntegrator* const fem_solver,
                 ReplicatedMesh& /*tube_mesh*/,
                 EquationSystems* tube_equation_systems,
                 const int /*iteration_num*/,
                 const double loop_time,
                 const string& /*data_dump_dirname*/)
{
    System& X_system = tube_equation_systems->get_system<System>(fem_solver->getCurrentCoordinatesSystemName());
    NumericVector<double>* X_vec = X_system.solution.get();
    std::unique_ptr<NumericVector<Number> > X_serial_vec = NumericVector<Number>::build(X_vec->comm());
    X_serial_vec->init(X_vec->size(), true, SERIAL);
    X_vec->localize(*X_serial_vec);
    DofMap& X_dof_map = X_system.get_dof_map();
    vector<unsigned int> vars(2);
    vars[0] = 0;
    vars[1] = 1;
    MeshFunction X_fcn(*tube_equation_systems, *X_serial_vec, X_dof_map, vars);
    X_fcn.init();
    DenseVector<double> X_A(2);
    X_fcn(libMesh::Point(0.0, 1.95 + H0, 0), 0.0, X_A);
    if (SAMRAI_MPI::getRank() == 0)
    {
        dx_posn_stream.precision(12);
        dx_posn_stream.setf(ios::fixed, ios::floatfield);
        dx_posn_stream << loop_time << " " << sqrt(X_A(0) * X_A(0) + (X_A(1) - 1.95 - H0) * (X_A(1) - 1.95 - H0))
                       << endl;
    }
    return;
} // postprocess_data
