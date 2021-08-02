// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/libmesh_utilities.h>

#include <tbox/SAMRAIManager.h>

#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/libmesh_version.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

#include <SAMRAI_config.h>

// Set up application namespace declarations
#include <fstream>
#include <iostream>

#include <ibamr/app_namespaces.h>

System&
setup_deformation_system(ReplicatedMesh& mesh, EquationSystems& equation_systems, const libMesh::Order order)
{
    // Set up the system
    auto& X_system = equation_systems.add_system<ExplicitSystem>("X");
    const auto X_sys_num = X_system.number();
    const unsigned int dim = mesh.mesh_dimension();
    for (unsigned int d = 0; d < dim; ++d)
    {
        X_system.add_variable("X_" + std::to_string(d), order, LAGRANGE);
    }
    equation_systems.init();

    // Set up the system solution to match the current mesh coordinates
    auto& X_solution = *X_system.solution;
    const auto el_begin = mesh.active_local_elements_begin();
    const auto el_end = mesh.active_local_elements_end();
    for (auto el_it = el_begin; el_it != el_end; ++el_it)
    {
        const libMesh::Elem* const elem = *el_it;
        const unsigned int n_nodes = elem->n_nodes();
        for (unsigned int k = 0; k < n_nodes; ++k)
        {
            const libMesh::Node* const node = elem->node_ptr(k);
            for (unsigned int d = 0; d < dim; ++d)
            {
                TBOX_ASSERT(node->n_dofs(X_sys_num, d) == 1);
                const auto dof = node->dof_number(X_sys_num, d, 0);
                X_solution.set(dof, (*node)(d));
            }
        }
    }
    X_solution.close();

    return X_system;
}

void
test(LibMeshInit& init,
     const unsigned int dim,
     const std::string geometry,
     const ElemType elem_type,
     const Order order,
     std::ofstream& out,
     const bool use_nodal_quadrature = true,
     const bool use_amr = false)
{
    ReplicatedMesh mesh(init.comm(), dim);
    if (geometry == "sphere")
    {
        MeshTools::Generation::build_sphere(mesh, 0.1, 1, elem_type);
    }
    else
    {
        TBOX_ASSERT(geometry == "cube");
        if (dim == 2)
            MeshTools::Generation::build_square(mesh, 2, 3, 0.0, 0.5, 0.0, 2.0, elem_type);
        else
            MeshTools::Generation::build_cube(mesh, 2, 3, 2, 0.0, 0.5, 0.0, 0.25, 0.0, 8.0, elem_type);
    }
    MeshRefinement mesh_refinement(mesh);
    if (use_amr)
    {
        const auto el_begin = mesh.local_elements_begin();
        const auto el_end = mesh.local_elements_end();
        std::size_t i = 0;
        for (auto el_it = el_begin; el_it != el_end; ++el_it, ++i)
        {
            if (i % 8 == 0) (*el_it)->set_refinement_flag(Elem::REFINE);
        }
        mesh_refinement.refine_and_coarsen_elements();
    }
    EquationSystems equation_systems(mesh);
    auto& X_system = setup_deformation_system(mesh, equation_systems, order == FIRST ? FIRST : SECOND);

    // Verify that the default settings are the same as both what we get with
    // nodal quadratures and what libMesh computes
    const auto quad_type = use_nodal_quadrature ? (order == FIRST ? QTRAP : QSIMPSON) : QGAUSS;
    std::vector<libMeshWrappers::BoundingBox> nodal_bboxes_1 =
        get_local_element_bounding_boxes(mesh, X_system, quad_type, order, false, 1.0, true, 1.0);
    std::vector<libMeshWrappers::BoundingBox> nodal_bboxes_2 = get_local_element_bounding_boxes(mesh, X_system);
    TBOX_ASSERT(nodal_bboxes_1.size() == nodal_bboxes_2.size());

    const auto el_begin = mesh.local_elements_begin();
    const auto el_end = mesh.local_elements_end();
    std::size_t i = 0;
    for (auto el_it = el_begin; el_it != el_end; ++el_it, ++i)
    {
        const Elem* elem = *el_it;
        if (!elem->active()) continue;
        const libMeshWrappers::BoundingBox& box_1 = nodal_bboxes_1[i];
        const libMeshWrappers::BoundingBox& box_2 = nodal_bboxes_2[i];

        if (use_nodal_quadrature)
        {
            // boxes should be the same
            TBOX_ASSERT((box_1.min() - box_2.min()).norm() < std::max(1.0, box_1.min().norm()) * 1e-16);
            TBOX_ASSERT((box_1.max() - box_2.max()).norm() < std::max(1.0, box_1.max().norm()) * 1e-16);

            // nodes should be in the box
            const unsigned int n_nodes = elem->n_nodes();
            for (unsigned int k = 0; k < n_nodes; ++k)
            {
                TBOX_ASSERT(box_1.contains_point(*elem->node_ptr(k)));
            }
        }
        else
        {
            // the non-nodal box should be a subset of the nodal one
#if !LIBMESH_VERSION_LESS_THAN(1, 2, 0)
            libMeshWrappers::BoundingBox box_union(box_2);
            box_union.union_with(box_1);
            TBOX_ASSERT(box_union.min() == box_2.min());
            TBOX_ASSERT(box_union.max() == box_2.max());
#endif

            // box 2 should be a superset of box 1: i.e., the corners of box 1
            // should be in box 2
#if !LIBMESH_VERSION_LESS_THAN(1, 4, 0)
            TBOX_ASSERT(box_2.signed_distance(box_1.min()) <= 0.0);
            TBOX_ASSERT(box_2.signed_distance(box_1.max()) <= 0.0);
#endif
        }
        out << box_2.first << ", " << box_2.second << std::endl;

        // Since X is just the initial coordinates of the mesh, we should
        // match the bounding box which is computed directly from the
        // element too.
#if !LIBMESH_VERSION_LESS_THAN(1, 2, 0)
        libMeshWrappers::BoundingBox box_3 = (*el_it)->loose_bounding_box();
        if (order == FIRST)
        {
            TBOX_ASSERT((box_2.min() - box_3.min()).norm() < std::max(1.0, box_2.min().norm()) * 1e-16);
            TBOX_ASSERT((box_2.max() - box_3.max()).norm() < std::max(1.0, box_2.max().norm()) * 1e-16);
        }

        // in both cases check that we are inside the libMesh box: it claims
        // to be a loose fit
        libMeshWrappers::BoundingBox box_union(box_3);
        box_union.union_with(box_2);
        TBOX_ASSERT((box_union.min() - box_3.min()).norm() < std::max(1.0, box_union.min().norm()) * 1e-16);
        TBOX_ASSERT((box_union.max() - box_3.max()).norm() < std::max(1.0, box_union.max().norm()) * 1e-16);
#endif

        // box 3 should be a superset of box 2: i.e., the corners of box 2
        // should be in box 3
#if !LIBMESH_VERSION_LESS_THAN(1, 4, 0)
        TBOX_ASSERT(box_2.signed_distance(box_1.min()) <= 0.0);
        TBOX_ASSERT(box_2.signed_distance(box_1.max()) <= 0.0);
#endif
    }
}

int
main(int argc, char** argv)
{
    std::ofstream out("output");

    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    LibMeshInit& init = ibtk_init.getLibMeshInit();

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    const bool use_amr = input_db->getBoolWithDefault("use_amr", false);

    if (NDIM == 2)
    {
        // test 2D
        out << "2D:" << std::endl;
        test(init, 2, "sphere", TRI3, FIRST, out, true, use_amr);
        test(init, 2, "sphere", TRI3, FIRST, out, false, use_amr);
        out << "TRI3 OK" << std::endl;
        test(init, 2, "sphere", TRI6, THIRD, out, true, use_amr);
        test(init, 2, "sphere", TRI6, THIRD, out, false, use_amr);
        out << "TRI6 OK" << std::endl;
        test(init, 2, "sphere", QUAD4, FIRST, out, true, use_amr);
        test(init, 2, "sphere", QUAD4, FIRST, out, false, use_amr);
        out << "QUAD4 OK" << std::endl;
        test(init, 2, "sphere", QUAD9, THIRD, out, true, use_amr);
        test(init, 2, "sphere", QUAD9, THIRD, out, false, use_amr);
        out << "QUAD9 OK" << std::endl;
    }

    if (NDIM == 3)
    {
        // test 3D
        out << "3D:" << std::endl;
        test(init, 3, "cube", TET4, FIRST, out, true, use_amr);
        test(init, 3, "cube", TET4, FIRST, out, false, use_amr);
        out << "TET4 OK" << std::endl;
        test(init, 3, "cube", TET10, THIRD, out, true, use_amr);
        test(init, 3, "cube", TET10, THIRD, out, false, use_amr);
        out << "TET10 OK" << std::endl;
        test(init, 3, "sphere", HEX8, FIRST, out, true, use_amr);
        test(init, 3, "sphere", HEX8, FIRST, out, false, use_amr);
        out << "HEX8 OK" << std::endl;
        test(init, 3, "sphere", HEX27, THIRD, out, true, use_amr);
        test(init, 3, "sphere", HEX27, THIRD, out, false, use_amr);
        out << "HEX27 OK" << std::endl;
    }
}
