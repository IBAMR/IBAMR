// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

// Headers for basic libMesh objects
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <boost/multi_array.hpp>

#include <ibamr/app_namespaces.h>

// Verify that the new function, IBTK::get_max_edge_length, prints out the
// same result as the old get_elem_hmax function (that was internal to
// FEDataManager). Test with several 3D meshes.

void
log_max_edge_length(const ReplicatedMesh& mesh)
{
    const unsigned int dim = mesh.mesh_dimension();
    boost::multi_array<double, 2> X_node;
    for (auto elem_iter = mesh.active_local_elements_begin(); elem_iter != mesh.active_local_elements_end();
         ++elem_iter)
    {
        auto elem = *elem_iter;
        // Don't bother with deformation: just use the material coordinates
        boost::multi_array<double, 2>::extent_gen extent;
        const unsigned int n_nodes = elem->n_nodes();
        X_node.resize(extent[n_nodes][dim]);

        const Node* const* nodes = elem->get_nodes();
        for (unsigned int node_n = 0; node_n < n_nodes; ++node_n)
            for (unsigned int d = 0; d < dim; ++d) X_node[node_n][d] = (*nodes[node_n])(d);

        plog << std::setprecision(12) << get_max_edge_length(elem, X_node) << std::endl;
    }
}

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const double radius = 2.0;
        const unsigned int n_refinements = 1;

        {
            plog << "Test 1: tet4" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_cube(mesh, 2, 3, 3, 0.0, 0.7, 0.1, 0.9, 1.0, 2.0, TET4);
            log_max_edge_length(mesh);
        }

        {
            plog << std::endl << "Test 1: tet10" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_cube(mesh, 2, 3, 3, 0.0, 0.7, 0.1, 0.9, 1.0, 2.0, TET10);
            log_max_edge_length(mesh);
        }

        // like the last test, but wobble the triangulation a bit so that the
        // edges are not all straight lines
        {
            plog << std::endl << "Test 2: tet10" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_cube(mesh, 2, 3, 3, 0.0, 0.7, 0.1, 0.9, 1.0, 2.0, TET10);
            for (std::size_t node_n = 0; node_n < mesh.n_nodes(); ++node_n)
            {
                const double y = (*mesh.node_ptr(node_n))(1);
                (*mesh.node_ptr(node_n))(0) += 2 * y * (1 - y);
            }
            log_max_edge_length(mesh);
        }

        {
            plog << std::endl << "Test 3: pyramid5" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_cube(mesh, 2, 3, 3, 0.0, 0.7, 0.1, 0.9, 1.0, 2.0, PYRAMID5);
            log_max_edge_length(mesh);
        }

        {
            plog << std::endl << "Test 4: hex8" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, HEX8);
            log_max_edge_length(mesh);
        }

        {
            plog << std::endl << "Test 5: hex27" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, HEX27);
            log_max_edge_length(mesh);
        }
    }
} // main
