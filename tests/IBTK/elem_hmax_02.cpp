// Copyright (c) 2019, Boyce Griffith
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

// Headers for basic libMesh objects
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include <boost/multi_array.hpp>

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
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

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

        {
            plog << std::endl << "Test 2: pyramid5" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_cube(mesh, 2, 3, 3, 0.0, 0.7, 0.1, 0.9, 1.0, 2.0, PYRAMID5);
            log_max_edge_length(mesh);
        }

        {
            plog << std::endl << "Test 3: hex8" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, HEX8);
            log_max_edge_length(mesh);
        }

        {
            plog << std::endl << "Test 4: hex27" << std::endl;
            ReplicatedMesh mesh(init.comm(), 3);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, HEX27);
            log_max_edge_length(mesh);
        }
    }

    SAMRAIManager::shutdown();
} // main
