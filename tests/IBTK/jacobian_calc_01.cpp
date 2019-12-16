// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
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
#include <IBAMR_config.h>
#include <IBTK_config.h>

#include <SAMRAI_config.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/FEMapCache.h>
#include <ibtk/JacobianCalculator.h>
#include <ibtk/QuadratureCache.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include <boost/multi_array.hpp>

// Verify that JacobianCalc and descendants output the same values as libMesh::FEMap.
using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

void
test_cube(LibMeshInit& init, JacobianCalculator& jc_1, JacobianCalculator& jc_2, JacobianCalculator& jc_boundary, const int dim, const key_type key)
{
    const auto elem_type = std::get<0>(key);

    FEMapCache map_cache(dim);
    QuadratureCache quad_cache(dim);

    ReplicatedMesh mesh(init.comm(), dim);
    if (dim == 2)
        MeshTools::Generation::build_square(mesh, 10, 10, 0.0, 0.5, 0.0, 2.0, elem_type);
    else if (dim == 3)
        MeshTools::Generation::build_cube(mesh, 10, 10, 10, 0.0, 0.5, 0.0, 0.25, 0.0, 8.0, elem_type);
    else
        TBOX_ASSERT(false);

    // check that we get the same thing with both calculators and libMesh's general code:
    const std::vector<double>& JxW = jc_1.get_JxW(*mesh.active_local_elements_begin());
    for (const double jxw : JxW) plog << std::setprecision(12) << jxw << '\n';

    double volume = 0;
    double volume_2 = 0;
    for (auto elem_iter = mesh.active_local_elements_begin(); elem_iter != mesh.active_local_elements_end();
         ++elem_iter)
    {
        FEMap& fe_map = map_cache[key];
        QBase& quad = quad_cache[key];
        fe_map.compute_map(dim, quad.get_weights(), *elem_iter, false);
        // all computed JxW values should agree
        const std::vector<double>& JxW = jc_1.get_JxW(*elem_iter);
        const std::vector<double>& JxW_2 = jc_2.get_JxW(*elem_iter);
        const std::vector<double>& JxW_3 = fe_map.get_JxW();
        for (unsigned int i = 0; i < JxW.size(); ++i)
        {
            TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
            TBOX_ASSERT(std::abs(JxW[i] - JxW_3[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
        }
        volume += std::accumulate(JxW.begin(), JxW.end(), 0.0);
        volume_2 += std::accumulate(JxW_2.begin(), JxW_2.end(), 0.0);
    }
    TBOX_ASSERT(std::abs(volume - volume_2) < 1e-15 * volume);
    plog << "volume is " << volume << '\n';

    // also test the surface mesh:
    {
        BoundaryMesh boundary_mesh(mesh.comm(), mesh.mesh_dimension() - 1);
        mesh.boundary_info->sync(boundary_mesh);
        boundary_mesh.prepare_for_use();
        TBOX_ASSERT(boundary_mesh.spatial_dimension() == mesh.mesh_dimension());
        TBOX_ASSERT(boundary_mesh.mesh_dimension() == mesh.mesh_dimension() - 1);

        FEMapCache boundary_map_cache(dim - 1);
        QuadratureCache boundary_quad_cache(dim - 1);

        for (auto elem_iter = boundary_mesh.active_local_elements_begin();
             elem_iter != boundary_mesh.active_local_elements_end();
             ++elem_iter)
        {
            key_type boundary_key = key;
            std::get<0>(boundary_key) = (*elem_iter)->type();

            FEMap& fe_map = boundary_map_cache[boundary_key];
            QBase& quad = boundary_quad_cache[boundary_key];
            fe_map.compute_map(dim - 1, quad.get_weights(), *elem_iter, false);
            // all computed JxW values should agree
            const std::vector<double>& JxW = jc_boundary.get_JxW(*elem_iter);
            const std::vector<double>& JxW_2 = fe_map.get_JxW();
            for (unsigned int i = 0; i < JxW.size(); ++i)
            {
                TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
            }
        }
    }
}

void
test_circle(LibMeshInit& init,
            JacobianCalculator& jc_1,
            JacobianCalculator& jc_2,
            JacobianCalculator& jc_boundary,
            const int n_refines,
            const int dim,
            const key_type key)
{
    const auto elem_type = std::get<0>(key);

    FEMapCache map_cache(dim);
    QuadratureCache quad_cache(dim);

    const double radius = 1.0;
    ReplicatedMesh mesh(init.comm(), dim);
    MeshTools::Generation::build_sphere(mesh, radius, n_refines, elem_type);

    // Ensure nodes on the surface are on the analytic boundary.
    MeshBase::element_iterator el_end = mesh.elements_end();
    for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
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
                n = radius * n.unit();
            }
        }
    }
    mesh.prepare_for_use();

    // check that we get the same thing with both calculators and libMesh's general code:
    const std::vector<double>& JxW = jc_1.get_JxW(*mesh.active_local_elements_begin());
    for (const double jxw : JxW) plog << std::setprecision(12) << jxw << '\n';

    double volume = 0;
    double volume_2 = 0;
    for (auto elem_iter = mesh.active_local_elements_begin(); elem_iter != mesh.active_local_elements_end();
         ++elem_iter)
    {
        FEMap& fe_map = map_cache[key];
        QBase& quad = quad_cache[key];
        fe_map.compute_map(dim, quad.get_weights(), *elem_iter, false);
        // all computed JxW values should agree
        const std::vector<double>& JxW = jc_1.get_JxW(*elem_iter);
        const std::vector<double>& JxW_2 = jc_2.get_JxW(*elem_iter);
        const std::vector<double>& JxW_3 = fe_map.get_JxW();
        for (unsigned int i = 0; i < JxW.size(); ++i)
        {
            TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
            TBOX_ASSERT(std::abs(JxW[i] - JxW_3[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
        }
        volume += std::accumulate(JxW.begin(), JxW.end(), 0.0);
        volume_2 += std::accumulate(JxW_2.begin(), JxW_2.end(), 0.0);
    }
    TBOX_ASSERT(std::abs(volume - volume_2) < 1e-15 * volume);
    plog << "volume is " << volume << '\n';

    // also test the surface mesh:
    {
        BoundaryMesh boundary_mesh(mesh.comm(), mesh.mesh_dimension() - 1);
        mesh.boundary_info->sync(boundary_mesh);
        boundary_mesh.prepare_for_use();
        TBOX_ASSERT(boundary_mesh.spatial_dimension() == mesh.mesh_dimension());
        TBOX_ASSERT(boundary_mesh.mesh_dimension() == mesh.mesh_dimension() - 1);

        FEMapCache boundary_map_cache(dim - 1);
        QuadratureCache boundary_quad_cache(dim - 1);

        for (auto elem_iter = boundary_mesh.active_local_elements_begin();
             elem_iter != boundary_mesh.active_local_elements_end();
             ++elem_iter)
        {
            key_type boundary_key = key;
            std::get<0>(boundary_key) = (*elem_iter)->type();

            FEMap& fe_map = boundary_map_cache[boundary_key];
            QBase& quad = boundary_quad_cache[boundary_key];
            fe_map.compute_map(dim - 1, quad.get_weights(), *elem_iter, false);
            // all computed JxW values should agree
            const std::vector<double>& JxW = jc_boundary.get_JxW(*elem_iter);
            const std::vector<double>& JxW_2 = fe_map.get_JxW();
            for (unsigned int i = 0; i < JxW.size(); ++i)
            {
                TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
            }
        }
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

        unsigned int test_n = 1;
        {
            plog << "Test " << test_n << ": TRI3 square" << std::endl;
            const key_type key(TRI3, QGAUSS, THIRD);
            Tri3JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<2> jac_calc_2(key);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD);
            LagrangeJacobianCalculator<1, 2> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 2, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": TRI3 square" << std::endl;
            const key_type key(TRI3, QGAUSS, THIRD);
            Tri3JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<2> jac_calc_2(key);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD);
            LagrangeJacobianCalculator<1, 2> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 2, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad4 square" << std::endl;
            const key_type key(QUAD4, QGAUSS, THIRD);
            Quad4JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<2> jac_calc_2(key);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD);
            LagrangeJacobianCalculator<1, 2> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 2, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad4 circle" << std::endl;
            const key_type key(QUAD4, QGAUSS, THIRD);
            Quad4JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<2> jac_calc_2(key);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD);
            LagrangeJacobianCalculator<1, 2> jac_calc_b(boundary_key);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 5, 2, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad9 square" << std::endl;
            const key_type key(QUAD9, QGAUSS, FOURTH);
            Quad9JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<2> jac_calc_2(key);
            const key_type boundary_key(EDGE3, QGAUSS, FOURTH);
            LagrangeJacobianCalculator<1, 2> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 2, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad9 circle" << std::endl;
            const key_type key(QUAD9, QGAUSS, FOURTH);
            Quad9JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<2> jac_calc_2(key);
            const key_type boundary_key(EDGE3, QGAUSS, FOURTH);
            LagrangeJacobianCalculator<1, 2> jac_calc_b(boundary_key);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 4, 2, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": TET4 cube" << std::endl;
            const key_type key(TET4, QGAUSS, THIRD);
            Tet4JacobianCalculator jac_calc_1(key);
            LagrangeJacobianCalculator<3> jac_calc_2(key);
            const key_type boundary_key(TRI3, QGAUSS, THIRD);
            LagrangeJacobianCalculator<2, 3> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 3, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX8 square" << std::endl;
            const key_type key(HEX8, QGAUSS, THIRD);
            // HEX8 doesn't have a custom calculator yet
            LagrangeJacobianCalculator<3> jac_calc_1(key);
            LagrangeJacobianCalculator<3> jac_calc_2(key);
            const key_type boundary_key(QUAD4, QGAUSS, THIRD);
            LagrangeJacobianCalculator<2, 3> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 3, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX8 circle" << std::endl;
            const key_type key(HEX8, QGAUSS, THIRD);
            // HEX8 doesn't have a custom calculator yet
            LagrangeJacobianCalculator<3> jac_calc_1(key);
            LagrangeJacobianCalculator<3> jac_calc_2(key);
            const key_type boundary_key(QUAD4, QGAUSS, THIRD);
            LagrangeJacobianCalculator<2, 3> jac_calc_b(boundary_key);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 4, 3, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX27 square" << std::endl;
            const key_type key(HEX27, QGAUSS, FOURTH);
            // HEX27 doesn't have a custom calculator yet
            LagrangeJacobianCalculator<3> jac_calc_1(key);
            LagrangeJacobianCalculator<3> jac_calc_2(key);
            const key_type boundary_key(QUAD9, QGAUSS, FOURTH);
            LagrangeJacobianCalculator<2, 3> jac_calc_b(boundary_key);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, 3, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX27 circle" << std::endl;
            const key_type key(HEX27, QGAUSS, FOURTH);
            // HEX27 doesn't have a custom calculator yet
            LagrangeJacobianCalculator<3> jac_calc_1(key);
            LagrangeJacobianCalculator<3> jac_calc_2(key);
            const key_type boundary_key(QUAD9, QGAUSS, FOURTH);
            LagrangeJacobianCalculator<2, 3> jac_calc_b(boundary_key);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 2, 3, key);
            ++test_n;
        }
    }

    SAMRAIManager::shutdown();
} // main
