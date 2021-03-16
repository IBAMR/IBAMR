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

// Config files

#include <SAMRAI_config.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/FEMapping.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/QuadratureCache.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <boost/multi_array.hpp>

#include <ibamr/app_namespaces.h>

#include "../tests.h"

// Verify that FEMapping and descendants output the same values as libMesh::FEMap.
using key_type = quadrature_key_type;

template <int dim, int spacedim>
void
test_cube(LibMeshInit& init,
          FEMapping<dim, spacedim>& jc_1,
          FEMapping<dim, spacedim>& jc_2,
          FEMapping<dim - 1, spacedim>& jc_boundary,
          const key_type key)
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
    jc_1.reinit(*mesh.active_local_elements_begin());
    const std::vector<double>& JxW = jc_1.getJxW();
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
        jc_1.reinit(*elem_iter);
        jc_2.reinit(*elem_iter);
        const std::vector<double>& JxW = jc_1.getJxW();
        const std::vector<double>& JxW_2 = jc_2.getJxW();
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
        BoundaryInfo& boundary_info = mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
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
            jc_boundary.reinit(*elem_iter);
            const std::vector<double>& JxW = jc_boundary.getJxW();
            const std::vector<double>& JxW_2 = fe_map.get_JxW();
            for (unsigned int i = 0; i < JxW.size(); ++i)
            {
                TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
            }
        }
    }
}

template <int dim, int spacedim>
void
test_circle(LibMeshInit& init,
            FEMapping<dim, spacedim>& jc_1,
            FEMapping<dim, spacedim>& jc_2,
            FEMapping<dim - 1, spacedim>& jc_boundary,
            const int n_refines,
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
    jc_1.reinit(*mesh.active_local_elements_begin());
    const std::vector<double>& JxW = jc_1.getJxW();
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
        jc_1.reinit(*elem_iter);
        jc_2.reinit(*elem_iter);
        const std::vector<double>& JxW = jc_1.getJxW();
        const std::vector<double>& JxW_2 = jc_2.getJxW();
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
        BoundaryInfo& boundary_info = mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
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
            jc_boundary.reinit(*elem_iter);
            const std::vector<double>& JxW = jc_boundary.getJxW();
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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    LibMeshInit& init = ibtk_init.getLibMeshInit();

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();

        unsigned int test_n = 1;
        {
            plog << "Test " << test_n << ": TRI3 square" << std::endl;
            const key_type key(TRI3, QGAUSS, THIRD, true);
            Tri3Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": TRI3 square" << std::endl;
            const key_type key(TRI3, QGAUSS, THIRD, true);
            Tri3Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": TRI6 square" << std::endl;
            const key_type key(TRI6, QGAUSS, FIFTH, true);
            Tri6Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE2, QGAUSS, FIFTH, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": TRI6 square" << std::endl;
            const key_type key(TRI6, QGAUSS, FIFTH, true);
            Tri6Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE2, QGAUSS, FIFTH, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad4 square" << std::endl;
            const key_type key(QUAD4, QGAUSS, THIRD, true);
            Quad4Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad4 circle" << std::endl;
            const key_type key(QUAD4, QGAUSS, THIRD, true);
            Quad4Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE2, QGAUSS, THIRD, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 5, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad9 square" << std::endl;
            const key_type key(QUAD9, QGAUSS, FOURTH, true);
            Quad9Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE3, QGAUSS, FOURTH, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": Quad9 circle" << std::endl;
            const key_type key(QUAD9, QGAUSS, FOURTH, true);
            Quad9Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<2> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(EDGE3, QGAUSS, FOURTH, true);
            FELagrangeMapping<1, 2> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 4, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": TET4 cube" << std::endl;
            const key_type key(TET4, QGAUSS, THIRD, true);
            Tet4Mapping jac_calc_1(key, FEUpdateFlags::update_JxW);
            FELagrangeMapping<3> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(TRI3, QGAUSS, THIRD, true);
            FELagrangeMapping<2, 3> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX8 square" << std::endl;
            const key_type key(HEX8, QGAUSS, THIRD, true);
            // HEX8 doesn't have a custom calculator yet
            FELagrangeMapping<3> jac_calc_1(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            FELagrangeMapping<3> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(QUAD4, QGAUSS, THIRD, true);
            FELagrangeMapping<2, 3> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX8 circle" << std::endl;
            const key_type key(HEX8, QGAUSS, THIRD, true);
            // HEX8 doesn't have a custom calculator yet
            FELagrangeMapping<3> jac_calc_1(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            FELagrangeMapping<3> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(QUAD4, QGAUSS, THIRD, true);
            FELagrangeMapping<2, 3> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 4, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX27 square" << std::endl;
            const key_type key(HEX27, QGAUSS, FOURTH, true);
            // HEX27 doesn't have a custom calculator yet
            FELagrangeMapping<3> jac_calc_1(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            FELagrangeMapping<3> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(QUAD9, QGAUSS, FOURTH, true);
            FELagrangeMapping<2, 3> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_cube(init, jac_calc_1, jac_calc_2, jac_calc_b, key);
            ++test_n;
        }

        {
            plog << "Test " << test_n << ": HEX27 circle" << std::endl;
            const key_type key(HEX27, QGAUSS, FOURTH, true);
            // HEX27 doesn't have a custom calculator yet
            FELagrangeMapping<3> jac_calc_1(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            FELagrangeMapping<3> jac_calc_2(key, std::get<0>(key), FEUpdateFlags::update_JxW);
            const key_type boundary_key(QUAD9, QGAUSS, FOURTH, true);
            FELagrangeMapping<2, 3> jac_calc_b(boundary_key, std::get<0>(boundary_key), FEUpdateFlags::update_JxW);
            test_circle(init, jac_calc_1, jac_calc_2, jac_calc_b, 2, key);
            ++test_n;
        }
    }
} // main
