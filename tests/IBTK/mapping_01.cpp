// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
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
#include <libmesh/fe_base.h>
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
#include <ibamr/app_namespaces.h>

#include "../tests.h"

template <int dim, libMesh::Order order, libMesh::FEFamily fe_family, libMesh::ElemType elem_type>
void
test(LibMeshInit& init)
{
    ReplicatedMesh mesh(init.comm(), dim);
    switch (elem_type)
    {
    case QUAD4:
    case QUAD9:
    case HEX8:
    case HEX27:
        MeshTools::Generation::build_sphere(mesh, 1.0, 2, elem_type);
        break;
    case TRI3:
    case TRI6:
        MeshTools::Generation::build_square(mesh, 3, 3, 0.0, 0.5, 0.0, 2.0, elem_type);
        break;
    case TET4:
    case TET10:
        MeshTools::Generation::build_cube(mesh, 3, 3, 3, 0.0, 0.5, 0.0, 0.25, 0.0, 8.0, elem_type);
        break;
    default:
        TBOX_ASSERT(false);
    }

    std::vector<libMesh::QuadratureType> quad_types;
    quad_types.push_back(QGAUSS);
    quad_types.push_back(QGRID);
#if !LIBMESH_VERSION_LESS_THAN(1, 5, 0)
    quad_types.push_back(QNODAL);
#endif

    for (const libMesh::QuadratureType quad_type : quad_types)
    {
        std::unique_ptr<QBase> quad = QBase::build(quad_type, dim, THIRD);
        quad->init(elem_type);
        FEType fe_type(order, fe_family);
        std::unique_ptr<FEBase> libmesh_fe = FEBase::build(dim, fe_type);
        libmesh_fe->attach_quadrature_rule(quad.get());
        libmesh_fe->get_xyz();
        libmesh_fe->get_JxW();

        const quadrature_key_type key{ elem_type, quad_type, THIRD, true };
        std::unique_ptr<FEMapping<dim> > mapping =
            FEMapping<dim>::build(key, FEUpdateFlags::update_JxW | FEUpdateFlags::update_quadrature_points);

        for (auto elem_iter = mesh.active_local_elements_begin(); elem_iter != mesh.active_local_elements_end();
             ++elem_iter)
        {
            libmesh_fe->reinit(*elem_iter);
            mapping->reinit(*elem_iter);

            const std::vector<double>& JxW = libmesh_fe->get_JxW();
            const std::vector<double>& JxW_2 = mapping->getJxW();
            for (unsigned int i = 0; i < JxW.size(); ++i)
            {
                TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < 1e-14 * std::max(1.0, std::abs(JxW[i])));
            }

            const std::vector<libMesh::Point>& q = libmesh_fe->get_xyz();
            const std::vector<libMesh::Point>& q_2 = mapping->getQuadraturePoints();

            for (unsigned int i = 0; i < q.size(); ++i)
            {
                TBOX_ASSERT(q[i].relative_fuzzy_equals(q_2[i], 1e-14));
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

    // 2d
    test<2, FIRST, LAGRANGE, TRI3>(init);
    test<2, SECOND, LAGRANGE, TRI6>(init);
    test<2, FIRST, LAGRANGE, QUAD4>(init);
    test<2, SECOND, LAGRANGE, QUAD9>(init);

    // 3d
    test<3, FIRST, LAGRANGE, TET4>(init);
    test<3, SECOND, LAGRANGE, TET10>(init);
    test<3, FIRST, LAGRANGE, HEX8>(init);
    test<3, SECOND, LAGRANGE, HEX27>(init);

    std::ofstream output("output");
}
