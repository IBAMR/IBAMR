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
#include <libmesh/abaqus_io.h>
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/fe_base.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/FEValues.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Test that our own FEValues matches what libMesh calculates for surface meshes

enum class MeshType
{
    libmesh = 0,
    readin = 1
};

template <int dim, libMesh::Order order, libMesh::FEFamily fe_family, libMesh::ElemType elem_type>
void
test(LibMeshInit& init, const MeshType mesh_type = MeshType::libmesh)
{
    constexpr int spacedim = dim + 1;

    ReplicatedMesh solid_mesh(init.comm(), spacedim);
    switch (mesh_type)
    {
    case MeshType::libmesh:
        switch (elem_type)
        {
        case QUAD4:
        case QUAD9:
        case HEX8:
        case HEX27:
            MeshTools::Generation::build_sphere(solid_mesh, 1.0, 2, elem_type);
            break;
        case TRI3:
        case TRI6:
            MeshTools::Generation::build_square(solid_mesh, 3, 3, 0.0, 0.5, 0.0, 2.0, elem_type);
            break;
        case TET4:
        case TET10:
            MeshTools::Generation::build_cube(solid_mesh, 3, 3, 3, 0.0, 0.5, 0.0, 0.25, 0.0, 8.0, elem_type);
            break;
        default:
            TBOX_ASSERT(false);
        }
        break;
    case MeshType::readin:
        switch (elem_type)
        {
        case TET4:
        case TET10:
        {
            AbaqusIO abaqus_io(solid_mesh);
            abaqus_io.read(SOURCE_DIR "/klein.inp");
            if (elem_type == TET10) solid_mesh.all_second_order();
            // this seems like a bug in libMesh, but nevertheless is
            // required (at least for 1.5.1)
            solid_mesh.prepare_for_use();
        }
        break;
        default:
            TBOX_ASSERT(false);
        }
        break;
    default:
        TBOX_ASSERT(false);
    }

    // set up boundary mesh:
    BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
    BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
    boundary_info.sync(boundary_mesh);
    boundary_mesh.prepare_for_use();
    Mesh& mesh = boundary_mesh;
    // the rest is the same as fe_values_01 (with the small exception of using
    // the face element type)
    // TODO - surely there is a better way to do this
    ElemType face_elem_type{};
    switch (elem_type)
    {
    case TRI3:
    case QUAD4:
        face_elem_type = EDGE2;
        break;
    case TRI6:
    case QUAD9:
        face_elem_type = EDGE3;
        break;
    case HEX8:
        face_elem_type = QUAD4;
        break;
    case HEX27:
        face_elem_type = QUAD9;
        break;
    case TET4:
        face_elem_type = TRI3;
        break;
    case TET10:
        face_elem_type = TRI6;
        break;
    default:
        TBOX_ASSERT(false);
    }

    // libMesh values:
    std::unique_ptr<QBase> quad = QBase::build(QGAUSS, dim, THIRD);
    quad->init(face_elem_type);
    FEType fe_type(order, fe_family);
    std::unique_ptr<FEBase> libmesh_fe = FEBase::build(dim, fe_type);
    libmesh_fe->attach_quadrature_rule(quad.get());
    libmesh_fe->get_xyz();
    libmesh_fe->get_JxW();
    libmesh_fe->get_phi();
    libmesh_fe->get_dphi();

    // IBTK values:
    std::unique_ptr<QBase> quad_2 = QBase::build(QGAUSS, dim, THIRD);
    quad_2->init(face_elem_type);
    IBTK::FEValues<dim, spacedim> ibtk_fe(quad_2.get(),
                                          fe_type,
                                          IBTK::update_quadrature_points | IBTK::update_JxW | IBTK::update_phi |
                                              IBTK::update_dphi);

    for (auto elem_iter = mesh.active_local_elements_begin(); elem_iter != mesh.active_local_elements_end();
         ++elem_iter)
    {
        libmesh_fe->reinit(*elem_iter);
        ibtk_fe.reinit(*elem_iter);

        // check mapping values:
        const std::vector<double>& JxW = libmesh_fe->get_JxW();
        const std::vector<double>& JxW_2 = ibtk_fe.getJxW();
        for (unsigned int i = 0; i < JxW.size(); ++i)
        {
            // the read-in test requires a lighter tolerance
            const double tol = mesh_type == MeshType::readin ? 1e-13 : 1e-14;
            TBOX_ASSERT(std::abs(JxW[i] - JxW_2[i]) < tol * std::max(1.0, std::abs(JxW[i])));
        }

        const std::vector<libMesh::Point>& q = libmesh_fe->get_xyz();
        const std::vector<libMesh::Point>& q_2 = ibtk_fe.getQuadraturePoints();

        for (unsigned int i = 0; i < q.size(); ++i)
        {
            TBOX_ASSERT(q[i].relative_fuzzy_equals(q_2[i], 1e-15));
        }

        // Shape functions are tabulated in the same way so they should be
        // identical, down to the bit:
        const std::vector<std::vector<double> >& phi = libmesh_fe->get_phi();
        const std::vector<std::vector<double> >& phi_2 = ibtk_fe.getShapeValues();
        TBOX_ASSERT(phi == phi_2);

        // gradients are calculated with our own mapping classes so these will
        // be slightly different from what libMesh gets:
        const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi = libmesh_fe->get_dphi();
        const std::vector<std::vector<libMesh::VectorValue<double> > >& dphi_2 = ibtk_fe.getShapeGradients();
        for (unsigned int i = 0; i < dphi.size(); ++i)
        {
            for (unsigned int q = 0; q < dphi[i].size(); ++q)
            {
                // the read-in test requires a lighter tolerance
                const double tol = mesh_type == MeshType::readin ? 1e-11 : 1e-12;
                TBOX_ASSERT(dphi[i][q].relative_fuzzy_equals(dphi_2[i][q], tol));
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

    // 2d, libMesh meshes:
    test<1, FIRST, LAGRANGE, TRI3>(init);
    test<1, SECOND, LAGRANGE, TRI6>(init);
    test<1, FIRST, LAGRANGE, QUAD4>(init);
    test<1, SECOND, LAGRANGE, QUAD9>(init);

    test<1, FIRST, MONOMIAL, TRI3>(init);
    test<1, SECOND, MONOMIAL, TRI6>(init);
    test<1, FIRST, MONOMIAL, QUAD4>(init);
    test<1, SECOND, MONOMIAL, QUAD9>(init);

    test<1, FIRST, L2_LAGRANGE, TRI3>(init);
    test<1, SECOND, L2_LAGRANGE, TRI6>(init);
    test<1, FIRST, L2_LAGRANGE, QUAD4>(init);
    test<1, SECOND, L2_LAGRANGE, QUAD9>(init);

    test<1, FIRST, SCALAR, TRI3>(init);
    test<1, SECOND, SCALAR, TRI6>(init);
    test<1, FIRST, SCALAR, QUAD4>(init);
    test<1, SECOND, SCALAR, QUAD9>(init);

    // 3d, libMesh meshes:
    test<2, FIRST, LAGRANGE, TET4>(init);
    test<2, SECOND, LAGRANGE, TET10>(init);
    test<2, FIRST, LAGRANGE, HEX8>(init);
    test<2, SECOND, LAGRANGE, HEX27>(init);

    test<2, FIRST, MONOMIAL, TET4>(init);
    test<2, SECOND, MONOMIAL, TET10>(init);
    test<2, FIRST, MONOMIAL, HEX8>(init);
    test<2, SECOND, MONOMIAL, HEX27>(init);

    test<2, FIRST, L2_LAGRANGE, TET4>(init);
    test<2, SECOND, L2_LAGRANGE, TET10>(init);
    test<2, FIRST, L2_LAGRANGE, HEX8>(init);
    test<2, SECOND, L2_LAGRANGE, HEX27>(init);

    test<2, FIRST, SCALAR, TET4>(init);
    test<2, SECOND, SCALAR, TET10>(init);
    test<2, FIRST, SCALAR, HEX8>(init);
    test<2, SECOND, SCALAR, HEX27>(init);

    // 3d but with read meshes:
    test<2, FIRST, LAGRANGE, TET4>(init, MeshType::readin);
    test<2, SECOND, LAGRANGE, TET10>(init, MeshType::readin);

    test<2, FIRST, MONOMIAL, TET4>(init, MeshType::readin);
    test<2, SECOND, MONOMIAL, TET10>(init, MeshType::readin);

    test<2, FIRST, L2_LAGRANGE, TET4>(init, MeshType::readin);
    test<2, SECOND, L2_LAGRANGE, TET10>(init, MeshType::readin);

    test<2, FIRST, SCALAR, TET4>(init, MeshType::readin);
    test<2, SECOND, SCALAR, TET10>(init, MeshType::readin);

    std::ofstream output("output");
}
