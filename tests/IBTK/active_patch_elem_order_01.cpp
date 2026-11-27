// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Verify that FEDataManager stores the active elements of each patch in
// increasing id order, even when the allocation order of the elements (and so
// their memory addresses) is the reverse of their id order.

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/FEDataManager.h>
#include <ibtk/IBTKInit.h>

// Set up application namespace declarations
#include "../tests.h"

#include <ibamr/app_namespaces.h>

class TestInit : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    void initializeLevelData(const Pointer<BasePatchHierarchy<NDIM>> /*hierarchy*/,
                             const int /*level_number*/,
                             const double /*init_data_time*/,
                             const bool /*can_be_refined*/,
                             const bool /*initial_time*/,
                             const tbox::Pointer<hier::BasePatchLevel<NDIM>> /*old_level*/ = nullptr,
                             const bool /*allocate_data*/ = true) override
    {
    }

    void resetHierarchyConfiguration(const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                     const int /*coarsest_level*/,
                                     const int /*finest_level*/) override
    {
    }

    void applyGradientDetector(const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> /*hierarchy*/,
                               const int /*level_number*/,
                               const double /*error_data_time*/,
                               const int /*tag_index*/,
                               const bool /*initial_time*/,
                               const bool /*uses_richardson_extrapolation_too*/) override
    {
    }
};

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "active_patch_elem_order_01.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Build an n x n mesh of QUAD4 elements on [-0.55, 0.45] x [-0.45, 0.55].
    // No element edge lies on a patch boundary grown by one cell, so each patch
    // gets a different subset of the elements. The elements are created in
    // decreasing id order, so their addresses (usually) decrease as their ids
    // increase.
    const unsigned int n = 4;
    ReplicatedMesh mesh(init.comm(), NDIM);
    mesh.allow_renumbering(false);
    for (unsigned int j = 0; j <= n; ++j)
    {
        for (unsigned int i = 0; i <= n; ++i)
        {
            mesh.add_point(libMesh::Point(-0.55 + static_cast<double>(i) / n, -0.45 + static_cast<double>(j) / n),
                           j * (n + 1) + i);
        }
    }
    for (int e = static_cast<int>(n * n) - 1; e >= 0; --e)
    {
        const unsigned int i = e % n;
        const unsigned int j = e / n;
        std::unique_ptr<Elem> elem = Elem::build(QUAD4);
        elem->set_id(e);
        elem->set_node(0) = mesh.node_ptr(j * (n + 1) + i);
        elem->set_node(1) = mesh.node_ptr(j * (n + 1) + i + 1);
        elem->set_node(2) = mesh.node_ptr((j + 1) * (n + 1) + i + 1);
        elem->set_node(3) = mesh.node_ptr((j + 1) * (n + 1) + i);
        mesh.add_elem(std::move(elem));
    }
    mesh.prepare_for_use();

    EquationSystems equation_systems(mesh);
    auto& X_system = equation_systems.add_system<ExplicitSystem>("X");
    for (unsigned int d = 0; d < NDIM; ++d) X_system.add_variable("X_" + std::to_string(d), FIRST, LAGRANGE);
    equation_systems.init();
    for (const Node* const node : mesh.local_node_ptr_range())
    {
        for (unsigned int d = 0; d < NDIM; ++d)
            X_system.solution->set(node->dof_number(X_system.number(), d, 0), (*node)(d));
    }
    X_system.solution->close();

    auto fe_data = std::make_shared<FEData>("fe_data", equation_systems, true);
    FEDataManager::InterpSpec interp_spec("IB_4", QGAUSS, THIRD, true, 0.9, true, false);
    FEDataManager::SpreadSpec spread_spec("IB_4", QGAUSS, THIRD, true, 0.9, false);
    FEDataManager::WorkloadSpec workload_spec;
    FEDataManager* fe_data_manager = FEDataManager::getManager(fe_data,
                                                               "fe_data_manager",
                                                               input_db->getDatabase("FEDataManager"),
                                                               input_db->getInteger("MAX_LEVELS"),
                                                               interp_spec,
                                                               spread_spec,
                                                               workload_spec);
    fe_data_manager->setCurrentCoordinatesSystemName(X_system.name());

    // Set up a single-level hierarchy with several patches.
    TestInit test_init;
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", &test_init, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);
    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

    fe_data_manager->setPatchHierarchy(patch_hierarchy);
    fe_data_manager->reinitElementMappings();

    // Print the element ids of each patch in the order that the map stores them.
    const std::vector<std::vector<Elem*>>& patch_elems = fe_data_manager->getActivePatchElementMap();
    for (std::size_t patch_n = 0; patch_n < patch_elems.size(); ++patch_n)
    {
        plog << "patch " << patch_n << " element ids:";
        for (const Elem* const elem : patch_elems[patch_n]) plog << " " << elem->id();
        plog << "\n";
    }
}
