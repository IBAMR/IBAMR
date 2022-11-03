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

// Verify that we can correctly tag cells based on input information to
// FEDataManager so that a single FE mesh ends up on multiple patch levels.

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchyCellDataOpsInteger.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/gmv_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include "../tests.h"

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

class TestTag : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    TestTag(const Parallel::Communicator& comm, Pointer<Database>& input_db) : d_mesh(comm), d_es(d_mesh)
    {
        MeshTools::Generation::build_sphere(d_mesh, 0.5, 4, NDIM == 2 ? TRI3 : HEX8);
        MeshBase::element_iterator el_end = d_mesh.elements_end();
        for (MeshBase::element_iterator el = d_mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* elem = *el;
#if LIBMESH_VERSION_LESS_THAN(1, 7, 0)
            const libMesh::Point centroid = elem->centroid();
#else
            const libMesh::Point centroid = elem->vertex_average();
#endif
#if NDIM == 2
            if (centroid(1) > 0.0)
                elem->subdomain_id() = 1;
            else
                elem->subdomain_id() = 2;
#else

            if (centroid(1) > 0.0 && centroid(2) > 0.0)
                elem->subdomain_id() = 1;
            else
                elem->subdomain_id() = 2;
#endif
        }

        GMVIO gmv_io(d_mesh);
        gmv_io.write("mesh.gmv");

        auto& X_system = setup_deformation_system(d_mesh, d_es, FIRST);

        auto fe_data = std::make_shared<FEData>("fe_data", d_es, true);

        FEDataManager::InterpSpec interp_spec("IB_4", QGAUSS, THIRD, true, 0.9, true, false);
        FEDataManager::SpreadSpec spread_spec("IB_4", QGAUSS, THIRD, true, 0.9, false);
        FEDataManager::WorkloadSpec workload_spec;

        d_fe_data_manager = FEDataManager::getManager(fe_data,
                                                      "fe_data_manager",
                                                      input_db->getDatabase("FEDataManager"),
                                                      input_db->getInteger("MAX_LEVELS"),
                                                      interp_spec,
                                                      spread_spec,
                                                      workload_spec);
        d_fe_data_manager->setCurrentCoordinatesSystemName(X_system.name());
    }

    void initializeLevelData(const Pointer<BasePatchHierarchy<NDIM> > /*hierarchy*/,
                             const int /*level_number*/,
                             const double /*init_data_time*/,
                             const bool /*can_be_refined*/,
                             const bool /*initial_time*/,
                             const tbox::Pointer<hier::BasePatchLevel<NDIM> > /*old_level*/ = nullptr,
                             const bool /*allocate_data*/ = true) override
    {
    }

    void resetHierarchyConfiguration(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > /*hierarchy*/,
                                     const int /*coarsest_level*/,
                                     const int /*finest_level*/) override
    {
    }

    void applyGradientDetector(const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
                               const int level_number,
                               const double error_data_time,
                               const int tag_index,
                               const bool initial_time,
                               const bool uses_richardson_extrapolation_too) override
    {
        d_fe_data_manager->setPatchHierarchy(hierarchy);
        d_fe_data_manager->applyGradientDetector(
            hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    }

protected:
    FEDataManager* d_fe_data_manager;

    ReplicatedMesh d_mesh;

    EquationSystems d_es;
};

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // Parse command line options, set some standard options from the input
    // file, and enable file logging.
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "multilevel_fe_01.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    TestTag tag(init.comm(), input_db);

    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", &tag, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM> > load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    int level_number = 0;
    while (gridding_algorithm->levelCanBeRefined(level_number))
    {
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, true, 0);
        ++level_number;
    }
#if 0
    tbox::Array<int> tag_buffer;
    const int finest_hier_ln = patch_hierarchy->getFinestLevelNumber();
    tag_buffer.resizeArray(finest_hier_ln + 1);
    for (int i = 0; i < tag_buffer.getSize(); ++i)
        tag_buffer[i] = 0;
    gridding_algorithm->regridAllFinerLevels(patch_hierarchy, 0, 0.0, tag_buffer);
#endif

    // plot stuff
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("context");
    Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc");
    const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(1));
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        TBOX_ASSERT(level);
        level->allocatePatchData(u_cc_idx, 0.0);
    }

    // Register variables for plotting.
    Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
    TBOX_ASSERT(visit_data_writer);
    visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "SCALAR", u_cc_idx);
    visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    print_partitioning_on_plog_0(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
}
