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

// Collection of utility functions that are useful in tests.

#ifndef included_ibamr_tests_h
#define included_ibamr_tests_h

#include <ibamr/config.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartSideDoubleSpecializedLinearRefine.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>
#include <tbox/PIO.h>
#include <tbox/SAMRAI_MPI.h>

#include <BergerRigoutsos.h>
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
IBTK_DISABLE_EXTRA_WARNINGS
#include <StandardTagAndInitialize.h>
IBTK_ENABLE_EXTRA_WARNINGS

#ifdef IBAMR_HAVE_LIBMESH
#include <ibtk/QuadratureCache.h>
#include <ibtk/libmesh_utilities.h>

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/fe_map.h>
#include <libmesh/quadrature.h>
#endif

#include <mpi.h>

#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// A utility function that prints @p out to plog by sending each string to
// rank 0.
inline void
print_strings_on_plog_0(const std::string& out)
{
    using namespace SAMRAI::tbox;
    using namespace IBTK;
    const int n_nodes = IBTK_MPI::getNodes();
    std::vector<unsigned long> string_sizes(n_nodes);

    const unsigned long size = out.size();
    int ierr = MPI_Gather(
        &size, 1, MPI_UNSIGNED_LONG, string_sizes.data(), 1, MPI_UNSIGNED_LONG, 0, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);

    // MPI_Gatherv would be more efficient, but this just a test so its
    // not too important
    if (IBTK_MPI::getRank() == 0)
    {
        plog << out;
        for (int r = 1; r < n_nodes; ++r)
        {
            std::string input;
            input.resize(string_sizes[r]);
            ierr = MPI_Recv(&input[0], string_sizes[r], MPI_CHAR, r, 0, IBTK_MPI::getCommunicator(), MPI_STATUS_IGNORE);
            TBOX_ASSERT(ierr == 0);
            plog << input;
        }
    }
    else
        MPI_Send(out.data(), size, MPI_CHAR, 0, 0, IBTK_MPI::getCommunicator());
}

/**
 * Print the parallel partitioning (i.e., the boxes) on all processes to plog
 * stream on processor 0.
 */
template <int DIM>
inline void
print_partitioning_on_plog_0(const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<DIM> >& patch_hierarchy,
                             const int coarsest_ln,
                             const int finest_ln)
{
    using namespace SAMRAI;
    using namespace IBTK;
    std::ostringstream out;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        tbox::Pointer<hier::PatchLevel<DIM> > patch_level = patch_hierarchy->getPatchLevel(ln);
        out << "rank: " << IBTK_MPI::getRank() << " level: " << ln << " boxes:\n";
        for (typename hier::PatchLevel<DIM>::Iterator p(patch_level); p; p++)
        {
            const hier::Box<DIM> box = patch_level->getPatch(p())->getBox();
            out << box << '\n';
        }
    }
    print_strings_on_plog_0(out.str());
}

#ifdef IBAMR_HAVE_LIBMESH
/**
 * \brief Class storing multiple libMesh::FEMap objects, each corresponding to
 * a different quadrature rule. Each FEMap object is configured with a
 * quadrature rule corresponding to the provided <code>quad_key</code>
 * parameter.
 *
 * This class essentially provides a wrapper around std::map to manage FEMap
 * objects and the quadrature rules they use. The keys are descriptions of
 * quadrature rules.
 *
 * This class used to be part of IBTK but it is no longer needed since we
 * reimplemented all the mappings ourselves. However, it is still useful in the
 * test suite so that we can easily compare results against what libMesh
 * calculates for a given quadrature key.
 */
class FEMapCache
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = IBTK::quadrature_key_type;

    /**
     * Type of values stored by this class that are accessible through
     * <code>operator[]</code>.
     */
    using value_type = libMesh::FEMap;

    /**
     * Constructor. Sets up a cache of FE objects calculating values for the
     * given FEType argument. All cached FE objects have the same FEType.
     *
     * @param dim The topological dimension of the relevant libMesh::Mesh: see
     * libMesh::MeshBase::mesh_dimension() for more information.
     */
    FEMapCache(const unsigned int dim);

    /**
     * Return a reference to an FEMap object that matches the specified
     * quadrature rule type and order.
     *
     * @param quad_key a tuple of enums that completely describes
     * a libMesh quadrature rule.
     */
    libMesh::FEMap& operator[](const key_type& quad_key);

protected:
    /**
     * Topological dimension of the FE mesh.
     */
    const unsigned int d_dim;

    /**
     * Managed libMesh::Quadrature objects. These are used to partially
     * initialize (i.e., points but not weights are stored) the FEMap objects.
     */
    IBTK::QuadratureCache d_quadrature_cache;

    /**
     * Managed libMesh::FEMap objects of specified dimension and family.
     */
    std::map<key_type, libMesh::FEMap> d_fe_maps;
};

inline FEMapCache::FEMapCache(const unsigned int dim) : d_dim(dim), d_quadrature_cache(d_dim)
{
}

inline libMesh::FEMap&
FEMapCache::operator[](const FEMapCache::key_type& quad_key)
{
    TBOX_ASSERT(static_cast<unsigned int>(IBTK::get_dim(std::get<0>(quad_key))) == d_dim);
    auto it = d_fe_maps.find(quad_key);
    if (it == d_fe_maps.end())
    {
        libMesh::QBase& quad = d_quadrature_cache[quad_key];
        libMesh::FEMap& fe_map = d_fe_maps[quad_key];
        // Calling this function enables JxW calculations
        fe_map.get_JxW();

        // Doing this may not work with future (1.4.0 or up) versions of
        // libMesh. In particular; init_reference_to_physical_map is
        // undocumented (and almost surely is not intended for use by anyone
        // but libMesh developers) and *happens* to not read any geometric or
        // topological information from the Elem argument (just the default
        // order and type).
        const libMesh::ElemType elem_type = std::get<0>(quad_key);
        std::unique_ptr<libMesh::Elem> exemplar_elem(libMesh::Elem::build(elem_type));

        // This is one of very few functions in libMesh that is templated on
        // the dimension (not spatial dimension) of the mesh
        switch (d_dim)
        {
        case 1:
            fe_map.init_reference_to_physical_map<1>(quad.get_points(), exemplar_elem.get());
            break;
        case 2:
            fe_map.init_reference_to_physical_map<2>(quad.get_points(), exemplar_elem.get());
            break;
        case 3:
            fe_map.init_reference_to_physical_map<3>(quad.get_points(), exemplar_elem.get());
            break;
        default:
            TBOX_ASSERT(false);
        }
        return fe_map;
    }
    else
    {
        return it->second;
    }
}
#endif // ifdef IBAMR_HAVE_LIBMESH

// utility function for getting the number of components in the test subdatabase
// of the input database.
int
get_n_f_components(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    auto test_db = input_db->getDatabase("test");
    int n_f_components = 0;
    if (test_db->keyExists("f"))
        n_f_components = test_db->getDatabase("f")->getAllKeys().getSize();
    else
        n_f_components = test_db->getIntegerWithDefault("n_components", 1);

    return n_f_components;
}

// A utility function that does the normal SAMRAI initialization stuff.
template <int spacedim>
std::tuple<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<spacedim> >,
           SAMRAI::tbox::Pointer<SAMRAI::mesh::StandardTagAndInitialize<spacedim> >,
           SAMRAI::tbox::Pointer<SAMRAI::mesh::BergerRigoutsos<spacedim> >,
           SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<spacedim> >,
           SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<spacedim> >,
           int>
setup_hierarchy(SAMRAI::tbox::Pointer<IBTK::AppInitializer> app_initializer,
                SAMRAI::mesh::StandardTagAndInitStrategy<spacedim>* tag = nullptr)
{
    using namespace SAMRAI;

    auto input_db = app_initializer->getInputDatabase();

    // database specific to tests, should it exist
    auto test_db = input_db;
    if (input_db->keyExists("test")) test_db = input_db->getDatabase("test");

    // Set up basic SAMRAI stuff:
    tbox::Pointer<geom::CartesianGridGeometry<spacedim> > grid_geometry =
        new geom::CartesianGridGeometry<spacedim>("CartesianGeometry", input_db->getDatabase("CartesianGeometry"));

    // More questionable SAMRAI design decisions - we have to register refine
    // operations associated with different variables with the grid geometry class
    grid_geometry->addSpatialRefineOperator(new geom::CartesianCellDoubleLinearRefine<spacedim>());
    grid_geometry->addSpatialRefineOperator(new IBTK::CartSideDoubleSpecializedLinearRefine());

    tbox::Pointer<hier::PatchHierarchy<spacedim> > patch_hierarchy =
        new hier::PatchHierarchy<spacedim>("PatchHierarchy", grid_geometry);
    tbox::Pointer<mesh::StandardTagAndInitialize<spacedim> > error_detector =
        new mesh::StandardTagAndInitialize<spacedim>(
            "StandardTagAndInitialize", tag, input_db->getDatabase("StandardTagAndInitialize"));

    tbox::Pointer<mesh::BergerRigoutsos<spacedim> > box_generator = new mesh::BergerRigoutsos<spacedim>();
    tbox::Pointer<mesh::LoadBalancer<spacedim> > load_balancer =
        new mesh::LoadBalancer<spacedim>("LoadBalancer", input_db->getDatabase("LoadBalancer"));
    tbox::Pointer<mesh::GriddingAlgorithm<spacedim> > gridding_algorithm = new mesh::GriddingAlgorithm<spacedim>(
        "GriddingAlgorithm", input_db->getDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

    // Set up a variable so that we can actually output the grid. Note that this
    // has to happen before we make any levels since this is where we set the
    // maximum ghost width.
    auto* var_db = hier::VariableDatabase<spacedim>::getDatabase();
    tbox::Pointer<hier::VariableContext> ctx = var_db->getContext("context");

    const int n_f_components = get_n_f_components(input_db);

    tbox::Pointer<hier::Variable<spacedim> > f_var;
    const std::string f_data_type = test_db->getStringWithDefault("f_data_type", "CELL");
    if (f_data_type == "CELL")
        f_var = new pdat::CellVariable<spacedim, double>("f_cc", n_f_components);
    else if (f_data_type == "SIDE")
    // This one is different since side-centered data already has implicitly
    // spacedim 'depth' in a different sense
    //
    // TODO this should probably be clarified somehow
    {
        TBOX_ASSERT(n_f_components == spacedim);
        f_var = new pdat::SideVariable<spacedim, double>("f_sc", 1);
    }
    else if (f_data_type == "NODE")
        f_var = new pdat::SideVariable<spacedim, double>("f_nc", n_f_components);
    else
        TBOX_ASSERT(false);

    // BSPLINE-3 needs 3 ghost cells
    const int ghost_width = test_db->getIntegerWithDefault("ghost_width", 3);
    const int f_idx = var_db->registerVariableAndContext(f_var, ctx, hier::IntVector<spacedim>(ghost_width));

    // set up the SAMRAI grid:
    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    int level_number = 0;
    const int tag_buffer = 1;
    while (gridding_algorithm->levelCanBeRefined(level_number))
    {
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
        ++level_number;
    }

    const int finest_level = patch_hierarchy->getFinestLevelNumber();
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        tbox::Pointer<hier::PatchLevel<spacedim> > level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(f_idx, 0.0);
    }

    auto visit_data_writer = app_initializer->getVisItDataWriter();
    TBOX_ASSERT(visit_data_writer);
    // non-cell data cannot be plotted
    int plot_cc_idx = 0;
    tbox::Pointer<hier::Variable<spacedim> > plot_cc_var;
    if (f_data_type == "CELL")
        for (int d = 0; d < n_f_components; ++d)
            visit_data_writer->registerPlotQuantity(f_var->getName() + std::to_string(d), "SCALAR", f_idx, d);
    else
    {
        plot_cc_var = new pdat::CellVariable<spacedim, double>("f_cc", n_f_components);
        plot_cc_idx = var_db->registerVariableAndContext(plot_cc_var, ctx, hier::IntVector<spacedim>(0));

        for (int d = 0; d < n_f_components; ++d)
        {
            visit_data_writer->registerPlotQuantity(
                plot_cc_var->getName() + std::to_string(d), "SCALAR", plot_cc_idx, d);
        }
    }

    // If it exists, set up an initial condition
    tbox::Pointer<tbox::Database> f_db;
    if (test_db->keyExists("f"))
        f_db = test_db->getDatabase("f");
    else
    {
        f_db = new tbox::MemoryDatabase("f_db");
        if (n_f_components == 1)
        {
            f_db->putString("function", "0.0");
        }
        else
        {
            for (int i = 0; i < n_f_components; ++i)
            {
                f_db->putString("function" + std::to_string(i), "0.0");
            }
        }
    }

    {
        IBTK::muParserCartGridFunction f_fcn("f", f_db, patch_hierarchy->getGridGeometry());
        f_fcn.setDataOnPatchHierarchy(f_idx, f_var, patch_hierarchy, 0.0);

        using ITC = IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_components(1);
        const std::string refine_type = f_data_type == "CELL" ? "LINEAR_REFINE" : "SPECIALIZED_LINEAR_REFINE";
        ghost_cell_components[0] = ITC(f_idx,
                                       refine_type,
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // f_bc_coefs
                                       nullptr);
        IBTK::HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);
    }

    return std::make_tuple(patch_hierarchy, error_detector, box_generator, load_balancer, gridding_algorithm, f_idx);
}

/**
 * SAMRAI's default abort error message printer prints some things that are not
 * compatible with the test suite (the line number and full path to the file).
 * Get rid of them as needed in the test suite by replacing its abort appender
 * with our own (see SAMRAI::tbox::Logger::setAbortAppender)
 */
class TestAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    virtual void logMessage(const std::string& message, const std::string&, const int) override
    {
        // SAMRAI appends an extra NUL to the message string so we have to
        // convert it to a C string to prevent it from being printed
        SAMRAI::tbox::pout << "Program abort called with error message\n\n    " << message.c_str() << std::endl
                           << std::flush;
    }
};

#endif // define included_ibamr_tests_h
