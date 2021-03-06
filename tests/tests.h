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

#include <ibtk/IBTK_MPI.h>

#include <tbox/Logger.h>
#include <tbox/PIO.h>

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

#include <string>
#include <vector>

// A utility function that prints @p out to plog by sending each string to
// rank 0.
inline void
print_strings_on_plog_0(const std::string& out)
{
    using namespace SAMRAI::tbox;
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
print_partitioning_on_plog_0(const Pointer<PatchHierarchy<DIM> >& patch_hierarchy,
                             const int coarsest_ln,
                             const int finest_ln)
{
    std::ostringstream out;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<DIM> > patch_level = patch_hierarchy->getPatchLevel(ln);
        out << "rank: " << IBTK_MPI::getRank() << " level: " << ln << " boxes:\n";
        for (typename PatchLevel<DIM>::Iterator p(patch_level); p; p++)
        {
            const Box<DIM> box = patch_level->getPatch(p())->getBox();
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
    using key_type = quadrature_key_type;

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
    QuadratureCache d_quadrature_cache;

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
    TBOX_ASSERT(static_cast<unsigned int>(get_dim(std::get<0>(quad_key))) == d_dim);
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
        pout << "Program abort called with error message\n\n    " << message.c_str() << std::endl << std::flush;
    }
};

#endif // define included_ibamr_tests_h
