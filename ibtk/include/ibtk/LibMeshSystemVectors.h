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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_LibMeshSystemVectors
#define included_IBTK_LibMeshSystemVectors

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include <libmesh/equation_systems.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/system.h>

#include <map>
#include <memory>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LibMeshSystemVectors is a convenience class that manages
 * access to libMesh vectors for the same system defined on multiple
 * parts. This class only supports access to vectors ghosted with the
 * Lagrangian partitioning (i.e., libMesh's computed partitioning).  The
 * subclass LibMeshSystemIBVectors provides access to the IB partitioning
 * (i.e., the partitioning based on the distribution of SAMRAI data).
 *
 * A libMesh::System stores vectors in two different ways: the
 * <tt>solution</tt> and <tt>current_local_solution</tt> vectors are stored
 * explicitly by the System object (see libMesh::System::solution and
 * libMesh::System::current_local_solution). They are <em>not</em> accessible
 * by libMesh::System::request_vector() and related functions. All other
 * vectors (including the right-hand side) are stored inside a map and are
 * accessible by <tt>request_vector()</tt>. IBTK::LibMeshSystemVectors::get(),
 * for the sake of a uniform interface, is inconsistent with libMesh in this
 * regard: the vector name <tt>"solution"</tt> refers to the solution vector
 * (which is not stored in the map). The same is true for
 * <tt>"current_local_solution"</tt>. For convenience, this class also accepts
 * <tt>"current"</tt> as an alias for <tt>"current_local_solution"</tt>. All
 * other vector names are stored in the libMesh::System vector map.
 *
 * All vectors added to a libMesh::System by this class are ghosted (with
 * ghost data computed by the System) and registered for projection (i.e., if
 * the mesh changes, these vectors will be automatically updated).
 *
 * @note This class is intended for internal use in IBAMR. It is used in
 * IBFEMethod and IBFESurfaceMethod.
 */
class LibMeshSystemVectors
{
public:
    /*!
     * Constructor.
     *
     * @param[in] equation_systems libMesh::EquationSystems objects for each
     * part.
     *
     * @param[in] part_mask A vector indicating on which parts the given
     * system actually exists (for example, only some parts have fluid
     * sources)
     *
     * @param[in] system_name Name of the libMesh::System whose vectors we are accessing.
     */
    LibMeshSystemVectors(const std::vector<libMesh::EquationSystems*>& equation_systems,
                         const std::vector<bool>& part_mask,
                         std::string system_name);

    /*!
     * Constructor.
     *
     * @param[in] equation_systems libMesh::EquationSystems objects for each
     * part.
     *
     * @param[in] system_name Name of the libMesh::System whose vectors we are accessing.
     */
    LibMeshSystemVectors(const std::vector<libMesh::EquationSystems*>& equation_systems, std::string system_name);

    /*!
     * Destructor. Made virtual to aid inheritance.
     */
    virtual ~LibMeshSystemVectors() = default;

    /*!
     * Get a specific vector for a specific part. These vectors are managed by
     * libMesh::System objects.
     *
     * @param[in] vec_name Name of the vector - see the note in the general
     * class description for an explanation of the inconsistency of this term
     * versus what libMesh expects.
     *
     * @param[in] part Part number.
     */
    libMesh::PetscVector<double>& get(const std::string& vec_name, const unsigned int part);

    /*!
     * Get, for each part, the vector corresponding to the name @p
     * vec_name. These vectors are managed by libMesh::System objects.
     *
     * @param[in] vec_name Name of the vector - see the note in the general
     * class description for an explanation of the inconsistency of this term
     * versus what libMesh expects.
     */
    std::vector<libMesh::PetscVector<double>*> get(const std::string& vec_name);

    /*!
     * Reinitialize the object.
     */
    virtual void reinit();

    /*!
     * Convenience function for copying libMesh vectors.
     *
     * @param[in] source Name of the input vectors. These should already be
     * present in the libMesh::System.
     *
     * @param[in] dests Names of the output vectors.
     */
    void copy(const std::string& source, const std::vector<std::string>& dests);

    /*!
     * Zero vectors owned by a libMesh::System.
     */
    void zero(const std::string& vec_name);

protected:
    /// Returns false if a standard libMesh vector; otherwise returns true.
    inline static bool vec_stored_in_map(const std::string& vec_name)
    {
        if (vec_name == "solution" || vec_name == "current" || vec_name == "current_local_solution") return false;
        return true;
    }

    /// Check if a vector exists and, if not, add it.
    void maybeAdd(const std::string& vec_name);

    /// Mask indicating which parts actually have a system with the given name.
    std::vector<bool> d_part_mask;

    /// Name of the system.
    const std::string d_system_name;

    /// libMesh::System objects.
    std::vector<libMesh::System*> d_systems;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_AppInitializer
