// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_LibMeshSystemIBVectors
#define included_IBTK_LibMeshSystemIBVectors

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/FEDataManager.h>
#include <ibtk/LibMeshSystemVectors.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LibMeshSystemIBVectors is a convenience class that manages
 * access to libMesh vectors for the same system defined on multiple
 * parts. This class supports access to vectors ghosted with the Lagrangian
 * partitioning (i.e., libMesh's computed partitioning) as well as the IB
 * partitioning (i.e., the partitioning based on the distribution of SAMRAI
 * data).
 *
 * A libMesh::System stores vectors in two different ways: the
 * <tt>solution</tt> and <tt>current_local_solution</tt> vectors are stored
 * explicitly by the System object (see libMesh::System::solution and
 * libMesh::System::current_local_solution). They are <em>not</em> accessible
 * by libMesh::System::request_vector() and related functions. All other
 * vectors (including the right-hand side) are stored inside a map and are
 * accessible by <tt>request_vector()</tt>. IBTK::LibMeshSystemIBVectors::get(),
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
 * This class stores information that depends on SAMRAI's parallel
 * partitioning: LibMeshSystemIBVectors::reinit() should be called after
 * regridding to clear this data.
 *
 * @note This class is intended for internal use in IBAMR. It is used in
 * IBFEMethod and IBFESurfaceMethod.
 */
class LibMeshSystemIBVectors : public LibMeshSystemVectors
{
public:
    /*!
     * Constructor.
     *
     * @param[in] fe_data_managers IBTK::FEDataManager objects for each
     * part. These are used to get libMesh vectors with ghost regions
     * corresponding to the Eulerian partitioning.
     *
     * @param[in] system_name Name of the libMesh::System whose vectors we are accessing.
     */
    LibMeshSystemIBVectors(const std::vector<IBTK::FEDataManager*>& fe_data_managers, std::string system_name);

    /*!
     * Constructor, taking a mask argument indicating which systems actually exist.
     *
     * @param[in] fe_data_managers IBTK::FEDataManager objects for each
     * part. These are used to get libMesh vectors with ghost regions
     * corresponding to the Eulerian partitioning.
     *
     * @param[in] part_mask A vector indicating on which parts the given
     * system actually exists (for example, only some parts have fluid
     * sources)
     *
     * @param[in] system_name Name of the libMesh::System whose vectors we are accessing.
     */
    LibMeshSystemIBVectors(const std::vector<IBTK::FEDataManager*>& fe_data_managers,
                           const std::vector<bool>& part_mask,
                           std::string system_name);

    /*!
     * Get an IB-ghosted vector for a specific part.
     *
     * @param[in] vec_name Name of the vector.
     *
     * @param[in] part Part number.
     */
    libMesh::PetscVector<double>& getIBGhosted(const std::string& vec_name, const unsigned int part);

    /**
     * Get, for each part, the IB-ghosted vectors corresponding to the name @p
     * vec_name These vectors are ghosted with values determined by
     * IBTK::FEDataManager (i.e., based on SAMRAI's parallel partitioning).
     *
     * @note These vectors are managed by this object: the ghosted information
     * will not be valid after a regrid.
     */
    std::vector<libMesh::PetscVector<double>*> getIBGhosted(const std::string& vec_name);

    /*!
     * Reinitialize the object. This method must be called after the Eulerian
     * data is updated.
     */
    void reinit() override;

protected:
    /// Add an IB ghosted vector.
    std::vector<std::unique_ptr<libMesh::PetscVector<double> > >& maybeAddIBGhosted(const std::string& vec_name);

    /// Pointers to IBTK::FEDataManager objects. These are needed to get IB ghosted vectors.
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;

    /// Previously computed IB ghosted vectors.
    std::map<std::string, std::vector<std::unique_ptr<libMesh::PetscVector<double> > > > d_ib_ghosted_vectors;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_AppInitializer
