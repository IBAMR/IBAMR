// Filename: LMeshData.h
// Created on 08 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_LMeshData
#define included_LMeshData

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// PETSc INCLUDES
#include <petscvec.h>

// C++ STDLIB INCLUDES
#include <vector>

// BLITZ++ INCLUDES
#include <blitz/array.h>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LDataManager;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMeshData provides storage for a single scalar- or
 * vector-valued quantity defined at the nodes of the Lagrangian mesh.
 */
class LMeshData
    : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Destructor.
     */
    ~LMeshData();

    /*!
     * \brief Returns a const reference to the name of this LMeshData
     * object.
     */
    const std::string&
    getName() const;

    /*!
     * \brief Returns the total number of nodes of the Lagrangian data.
     */
    int
    getGlobalNodeCount() const;

    /*!
     * \brief Returns the number of local (i.e., on processor) nodes of the
     * Lagrangian data.
     *
     * \note This count does not include ghost nodes.
     */
    int
    getLocalNodeCount() const;

    /*!
     * \brief Returns the number of local ghost nodes of the Lagrangian data.
     */
    int
    getLocalGhostNodeCount() const;

    /*!
     * \brief Returns the depth (i.e., number of components per node) of the
     * Lagrangian data.
     */
    int
    getDepth() const;

    /*!
     * \brief Returns the PETSc Vec object that contains the mesh data.
     *
     * \note getVec() automatically calls restoreArrays().
     *
     * \see restoreArrays()
     */
    Vec&
    getVec();

    /*!
     * \brief Returns a \em reference to a blitz::Array object that wraps the
     * array corresponding to the PETSc Vec object.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em global PETSc indexing
     * scheme.
     */
    blitz::Array<double,1>*
    getArray();

    /*!
     * \brief Returns a \em reference to a blitz::Array object that wraps the
     * array corresponding to the PETSc Vec object.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em local PETSc indexing
     * scheme.
     */
    blitz::Array<double,1>*
    getLocalFormArray();

    /*!
     * \brief Returns a reference to a blitz::Array object that wraps the array
     * corresponding to the ghosted local part of the PETSc Vec object.
     *
     * \note Only local (i.e., on processor) and ghost node data are accessible
     * via this array.  All other nonlocal data must be accessed via appropriate
     * PETSc function calls.
     *
     * \note The returned array object is indexed using the \em local PETSc
     * indexing scheme.
     */
    blitz::Array<double,1>*
    getGhostedLocalFormArray();

    /*!
     * \brief Returns a \em reference to a blitz::Array object that wraps the
     * array corresponding to the PETSc Vec object.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em global PETSc indexing
     * scheme.
     */
    blitz::Array<double,2>*
    getVecArray();

    /*!
     * \brief Returns a \em reference to a blitz::Array object that wraps the
     * array corresponding to the PETSc Vec object.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em local PETSc indexing
     * scheme.
     */
    blitz::Array<double,2>*
    getLocalFormVecArray();

    /*!
     * \brief Returns a reference to a blitz::Array object that wraps the array
     * corresponding to the ghosted local part of the PETSc Vec object.
     *
     * \note Only local (i.e., on processor) and ghost node data are accessible
     * via this array.  All other nonlocal data must be accessed via appropriate
     * PETSc function calls.
     *
     * \note The returned array object is indexed using the \em local PETSc
     * indexing scheme.
     */
    blitz::Array<double,2>*
    getGhostedLocalFormVecArray();

    /*!
     * \brief Restore any arrays extracted via calls to getArray(),
     * getLocalFormArray(), and getGhostedLocalFormArray().
     */
    void
    restoreArrays();

    /*!
     * \brief Begin updating ghost values.
     */
    void
    beginGhostUpdate();

    /*!
     * \brief End updating ghost values.
     */
    void
    endGhostUpdate();

    /*!
     * \brief Write out object state to the given database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    friend class LDataManager;

    /*!
     * \brief Constructor.
     */
    LMeshData(
        const std::string& name,
        const int num_local_nodes,
        const int depth,
        const std::vector<int>& nonlocal_petsc_indices=std::vector<int>(0));

    /*!
     * \brief Constructor.
     */
    LMeshData(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Reset the PETSc Vec and related data.
     */
    void
    resetData(
        Vec& new_global_vec,
        const std::vector<int>& new_nonlocal_petsc_indices);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LMeshData();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LMeshData(
        const LMeshData& from);

    /*
     * Extract the array data.
     */
    void getArrayCommon();
    void getGhostedLocalFormArrayCommon();

    /*
     * The name of the LMeshData object.
     */
    std::string d_name;

    /*
     * The total number of nodes in the mesh.
     */
    int d_global_node_count;

    /*
     * The number of local nodes in the mesh.
     */
    int d_local_node_count;

    /*
     * The number of local ghost nodes.
     */
    int d_local_ghost_node_count;

    /*
     * The depth (i.e., number of components per node) of the LMeshData.
     */
    int d_depth;

    /*
     * The nonlocal PETSc indices available in the ghosted local form.
     */
    std::vector<int> d_nonlocal_petsc_indices;

    /*
     * The global PETSc Vec object that contains the mesh data, its underlying
     * array, and a blitz::Array object that wraps that array.
     */
    Vec d_global_vec;
    double* d_array;
    blitz::Array<double,1> d_blitz_array, d_blitz_local_array;
    blitz::Array<double,2> d_blitz_vec_array, d_blitz_vec_local_array;

    /*
     * The array corresponding to the PETSc Vec object in local form, its
     * underyling array, and a blitz::Array object that wraps that array.
     */
    Vec d_ghosted_local_vec;
    double* d_ghosted_local_array;
    blitz::Array<double,1> d_blitz_ghosted_local_array;
    blitz::Array<double,2> d_blitz_vec_ghosted_local_array;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LMeshData.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LMeshData
