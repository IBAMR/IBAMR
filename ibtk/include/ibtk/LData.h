// Filename: LData.h
// Created on 08 Mar 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_LData
#define included_LData

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "petscvec.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
namespace boost
{
template <typename T, std::size_t NumDims>
class multi_array_ref;
} // namespace boost

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LData provides storage for a single scalar- or vector-valued
 * Lagrangian quantity.
 */
class LData : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Constructor.
     *
     * \note This constructor will allocate an appropriately sized PETSc Vec
     * object.  Data management for this PETSc Vec object is handled by the
     * LData object.
     */
    LData(const std::string& name,
          unsigned int num_local_nodes,
          unsigned int depth,
          const std::vector<int>& nonlocal_petsc_indices = std::vector<int>(0));

    /*!
     * \brief Constructor.
     *
     * \note This constructor \em does \em not allocate a PETSc Vec object.
     * Instead, it can assume responsibilities for data management for the supplied
     * PETSc Vec object.  In particular, the caller \em must \em not destroy the
     * PETSc Vec object provided to the class constructor, if it is instructed to
     * manage the PETSc Vec.
     *
     * \note The blocksize of the supplied PETSc Vec object \em must be set
     * appropriately.  Its value is used to determine the data depth (i.e., the
     * number of data components per node).
     */
    LData(const std::string& name,
          Vec vec,
          const std::vector<int>& nonlocal_petsc_indices = std::vector<int>(0),
          const bool manage_petsc_vec = true);

    /*!
     * \brief Constructor.
     */
    LData(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~LData();

    /*!
     * \brief Reset data items.
     *
     * \note The LData object assumes responsibilities for data management for
     * the supplied PETSc Vec object.  In particular, the caller \em must \em
     * not destroy the PETSc Vec object provided to resetData().
     *
     * \note The blocksize of the supplied PETSc Vec object \em must be set
     * appropriately.  Its value is used to determine the data depth (i.e., the
     * number of data components per node).
     */
    virtual void resetData(Vec vec,
                           const std::vector<int>& nonlocal_petsc_indices = std::vector<int>(0),
                           const bool manage_petsc_vec = true);

    /*!
     * \brief Returns a const reference to the name of this LData object.
     */
    const std::string& getName() const;

    /*!
     * \brief Returns the total number of Lagrangian nodes.
     */
    unsigned int getGlobalNodeCount() const;

    /*!
     * \brief Returns the number of local (i.e., on processor) Lagrangian nodes.
     *
     * \note This count does not include ghost nodes (if any).
     */
    unsigned int getLocalNodeCount() const;

    /*!
     * \brief Returns the number of local ghost nodes associated with the LData
     * object.
     */
    unsigned int getGhostNodeCount() const;

    /*!
     * \brief Returns the depth (i.e., the number of components per node) of the
     * Lagrangian data.
     */
    unsigned int getDepth() const;

    /*!
     * \brief Returns the PETSc Vec object that stores the data.
     *
     * \note getVec() calls restoreArrays(), which invalidates any oustanding
     * references to the underlying array data.
     *
     * \see restoreArrays()
     */
    Vec getVec();

    /*!
     * \brief Returns a \em pointer to a boost::multi_array_ref object that
     * wraps the array corresponding to the PETSc Vec object.  This method is
     * appropriate only for \em scalar-valued quantities.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em global PETSc indexing
     * scheme.
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     *
     * \see restoreArrays()
     */
    boost::multi_array_ref<double, 1>* getArray();

    /*!
     * \brief Returns a \em pointer to a boost::multi_array_ref object that
     * wraps the array corresponding to the PETSc Vec object.  This method is
     * appropriate only for \em scalar-valued quantities.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em local PETSc indexing
     * scheme.
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     *
     * \see restoreArrays()
     */
    boost::multi_array_ref<double, 1>* getLocalFormArray();

    /*!
     * \brief Returns a \em pointer to a boost::multi_array_ref object that
     * wraps the array corresponding to the \em ghosted local part of the PETSc
     * Vec object.  This method is appropriate only for \em scalar-valued
     * quantities.
     *
     * \note Only local (i.e., on processor) and ghost node data are accessible
     * via this array.  All other nonlocal data must be accessed via appropriate
     * PETSc function calls.
     *
     * \note The returned array object is indexed using the \em local PETSc
     * indexing scheme.
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     *
     * \see restoreArrays()
     */
    boost::multi_array_ref<double, 1>* getGhostedLocalFormArray();

    /*!
     * \brief Returns a \em pointer to a boost::multi_array_ref object that
     * wraps the array corresponding to the PETSc Vec object.  This method is
     * appropriate for \em either scalar- or vector-valued quantities.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em global PETSc indexing
     * scheme.
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     *
     * \see restoreArrays()
     */
    boost::multi_array_ref<double, 2>* getVecArray();

    /*!
     * \brief Returns a \em pointer to a boost::multi_array_ref object that
     * wraps the array corresponding to the PETSc Vec object.  This method is
     * appropriate for \em either scalar- or vector-valued quantities.
     *
     * \note Only local data are accessible via the returned array.  Nonlocal
     * data must be accessed via appropriate PETSc function calls.
     *
     * \note The returned array is indexed using the \em local PETSc indexing
     * scheme.
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     *
     * \see restoreArrays()
     */
    boost::multi_array_ref<double, 2>* getLocalFormVecArray();

    /*!
     * \brief Returns a \em pointer to a boost::multi_array_ref object that
     * wraps the array corresponding to the \em ghosted local part of the PETSc
     * Vec object.  This method is appropriate for \em either scalar- or
     * vector-valued quantities.
     *
     * \note Only local (i.e., on processor) and ghost node data are accessible
     * via this array.  All other nonlocal data must be accessed via appropriate
     * PETSc function calls.
     *
     * \note The returned array object is indexed using the \em local PETSc
     * indexing scheme.
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     *
     * \see restoreArrays()
     */
    boost::multi_array_ref<double, 2>* getGhostedLocalFormVecArray();

    /*!
     * \brief Restore any arrays extracted via calls to getArray(),
     * getLocalFormArray(), and getGhostedLocalFormArray().
     *
     * \note Any outstanding references to the underlying array data are
     * invalidated by restoreArrays().
     */
    void restoreArrays();

    /*!
     * \brief Begin updating ghost values.
     */
    void beginGhostUpdate();

    /*!
     * \brief End updating ghost values.
     */
    void endGhostUpdate();

    /*!
     * \brief Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LData();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LData(const LData& from);

    /*
     * Extract the array data.
     */
    void getArrayCommon();
    void getGhostedLocalFormArrayCommon();

    /*
     * The name of the LData object.
     */
    std::string d_name;

    /*
     * The total number of nodes in the mesh.
     */
    unsigned int d_global_node_count;

    /*
     * The number of local nodes in the mesh.
     */
    unsigned int d_local_node_count;

    /*
     * The number of local ghost nodes.
     */
    unsigned int d_ghost_node_count;

    /*
     * The depth (i.e., number of components per node) of the LData.
     */
    unsigned int d_depth;

    /*
     * The nonlocal PETSc indices available in the ghosted local form.
     */
    std::vector<int> d_nonlocal_petsc_indices;

    /*
     * The global PETSc Vec object that contains the mesh data, its underlying
     * array, and a boost::multi_array_ref object that wraps that array.
     */
    Vec d_global_vec;
    bool d_managing_petsc_vec;
    double* d_array;
    boost::multi_array_ref<double, 1> *d_boost_array, *d_boost_local_array;
    boost::multi_array_ref<double, 2> *d_boost_vec_array, *d_boost_local_vec_array;

    /*
     * The array corresponding to the PETSc Vec object in local form, its
     * underyling array, and a boost::multi_array_ref object that wraps that
     * array.
     */
    Vec d_ghosted_local_vec;
    double* d_ghosted_local_array;
    boost::multi_array_ref<double, 1>* d_boost_ghosted_local_array;
    boost::multi_array_ref<double, 2>* d_boost_vec_ghosted_local_array;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LData-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LData
