// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_LData
#define included_IBTK_LData

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include "petscvec.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <string>
#include <vector>

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

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
    LData(std::string name,
          unsigned int num_local_nodes,
          unsigned int depth,
          std::vector<int> nonlocal_petsc_indices = std::vector<int>(0));

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
    LData(std::string name,
          Vec vec,
          std::vector<int> nonlocal_petsc_indices = std::vector<int>(0),
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
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LData() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LData(const LData& from) = delete;

    /*!
     * Convenience method for doing in-place destruction of array references.
     */
    template <std::size_t dim>
    void destroy_ref(boost::multi_array_ref<double, dim>& ref);

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
    unsigned int d_global_node_count = 0;

    /*
     * The number of local nodes in the mesh.
     */
    unsigned int d_local_node_count = 0;

    /*
     * The number of local ghost nodes.
     */
    unsigned int d_ghost_node_count = 0;

    /*
     * The depth (i.e., number of components per node) of the LData.
     */
    unsigned int d_depth = 0;

    /*
     * The nonlocal PETSc indices available in the ghosted local form.
     */
    std::vector<int> d_nonlocal_petsc_indices;

    /*
     * The global PETSc Vec object that contains the mesh data, its underlying
     * array, and a boost::multi_array_ref object that wraps that array.
     *
     * @note We initialize these with std::vectors here since, with GCC 5 and
     * glibcxx debugging mode enabled, std::array fails to pass a boost
     * concept check.
     */
    Vec d_global_vec = nullptr;
    bool d_managing_petsc_vec = true;
    double* d_array = nullptr;
    boost::multi_array_ref<double, 1> d_boost_array{ nullptr, std::vector<int>{ 0 } };
    boost::multi_array_ref<double, 1> d_boost_local_array{ nullptr, std::vector<int>{ 0 } };
    boost::multi_array_ref<double, 2> d_boost_vec_array{ nullptr, std::vector<int>{ 0, 0 } };
    boost::multi_array_ref<double, 2> d_boost_local_vec_array{ nullptr, std::vector<int>{ 0, 0 } };

    /*
     * The array corresponding to the PETSc Vec object in local form, its
     * underyling array, and a boost::multi_array_ref object that wraps that
     * array.
     */
    Vec d_ghosted_local_vec = nullptr;
    double* d_ghosted_local_array = nullptr;
    boost::multi_array_ref<double, 1> d_boost_ghosted_local_array{ nullptr, std::vector<int>{ 0 } };
    boost::multi_array_ref<double, 2> d_boost_vec_ghosted_local_array{ nullptr, std::vector<int>{ 0, 0 } };
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LData-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LData
