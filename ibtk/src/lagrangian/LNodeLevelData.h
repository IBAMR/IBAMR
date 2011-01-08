// Filename: LNodeLevelData.h
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

#ifndef included_LNodeLevelData
#define included_LNodeLevelData

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// PETSc INCLUDES
#include <petscvec.h>
#include <petscao.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class LDataManager;
}// namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNodeLevelData provides storage for a single variable defined at
 * the nodes of the Lagrangian mesh.
 */
class LNodeLevelData
    : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Destructor.
     */
    ~LNodeLevelData();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNodeLevelData&
    operator=(
        const LNodeLevelData& that);

    /*!
     * \brief This operator() method provides a non-const reference to the
     * Lagrangian nodal quantity at local index i at the specified depth.
     *
     * \return A non-const reference to the Lagrangian nodal quantity at local
     * index i at the specified depth.
     */
    double&
    operator()(
        const int i,
        const int depth=0);

    /*!
     * \return The total number of nodes of the Lagrangian data.
     */
    int
    getGlobalNodeCount();

    /*!
     * \return The number of local (i.e., on processor) nodes of the Lagrangian
     * data.
     *
     * \note This count does not include ghost nodes.
     */
    int
    getLocalNodeCount();

    /*!
     * \return The number of local ghost nodes of the Lagrangian data.
     */
    int
    getLocalGhostNodeCount();

    /*!
     * \return The depth (i.e., number of components per node) of the Lagrangian
     * data.
     */
    int
    getDepth() const;

    /*!
     * \brief Begin updating the ghost nodes.
     */
    void
    beginGhostUpdate();

    /*!
     * \brief End updating the ghost nodes.
     */
    void
    endGhostUpdate();

    /*!
     * \return The global PETSc Vec object which contains the level data.
     */
    Vec&
    getGlobalVec();

    /*!
     * \return The ghosted local form of the global PETSc Vec object which
     * contains the level data.
     */
    Vec&
    getLocalFormVec();

    /*!
     * \brief Restore the ghosted local form of the global PETSc Vec object.
     */
    void
    restoreLocalFormVec();

    /*!
     * \return The array corresponding to the PETSc Vec object in ghosted local
     * form.
     */
    PetscScalar*
    getLocalFormArray();

    /*!
     * \brief Restore the array to the ghosted local form of the global PETSc
     * Vec object.
     */
    void
    restoreLocalFormArray();

    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///

    /*!
     * Store the contents of the class in a database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    friend class LDataManager;

    /*!
     * \brief Constructor.
     */
    LNodeLevelData(
        const std::string& name,
        const int num_local_nodes,
        const int depth,
        const std::vector<int>& nonlocal_petsc_indices=std::vector<int>(0));

    /*!
     * \brief Constructor.
     */
    LNodeLevelData(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNodeLevelData(
        const LNodeLevelData& from);

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
    LNodeLevelData();

    /*
     * The name of the LNodeLevelData.
     */
    std::string d_name;

    /*
     * The depth (i.e. number of components per node) of the LNodeLevelData.
     */
    int d_depth;

    /*
     * The nonlocal PETSc indices available in the ghosted local form.
     */
    std::vector<int> d_nonlocal_petsc_indices;

    /*
     * The global PETSc Vec object which contains the level data and its ghosted
     * local form.
     */
    Vec d_global_vec, d_local_vec;
    bool d_in_local_form;

    /*
     * The array corresponding to the PETSc Vec object in local form.
     */
    PetscScalar* d_local_vec_array;
    bool d_extracted_local_array;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LNodeLevelData.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNodeLevelData
