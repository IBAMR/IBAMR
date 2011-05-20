// Filename: LMesh.h
// Created on 05 May 2011 by Boyce Griffith
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

#ifndef included_LMesh
#define included_LMesh

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNode.h>

// SAMRAI INCLUDES
#include <tbox/DescribedClass.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMesh is a collection of LNode objects.
 */
class LMesh
    : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    LMesh(
        const std::string& object_name);

    /*!
     * \brief Destructor.
     */
    ~LMesh();

    /*!
     * \brief Set the local LNode objects.
     */
    void
    setNodes(
        const std::vector<LNode*>& nodes,
        const bool sorted=false);

    /*!
     * \brief Return a const reference to the set of local LNode objects.
     */
    const std::vector<LNode*>&
    getNodes() const;

    /*!
     * \brief Return a const reference to the set of Lagrangian indices for the
     * local LNode objects.
     */
    const std::vector<int>&
    getLagrangianIndices() const;

    /*!
     * \brief Return a const reference to the set of global PETSc indices for
     * the local LNode objects.
     */
    const std::vector<int>&
    getGlobalPETScIndices() const;

    /*!
     * \brief Return a const reference to the set of local PETSc indices for the
     * local LNode objects.
     */
    const std::vector<int>&
    getLocalPETScIndices() const;

    /*!
     * \brief Set the local ghost LNode objects.
     */
    void
    setGhostNodes(
        const std::vector<LNode*>& ghost_nodes,
        const bool sorted=false);

    /*!
     * \brief Return a const reference to the set of local ghost LNode objects.
     */
    const std::vector<LNode*>&
    getGhostNodes() const;

    /*!
     * \brief Return a const reference to the set of Lagrangian indices for the
     * local ghost LNode objects.
     */
    const std::vector<int>&
    getGhostLagrangianIndices() const;

    /*!
     * \brief Return a const reference to the set of global PETSc indices for
     * the local ghost LNode objects.
     */
    const std::vector<int>&
    getGhostGlobalPETScIndices() const;

    /*!
     * \brief Return a const reference to the set of local PETSc indices for the
     * local ghost LNode objects.
     */
    const std::vector<int>&
    getGhostLocalPETScIndices() const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LMesh(
        const LMesh& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMesh&
    operator=(
        const LMesh& that);

    const std::string& d_object_name;
    std::vector<LNode*> d_nodes;
    std::vector<int> d_lag_idxs, d_global_petsc_idxs, d_local_petsc_idxs;
    std::vector<LNode*> d_ghost_nodes;
    std::vector<int> d_ghost_lag_idxs, d_ghost_global_petsc_idxs, d_ghost_local_petsc_idxs;
};

}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LMesh.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LMesh
