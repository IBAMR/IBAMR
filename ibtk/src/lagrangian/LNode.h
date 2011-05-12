// Filename: LNode.h
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

#ifndef included_LNode
#define included_LNode

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeIndex.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LNode is the basic element of an LMesh.
 *
 * Class LNode provides Lagrangian and <A
 * HREF="http://www-unix.mcs.anl.gov/petsc">PETSc</A> indexing information and
 * data storage for a single node of a Lagrangian mesh.
 */
class LNode
    : public LNodeIndex
{
public:
    /*!
     * \brief Default constructor.
     */
    LNode(
        const int lagrangian_nidx=-1,
        const int global_petsc_nidx=-1,
        const int local_petsc_nidx=-1,
        const double* const X_ptr=NULL,
        const SAMRAI::hier::IntVector<NDIM>& periodic_offset=SAMRAI::hier::IntVector<NDIM>(0),
        const std::vector<double>& periodic_displacement=std::vector<double>(NDIM,0.0),
        const std::vector<SAMRAI::tbox::Pointer<Streamable> >& node_data=std::vector<SAMRAI::tbox::Pointer<Streamable> >());

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     */
    LNode(
        const LNode& from);

    /*!
     * \brief Virtuald destructor.
     *
     * The LNode destructor does nothing interesting.
     */
    virtual ~LNode();

    /*!
     * \brief Assignment operator.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LNode&
    operator=(
        const LNode& that);

private:
};

/*!
 * \brief Less-than comparison operator.
 *
 * \return Whether lhs < rhs.
 *
 * The ordering is determined on the PETSc indexing of the nodes.
 */
bool
operator<(
    const LNode& lhs,
    const LNode& rhs);

}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/LNode.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LNode
