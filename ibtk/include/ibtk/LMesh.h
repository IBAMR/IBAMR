// Filename: LMesh.h
// Created on 05 May 2011 by Boyce Griffith
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

#ifndef included_LMesh
#define included_LMesh

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "tbox/DescribedClass.h"

namespace IBTK
{
class LNode;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMesh is a collection of LNode objects.
 */
class LMesh : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    LMesh(const std::string& object_name,
          const std::vector<LNode*>& local_nodes,
          const std::vector<LNode*>& ghost_nodes);

    /*!
     * \brief Destructor.
     */
    ~LMesh();

    /*!
     * \brief Return a const reference to the set of local LNode objects.
     */
    const std::vector<LNode*>& getLocalNodes() const;

    /*!
     * \brief Return a const reference to the set of local ghost LNode objects.
     */
    const std::vector<LNode*>& getGhostNodes() const;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LMesh(const LMesh& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMesh& operator=(const LMesh& that);

    const std::string& d_object_name;
    const std::vector<LNode*> d_local_nodes;
    const std::vector<LNode*> d_ghost_nodes;
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LMesh-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LMesh
