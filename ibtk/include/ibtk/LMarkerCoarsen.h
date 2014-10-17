// Filename: LMarkerCoarsen.h
// Created on 13 Sep 2007 by Boyce Griffith
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

#ifndef included_LMarkerCoarsen
#define included_LMarkerCoarsen

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "Box.h"
#include "CoarsenOperator.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LMarkerCoarsen is a concrete SAMRAI::xfer::CoarsenOperator for
 * restricting IB marker data from finer levels to coarser levels in the patch
 * hierarchy.
 */
class LMarkerCoarsen : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    LMarkerCoarsen();

    /*!
     * \brief Destructor.
     */
    ~LMarkerCoarsen();

    /*!
     * \name Implementation of SAMRAI::xfer::CoarsenOperator interface.
     */
    //\{

    /*!
     * Return true if the coarsening operation matches the variable and name
     * string identifier request; false, otherwise.
     */
    bool findCoarsenOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
                             const std::string& op_name) const;

    /*!
     * Return name string identifier of the coarsening operation.
     */
    const std::string& getOperatorName() const;

    /*!
     * Return the priority of this operator relative to other coarsening
     * operators.  The SAMRAI transfer routines guarantee that coarsening using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    int getOperatorPriority() const;

    /*!
     * Return the stencil width associated with the coarsening operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each coarsening operator.
     */
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const;

    /*!
     * Coarsen the source component on the fine patch to the destination
     * component on the coarse patch. The coarsening operation is performed on
     * the intersection of the destination patch and the coarse box.  The fine
     * patch is guaranteed to contain sufficient data for the stencil width of
     * the coarsening operator.
     */
    void coarsen(SAMRAI::hier::Patch<NDIM>& coarse,
                 const SAMRAI::hier::Patch<NDIM>& fine,
                 int dst_component,
                 int src_component,
                 const SAMRAI::hier::Box<NDIM>& coarse_box,
                 const SAMRAI::hier::IntVector<NDIM>& ratio) const;

    //\}

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LMarkerCoarsen(const LMarkerCoarsen& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LMarkerCoarsen& operator=(const LMarkerCoarsen& that);

    /*!
     * The operator name.
     */
    static const std::string s_op_name;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LMarkerCoarsen
