// Filename: CartSideDoubleCubicCoarsen.h
// Created on 02 May 2008 by Boyce Griffith
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

#ifndef included_CartSideDoubleCubicCoarsen
#define included_CartSideDoubleCubicCoarsen

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "Box.h"
#include "CartesianSideDoubleWeightedAverage.h"
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
 * \brief Class CartSideDoubleCubicCoarsen is a concrete
 * SAMRAI::xfer::CoarsenOperator for restricting side-centered double precision
 * patch data via cubic interpolation.
 *
 * \note This coarsen operator requires that the refinement ratio between coarse
 * and fine index spaces is even.
 *
 * \note To avoid complicating the interpolation scheme at coarse-fine
 * interfaces, the use of cubic coarsening is restricted to cases where the
 * refinement ratio is at least 4.  For refinement ratios less than 4, a warning
 * is emitted and simple weighted averaging is used instead.
 */
class CartSideDoubleCubicCoarsen : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    CartSideDoubleCubicCoarsen();

    /*!
     * \brief Destructor.
     */
    ~CartSideDoubleCubicCoarsen();

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
    CartSideDoubleCubicCoarsen(const CartSideDoubleCubicCoarsen& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartSideDoubleCubicCoarsen& operator=(const CartSideDoubleCubicCoarsen& that);

    /*!
     * The operator name.
     */
    static const std::string s_op_name;

    /*!
     * "Backup" coarsen operator for even refinement ratios less than 4.
     */
    SAMRAI::geom::CartesianSideDoubleWeightedAverage<NDIM> d_weighted_average_coarsen_op;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartSideDoubleCubicCoarsen
