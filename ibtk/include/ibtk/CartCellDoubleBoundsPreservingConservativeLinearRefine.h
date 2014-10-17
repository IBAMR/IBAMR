// Filename: CartCellDoubleBoundsPreservingConservativeLinearRefine.h
// Created on 06 Jul 2010 by Boyce Griffith
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

#ifndef included_CartCellDoubleBoundsPreservingConservativeLinearRefine
#define included_CartCellDoubleBoundsPreservingConservativeLinearRefine

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "Box.h"
#include "CartesianCellDoubleConservativeLinearRefine.h"
#include "CellDoubleConstantRefine.h"
#include "IntVector.h"
#include "RefineOperator.h"
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
 * \brief Class CartCellDoubleBoundsPreservingConservativeLinearRefine is a
 * concrete SAMRAI::xfer::RefineOperator object which prolongs cell-centered
 * double precision patch data via conservative linear interpolation with an
 * additional bounds preservation repair step.
 */
class CartCellDoubleBoundsPreservingConservativeLinearRefine : public SAMRAI::xfer::RefineOperator<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     */
    CartCellDoubleBoundsPreservingConservativeLinearRefine();

    /*!
     * \brief Destructor.
     */
    ~CartCellDoubleBoundsPreservingConservativeLinearRefine();

    /*!
     * \name Implementation of SAMRAI::xfer::RefineOperator interface.
     */
    //\{

    /*!
     * Return true if the refining operation matches the variable and name
     * string identifier request; false, otherwise.
     */
    bool findRefineOperator(const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
                            const std::string& op_name) const;

    /*!
     * Return name string identifier of the refining operation.
     */
    const std::string& getOperatorName() const;

    /*!
     * Return the priority of this operator relative to other refining
     * operators.  The SAMRAI transfer routines guarantee that refining using
     * operators with lower priority will be performed before those with higher
     * priority.
     */
    int getOperatorPriority() const;

    /*!
     * Return the stencil width associated with the refining operator.  The
     * SAMRAI transfer routines guarantee that the source patch will contain
     * sufficient ghost cell data surrounding the interior to satisfy the
     * stencil width requirements for each refining operator.
     */
    SAMRAI::hier::IntVector<NDIM> getStencilWidth() const;

    /*!
     * Refine the source component on the fine patch to the destination
     * component on the coarse patch. The refining operation is performed on the
     * intersection of the destination patch and the coarse box.  The fine patch
     * is guaranteed to contain sufficient data for the stencil width of the
     * refining operator.
     */
    void refine(SAMRAI::hier::Patch<NDIM>& fine,
                const SAMRAI::hier::Patch<NDIM>& coarse,
                int dst_component,
                int src_component,
                const SAMRAI::hier::Box<NDIM>& fine_box,
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
    CartCellDoubleBoundsPreservingConservativeLinearRefine(
        const CartCellDoubleBoundsPreservingConservativeLinearRefine& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartCellDoubleBoundsPreservingConservativeLinearRefine&
    operator=(const CartCellDoubleBoundsPreservingConservativeLinearRefine& that);

    /*!
     * The operator name.
     */
    static const std::string s_op_name;

    /*!
     * The basic, non-bounds preserving conservative linear refine operator.
     */
    SAMRAI::geom::CartesianCellDoubleConservativeLinearRefine<NDIM> d_conservative_linear_refine_op;

    /*!
     * The constant refine operator.
     */
    SAMRAI::pdat::CellDoubleConstantRefine<NDIM> d_constant_refine_op;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartCellDoubleBoundsPreservingConservativeLinearRefine
