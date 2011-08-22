// Filename: VecCellRefineAdapter.h
// Created on 09 Apr 2010 by Boyce Griffith
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

#ifndef included_VecCellRefineAdapter
#define included_VecCellRefineAdapter

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <RefineOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class VecCellRefineAdapter is used to convert a refine operator which is
 * defined for regular SAMRAI::pdat::CellData patch data into a refine operator
 * which is defined for vector-valued VecCellData patch data.
 *
 * \see SAMRAI::xfer::RefineOperator
 */
class VecCellRefineAdapter
    : public SAMRAI::xfer::RefineOperator<NDIM>
{
public:
    /*!
     * Uninteresting constructor.
     */
    VecCellRefineAdapter(
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > cell_refine_op);

    /*!
     * Destructor.
     */
    ~VecCellRefineAdapter();

    /*!
     * Return true if the variable and name string match their expected values;
     * otherwise, return false.
     */
    bool
    findRefineOperator(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
        const std::string& op_name) const;

    /*!
     * Return name string identifier of this refinement operator.
     */
    const std::string&
    getOperatorName() const;

    /*!
     * The priority of the operator is the same as the cell-centered refine
     * operator encapsulated by this object.
     */
    int
    getOperatorPriority() const;

    /*!
     * The stencil width of the operator is the same as the cell-centered refine
     * operator encapsulated by this object.
     */
    SAMRAI::hier::IntVector<NDIM>
    getStencilWidth() const;

    /*!
     * Refine the source component on the coarse patch to the destination
     * component on the fine patch using the encapsulated cell-centered
     * interpolation operator.  Interpolation is performed on the intersection
     * of the destination patch and the fine box.  It is assumed that the coarse
     * patch contains sufficient data for the stencil width of the refinement
     * operator.
     */
    void
    refine(
        SAMRAI::hier::Patch<NDIM>& fine,
        const SAMRAI::hier::Patch<NDIM>& coarse,
        int dst_component,
        int src_component,
        const SAMRAI::hier::Box<NDIM>& fine_box,
        const SAMRAI::hier::IntVector<NDIM>& ratio) const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VecCellRefineAdapter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VecCellRefineAdapter(
        const VecCellRefineAdapter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VecCellRefineAdapter&
    operator=(
        const VecCellRefineAdapter& that);

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_cell_refine_op;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/VecCellRefineAdapter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VecCellRefineAdapter
