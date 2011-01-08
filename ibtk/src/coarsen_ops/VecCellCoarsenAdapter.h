// Filename: VecCellCoarsenAdapter.h
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

#ifndef included_VecCellCoarsenAdapter
#define included_VecCellCoarsenAdapter

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CoarsenOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class VecCellCoarsenAdapter is used to convert a coarsen operator which is
 * defined for regular SAMRAI::pdat::CellData patch data into a coarsen operator
 * which is defined for vector-valued VecCellData patch data.
 *
 * \see SAMRAI::xfer::CoarsenOperator
 */
class VecCellCoarsenAdapter
    : public SAMRAI::xfer::CoarsenOperator<NDIM>
{
public:
    /*!
     * Uninteresting constructor.
     */
    VecCellCoarsenAdapter(
        SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > cell_coarsen_op);

    /*!
     * Uninteresting virtual destructor.
     */
    virtual
    ~VecCellCoarsenAdapter();

    /*!
     * Return true if the variable and name string match their expected values;
     * otherwise, return false.
     */
    virtual bool
    findCoarsenOperator(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& var,
        const std::string& op_name) const;

    /*!
     * Return name string identifier of this coarsening operation.
     */
    virtual const std::string&
    getOperatorName() const;

    /*!
     * The priority of the operator is the same as the cell-centered coarsen
     * operator encapsulated by this object.
     */
    virtual int
    getOperatorPriority() const;

    /*!
     * The stencil width of the operator is the same as the cell-centered
     * coarsen operator encapsulated by this object.
     */
    virtual SAMRAI::hier::IntVector<NDIM>
    getStencilWidth() const;

    /*!
     * Coarsen the source component on the fine patch to the destination
     * component on the coarse patch using the encapsulated cell-centered
     * coarsening operator.  Coarsening is performed on the intersection of the
     * destination patch and the coarse box.  It is assumed that the fine patch
     * contains sufficient data for the stencil width of the coarsening
     * operator.
     */
    virtual void
    coarsen(
        SAMRAI::hier::Patch<NDIM>& coarse,
        const SAMRAI::hier::Patch<NDIM>& fine,
        const int dst_component,
        const int src_component,
        const SAMRAI::hier::Box<NDIM>& coarse_box,
        const SAMRAI::hier::IntVector<NDIM>& ratio) const;

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VecCellCoarsenAdapter();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VecCellCoarsenAdapter(
        const VecCellCoarsenAdapter& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VecCellCoarsenAdapter&
    operator=(
        const VecCellCoarsenAdapter& that);

    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_cell_coarsen_op;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/VecCellCoarsenAdapter.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_VecCellCoarsenAdapter
