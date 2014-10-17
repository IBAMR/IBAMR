// Filename: CartSideDoubleDivPreservingRefine.h
// Created on 09 Nov 2008 by Boyce Griffith
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

#ifndef included_CartSideDoubleDivPreservingRefine
#define included_CartSideDoubleDivPreservingRefine

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "Box.h"
#include "CoarsenOperator.h"
#include "IntVector.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartSideDoubleDivPreservingRefine is a concrete
 * SAMRAI::xfer::RefinePatchStrategy which prolongs side-centered double
 * precision patch data via conservative linear interpolation with divergence-
 * and curl-preserving corrections.
 */
class CartSideDoubleDivPreservingRefine : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    CartSideDoubleDivPreservingRefine(int u_dst_idx,
                                      int u_src_idx,
                                      int indicator_idx,
                                      SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_op,
                                      SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_op,
                                      double fill_time,
                                      SAMRAI::xfer::RefinePatchStrategy<NDIM>* phys_bdry_op);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~CartSideDoubleDivPreservingRefine();

    /*!
     * \brief The number of required ghost cells.
     *
     * \note This value is chosen to allow refinement ratios up to 4.  A larger
     * value would be necessary for refinement ratios greater than 4.
     */
    static const int REFINE_OP_STENCIL_WIDTH = 4;

    /*!
     * \name Implementation of SAMRAI::xfer::RefinePatchStrategy interface.
     */
    //\{

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * specific boundary conditions are determined by the user.  The patch data
     * components set in this routine correspond to the scratch components
     * specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     *all
     *registered scratch components.
     */
    virtual void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                               double fill_time,
                                               const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    virtual SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const;

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called before standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The preprocess function refines data from the scratch components
     * of the coarse patch into the scratch components of the fine patch on the
     * specified fine box region.  Recall that the scratch components are
     * specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    virtual void preprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                                  const SAMRAI::hier::Patch<NDIM>& coarse,
                                  const SAMRAI::hier::Box<NDIM>& fine_box,
                                  const SAMRAI::hier::IntVector<NDIM>& ratio);

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called after standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The postprocess function refines data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    virtual void postprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                                   const SAMRAI::hier::Patch<NDIM>& coarse,
                                   const SAMRAI::hier::Box<NDIM>& fine_box,
                                   const SAMRAI::hier::IntVector<NDIM>& ratio);

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CartSideDoubleDivPreservingRefine();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartSideDoubleDivPreservingRefine(const CartSideDoubleDivPreservingRefine& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartSideDoubleDivPreservingRefine& operator=(const CartSideDoubleDivPreservingRefine& that);

    /*!
     * Patch data indices.
     */
    const int d_u_dst_idx;
    const int d_u_src_idx;
    const int d_indicator_idx;

    /*!
     * Routines for setting physical boundary conditions.
     */
    const double d_fill_time;
    SAMRAI::xfer::RefinePatchStrategy<NDIM>* const d_phys_bdry_op;

    /*!
     * The basic linear refine operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_refine_op;

    /*!
     * The basic coarsening operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_coarsen_op;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartSideDoubleDivPreservingRefine
