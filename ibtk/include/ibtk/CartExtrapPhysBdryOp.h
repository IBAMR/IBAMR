// Filename: CartExtrapPhysBdryOp.h
// Created on 30 Sep 2006 by Boyce Griffith
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

#ifndef included_CartExtrapPhysBdryOp
#define included_CartExtrapPhysBdryOp

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Box.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "RefinePatchStrategy.h"

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
 * \brief Class CartExtrapPhysBdryOp is a concrete
 * SAMRAI::xfer::RefinePatchStrategy for setting ghost cell values at physical
 * boundaries via constant, linear, or quadratic extrapolation from interior
 * values.
 */
class CartExtrapPhysBdryOp : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor leaves the class in an undefined state.  In
     * particular, the patch data indices requiring cell filling and
     * extrapolation type must each be set prior to use of method
     * setPhysicalBoundaryConditions().
     *
     * \see setPatchDataIndex
     * \see setPatchDataIndices
     * \see setExtrapolationType
     */
    CartExtrapPhysBdryOp();

    /*!
     * \brief Constructor.
     *
     * \param patch_data_index  Patch data index requiring ghost cell filling.
     * \param extrap_type       Type of extrapolation to perform: \p "CONSTANT" specifies
     *constant
     *extrapolation, \p "LINEAR" specifies linear extrapolation, and \p "QUADRATIC" specifies
     *quadratic extrapolation.
     */
    CartExtrapPhysBdryOp(int patch_data_index, const std::string& extrap_type = "CONSTANT");

    /*!
     * \brief Constructor.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell
     *filling.
     * \param extrap_type         Type of extrapolation to perform: \p "CONSTANT" specifies
     *constant
     *extrapolation, \p "LINEAR" specifies linear extrapolation, and \p "QUADRATIC" specifies
     *quadratic extrapolation.
     */
    CartExtrapPhysBdryOp(const std::set<int>& patch_data_indices, const std::string& extrap_type = "CONSTANT");

    /*!
     * \brief Constructor.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell
     *filling.
     * \param extrap_type         Type of extrapolation to perform: \p "CONSTANT" specifies
     *constant
     *extrapolation, \p "LINEAR" specifies linear extrapolation, and \p "QUADRATIC" specifies
     *quadratic extrapolation.
     */
    CartExtrapPhysBdryOp(const SAMRAI::hier::ComponentSelector& patch_data_indices,
                         const std::string& extrap_type = "CONSTANT");

    /*!
     * \brief Destructor.
     */
    ~CartExtrapPhysBdryOp();

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    void setPatchDataIndex(int patch_data_index);

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const std::set<int>& patch_data_indices);

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void setPatchDataIndices(const SAMRAI::hier::ComponentSelector& patch_data_indices);

    /*!
     * \brief Set the extrapolation type employed by this class.
     *
     * \param extrap_type  Type of extrapolation to perform: \p "CONSTANT" specifies constant
     *extrapolation, \p "LINEAR" specifies linear extrapolation, and \p "QUADRATIC" specifies
     *quadratic extrapolation.
     */
    void setExtrapolationType(const std::string& extrap_type);

    /*!
     * \name Implementation of SAMRAI::xfer::RefinePatchStrategy interface.
     */
    //\{

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * patch data components set in this routine correspond to the "scratch"
     * components specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     *all
     *registered scratch components.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     *
     * Presently, the refine operator stencil width is zero.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const;

    /*!
     * Function to perform user-defined preprocess data refine operations.  This
     * member function is called before standard refine operations (expressed
     * using concrete subclasses of the SAMRAI::xfer::RefineOperator base
     * class).  The preprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void preprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                          const SAMRAI::hier::Patch<NDIM>& coarse,
                          const SAMRAI::hier::Box<NDIM>& fine_box,
                          const SAMRAI::hier::IntVector<NDIM>& ratio);

    /*!
     * Function to perform user-defined postprocess data refine operations.
     * This member function is called after standard refine operations
     * (expressed using concrete subclasses of the SAMRAI::xfer::RefineOperator
     * base class).  The postprocess function must refine data from the scratch
     * components of the coarse patch into the scratch components of the fine
     * patch on the specified fine box region.  Recall that the scratch
     * components are specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
     *
     * \param fine      Fine patch containing destination data.
     * \param coarse    Coarse patch containing source data.
     * \param fine_box  Box region on fine patch into which data is refined.
     * \param ratio     Integer vector containing ratio relating index space between coarse and
     *fine
     *patches.
     */
    void postprocessRefine(SAMRAI::hier::Patch<NDIM>& fine,
                           const SAMRAI::hier::Patch<NDIM>& coarse,
                           const SAMRAI::hier::Box<NDIM>& fine_box,
                           const SAMRAI::hier::IntVector<NDIM>& ratio);

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
    CartExtrapPhysBdryOp(const CartExtrapPhysBdryOp& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartExtrapPhysBdryOp& operator=(const CartExtrapPhysBdryOp& that);

    /*!
     * \brief The implementation of setPhysicalBoundaryConditions() for
     * cell-centered quantities.
     */
    void setPhysicalBoundaryConditions_cell(
        SAMRAI::hier::Patch<NDIM>& patch,
        const std::vector<std::pair<SAMRAI::hier::Box<NDIM>, std::pair<int, int> > >& bdry_fill_boxes);

    /*!
     * \brief The implementation of setPhysicalBoundaryConditions() for
     * face-centered quantities.
     */
    void setPhysicalBoundaryConditions_face(
        SAMRAI::hier::Patch<NDIM>& patch,
        const std::vector<std::pair<SAMRAI::hier::Box<NDIM>, std::pair<int, int> > >& bdry_fill_boxes);

    /*!
     * \brief The implementation of setPhysicalBoundaryConditions() for
     * node-centered quantities.
     */
    void setPhysicalBoundaryConditions_node(
        SAMRAI::hier::Patch<NDIM>& patch,
        const std::vector<std::pair<SAMRAI::hier::Box<NDIM>, std::pair<int, int> > >& bdry_fill_boxes);

    /*!
     * \brief The implementation of setPhysicalBoundaryConditions() for
     * side-centered quantities.
     */
    void setPhysicalBoundaryConditions_side(
        SAMRAI::hier::Patch<NDIM>& patch,
        const std::vector<std::pair<SAMRAI::hier::Box<NDIM>, std::pair<int, int> > >& bdry_fill_boxes);

    /*
     * The patch data indices corresponding to the "scratch" patch data that
     * requires extrapolation of ghost cell values at physical boundaries.
     */
    std::set<int> d_patch_data_indices;

    /*
     * The type of extrapolation to perform.  Choices are presently "CONSTANT",
     * "LINEAR", or "QUADRATIC".
     */
    std::string d_extrap_type;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartExtrapPhysBdryOp
