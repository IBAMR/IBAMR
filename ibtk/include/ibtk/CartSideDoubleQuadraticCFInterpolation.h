// Filename: CartSideDoubleQuadraticCFInterpolation.h
// Created on 30 Apr 2008 by Boyce Griffith
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

#ifndef included_CartSideDoubleQuadraticCFInterpolation
#define included_CartSideDoubleQuadraticCFInterpolation

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <set>
#include <vector>

#include "Box.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineOperator.h"
#include "SideVariable.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class CoarseFineBoundary;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartSideDoubleQuadraticCFInterpolation is a concrete
 * SAMRAI::xfer::RefinePatchStrategy which sets coarse-fine interface ghost cell
 * values for side-centered double precision patch data via quadratic
 * interpolation in the normal and tangential directions at coarse-fine
 * interfaces.
 */
class CartSideDoubleQuadraticCFInterpolation : public CoarseFineBoundaryRefinePatchStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    CartSideDoubleQuadraticCFInterpolation();

    /*!
     * \brief Destructor.
     */
    ~CartSideDoubleQuadraticCFInterpolation();

    /*!
     * \name SAMRAI::xfer::RefinePatchStrategy interface.
     */
    //\{

    /*!
     * Function to set data associated with the given list of patch data indices
     * at patch boundaries that intersect the physical domain boundary.  The
     * patch data components set in this routine correspond to the "scratch"
     * components specified in calls to the registerRefine() function in the
     * SAMRAI::xfer::RefineAlgorithm class.
     *
     * Presently, the implementation does nothing.
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

    /*!
     * \name Extension of SAMRAI::xfer::RefinePatchStrategy interface to support more
     * complex coarse-fine interface discretizations.
     */
    //\{

    /*!
     * Whether or not to employ a consistent interpolation scheme at "Type 2"
     * coarse-fine interface ghost cells.
     */
    void setConsistentInterpolationScheme(bool consistent_type_2_bdry);

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
     * Set the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Clear the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    void clearPatchHierarchy();

    /*!
     * Compute the normal extension of fine data at coarse-fine interfaces.
     */
    void computeNormalExtension(SAMRAI::hier::Patch<NDIM>& patch,
                                const SAMRAI::hier::IntVector<NDIM>& ratio,
                                const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

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
    CartSideDoubleQuadraticCFInterpolation(const CartSideDoubleQuadraticCFInterpolation& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartSideDoubleQuadraticCFInterpolation& operator=(const CartSideDoubleQuadraticCFInterpolation& that);

    /*!
     * The patch data indices corresponding to the "scratch" patch data that is
     * operated on by this class.
     */
    std::set<int> d_patch_data_indices;

    /*!
     * Boolean value indicating whether we are enforcing a consistent
     * interpolation scheme at "Type 2" coarse-fine interface ghost cells.
     */
    bool d_consistent_type_2_bdry;

    /*!
     * Refine operator employed to fill coarse grid ghost cell values.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_refine_op;

    /*!
     * Cached hierarchy-related information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    std::vector<SAMRAI::hier::CoarseFineBoundary<NDIM>*> d_cf_boundary;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int> > d_sc_indicator_var;
    int d_sc_indicator_idx;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartSideDoubleQuadraticCFInterpolation
