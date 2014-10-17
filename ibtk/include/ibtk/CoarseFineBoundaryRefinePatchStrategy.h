// Filename: CoarseFineBoundaryRefinePatchStrategy.h
// Created on 06 Nov 2007 by Boyce Griffith
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

#ifndef included_CoarseFineBoundaryRefinePatchStrategy
#define included_CoarseFineBoundaryRefinePatchStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <set>

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
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CoarseFineBoundaryRefinePatchStrategy is a subclass of the
 * abstract base class SAMRAI::xfer::RefinePatchStrategy that extends the
 * functionality of SAMRAI::xfer::RefinePatchStrategy to facilitate the
 * implementation of coarse-fine interface discretizations.
 */
class CoarseFineBoundaryRefinePatchStrategy : public SAMRAI::xfer::RefinePatchStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    CoarseFineBoundaryRefinePatchStrategy();

    /*!
     * \brief Destructor.
     */
    virtual ~CoarseFineBoundaryRefinePatchStrategy();

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
     * \param patch                Patch on which to fill boundary data.
     * \param fill_time            Double simulation time for boundary filling.
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over
     *all
     *registered scratch components.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) = 0;

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const = 0;

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
                          const SAMRAI::hier::IntVector<NDIM>& ratio) = 0;

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
                           const SAMRAI::hier::IntVector<NDIM>& ratio) = 0;

    //\}

    /*!
     * \name Extension of SAMRAI::xfer::RefinePatchStrategy interface to support more
     * complex coarse-fine interface discretizations.
     */
    //\{

    /*!
     * Whether or not to employ a consistent interpolation scheme at "Type 2"
     * coarse-fine interface ghost cells.
     *
     * \note This subclasses may choose not to support a consistent "Type 2"
     * coarse-fine interface ghost cell interpolation scheme.
     */
    virtual void setConsistentInterpolationScheme(bool consistent_type_2_bdry) = 0;

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    virtual void setPatchDataIndex(int patch_data_index) = 0;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    virtual void setPatchDataIndices(const std::set<int>& patch_data_indices) = 0;

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    virtual void setPatchDataIndices(const SAMRAI::hier::ComponentSelector& patch_data_indices) = 0;

    /*!
     * Set the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    virtual void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy) = 0;

    /*!
     * Clear the patch hierarchy used in constructing coarse-fine interface
     * boundary boxes.
     */
    virtual void clearPatchHierarchy() = 0;

    /*!
     * Compute the normal extension of fine data at coarse-fine interfaces.
     */
    virtual void computeNormalExtension(SAMRAI::hier::Patch<NDIM>& patch,
                                        const SAMRAI::hier::IntVector<NDIM>& ratio,
                                        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) = 0;

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
    CoarseFineBoundaryRefinePatchStrategy(const CoarseFineBoundaryRefinePatchStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CoarseFineBoundaryRefinePatchStrategy& operator=(const CoarseFineBoundaryRefinePatchStrategy& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CoarseFineBoundaryRefinePatchStrategy
