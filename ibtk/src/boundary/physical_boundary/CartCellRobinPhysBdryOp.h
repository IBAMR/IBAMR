// Filename: CartCellRobinPhysBdryOp.h
// Created on 10 Feb 2007 by Boyce Griffith
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

#ifndef included_CartCellRobinPhysBdryOp
#define included_CartCellRobinPhysBdryOp

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/CartExtrapPhysBdryOp.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <ComponentSelector.h>
#include <IntVector.h>
#include <Patch.h>
#include <RefinePatchStrategy.h>
#include <RobinBcCoefStrategy.h>
#include <tbox/DescribedClass.h>

// C++ STDLIB INCLUDES
#include <set>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartCellRobinPhysBdryOp is a concrete
 * SAMRAI::xfer::RefinePatchStrategy for setting Robin boundary conditions at
 * physical boundaries for cell-centered scalar- and vector-valued quantities.
 *
 * This implementation works similarly to the
 * SAMRAI::solv::CartesianRobinBcHelper class.  An important difference between
 * the two classes is that class CartCellRobinPhysBdryOp allows for the
 * specification of boundary conditions for vector-valued quantities.
 */
class CartCellRobinPhysBdryOp
    : public SAMRAI::xfer::RefinePatchStrategy<NDIM>,
      public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor leaves the class in an undefined state.  In
     * particular, the patch data indices requiring cell filling and the
     * boundary condition specifications must be set prior to use of method
     * setPhysicalBoundaryConditions().
     *
     * \see setPatchDataIndex
     * \see setPatchDataIndices
     * \see setPhysicalBcCoef
     * \see setPhysicalBcCoefs
     * \see setHomogeneousBc
     */
    CartCellRobinPhysBdryOp();

    /*!
     * \brief Constructor to fill boundary conditions for scalar-valued
     * quantities.
     *
     * \param patch_data_index  Patch data index requiring ghost cell filling.
     * \param bc_coef           Robin boundary conditions to use with this class.
     * \param homogeneous_bc    Whether to employ the homogeneous form of the boundary conditions.
     */
    CartCellRobinPhysBdryOp(
        int patch_data_index,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
        bool homogeneous_bc=false);

    /*!
     * \brief Constructor to fill boundary conditions for scalar-valued
     * quantities.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell filling.
     * \param bc_coef             Robin boundary conditions to use with this class.
     * \param homogeneous_bc      Whether to employ the homogeneous form of the boundary conditions.
     */
    CartCellRobinPhysBdryOp(
        const std::set<int>& patch_data_indices,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
        bool homogeneous_bc=false);

    /*!
     * \brief Constructor to fill boundary conditions for scalar-valued
     * quantities.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell filling.
     * \param bc_coef             Robin boundary conditions to use with this class.
     * \param homogeneous_bc      Whether to employ the homogeneous form of the boundary conditions.
     */
    CartCellRobinPhysBdryOp(
        const SAMRAI::hier::ComponentSelector& patch_data_indices,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
        bool homogeneous_bc=false);

    /*!
     * \brief Constructor to fill boundary conditions for vector-valued
     * quantities.  In this case, we require that distinct boundary condition
     * objects are provided for each data depth.
     *
     * \param patch_data_index  Patch data index requiring ghost cell filling.
     * \param bc_coefs          Robin boundary conditions to use with this class.
     * \param homogeneous_bc    Whether to employ the homogeneous form of the boundary conditions.
     */
    CartCellRobinPhysBdryOp(
        int patch_data_index,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        bool homogeneous_bc=false);

    /*!
     * \brief Constructor to fill boundary conditions for vector-valued
     * quantities.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell filling.
     * \param bc_coefs            Robin boundary conditions to use with this class.
     * \param homogeneous_bc      Whether to employ the homogeneous form of the boundary conditions.
     */
    CartCellRobinPhysBdryOp(
        const std::set<int>& patch_data_indices,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        bool homogeneous_bc=false);

    /*!
     * \brief Constructor to fill boundary conditions for vector-valued
     * quantities.
     *
     * \param patch_data_indices  Collection of patch data indices which require ghost cell filling.
     * \param bc_coefs            Robin boundary conditions to use with this class.
     * \param homogeneous_bc      Whether to employ the homogeneous form of the boundary conditions.
     */
    CartCellRobinPhysBdryOp(
        const SAMRAI::hier::ComponentSelector& patch_data_indices,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        bool homogeneous_bc=false);

    /*!
     * \brief Destructor.
     */
    ~CartCellRobinPhysBdryOp();

    /*!
     * \brief Reset the patch data index operated upon by this class.
     */
    void
    setPatchDataIndex(
        int patch_data_index);

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void
    setPatchDataIndices(
        const std::set<int>& patch_data_indices);

    /*!
     * \brief Reset the patch data indices operated upon by this class.
     */
    void
    setPatchDataIndices(
        const SAMRAI::hier::ComponentSelector& patch_data_indices);

    /*!
     * \brief Reset the Robin boundary condition specification object employed
     * by this class to set physical boundary conditions.
     *
     * \note \a bc_coef cannot be NULL.
     */
    void
    setPhysicalBcCoef(
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef);

    /*!
     * \brief Reset the Robin boundary condition specification object employed
     * by this class to set physical boundary conditions.
     *
     * \note None of the elements of \a bc_coefs can be NULL.
     */
    void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Set whether boundary filling should employ homogeneous boundary
     * conditions.
     *
     * \note By default, inhomogeneous boundary conditions are assumed.
     */
    void
    setHomogeneousBc(
        bool homogeneous_bc);

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
     * \param ghost_width_to_fill  Integer vector describing maximum ghost width to fill over all registered scratch components.
     */
    void
    setPhysicalBoundaryConditions(
        SAMRAI::hier::Patch<NDIM>& patch,
        double fill_time,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill);

    /*!
     * Function to return maximum stencil width needed over user-defined data
     * interpolation operations.  This is needed to determine the correct
     * interpolation data dependencies.
     *
     * Presently, the refine operator stencil width is zero.
     */
    SAMRAI::hier::IntVector<NDIM>
    getRefineOpStencilWidth() const;

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
     * \param ratio     Integer vector containing ratio relating index space between coarse and fine patches.
     */
    void
    preprocessRefine(
        SAMRAI::hier::Patch<NDIM>& fine,
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
     * \param ratio     Integer vector containing ratio relating index space between coarse and fine patches.
     */
    void
    postprocessRefine(
        SAMRAI::hier::Patch<NDIM>& fine,
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
    CartCellRobinPhysBdryOp(
        const CartCellRobinPhysBdryOp& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartCellRobinPhysBdryOp&
    operator=(
        const CartCellRobinPhysBdryOp& that);

    /*
     * The patch data indices corresponding to the "scratch" patch data that
     * requires extrapolation of ghost cell values at physical boundaries.
     */
    std::set<int> d_patch_data_indices;

    /*
     * The RobinBcCoefStrategy objects used to specify Robin boundary conditions
     * for each data depth.
     *
     * The boolean value indicates whether homogeneous boundary conditions
     * should be used.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    bool d_homogeneous_bc;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/CartCellRobinPhysBdryOp.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartCellRobinPhysBdryOp
