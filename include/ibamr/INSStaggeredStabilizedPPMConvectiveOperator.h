// Filename: INSStaggeredStabilizedPPMConvectiveOperator.h
// Created on 07 Sep 2012 by Boyce Griffith
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

#ifndef included_INSStaggeredStabilizedPPMConvectiveOperator
#define included_INSStaggeredStabilizedPPMConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "boost/array.hpp"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredStabilizedPPMConvectiveOperator is a concrete
 * ConvectiveOperator that implements a upwind convective differencing operator
 * based on the piecewise parabolic method (PPM).  This discretization is
 * "stabilized" at open boundaries by blending the higher-order discretization
 * with first-order upwind.
 *
 * Class INSStaggeredStabilizedPPMConvectiveOperator computes the convective
 * derivative of a side-centered velocity field using the xsPPM7 method of
 * Rider, Greenough, and Kamm.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredStabilizedPPMConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSStaggeredStabilizedPPMConvectiveOperator(const std::string& object_name,
                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                ConvectiveDifferencingType difference_form,
                                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredStabilizedPPMConvectiveOperator();

    /*!
     * \brief Static function to construct an
     * INSStaggeredStabilizedPPMConvectiveOperator.
     */
    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        return new INSStaggeredStabilizedPPMConvectiveOperator(object_name, input_db, difference_form, bc_coefs);
    } // allocate_operator

    /*!
     * \brief Compute the action of the convective operator.
     */
    void applyConvectiveOperator(int U_idx, int N_idx);

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for computing y=F[x] and
     * z=F[x]+y.
     *
     * The vector arguments for apply(), applyAdjoint(), etc, need not match
     * those for initializeOperatorState().  However, there must be a certain
     * degree of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the input and output vectors
     *
     * \note It is generally necessary to reinitialize the operator state when
     * the hierarchy configuration changes.
     *
     * It is safe to call initializeOperatorState() when the state is already
     * initialized.  In this case, the operator state is first deallocated and
     * then reinitialized.
     *
     * Conditions on arguments:
     * - input and output vectors must have same hierarchy
     * - input and output vectors must have same structure, depth, etc.
     *
     * Call deallocateOperatorState() to remove any data allocated by this
     * method.
     *
     * \see deallocateOperatorState
     *
     * \param in input vector
     * \param out output vector
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * \note It is safe to call deallocateOperatorState() when the operator
     * state is already deallocated.
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredStabilizedPPMConvectiveOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredStabilizedPPMConvectiveOperator(const INSStaggeredStabilizedPPMConvectiveOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredStabilizedPPMConvectiveOperator& operator=(const INSStaggeredStabilizedPPMConvectiveOperator& that);

    // Operator configuration.
    boost::array<bool, 2 * NDIM> d_open_bdry;
    boost::array<double, 2 * NDIM> d_width;

    // Boundary condition helper object.
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Cached communications operators.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    std::string d_bdry_extrap_type;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    int d_U_scratch_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredStabilizedPPMConvectiveOperator
