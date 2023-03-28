// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_INSStaggeredWavePropConvectiveOperator
#define included_IBAMR_INSStaggeredWavePropConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ConvectiveOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/ibtk_utilities.h"

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

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
 * \brief Class INSStaggeredWavePropConvectiveOperator is a concrete
 * ConvectiveOperator that implements a convective operator based on
 * the wave propagation method.
 *
 * Class INSStaggeredWavePropConvectiveOperator computes the
 * convective derivative of a cell-centered velocity field using the
 * wave propagation method of Ketcheson and LeVeque.
 *
 * \note INSStaggeredWavePropConvectiveOperator currently only works
 * in ADVECTIVE form.
 *
 * \todo Implement a CONSERVATIVE form for the operator. This will
 * involve writing a new Riemann solver for the flux.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredWavePropConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSStaggeredWavePropConvectiveOperator(std::string object_name,
                                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                           ConvectiveDifferencingType difference_form,
                                           std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> bc_coefs);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredWavePropConvectiveOperator();

    /*!
     * \brief Static function to construct an INSStaggeredWavePropConvectiveOperator.
     */
    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        return new INSStaggeredWavePropConvectiveOperator(object_name, input_db, difference_form, bc_coefs);
    } // allocate_operator

    /*!
     * \brief Compute the action of the convective operator.
     */
    void applyConvectiveOperator(int U_idx, int N_idx) override;

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
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * \note It is safe to call deallocateOperatorState() when the operator
     * state is already deallocated.
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState() override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredWavePropConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredWavePropConvectiveOperator(const INSStaggeredWavePropConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredWavePropConvectiveOperator& operator=(const INSStaggeredWavePropConvectiveOperator& that) = delete;

    // Boundary condition helper object.
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Cached communications operators.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    std::string d_bdry_extrap_type = "CONSTANT";
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    int d_U_scratch_idx = IBTK::invalid_index;

    // Reconstruction order 2*k-1.
    // Currently only available for k=3
    int d_k = 3;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSStaggeredWavePropConvectiveOperator
