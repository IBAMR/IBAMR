// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_INSCollocatedCenteredConvectiveOperator
#define included_IBAMR_INSCollocatedCenteredConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "FaceVariable.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefinePatchStrategy.h"
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
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSCollocatedCenteredConvectiveOperator is a concrete
 * ConvectiveOperator which implements a centered convective differencing
 * operator.
 *
 * \see INSCollocatedHierarchyIntegrator
 */
class INSCollocatedCenteredConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSCollocatedCenteredConvectiveOperator(std::string object_name,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                            ConvectiveDifferencingType difference_form,
                                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Destructor.
     */
    ~INSCollocatedCenteredConvectiveOperator();

    /*!
     * \brief Static function to construct an
     * INSCollocatedCenteredConvectiveOperator.
     */
    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        return new INSCollocatedCenteredConvectiveOperator(object_name, input_db, difference_form, bc_coefs);
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
    INSCollocatedCenteredConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSCollocatedCenteredConvectiveOperator(const INSCollocatedCenteredConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSCollocatedCenteredConvectiveOperator& operator=(const INSCollocatedCenteredConvectiveOperator& that) = delete;

    // Data communication algorithms, operators, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_coarsen_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_ghostfill_strategy;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_scheds;
    std::string d_bdry_extrap_type = "CONSTANT";

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_U_var;
    int d_U_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_extrap_var, d_u_flux_var;
    int d_u_extrap_idx = IBTK::invalid_index, d_u_flux_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSCollocatedCenteredConvectiveOperator
