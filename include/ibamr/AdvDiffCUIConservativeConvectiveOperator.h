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

#ifndef included_IBAMR_AdvDiffCUIConservativeConvectiveOperator
#define included_IBAMR_AdvDiffCUIConservativeConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/AdvDiffCUIConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"

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
 * \brief Class AdvDiffCUIConservativeConvectiveOperator is a concrete ConvectiveOperator
 * that implements a upwind convective differencing operator based on the
 * cubic upwind interpolation (CUI).
 *
 * Class AdvDiffCUIConservativeConvectiveOperator computes the convective derivative which is in
 * the form of div(Q P u) using the CUI method described by Waterson and Deconinck,
 * and Patel and Natarajan.
 *
 *
 * References
 * Waterson, NP. and Deconinck, H., <A HREF="https://www.sciencedirect.com/science/article/pii/S002199910700040X">
 * Design principles for bounded higher-order convection schemes – a unified approach</A>
 *
 * Patel, JK. and Natarajan, G., <A HREF="https://www.sciencedirect.com/science/article/pii/S0045793014004009">
 * A generic framework for design of interface capturing schemes for multi-fluid flows</A>
 *
 * \see AdvDiffSemiImplicitHierarchyIntegrator
 */
class AdvDiffCUIConservativeConvectiveOperator : public AdvDiffCUIConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    AdvDiffCUIConservativeConvectiveOperator(std::string object_name,
                                             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > P_var,
                                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                             std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> Q_bc_coefs,
                                             std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> P_bc_coefs);

    /*!
     * \brief Destructor.
     */
    ~AdvDiffCUIConservativeConvectiveOperator();

    /*!
     * \brief Compute the action of the convective operator.
     */
    void applyConvectiveOperator(int Q_idx, int N_idx) override;

    /*!
     * \brief Compute the action of the convective operator.
     */
    void applyConvectiveOperator(int Q_idx, int P_idx, int N_idx) override;

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
    AdvDiffCUIConservativeConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffCUIConservativeConvectiveOperator(const AdvDiffCUIConservativeConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffCUIConservativeConvectiveOperator& operator=(const AdvDiffCUIConservativeConvectiveOperator& that) = delete;

    // Data communication algorithms, operators, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_alg1;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_ghostfill_strategy1;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_scheds1;
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_P_bc_coefs;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P_var;
    unsigned int d_P_data_depth = 0;
    int d_P_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_p_extrap_var;
    int d_p_extrap_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AdvDiffCUIConservativeConvectiveOperator
