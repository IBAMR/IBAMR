// Filename: AdvDiffCUIConservativeConvectiveOperator.h
// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_AdvDiffCUIConservativeConvectiveOperator
#define included_IBAMR_AdvDiffCUIConservativeConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/AdvDiffCUIConvectiveOperator.h"
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
 * \brief Class AdvDiffCUIConservativeConvectiveOperator is a concrete ConvectiveOperator
 * that implements a upwind convective differencing operator based on the
 * cubic upwind interpolation (CUI).
 *
 * Class AdvDiffCUIConservativeConvectiveOperator computes the convective derivative of a
 * cell-centered field (in conservative form) using the CUI method described by Waterson and Deconinck,
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
 *
 */

class AdvDiffCUIConservativeConvectiveOperator : public AdvDiffCUIConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    AdvDiffCUIConservativeConvectiveOperator(std::string object_name,
                                             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                             ConvectiveDifferencingType difference_form,
                                             std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> bc_coefs);
    /*!
     * \brief Destructor.
     */
    ~AdvDiffCUIConservativeConvectiveOperator();

    /*!
     * \brief Static function to construct an AdvDiffCUIConservativeConvectiveOperator.
     */
    static SAMRAI::tbox::Pointer<ConvectiveOperator>
    allocate_operator(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_var,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      ConvectiveDifferencingType difference_form,
                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs)
    {
        return new AdvDiffCUIConservativeConvectiveOperator(object_name, Q_var, input_db, difference_form, bc_coefs);
    } // allocate_operator

    /*!
     * \brief Compute the action of the convective operator.
     */
    void applyConvectiveOperator(int Q_idx, int N_idx) override;

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
    /*!
     * \brief Set the mass density density variable
     * to be used when computing the conservative convective derivative.
     */
    void setMassDensityVariable(SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var);

    /*!
     * \brief Set the patch data index corresponding to the mass density
     * to be used when computing the conservative convective derivative.
     */
    void setMassDensity(int rho_idx);

    /*!
     * \brief Set the boundary condition object of the mass density
     * to be used when computing the conservative convective derivative.
     */
    void setMassDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* rho_bc_coef);

private:
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_rho_bc_coefs;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_rho_var;
    unsigned int d_rho_data_depth = 0;
    int d_rho_idx;
    int d_rho_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_rho_extrap_var, d_rho_flux_var;
    int d_rho_extrap_idx = IBTK::invalid_index, d_rho_flux_idx = IBTK::invalid_index;
};
} // namespace IBAMR

#endif //#ifndef included_IBAMR_AdvDiffCUIConservativeConvectiveOperator
