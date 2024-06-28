// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_CellConvectiveOperator
#define included_IBAMR_CellConvectiveOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/ConvectiveOperator.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CellConvectiveOperator is an abstract class for an implementation of
 * a convective differencing operator.
 */
class CellConvectiveOperator : public ConvectiveOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    CellConvectiveOperator(std::string object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > Q_cell_var,
                           int Q_min_ghost_cell_width,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           ConvectiveDifferencingType difference_form,
                           std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> bc_coefs);

    /*!
     * \brief Default destructor.
     */
    ~CellConvectiveOperator() = default;

    /*!
     * \brief Interpolate a cell-centered field Q to a face-centered field q, possibly using the provided advection
     * velocity field.
     */
    virtual void interpolateToFaceOnHierarchy(int q_interp_idx, int Q_cell_idx, int u_idx, bool synch_cf_bdry = true);

    /*!
     * \brief Evaluate the face-centered flux Q to a face-centered field q using the provided advection velocity field.
     */
    virtual void
    evaluateAdvectiveFluxOnHierarchy(int q_flux_idx, int Q_cell_idx, int u_idx, const bool synch_cf_bdry = true);

    /*!
     * \brief Compute the advective derivative N = u * grad Q on the patch hierarchy.
     */
    virtual void
    computeAdvectiveDerivativeOnHierarchy(int N_cell_idx, int q_interp_idx, int u_idx, bool synch_cf_bdry = true);

    /*!
     * \brief Compute the conservative derivative N = div(Q u) on the patch hierarchy.
     */
    virtual void computeConservativeDerivativeOnHierarchy(int N_cell_idx, int q_flux_idx, bool synch_cf_bdry = true);

    /*!
     * \brief Compute the skew-symmetric derivative N = 0.5[u * grad Q + div(Q u)] on the patch hierarchy.
     */
    virtual void computeSkewSymmetricDerivativeOnHierarchy(int N_cell_idx,
                                                           int q_interp_idx,
                                                           int q_flux_idx,
                                                           int u_idx,
                                                           bool synch_cf_bdry = true);

    /*!
     * \brief Interpolate a cell-centered field Q to a face-centered field q on a single grid patch.
     */
    virtual void interpolateToFaceOnPatch(SAMRAI::pdat::FaceData<NDIM, double>& q_interp_data,
                                          const SAMRAI::pdat::CellData<NDIM, double>& Q_cell_data,
                                          const SAMRAI::pdat::FaceData<NDIM, double>& u_data,
                                          const SAMRAI::hier::Patch<NDIM>& patch) = 0;

    /*!
     * \brief Evaluate the face-centered flux Q to a face-centered field q using the provided advection velocity field.
     *
     * A default implementation is provided that uses interpolateToFaceOnPatch to determine the advective fluxes.
     */
    virtual void evaluateAdvectiveFluxOnPatch(SAMRAI::pdat::FaceData<NDIM, double>& q_flux_data,
                                              const SAMRAI::pdat::CellData<NDIM, double>& Q_cell_data,
                                              const SAMRAI::pdat::FaceData<NDIM, double>& u_data,
                                              const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Compute N = u * grad Q or N = div(Q u).
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
     * The vector arguments for apply(), applyAdd(), etc, need not match those
     * for initializeOperatorState().  However, there must be a certain degree
     * of similarity, including
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
     */
    void deallocateOperatorState() override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CellConvectiveOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CellConvectiveOperator(const CellConvectiveOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CellConvectiveOperator& operator=(const CellConvectiveOperator& that) = delete;

    // Data communication algorithms, operators, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_Q_cell_refine_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_Q_cell_refine_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_Q_cell_refine_bdry_op;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_Q_cell_refine_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_q_flux_coarsen_op, d_q_interp_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_q_flux_coarsen_alg, d_q_interp_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_q_flux_coarsen_scheds,
        d_q_interp_coarsen_scheds;
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    std::string d_outflow_bdry_extrap_type = "CONSTANT";

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_cell_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_q_flux_var, d_q_interp_var, d_u_var;
    int d_Q_scratch_idx, d_Q_ghost_idx, d_q_flux_idx, d_q_interp_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CellConvectiveOperator
