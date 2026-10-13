// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#ifndef included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp
#define included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <BoxArray.h>
#include <IntVector.h>

#include <array>
#include <map>
#include <memory>

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredDivergenceFreePhysBdryOp fills the ghost values of a side-centered velocity outside the
 * physical boundaries of a rectangular domain so that the discrete divergence vanishes in every ghost cell there.
 *
 * The operator reads the boundary conditions that the fluid solver registered: \f$ a = 1 \f$ prescribes the velocity
 * \f$ g \f$ and \f$ b = 1 \f$ prescribes the traction (or pseudo-traction) \f$ g \f$. It emits an error for other
 * values, and it ignores the conditions on the normal component, setPhysicalBcCoef() and setPhysicalBcCoefs(). To fill
 * ghost width \f$ w \f$, the patch data need ghost width \f$ w + 1 \f$ and the domain must be at least that wide along
 * each non-periodic direction; the operator sets a ghost value that depends on a face outside the ghost box of the
 * patch data to quiet NaN. The caller fills the values inside the domain that a patch stores in its ghost cells before
 * it calls setPhysicalBoundaryConditions(), as a SAMRAI::xfer::RefineSchedule does. The physical domain must be a
 * single rectangular box, and the patch data must have depth one.
 */
class INSStaggeredDivergenceFreePhysBdryOp : public IBTK::RobinPhysBdryPatchStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param fluid_solver    Integrator that provides the boundary conditions. It must outlive this object.
     * \param homogeneous_bc  Whether to use homogeneous boundary conditions.
     */
    INSStaggeredDivergenceFreePhysBdryOp(const INSStaggeredHierarchyIntegrator* fluid_solver, bool homogeneous_bc);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredDivergenceFreePhysBdryOp();

    /*!
     * \brief Discard cached data; call this function when the types of the boundary conditions change, and after
     * regridding.
     */
    void clearCache();

    /*!
     * \name Implementation of SAMRAI::xfer::RefinePatchStrategy and
     * IBTK::RobinPhysBdryPatchStrategy interfaces.
     */
    //\{

    /*!
     * \brief Fill the velocity ghost values outside the physical boundaries.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    /*!
     * \brief Return the number of ghost layers beyond those that are filled that the extension reads.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override;

    /*!
     * \brief Accumulate the values outside the physical boundaries into the
     * values they depend on using the transpose of the extension with
     * homogeneous boundary conditions, and set them to zero.
     *
     * The caller combines the contributions to ghost cells inside the domain with those of the patches that own them.
     * A ghost value that depends on a face outside the ghost box of the patch data cannot be transposed: the function
     * emits an error if such a value is not zero.
     */
    void accumulateFromPhysicalBoundaryData(SAMRAI::hier::Patch<NDIM>& patch,
                                            double fill_time,
                                            const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    //\}

private:
    /*!
     * \brief A Tape is the list of linear assignments that computes the ghost values of one patch, built once and
     * replayed.
     */
    struct Tape;

    /*!
     * \brief Construction of a Tape.
     */
    class TapeBuilder;

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredDivergenceFreePhysBdryOp() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredDivergenceFreePhysBdryOp(const INSStaggeredDivergenceFreePhysBdryOp& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     */
    INSStaggeredDivergenceFreePhysBdryOp& operator=(const INSStaggeredDivergenceFreePhysBdryOp& that) = delete;

    /*!
     * \brief Return the tape for the velocity patch data patch_data_idx of
     * patch, which is computed once and stored if the patch is in the
     * hierarchy.
     */
    std::shared_ptr<const Tape> getTape(SAMRAI::hier::Patch<NDIM>& patch,
                                        int patch_data_idx,
                                        const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
                                        double fill_time);

    /*!
     * The integrator that provides the boundary conditions.
     */
    const INSStaggeredHierarchyIntegrator* const d_fluid_solver;

    /*!
     * The tapes of the patches in the hierarchy, by level number, patch box,
     * and the ghost widths of the patch data and of the fill.
     */
    std::map<std::array<int, 4 * NDIM + 1>, std::shared_ptr<const Tape>> d_tape_cache;

    /*!
     * The physical domain, with its periodic images, by refinement ratio.
     */
    std::map<std::array<int, NDIM>, SAMRAI::hier::BoxArray<NDIM>> d_physical_domain;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp
