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
 * \brief Class INSStaggeredDivergenceFreePhysBdryOp is a
 * SAMRAI::xfer::RefinePatchStrategy that fills the ghost values of a
 * side-centered velocity outside the physical boundaries of a rectangular
 * domain so that the discrete divergence
 * \f$ \sum_d (u_d(\mathrm{hi}) - u_d(\mathrm{lo}))/\Delta x_d \f$ vanishes in
 * every ghost cell outside the domain, for any values of the velocity in the
 * domain.
 *
 * The velocity boundary conditions are those that INSStaggeredVelocityBcCoef
 * interprets, specified by the user-supplied coefficients \f$ a \f$, \f$ b
 * \f$, and \f$ g \f$ that INSHierarchyIntegrator::getPhysicalBoundaryConditions()
 * returns: \f$ a = 1 \f$ prescribes the velocity \f$ g \f$, and \f$ b = 1 \f$
 * prescribes the traction (or pseudo-traction) \f$ g \f$ according to the
 * traction boundary condition type of the integrator. Any other values are an
 * error. Where the condition is of the second kind, the corner treatment
 * TractionBcCornerType of the integrator is used at the ends of a boundary.
 *
 * <h2>Definition</h2>
 *
 * A boundary is a lower or upper side of a non-periodic axis. A ghost face
 * lies beyond a set \f$ S \f$ of boundaries, which is its region; faces in the
 * domain, including the faces on its boundaries, have an empty region and keep
 * their values. The value of a ghost face is the mean over the boundaries
 * \f$ b \in S \f$ of \f$ W_b \f$, which treats \f$ b \f$ last. It takes the
 * values in the smaller regions as its interior, and
 * - a velocity component tangential to \f$ b \f$ is set by the one-dimensional
 *   rule for \f$ b \f$ from the value in the mirror cell across \f$ b \f$:
 *   twice the prescribed velocity minus the mirror value for a prescribed
 *   velocity, and the mirror value plus
 *   \f$ (2k+1)\Delta x_b \gamma \f$ at depth \f$ k \f$ (\f$ k = 0 \f$ next to
 *   the boundary) for a traction condition, where \f$ \gamma = \pm (g/\mu -
 *   D_t u_b) \f$ with the sign + for an upper boundary, \f$ D_t u_b \f$ the
 *   tangential difference of the normal velocity on \f$ b \f$ for TRACTION, and
 *   \f$ D_t u_b = 0 \f$ for PSEUDO_TRACTION;
 * - the component normal to \f$ b \f$ is set by marching
 *   \f$ \mathrm{Div}_h u = 0 \f$ outward from the boundary one ghost cell at a
 *   time, using \f$ W_b \f$ for the faces of the marched cell in region
 *   \f$ S \f$ and the stored values for the faces in smaller regions. The
 *   user's condition on the normal component is not used: the boundary face
 *   keeps whatever value the patch data holds.
 *
 * For a TRACTION condition, a face of the difference \f$ D_t u_b \f$ that lies
 * one cell beyond the end of \f$ b \f$ across a boundary that is not in \f$ S
 * \f$ is determined by the corner treatment. All other faces use the stored
 * values. Because the order average is symmetric, the values do not depend on
 * the order in which the boundaries are processed.
 *
 * <h2>Use</h2>
 *
 * The values depend only on the velocity in the domain and on the boundary
 * data, so every patch that stores a ghost face outside the domain computes the
 * same value for it. The values within the domain that a patch stores in its
 * ghost cells, such as copies from other patches, must be filled before
 * setPhysicalBoundaryConditions() is called, as a SAMRAI::xfer::RefineSchedule
 * does. The ghost faces filled are those outside the domain within
 * the ghost width of the patch data and within ghost_width_to_fill of the
 * patch.
 *
 * <h2>Ghost width</h2>
 *
 * A ghost value that depends on a face outside the ghost box of the patch data
 * is not computable. setPhysicalBoundaryConditions() sets such a value to quiet
 * NaN, and accumulateFromPhysicalBoundaryData() leaves it unchanged. A value at depth \f$ k \f$ depends on
 * values in the domain up to depth \f$ k \f$ from the boundary (the mirror
 * cell), and, near the ends of a TRACTION boundary, on one more layer
 * tangentially, so all values within a ghost width g are computable if the
 * patch data has ghost width g + 1 and the width of the domain along each
 * non-periodic direction is at least that ghost width.
 *
 * <h2>Transpose</h2>
 *
 * accumulateFromPhysicalBoundaryData() applies the exact transpose of the
 * extension by the homogeneous boundary conditions. It accumulates the values
 * in the ghost cells outside the domain into the values of the faces they
 * depend on, which are faces in the domain anywhere in the ghost box of the
 * patch data, including its ghost cells inside the domain, and then sets to
 * zero the values outside the domain that setPhysicalBoundaryConditions()
 * fills, which are the computable faces within the ghost width of the patch
 * data and within ghost_width_to_fill of the patch. The faces that are not
 * computable and the faces outside ghost_width_to_fill are left unchanged. The
 * caller is responsible for combining contributions to the ghost cells inside
 * the domain with those of the patches that own them.
 *
 * \note Only a single rectangular physical domain box is supported, the
 * coefficient types are read when the values for a patch are first computed and
 * are assumed not to change, and the patch data must have depth one.
 */
class INSStaggeredDivergenceFreePhysBdryOp : public IBTK::RobinPhysBdryPatchStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param patch_data_index  Side-centered patch data index of the velocity.
     * \param fluid_solver      Integrator that provides the boundary conditions, the type
     *                          of traction conditions, the corner treatment, the viscosity,
     *                          and the patch hierarchy. It must outlive this object.
     * \param homogeneous_bc    Whether to use homogeneous boundary conditions.
     */
    INSStaggeredDivergenceFreePhysBdryOp(int patch_data_index,
                                         const INSStaggeredHierarchyIntegrator* fluid_solver,
                                         bool homogeneous_bc);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredDivergenceFreePhysBdryOp();

    /*!
     * \brief Discard the stored descriptions of the ghost values of the
     * patches. A description depends on the level number, the patch box, the
     * ghost widths of the patch data and of the fill, the types of the
     * velocity boundary conditions at the boundary faces, the traction boundary
     * condition type, and the corner treatment. A description is stored under the
     * first four of these, so this function must be called after the boundary
     * condition types, the traction boundary condition type, or the corner
     * treatment change. After the hierarchy is regridded, a patch with the same level
     * number and box has the same description, so calling this function only
     * frees the memory of the descriptions of patches that no longer exist.
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
     * \brief Return the stencil width of the operator.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override;

    /*!
     * \brief Accumulate the values outside the physical boundaries into the
     * values they depend on using the transpose of the extension with
     * homogeneous boundary conditions, and set them to zero.
     */
    void accumulateFromPhysicalBoundaryData(SAMRAI::hier::Patch<NDIM>& patch,
                                            double fill_time,
                                            const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    //\}

private:
    /*!
     * \brief The assignments that compute the ghost values of a patch.
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
     * hierarchy. The coefficient types are determined at fill_time when the
     * tape is computed. The geometry of the patch is replaced temporarily
     * while the coefficients are evaluated.
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
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp
