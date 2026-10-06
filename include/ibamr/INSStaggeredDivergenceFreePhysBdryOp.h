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
 * physical boundaries of a rectangular domain so that the discrete divergence
 * \f$ \sum_d (u_d(\mathrm{hi}) - u_d(\mathrm{lo}))/\Delta x_d \f$ vanishes in every ghost cell outside the domain, for
 * any values of the velocity in the domain.
 *
 * The boundary conditions are the coefficients \f$ a \f$, \f$ b \f$, and \f$ g \f$ registered with the fluid solver:
 * \f$ a = 1 \f$ prescribes the velocity \f$ g \f$ and \f$ b = 1 \f$ prescribes the traction (or pseudo-traction)
 * \f$ g \f$, according to the traction boundary condition type of the fluid solver. Other values are an error. The
 * conditions on the normal component are not used: faces on a boundary keep their values. Because the boundary
 * conditions are read from the fluid solver, setPhysicalBcCoef() and setPhysicalBcCoefs() have no effect.
 *
 * A boundary \f$ B \f$ is the lower or upper side of a non-periodic direction of the domain, and \f$ \Delta x_B \f$ is
 * the cell width normal to it. Faces in the domain keep their values. The value of a ghost face that lies beyond the
 * set \f$ S \f$ of boundaries is
 * -# for \f$ S = \{B\} \f$, working outward from \f$ B \f$: for a component tangential to \f$ B \f$, twice the
 *    prescribed velocity minus the value at the mirror face across \f$ B \f$, or, for a traction condition, the
 *    mirror value plus \f$ (2k+1) \Delta x_B \gamma \f$ at depth \f$ k \f$ (\f$ k = 0 \f$ next to \f$ B \f$), where
 *    \f$ \gamma = \pm (g/\mu - D_t u_B) \f$ with the sign + at an upper boundary and \f$ D_t u_B \f$ the difference
 *    along the component of the normal velocity on \f$ B \f$ (zero for PSEUDO_TRACTION); for the component normal to
 *    \f$ B \f$, the value that makes the divergence of the adjacent ghost cell vanish. The normal velocity on \f$ B \f$
 *    one cell beyond a corner of \f$ B \f$ is, if the adjacent boundary prescribes the normal velocity \f$ u_b \f$,
 *    \f$ 2 u_b - u \f$ with \f$ u \f$ the value at the nearest face on \f$ B \f$, and is otherwise linearly
 *    extrapolated from the two nearest faces on \f$ B \f$ (the nearest face if \f$ B \f$ has no second one);
 * -# otherwise, the average over \f$ B \in S \f$ of the field obtained by applying (1) across \f$ B \f$ to the field
 *    already extended across the other boundaries in \f$ S \f$: the mirror values and the normal velocity on \f$ B \f$
 *    are those extended values, the coefficients of \f$ B \f$ are evaluated on its extension beyond the domain, and the
 *    divergence condition uses the values that the same \f$ B \f$ gives to the other faces of the cell. Each averaged
 *    field is divergence-free, so the average is.
 *
 * A ghost value at depth \f$ k \f$ depends on values up to depth \f$ k \f$ inside \f$ B \f$. On a TRACTION boundary
 * it also depends, through \f$ D_t u_B \f$, on the normal velocity on \f$ B \f$ one cell further along \f$ B \f$ than
 * the ghost value itself, either in the domain or, beyond a corner, extended. The values within ghost width \f$ w \f$
 * are therefore computed if the patch data have ghost width \f$ w + 1 \f$ and the domain is at least that wide along
 * each non-periodic direction. A ghost value that depends on a face outside the ghost box of the patch data is set to
 * quiet NaN by setPhysicalBoundaryConditions(). Values inside the domain that a patch stores in its ghost cells must be
 * filled before setPhysicalBoundaryConditions() is called, as a SAMRAI::xfer::RefineSchedule does. A communication
 * schedule calls a physical boundary operator only on patches that touch the physical boundary, so a patch that lies
 * within the ghost cell width of the boundary without touching it is not filled. The data that describe the ghost
 * values are computed once for each patch of the fluid solver's hierarchy, and are rebuilt at every call for any other
 * patch.
 *
 * accumulateFromPhysicalBoundaryData() applies the transpose of the extension with homogeneous boundary conditions.
 * It adds each ghost value outside the domain to the values of the faces it depends on, which lie anywhere in the
 * ghost box of the patch data inside the domain, and sets it to zero. The caller combines the contributions to the
 * ghost cells inside the domain with those of the patches that own them. A ghost value that depends on a face outside
 * the ghost box of the patch data is not transposed and is left unchanged; a caller that can place nonzero values there
 * must transpose them on a patch whose data contain those faces.
 *
 * \note The physical domain must be a single rectangular box, and the patch data must have depth one.
 */
class INSStaggeredDivergenceFreePhysBdryOp : public IBTK::RobinPhysBdryPatchStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param fluid_solver    Integrator that provides the boundary conditions, the type of traction conditions, the
     *                        viscosity, and the patch hierarchy. It must outlive this object.
     * \param homogeneous_bc  Whether to use homogeneous boundary conditions.
     *
     * \note The side-centered velocity patch data are set by setPatchDataIndex().
     */
    INSStaggeredDivergenceFreePhysBdryOp(const INSStaggeredHierarchyIntegrator* fluid_solver, bool homogeneous_bc);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredDivergenceFreePhysBdryOp();

    /*!
     * \brief Discard cached data. Whether each boundary face has a velocity or a traction condition is determined when
     * a patch is first filled; call this function when that changes or when the traction boundary condition type
     * changes, and after regridding to release the data of patches that have been removed.
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
     * \brief Return the number of ghost layers of the velocity, beyond those that are filled, that the extension
     * reads.
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
     * \brief A Tape is the list of linear assignments that computes every ghost value of one patch from the values in
     * the domain and the boundary data, for patch data with a given ghost box. It is built once per patch and
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

    /*!
     * The physical domain, with its periodic images, refined by each refinement ratio for which it has been requested.
     */
    std::map<std::array<int, NDIM>, SAMRAI::hier::BoxArray<NDIM>> d_physical_domain;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp
