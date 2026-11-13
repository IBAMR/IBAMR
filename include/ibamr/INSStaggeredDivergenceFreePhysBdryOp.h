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
#include <vector>

namespace IBAMR
{
class INSHierarchyIntegrator;
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
 * \brief Class INSStaggeredDivergenceFreePhysBdryOp sets side-centered velocity ghost values outside the physical
 * boundaries so that each ghost cell has zero discrete divergence, or the divergence of the adjacent interior cell if
 * the fluid solver has a divergence function.
 *
 * The boundary conditions come from the fluid solver and must prescribe the velocity or the traction. Ghost values that
 * depend on a face outside the patch box grown by the filled ghost width are NaN. Requirements:
 * - The domain is a single box at least as wide as the filled ghost width along each non-periodic direction.
 * - Each patch that touches a physical boundary is at least as wide as the filled ghost width normal to it.
 * - The velocity patch data have depth one.
 * - With a variable viscosity, the interpolated viscosity data of the fluid solver have at least the filled ghost
 *   width.
 * - A change of condition type takes effect after clearCache().
 */
class INSStaggeredDivergenceFreePhysBdryOp : public IBTK::RobinPhysBdryPatchStrategy
{
public:
    /*!
     * \brief Constructor.
     *
     * \param fluid_solver    Integrator that provides the boundary conditions. The integrator must outlive this object.
     * \param homogeneous_bc  Whether to use homogeneous boundary conditions.
     */
    INSStaggeredDivergenceFreePhysBdryOp(const INSHierarchyIntegrator* fluid_solver, bool homogeneous_bc);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredDivergenceFreePhysBdryOp();

    /*!
     * \brief Discard the stored boundary extensions. Call this function when the types of the boundary conditions
     * change.
     */
    void clearCache();

    /*!
     * \name Implementation of SAMRAI::xfer::RefinePatchStrategy and
     * IBTK::RobinPhysBdryPatchStrategy interfaces.
     */
    //\{

    /*!
     * \brief Emit an error, because this class uses the velocity boundary conditions of the fluid solver.
     */
    void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef) override;

    /*!
     * \brief Emit an error, because this class uses the velocity boundary conditions of the fluid solver.
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs) override;

    /*!
     * \brief Set the velocity ghost values outside the physical boundaries within ghost_width_to_fill. The fill neither
     * uses nor changes faces beyond it. The caller sets the ghost values inside the domain first.
     */
    void setPhysicalBoundaryConditions(SAMRAI::hier::Patch<NDIM>& patch,
                                       double fill_time,
                                       const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    /*!
     * \brief Return zero: the fill uses no ghost layer beyond those that it fills.
     */
    SAMRAI::hier::IntVector<NDIM> getRefineOpStencilWidth() const override;

    /*!
     * \brief Apply the transpose of the homogeneous extension: add each value outside the physical boundaries to the
     * values that it depends on, and set it to zero. A nonzero value that depends on a face outside the patch box grown
     * by ghost_width_to_fill is an error. The caller combines the contributions to ghost cells inside the domain with
     * those of the patches that own them.
     */
    void accumulateFromPhysicalBoundaryData(SAMRAI::hier::Patch<NDIM>& patch,
                                            double fill_time,
                                            const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill) override;

    //\}

private:
    /*!
     * \brief The linear assignments that set the ghost values of one patch, built once and reused for each fill.
     */
    struct Extension;

    /*!
     * \brief Construction of an Extension.
     */
    class ExtensionBuilder;

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
     * \brief Return the Extension for patch and the velocity patch data patch_data_idx, building and storing it if
     * absent.
     */
    std::shared_ptr<const Extension> getExtension(SAMRAI::hier::Patch<NDIM>& patch,
                                                  int patch_data_idx,
                                                  const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
                                                  double fill_time);

    /*!
     * The integrator that provides the boundary conditions.
     */
    const INSHierarchyIntegrator* const d_fluid_solver;

    /*!
     * The Extensions, keyed by refinement ratio, patch box, and the ghost widths of the patch data and of the fill.
     */
    std::map<std::array<int, 5 * NDIM>, std::shared_ptr<const Extension>> d_extension_cache;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp
