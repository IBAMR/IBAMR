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
 * The operator reads the boundary conditions registered with the fluid solver: (a, b) = (1, 0) prescribes the velocity
 * g and (0, 1) the traction (or pseudo-traction) g; other values are an error. It ignores the conditions on the normal
 * component, and setPhysicalBcCoef() and setPhysicalBcCoefs() emit an error. The domain must be a single box and the
 * patch data must have depth one. To fill ghost width w, the patch data need ghost width w + 1 and the domain must be
 * at least that wide along each non-periodic direction. A ghost value that depends on a face outside the ghost box is
 * set to quiet NaN. The caller fills the values inside the domain that a patch stores in its ghost cells before
 * setPhysicalBoundaryConditions().
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
     * \brief Discard the stored Extensions; call this function when the types of the boundary conditions change, and
     * after regridding.
     */
    void clearCache();

    /*!
     * \name Implementation of SAMRAI::xfer::RefinePatchStrategy and
     * IBTK::RobinPhysBdryPatchStrategy interfaces.
     */
    //\{

    /*!
     * \brief Emit an error; this class uses the velocity boundary conditions of the fluid solver.
     */
    void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef) override;

    /*!
     * \brief Emit an error; this class uses the velocity boundary conditions of the fluid solver.
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs) override;

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
     * \brief Accumulate the values outside the physical boundaries into the values they depend on, with the transpose
     * of the extension with homogeneous boundary conditions, and set them to zero.
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
     * \brief The list of linear assignments that gives the ghost values of one patch, built once and applied again for
     * each fill.
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
     * \brief Return the Extension for the velocity patch data patch_data_idx of patch; it is stored if the patch is in
     * the hierarchy.
     */
    std::shared_ptr<const Extension> getExtension(SAMRAI::hier::Patch<NDIM>& patch,
                                                  int patch_data_idx,
                                                  const SAMRAI::hier::IntVector<NDIM>& ghost_width_to_fill,
                                                  double fill_time);

    /*!
     * The integrator that provides the boundary conditions.
     */
    const INSStaggeredHierarchyIntegrator* const d_fluid_solver;

    /*!
     * The Extensions of the patches in the hierarchy, by level number, patch box, and the ghost widths of the patch
     * data and of the fill.
     */
    std::map<std::array<int, 4 * NDIM + 1>, std::shared_ptr<const Extension>> d_extension_cache;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredDivergenceFreePhysBdryOp
