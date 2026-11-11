// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_StokesBcCoefStrategy
#define included_IBAMR_StokesBcCoefStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/ibamr_enums.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/ibtk_utilities.h>

#include <SideIndex.h>

#include <string>
#include <utility>
#include <vector>

namespace IBAMR
{
class StokesSpecifications;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class BoxArray;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StokesBcCoefStrategy is a subclass of the abstract base class
 * IBTK::ExtendedRobinBcCoefStrategy to allow for specialization needed to treat
 * physical boundary conditions for coupled Stokes solvers.
 */
class StokesBcCoefStrategy : public IBTK::ExtendedRobinBcCoefStrategy
{
public:
    /*!
     * \brief Empty default constructor.
     */
    StokesBcCoefStrategy() = default;

    /*!
     * \brief Empty destructor.
     */
    ~StokesBcCoefStrategy() = default;

    /*!
     * \brief Set the StokesSpecifications object used by this boundary condition
     * specification object.
     *
     * \param problem_coefs   Problem coefficients
     */
    virtual void setStokesSpecifications(const StokesSpecifications* problem_coefs);

    /*!
     * \brief Set the target velocity data index to use when setting physical
     * boundary conditions and the time at which it is defined.
     */
    virtual void setTargetVelocityPatchDataIndex(int u_target_data_idx);

    /*!
     * \brief Clear the target velocity data index used when setting physical
     * boundary conditions.
     */
    virtual void clearTargetVelocityPatchDataIndex();

    /*!
     * \brief Set the target pressure data index to use when setting physical
     * boundary conditions and the time at which it is defined.
     */
    virtual void setTargetPressurePatchDataIndex(int u_target_data_idx);

    /*!
     * \brief Clear the target pressure data index used when setting physical
     * boundary conditions.
     */
    virtual void clearTargetPressurePatchDataIndex();

    /*!
     * \brief Set the type of traction boundary conditions.  Supported options
     * are: TRACTION_BOUNDARY_CONDITIONS and
     * PSEUDO_TRACTION_BOUNDARY_CONDITIONS.
     *
     * The default is TRACTION_BOUNDARY_CONDITIONS.
     */
    virtual void setTractionBcType(TractionBcType bc_type);

    /*!
     * \brief Get the type of traction boundary conditions.
     */
    virtual TractionBcType getTractionBcType() const;

    /*!
     * \brief Return the side indices of the normal velocities on which the inhomogeneous coefficient \f$ \gamma \f$ of
     * the tangential TRACTION condition depends at \a i, each with the derivative of \f$ \gamma \f$ with respect to it.
     * An index may appear more than once: its derivative is the sum of its entries.
     *
     * \param i  Location, in the index space of the Robin coefficients, at which the traction is prescribed.
     * \param patch  Patch with the geometry that setBcCoefs() requires for this component.
     * \param bdry_box  Boundary box that contains \a i.
     * \param domain  Physical domain of the level of the patch with its periodic images.
     * \param ghost_box  Ghost box of the data that the returned side indices refer to.
     * \param bdry_interp_type  Interpolation type, "LINEAR" or "QUADRATIC", of the assembled matrix. The type
     * determines the extrapolation of the normal velocity beyond a corner.
     * \param fill_time  Time at which the boundary conditions are evaluated.
     *
     * \note The default implementation returns an empty vector for the PSEUDO_TRACTION type, whose conditions do not
     * depend on the normal velocity, and stops with an error for the TRACTION type.
     */
    virtual std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
    getTractionStencil(const SAMRAI::hier::Index<NDIM>& i,
                       const SAMRAI::hier::Patch<NDIM>& patch,
                       const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                       const SAMRAI::hier::BoxArray<NDIM>& domain,
                       const SAMRAI::hier::Box<NDIM>& ghost_box,
                       const std::string& bdry_interp_type,
                       double fill_time) const;

protected:
    /*
     * Problem coefficients.
     */
    const StokesSpecifications* d_problem_coefs;

    /*!
     * Patch data indices.
     */
    int d_u_target_data_idx = IBTK::invalid_index, d_p_target_data_idx = IBTK::invalid_index;

    /*
     * The type of traction boundary conditions.  Supported options are:
     * TRACTION_BOUNDARY_CONDITIONS and PSEUDO_TRACTION_BOUNDARY_CONDITIONS.
     */
    TractionBcType d_traction_bc_type = TRACTION;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StokesBcCoefStrategy(const StokesBcCoefStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StokesBcCoefStrategy& operator=(const StokesBcCoefStrategy& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StokesBcCoefStrategy
