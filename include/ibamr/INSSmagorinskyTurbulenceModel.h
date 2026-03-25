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

#ifndef included_IBAMR_INSSmagorinskyTurbulenceModel
#define included_IBAMR_INSSmagorinskyTurbulenceModel

#include <ibamr/config.h>

#include <ibamr/INSSGSKinematics.h>
#include <ibamr/INSSGSStressData.h>
#include <ibamr/INSTurbulenceModel.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <limits>

namespace IBAMR
{
/*!
 * \brief A Smagorinsky-type LES closure for INSStaggeredHierarchyIntegrator.
 *
 * The model adds the explicit forcing \f$\nabla \cdot (2 \mu_t S)\f$ with
 * \f$\mu_t = \rho (C_s \Delta)^2 |S|\f$.
 */
class INSSmagorinskyTurbulenceModel : public INSTurbulenceModel
{
public:
    /*!
     * \brief Constructor.
     */
    INSSmagorinskyTurbulenceModel(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Compute the explicit Smagorinsky forcing contribution.
     */
    void computeTurbulenceForce(int F_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> F_var,
                                int U_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> U_var,
                                const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc_coefs,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                double data_time,
                                const StokesSpecifications& problem_coefs) override;

private:
    /*!
     * Helper object providing ghost-filled velocity data and resolved
     * strain-rate reconstruction.
     */
    SAMRAI::tbox::Pointer<INSSGSKinematics> d_kinematics;

    /*!
     * Staggered SGS stress workspace used to store diagonal and shear
     * components on solver-appropriate locations.
     */
    SAMRAI::tbox::Pointer<INSSGSStressData> d_tau_sgs_data;

    /*!
     * Smagorinsky coefficient \f$C_s\f$.
     */
    double d_smagorinsky_constant = 0.17;

    /*!
     * Multiplicative factor used to define the filter width from the local grid
     * spacing.
     */
    double d_filter_width_scale = 1.0;

    /*!
     * Upper bound for the modeled turbulent viscosity.
     */
    double d_max_turbulent_viscosity = std::numeric_limits<double>::max();
};
} // namespace IBAMR

#endif
