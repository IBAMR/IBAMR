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

#ifndef included_IBAMR_CarreauYasudaViscosity
#define included_IBAMR_CarreauYasudaViscosity

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/ViscosityModel.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <CellVariable.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CarreauYasudaViscosity computes the viscosity of a generalized Newtonian fluid from the Carreau-Yasuda
 * model.
 *
 * The viscosity is
 * \f$ \mu(\dot{\gamma}) = \mu_\infty + (\mu_0 - \mu_\infty) (1 + (\lambda \dot{\gamma})^a)^{(n-1)/a} \f$,
 * in which \f$ \dot{\gamma} = \sqrt{2 E : E} \f$ is the shear rate and
 * \f$ E = \frac{1}{2}(\nabla u + (\nabla u)^T) \f$ is the strain rate tensor.
 * The parameters are read from the required input keys
 *
 * - <tt>zero_shear_viscosity</tt> (\f$ \mu_0 > 0 \f$),
 * - <tt>infinite_shear_viscosity</tt> (\f$ \mu_\infty \ge 0 \f$),
 * - <tt>relaxation_time</tt> (\f$ \lambda \ge 0 \f$),
 * - <tt>power_law_index</tt> (\f$ n > 0 \f$), and
 * - <tt>yasuda_exponent</tt> (\f$ a > 0 \f$).
 *
 * Give the model to INSVCStaggeredHierarchyIntegrator::setViscosityModel(). Because the integrator evaluates the
 * viscosity at the new time at the start of each cycle of a time step from the most recent approximation of the
 * velocity at the new time, the cycles form a fixed-point iteration for the nonlinear viscous term.
 */
class CarreauYasudaViscosity : public ViscosityModel
{
public:
    /*!
     * Read the model parameters from \p input_db.
     */
    CarreauYasudaViscosity(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * The centered differences of the strain rate read one ghost cell of the velocity.
     */
    int getVelocityGhostWidth() const override;

    /*!
     * Set the viscosity from the shear rate of the velocity.
     */
    void setDataOnPatchHierarchy(int mu_idx, int U_idx, IBTK::HierarchyMathOps& hier_math_ops, double time) override;

private:
    std::string d_object_name;

    /*!
     * Carreau-Yasuda model parameters.
     */
    double d_zero_shear_viscosity, d_infinite_shear_viscosity, d_relaxation_time, d_power_law_index, d_yasuda_exponent;

    /*!
     * Scratch strain rate tensor, which stores the NDIM diagonal components followed by the off-diagonal components.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_EE_var;
    int d_EE_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CarreauYasudaViscosity
