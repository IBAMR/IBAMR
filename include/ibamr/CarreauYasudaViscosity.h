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

#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <SideVariable.h>
#include <Variable.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CarreauYasudaViscosity sets the viscosity of an INSVCStaggeredHierarchyIntegrator from the
 * Carreau-Yasuda model of a generalized Newtonian fluid.
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
 * The constructor registers a cell-centered viscosity variable and a viscosity reset function with the integrator,
 * which must be constructed with <tt>mu_is_const = FALSE</tt> and must not already have a viscosity variable. At the
 * start of each cycle of a time step, the integrator then evaluates the viscosity at the new time from the most recent
 * approximation of the velocity at the new time, so the cycles form a fixed-point iteration for the nonlinear viscous
 * term. The viscosity is first evaluated when the first time step begins, so it is zero before then unless viscosity
 * initial conditions are registered.
 *
 * The object keeps a reference to the integrator, which it does not own, and must exist for as long as the integrator
 * advances the solution.
 */
class CarreauYasudaViscosity
{
public:
    /*!
     * Read the model parameters from \p input_db and register the viscosity with \p ins_integrator.
     */
    CarreauYasudaViscosity(std::string object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           INSVCStaggeredHierarchyIntegrator& ins_integrator);

    /*!
     * Deleted default constructor.
     */
    CarreauYasudaViscosity() = delete;

    /*!
     * Deleted copy constructor: the integrator stores the address of this object.
     */
    CarreauYasudaViscosity(const CarreauYasudaViscosity& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    CarreauYasudaViscosity& operator=(const CarreauYasudaViscosity& that) = delete;

private:
    /*!
     * Set the interior values of the viscosity with patch data index \p mu_idx from the shear rate of the velocity
     * with patch data index \p U_idx at time \p time.
     */
    void computeViscosity(int mu_idx, int U_idx, double time);

    /*!
     * Reset the viscosity from the velocity in the same context; \p ctx is this object.
     *
     * \see INSVCStaggeredHierarchyIntegrator::ResetFluidPropertiesFcnPtr
     */
    static void reset_viscosity(int mu_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> mu_var,
                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                int cycle_num,
                                double time,
                                double current_time,
                                double new_time,
                                void* ctx);

    std::string d_object_name;

    /*!
     * Integrator that calls this object, which this object does not own.
     */
    INSVCStaggeredHierarchyIntegrator& d_ins_integrator;

    /*!
     * Carreau-Yasuda model parameters.
     */
    double d_zero_shear_viscosity, d_infinite_shear_viscosity, d_relaxation_time, d_power_law_index, d_yasuda_exponent;

    /*!
     * Viscosity variable registered with the integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_mu_var;

    /*!
     * Scratch velocity with ghost cells.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double>> d_U_scratch_var;
    int d_U_scratch_idx = IBTK::invalid_index;

    /*!
     * Scratch strain rate tensor, which stores the NDIM diagonal components followed by the off-diagonal components.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_EE_var;
    int d_EE_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_CarreauYasudaViscosity
