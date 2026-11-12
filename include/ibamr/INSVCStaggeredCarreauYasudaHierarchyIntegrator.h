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

#ifndef included_IBAMR_INSVCStaggeredCarreauYasudaHierarchyIntegrator
#define included_IBAMR_INSVCStaggeredCarreauYasudaHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <PatchHierarchy.h>
#include <Variable.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSVCStaggeredCarreauYasudaHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations of a generalized Newtonian fluid whose
 * viscosity follows the Carreau-Yasuda model.
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
 * The integrator owns the cell-centered viscosity variable, which getViscosityVariable() returns. At the start of
 * each cycle of a time step, it evaluates the viscosity at the new time from the most recent approximation of the
 * velocity at the new time, so the cycles form a fixed-point iteration for the nonlinear viscous term. The
 * viscosity cannot be constant (<tt>mu_is_const</tt> must be <tt>FALSE</tt>), and callers must not register another
 * viscosity variable, viscosity initial conditions, or a transported viscosity variable.
 */
class INSVCStaggeredCarreauYasudaHierarchyIntegrator : public INSVCStaggeredNonConservativeHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSVCStaggeredCarreauYasudaHierarchyIntegrator reads the viscosity model
     * parameters and registers the viscosity variable in addition to the configuration read by
     * INSVCStaggeredNonConservativeHierarchyIntegrator.
     */
    INSVCStaggeredCarreauYasudaHierarchyIntegrator(std::string object_name,
                                                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                   bool register_for_restart = true);

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

    /*!
     * Initialize the AMR patch hierarchy and data defined on the hierarchy at
     * the start of a computation, including the viscosity of the initial velocity.
     */
    void initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM>> gridding_alg) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredCarreauYasudaHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredCarreauYasudaHierarchyIntegrator(const INSVCStaggeredCarreauYasudaHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredCarreauYasudaHierarchyIntegrator&
    operator=(const INSVCStaggeredCarreauYasudaHierarchyIntegrator& that) = delete;

    /*!
     * Set the interior values of the viscosity with patch data index \p mu_idx from the shear rate of the velocity
     * with patch data index \p U_idx at time \p time.
     */
    void computeViscosity(int mu_idx, int U_idx, double time);

    /*!
     * Reset the viscosity from the velocity in the same context; \p ctx is the integrator.
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

    /*!
     * Carreau-Yasuda model parameters.
     */
    double d_zero_shear_viscosity, d_infinite_shear_viscosity, d_relaxation_time, d_power_law_index, d_yasuda_exponent;

    /*!
     * Scratch strain rate tensor, which stores the NDIM diagonal components followed by the off-diagonal components.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_strain_rate_var;
    int d_strain_rate_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSVCStaggeredCarreauYasudaHierarchyIntegrator
