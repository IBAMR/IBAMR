// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_VCCollocatedStokesProjectionPreconditioner
#define included_IBAMR_VCCollocatedStokesProjectionPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/LinearSolver.h"

#include "CellVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

namespace IBTK
{
class PoissonSolver;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class VCCollocatedStokesProjectionPreconditioner is a concrete StokesSolver
 * that implements a collocated-grid discretization projection solver for the
 * incompressible Stokes operator with variable coefficients.
 *
 */
class VCCollocatedStokesProjectionPreconditioner : public IBTK::LinearSolver
{
public:
    /*!
     * \brief Class constructor
     */
    VCCollocatedStokesProjectionPreconditioner(const std::string& object_name,
                                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr);

    /*!
     * \brief Destructor.
     */
    ~VCCollocatedStokesProjectionPreconditioner();

    /*!
     * \brief Provide a velocity subdomain solver.
     */
    virtual void setVelocitySubdomainSolver(SAMRAI::tbox::Pointer<IBTK::PoissonSolver> velocity_solver);

    /*!
     * \brief Provide a pressure subdomain solver.
     */
    virtual void setPressureSubdomainSolver(SAMRAI::tbox::Pointer<IBTK::PoissonSolver> pressure_solver);

    /*!
     * Set the pressure stabilization term coefficient (Brezzi–Pitkäranta–type)
     */
    virtual void setPressureStabilizationCoeffient(double eps);

    /*!
     * Set the velocity correction coefficient in the projection method
     */
    virtual void setVelocityCorrectionCoeffient(double alpha);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be nullptr.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be nullptr; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the velocity
     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    /*!
     * \brief Static function to construct a
     * VCCollocatedStokesProjectionPreconditioner.
     */
    static SAMRAI::tbox::Pointer<IBTK::LinearSolver>
    allocate_solver(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    {
        return new VCCollocatedStokesProjectionPreconditioner(object_name, input_db);
    } // allocate_solver

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Compute the action of the preconditioner.
     */
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note It is safe to call initializeSolverState() when the solver state is
     * already initialized.
     *
     * \see deallocateSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    void deallocateSolverState() override;

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    void setInitialGuessNonzero(bool initial_guess_nonzero = true) override;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    void setMaxIterations(int max_iterations) override;

    //\}
protected:
    /*!
     * \brief Remove components in operator null space.
     */
    void correctNullSpace(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > U_vec,
                          SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > P_vec);

    // Pressure stabilization coefficient
    double d_eps = std::numeric_limits<double>::signaling_NaN();

    // Velocity correction coefficient
    double d_alpha = 1.0;

    // Subdomain solvers.
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_pressure_solver;

    // Hierarchy data.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > d_velocity_data_ops, d_pressure_data_ops;
    int d_velocity_wgt_idx = IBTK::invalid_index, d_pressure_wgt_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_U_bc_coefs;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCCollocatedStokesProjectionPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCCollocatedStokesProjectionPreconditioner(const VCCollocatedStokesProjectionPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCCollocatedStokesProjectionPreconditioner&
    operator=(const VCCollocatedStokesProjectionPreconditioner& that) = delete;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_Phi_bdry_fill_op, d_no_fill_op;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Phi_var, d_F_Phi_var;
    int d_Phi_scratch_idx, d_F_Phi_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_VCCollocatedStokesProjectionPreconditioner
