// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_StaggeredStokesBlockPreconditioner
#define included_IBAMR_StaggeredStokesBlockPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StaggeredStokesSolver.h"

#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearSolver.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/ibtk_utilities.h"

#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesBlockPreconditioner is an abstract base class for
 * Stokes solvers that are implemented using block (subdomain) solvers.
 */
class StaggeredStokesBlockPreconditioner : public IBTK::LinearSolver, public StaggeredStokesSolver
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesBlockPreconditioner(bool needs_velocity_solver, bool needs_pressure_solver);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesBlockPreconditioner() = default;

    /*!
     * \brief Indicate whether the preconditioner needs a velocity subdomain
     * solver.
     */
    virtual bool needsVelocitySubdomainSolver() const;

    /*!
     * \brief Provide a velocity subdomain solver.
     */
    virtual void setVelocitySubdomainSolver(SAMRAI::tbox::Pointer<IBTK::PoissonSolver> velocity_solver);

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs) override;

    /*!
     * \brief Indicate whether the preconditioner needs a pressure subdomain
     * solver.
     */
    virtual bool needsPressureSubdomainSolver() const;

    /*!
     * \brief Provide a pressure subdomain solver.
     */
    virtual void setPressureSubdomainSolver(SAMRAI::tbox::Pointer<IBTK::PoissonSolver> pressure_solver);

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the pressure-Poisson problem.
     */
    virtual void setPressurePoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& P_problem_coefs);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the velocity
     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef) override;

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

protected:
    /*!
     * \brief Remove components in operator null space.
     */
    void correctNullspace(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > U_vec,
                          SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > P_vec);

    // Subdomain solvers.
    const bool d_needs_velocity_solver;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver;
    SAMRAI::solv::PoissonSpecifications d_P_problem_coefs;
    const bool d_needs_pressure_solver;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_pressure_solver;

    // Hierarchy data.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > d_velocity_data_ops, d_pressure_data_ops;
    int d_velocity_wgt_idx = IBTK::invalid_index, d_pressure_wgt_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesBlockPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesBlockPreconditioner(const StaggeredStokesBlockPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesBlockPreconditioner& operator=(const StaggeredStokesBlockPreconditioner& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StaggeredStokesBlockPreconditioner
