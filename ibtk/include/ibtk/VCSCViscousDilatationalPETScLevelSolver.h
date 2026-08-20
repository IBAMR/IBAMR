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

#ifndef included_IBTK_VCSCViscousDilatationalPETScLevelSolver
#define included_IBTK_VCSCViscousDilatationalPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/SCPoissonPETScLevelSolver.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class VCSCViscousDilatationalPETScLevelSolver is a subclass of SCPoissonPETScLevelSolver
 * class which solves vector-valued elliptic equation of the form
 * \f$ \mbox{$L \vec{u}$} = C \vec{u} + \nabla \cdot \mu (\nabla \vec{u} + (\nabla \vec{u})^T) + \nabla (\lambda \nabla
 \cdot \vec{u})= \vec{f} \f$
 * on a \em single SAMRAI::hier::PatchLevel
 * using <A HREF="http://www.mcs.anl.gov/petsc/petsc-as">PETSc</A>.
 *
 * This solver class uses the PETSc library to solve linear equations of the
 * form \f$ C u + \nabla \cdot \mu (\nabla \vec{u} + (\nabla \vec{u})^T) + \nabla (\lambda \nabla \cdot \vec{u}) = f
 \f$,
 * in which \f$ C \f$, \f$ \mu \f$ and \f$ \lambda \f$ are spacially varying coefficients,
 * and \f$\vec{u}\f$ and \f$\vec{f}\f$ are side-centered quantities. The discretization
 * is second-order accurate. For physical problems $\mu$ and $\lambda$ are negative and $C$ is positive.
 * The class does not modify the signs of these quantities, and uses whatever the user
 * sets in the IBTK::ProblemSpecification object.
 *
 * Robin boundary conditions may be specified through the interface class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * The user must perform the following steps to use class
 * VCSCViscousDilatationalPETScLevelSolver:
 *
 * -# Create a VCSCViscousDilatationalPETScLevelSolver object.
 * -# Set the problem specification via setProblemSpecification(),
 *    setPhysicalBcCoef(), and setHomogeneousBc().
 * -# Initialize VCSCViscousDilatationalPETScLevelSolver object using the function
 *    initializeSolverState().
 * -# Solve the linear system using the member function solveSystem(), passing
 *    in SAMRAI::solv::SAMRAIVectorReal objects corresponding to \f$\vec{u}\f$ and
 *    \f$\vec{f}\f$.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 max_iterations = 10            // see setMaxIterations()
 absolute_residual_tol = 0.0    // see setAbsoluteTolerance() (only used by hypre Krylov
 solvers)
 rel_residual_tol = 1.0e-6      // see setRelativeTolerance()
 enable_logging = FALSE         // see setLoggingEnabled()
 options_prefix = ""            // see setOptionsPrefix()
 \endverbatim
 *
 * PETSc is developed at the Argonne National Laboratory Mathematics and
 * Computer Science Division.  For more information about \em PETSc, see <A
 * HREF="http://www.mcs.anl.gov/petsc">http://www.mcs.anl.gov/petsc</A>.
 */
class VCSCViscousDilatationalPETScLevelSolver : public SCPoissonPETScLevelSolver
{
public:
    /*!
     * \brief Constructor.
     */
    VCSCViscousDilatationalPETScLevelSolver(std::string object_name,
                                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                            std::string default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~VCSCViscousDilatationalPETScLevelSolver();

    /*!
     * \brief Static function to construct a VCSCViscousDilatationalPETScLevelSolver.
     */
    static SAMRAI::tbox::Pointer<PoissonSolver> allocate_solver(const std::string& object_name,
                                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                                const std::string& default_options_prefix)
    {
        return new VCSCViscousDilatationalPETScLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \brief Set the interpolation type to be used in interpolating the
     * viscosity.
     */
    void setShearViscosityInterpolationType(IBTK::VCInterpType mu_interp_type);

protected:
    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     */
    void initializeSolverStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                          const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCSCViscousDilatationalPETScLevelSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCSCViscousDilatationalPETScLevelSolver(const VCSCViscousDilatationalPETScLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCSCViscousDilatationalPETScLevelSolver& operator=(const VCSCViscousDilatationalPETScLevelSolver& that) = delete;

    /*
     * The interpolation type to be used for viscosity
     */
    IBTK::VCInterpType d_mu_interp_type;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_VCSCViscousDilatationalPETScLevelSolver
