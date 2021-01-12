// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_VCSCViscousPETScLevelSolver
#define included_IBTK_VCSCViscousPETScLevelSolver

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
 * \brief Class VCSCViscousPETScLevelSolver is a subclass of SCPoissonPETScLevelSolver
 * class which solves vector-valued elliptic equation of the form
 * \f$ \mbox{$L u$} = C u + \nabla \cdot \mu (\nabla u + (\nabla u)^T) = f \f$
 * on a \em single SAMRAI::hier::PatchLevel
 * using <A HREF="http://www.mcs.anl.gov/petsc/petsc-as">PETSc</A>.
 *
 * This solver class uses the PETSc library to solve linear equations of the
 * form \f$ \beta C u  + \alpha \nabla \cdot \mu (\nabla u + (\nabla u)^T) = f \f$,
 * in which \f$ C \f$ and \f$ \mu \f$ are spacially varying coefficients,
 * and \f$u\f$ and \f$f\f$ are side-centered arrays. The scaling factors of
 * \f$ C \f$ and \f$ \mu \f$ are stored separately in the class and are
 * denoted by \f$ \beta \f$ and \f$ \alpha \f$, respectively. The discretization
 * is second-order accurate.
 *
 * Robin boundary conditions may be specified through the interface class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * The user must perform the following steps to use class
 * VCSCViscousPETScLevelSolver:
 *
 * -# Create a VCSCViscousPETScLevelSolver object.
 * -# Set the problem specification via setPoissonSpecifications(),
 *    setPhysicalBcCoef(), and setHomogeneousBc().
 * -# Initialize VCSCViscousPETScLevelSolver object using the function
 *    initializeSolverState().
 * -# Solve the linear system using the member function solveSystem(), passing
 *    in SAMRAI::solv::SAMRAIVectorReal objects corresponding to \f$u\f$ and
 *    \f$f\f$.
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
class VCSCViscousPETScLevelSolver : public SCPoissonPETScLevelSolver
{
public:
    /*!
     * \brief Constructor.
     */
    VCSCViscousPETScLevelSolver(std::string object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                std::string default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~VCSCViscousPETScLevelSolver();

    /*!
     * \brief Static function to construct a VCSCViscousPETScLevelSolver.
     */
    static SAMRAI::tbox::Pointer<PoissonSolver> allocate_solver(const std::string& object_name,
                                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                                const std::string& default_options_prefix)
    {
        return new VCSCViscousPETScLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \brief Set the interpolation type to be used in interpolating the
     * viscosity.
     */
    void setViscosityInterpolationType(IBTK::VCInterpType mu_interp_type);

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
    VCSCViscousPETScLevelSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCSCViscousPETScLevelSolver(const VCSCViscousPETScLevelSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCSCViscousPETScLevelSolver& operator=(const VCSCViscousPETScLevelSolver& that) = delete;

    /*
     * The interpolation type to be used for viscosity
     */
    IBTK::VCInterpType d_mu_interp_type;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_VCSCViscousPETScLevelSolver
