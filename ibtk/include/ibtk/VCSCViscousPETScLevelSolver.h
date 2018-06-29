// Filename: VCSCViscousPETScLevelSolver.h
// Created on 24 Aug 2017 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Nishant Nangia
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBTK_VCSCViscousPETScLevelSolver
#define included_IBTK_VCSCViscousPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

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
    VCSCViscousPETScLevelSolver(const std::string& object_name,
                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                const std::string& default_options_prefix);

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
                                          const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * \brief Copy solution and right-hand-side data to the PETSc
     * representation, including any modifications to account for boundary
     * conditions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCSCViscousPETScLevelSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCSCViscousPETScLevelSolver(const VCSCViscousPETScLevelSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCSCViscousPETScLevelSolver& operator=(const VCSCViscousPETScLevelSolver& that);

    /*
     * The interpolation type to be used for viscosity
     */
    IBTK::VCInterpType d_mu_interp_type;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_VCSCViscousPETScLevelSolver
