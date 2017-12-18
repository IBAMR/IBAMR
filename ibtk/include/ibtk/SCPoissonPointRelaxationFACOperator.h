// Filename: SCPoissonPointRelaxationFACOperator.h
// Created on 13 Nov 2008 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBTK_SCPoissonPointRelaxationFACOperator
#define included_IBTK_SCPoissonPointRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <map>
#include <string>
#include <vector>

#include "ibtk/PoissonFACPreconditioner.h"
#include "ibtk/PoissonFACPreconditionerStrategy.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace boost
{
template <class T, std::size_t N>
class array;
} // namespace boost

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class BoxList;
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
 * \brief Class SCPoissonPointRelaxationFACOperator is a concrete
 * PoissonFACPreconditionerStrategy for solving elliptic equations of the form
 * \f$ \mbox{$L u$} = \mbox{$(C I + \nabla \cdot D \nabla) u$} = f \f$ using a
 * globally second-order accurate side-centered finite-difference
 * discretization, with support for Robin and periodic boundary conditions.
 *
 * This class provides operators that are used by class FACPreconditioner to
 * solve scalar Poisson-type equations of the form \f[ (C I + \nabla \cdot D
 * \nabla) u = f \f] using a side-centered, globally second-order accurate
 * finite-difference discretization, where
 *
 * - \f$ C \f$, \f$ D \f$ and \f$ f \f$ are independent of \f$ u \f$,
 * - \f$ C \f$ is a constant damping factor,
 * - \f$ D \f$ is a constant diffusion coefficient, and
 * - \f$ f \f$ is a side-centered scalar function.
 *
 * Robin boundary conditions may be specified at physical boundaries; see class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * By default, the class is configured to solve the Poisson problem \f$
 * -\nabla^2 u = f \f$, subject to homogeneous Dirichlet boundary conditions.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_type = "PATCH_GAUSS_SEIDEL"         // see setSmootherType()
 prolongation_method = "CONSTANT_REFINE"      // see setProlongationMethod()
 restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethod()
 coarse_solver_type = "HYPRE_LEVEL_SOLVER"    // see setCoarseSolverType()
 coarse_solver_rel_residual_tol = 1.0e-5      // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50     // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 1             // see setCoarseSolverMaxIterations()
 coarse_solver_db = { ... }                   // SAMRAI::tbox::Database for initializing coarse
 level solver
 \endverbatim
*/
class SCPoissonPointRelaxationFACOperator : public PoissonFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    SCPoissonPointRelaxationFACOperator(const std::string& object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                        const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~SCPoissonPointRelaxationFACOperator();

    /*!
     * \brief Static function to construct a PoissonFACPreconditioner with a
     * SCPoissonPointRelaxationFACOperator FAC strategy.
     */
    static SAMRAI::tbox::Pointer<PoissonSolver> allocate_solver(const std::string& object_name,
                                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                                const std::string& default_options_prefix)
    {
        SAMRAI::tbox::Pointer<PoissonFACPreconditionerStrategy> fac_operator = new SCPoissonPointRelaxationFACOperator(
            object_name + "::SCPoissonPointRelaxationFACOperator", input_db, default_options_prefix);
        return new PoissonFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
    } // allocate

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Specify the smoother type.
     *
     * Select from:
     * - \c "PATCH_GAUSS_SEIDEL"
     * - \c "PROCESSOR_GAUSS_SEIDEL"
     * - \c "RED_BLACK_GAUSS_SEIDEL"
     */
    void setSmootherType(const std::string& smoother_type);

    /*!
     * \brief Specify the coarse level solver.
     */
    void setCoarseSolverType(const std::string& coarse_solver_type);

    //\}

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are
     *being
     *performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are
     *being
     *performed
     */
    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps);

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                            int coarsest_ln);

    /*!
     * \brief Compute composite grid residual on the specified range of levels.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param coarsest_level_num coarsest level number
     * \param finest_level_num finest level number
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num);

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                            int coarsest_reset_ln,
                                            int finest_reset_ln);

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln);

    /*
     * Coarse level solvers and solver parameters.
     */
    SAMRAI::tbox::Pointer<PoissonSolver> d_coarse_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_coarse_solver_db;

    /*
     * Patch overlap data.
     */
    std::vector<std::vector<boost::array<SAMRAI::hier::BoxList<NDIM>, NDIM> > > d_patch_bc_box_overlap;
    std::vector<std::vector<boost::array<std::map<int, SAMRAI::hier::Box<NDIM> >, NDIM> > > d_patch_neighbor_overlap;

    /*
     * Dirichlet boundary condition utilities.
     */
    SAMRAI::tbox::Pointer<StaggeredPhysicalBoundaryHelper> d_bc_helper;
    int d_mask_idx;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SCPoissonPointRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SCPoissonPointRelaxationFACOperator(const SCPoissonPointRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SCPoissonPointRelaxationFACOperator& operator=(const SCPoissonPointRelaxationFACOperator& that);

};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_SCPoissonPointRelaxationFACOperator
