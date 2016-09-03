// Filename: StaggeredStokesIBLevelRelaxationFACOperator.h
// Created on 22 Mar 2015 by Amneet Bhalla
//
// Copyright (c) 2002-2015, Amneet Bhalla and Boyce Griffith
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

#ifndef included_StaggeredStokesIBLevelRelaxationFACOperator
#define included_StaggeredStokesIBLevelRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <utility>
#include <vector>

#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "SAMRAIVectorReal.h"
#include "VariableContext.h"
#include "VariableFillPattern.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesFACPreconditionerStrategy.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "petscao.h"
#include "petscmat.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

namespace boost
{
template <class T, std::size_t N>
class array;
} // namespace boost
namespace IBAMR
{
class StaggeredStokesPETScLevelSolver;
} // namespace IBAMR
namespace IBTK
{
class HierarchyGhostCellInterpolation;
class HierarchyMathOps;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoxList;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesIBLevelRelaxationFACOperator is a concrete
 * FACPreconditionerStrategy class that implements the operations required by
 * smoothers for staggered-grid (MAC) discretizations of the implicit
 * incompressible Stokes-IB equations.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_type = "ADDITIVE"                     // see setSmootherType()
 U_prolongation_method = "CONSTANT_REFINE"      // see setProlongationMethods()
 P_prolongation_method = "LINEAR_REFINE"        // see setProlongationMethods()
 U_restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethods()
 P_restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethods()
 coarse_solver_type = "LEVEL_SMOOTHER"          // see setCoarseSolverType()
 coarse_solver_rel_residual_tol = 1.0e-5        // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50       // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 10              // see setCoarseSolverMaxIterations()
 coarse_solver_db = { ... }                     // SAMRAI::tbox::Database for initializing
 coarse
 level solver
 \endverbatim
*/
class StaggeredStokesIBLevelRelaxationFACOperator : public StaggeredStokesFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesIBLevelRelaxationFACOperator(const std::string& object_name,
                                                SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                                const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesIBLevelRelaxationFACOperator();

    /*!
     * \brief Static function to construct a StaggeredStokesFACPreconditioner with a
     * StaggeredStokesIBLevelRelaxationFACOperator FAC strategy.
     */
    static SAMRAI::tbox::Pointer<StaggeredStokesSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        SAMRAI::tbox::Pointer<IBTK::FACPreconditionerStrategy> fac_operator =
            new StaggeredStokesIBLevelRelaxationFACOperator(
                object_name + "::StaggeredStokesIBLevelRelaxationFACOperator", input_db, default_options_prefix);
        return new StaggeredStokesFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Set the IB time stepping type.
     */
    void setIBTimeSteppingType(TimeSteppingType time_stepping_type);

    /*!
     * \brief Set the IB-force Jacobian at the finest patch level (where the
     * structure resides).
     */
    void setIBForceJacobian(Mat& A);

    /*!
     * \brief Set the IB-interpolation operator at the finest patch level (where the
     * structure resides).
     */
    void setIBInterpOp(Mat& J);

    //\}

    /*!
     * \name Functions for accessing the level operators and solvers.
     */
    //\{

    /*!
     * \brief Get the Staggered Stokes IB level solver.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesPETScLevelSolver> getStaggeredStokesPETScLevelSolver(int ln) const;

    /*!
     * \brief Get the Eulerian elasticity level operator.
     */
    Mat getEulerianElasticityLevelOp(int ln) const;

    /*!
     * \brief Get the prolongation level operator. The prolongation
     * operator prolongs data from level \em ln to level \em ln + 1.
     */
    Mat getProlongationOp(int ln) const;

    /*!
     * \brief Get the scaling for level restriction operator. The restriction
     * operator restricts data from level \em ln + 1 to level \em ln.
     * Restriction op is defined to be the scaled adjoint of prolongation
     * operator, i.e., R = L P^T.
     */
    Vec getRestrictionScalingOp(int ln) const;

    //\}

    /*!
     * \name Partial implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Compute the composite-grid residual on the specified range of
     * levels of the patch hierarchy.
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num);

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

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesIBLevelRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesIBLevelRelaxationFACOperator(const StaggeredStokesIBLevelRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesIBLevelRelaxationFACOperator& operator=(const StaggeredStokesIBLevelRelaxationFACOperator& that);

    /*
     * Level solvers and solver parameters.
     */
    std::string d_level_solver_type, d_level_solver_default_options_prefix;
    double d_level_solver_abs_residual_tol, d_level_solver_rel_residual_tol;
    int d_level_solver_max_iterations;
    std::vector<SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesPETScLevelSolver> > d_level_solvers;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_level_solver_db;

    /*
     * Velocity prolongation type.
     */
    std::string d_u_petsc_prolongation_method;

    /*
     * Application ordering of u from MAC DOFs on various patch levels.
     */
    std::vector<AO> d_u_app_ordering;

    /*
     * Eulerian data for storing u and p DOFs indexing.
     */
    std::vector<std::vector<int> > d_num_dofs_per_proc;
    int d_u_dof_index_idx, d_p_dof_index_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int> > d_u_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int> > d_p_dof_index_var;

    /*
     * The time stepping type.
     */
    TimeSteppingType d_time_stepping_type;

    /*
     * Jacobian of the elasticity force at the finest patch level.
     */
    Mat d_A_mat;

    /*
     * IB interpolation operator J for the finest patch level.
     */
    Mat d_J_mat;

    /*
     * Data structures for elasticity and prolongation operator respresentation
     * on various patch levels.
     */
    double d_SAJ_fill;
    std::vector<Mat> d_SAJ_mat, d_prolongation_mat;
    std::vector<Vec> d_scale_restriction_mat;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::vector<boost::array<SAMRAI::hier::BoxList<NDIM>, NDIM> > > d_patch_side_bc_box_overlap;
    std::vector<std::vector<SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StaggeredStokesIBLevelRelaxationFACOperator
