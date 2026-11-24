// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_SCPoissonPointRelaxationFACOperator
#define included_IBTK_SCPoissonPointRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/PoissonFACPreconditioner.h>
#include <ibtk/PoissonFACPreconditionerStrategy.h>
#include <ibtk/PoissonSolver.h>
#include <ibtk/StaggeredPhysicalBoundaryHelper.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <map>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
template <int DIM>
class BoxList;
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

 smoother_type = "PATCH_GAUSS_SEIDEL"         // or "PROCESSOR_GAUSS_SEIDEL", "RED_BLACK_GAUSS_SEIDEL"
 prolongation_method = "CONSTANT_REFINE"      // see setProlongationMethod()
 restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethod()
 coarse_solver_type = "HYPRE_LEVEL_SOLVER"    // smoother or solver type for the coarsest level
 coarse_solver_rel_residual_tol = 1.0e-5      // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50     // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 1             // see setCoarseSolverMaxIterations()
 coarse_solver_db = { ... }                   // SAMRAI::tbox::Database for initializing coarse
 level solver
 \endverbatim
 *
 * With \c "PROCESSOR_GAUSS_SEIDEL", a sweep uses values on neighboring patches of the same process that were updated
 * earlier in that sweep. The results therefore depend on the patch layout.
*/
class SCPoissonPointRelaxationFACOperator : public PoissonFACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    SCPoissonPointRelaxationFACOperator(const std::string& object_name,
                                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                        std::string default_options_prefix);

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
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients of the operator. C must be zero or constant, and D must
     * be constant.
     */
    void setPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& poisson_spec) override;

    /*!
     * \brief Return the number of coarse-level solver reuses.
     */
    long getCoarseSolverReuseCount() const;

    /*!
     * \brief Return the number of coarse-level solver initializations when reuse caching is enabled.
     */
    long getCoarseSolverReinitCount() const;

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
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                            SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                            int coarsest_ln) override;

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
                         SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num) override;

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                                            int coarsest_reset_ln,
                                            int finest_reset_ln) override;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln) override;

    /*
     * Coarse level solvers and solver parameters.
     */
    SAMRAI::tbox::Pointer<PoissonSolver> d_coarse_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_coarse_solver_db;

    /*
     * Boxes of the ghost sides of each component along the coarse-fine interface of each local patch, empty for a
     * patch without one. smoothError() saves and restores the values on these sides between sweeps.
     */
    std::vector<std::vector<std::array<SAMRAI::hier::BoxList<NDIM>, NDIM>>> d_patch_cf_bdry_ghost_boxes;

    /*
     * Patch overlap data.
     */
    std::vector<std::vector<std::array<std::map<int, SAMRAI::hier::Box<NDIM>>, NDIM>>> d_patch_neighbor_overlap;

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
    SCPoissonPointRelaxationFACOperator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SCPoissonPointRelaxationFACOperator(const SCPoissonPointRelaxationFACOperator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SCPoissonPointRelaxationFACOperator& operator=(const SCPoissonPointRelaxationFACOperator& that) = delete;

    // Coarsest-level solver caching state for opt-in reuse across resets.
    enum class CachedCoarseSolverKind
    {
        UNSUPPORTED = 0,
        VCSC_PETSC = 1
    };

    void clearCachedCoarseState(bool deallocate_coarse_solver);

    bool d_reuse_unchanged_coarse_solver = false;
    double d_coarse_solver_reuse_rel_tol = 1.0e-6;
    bool d_has_cached_coarse_state = false;
    bool d_in_destructor = false;
    CachedCoarseSolverKind d_cached_l0_solver_kind = CachedCoarseSolverKind::UNSUPPORTED;
    IBTK::VCInterpType d_cached_l0_mu_interp_type = IBTK::VC_HARMONIC_INTERP;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_cached_l0_level;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_cached_l0_solution;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_cached_l0_rhs;
    std::vector<SAMRAI::hier::Box<NDIM>> d_cached_l0_boxes;
    std::vector<int> d_cached_l0_procs;
    double d_cached_l0_dt = 0.0;
    bool d_cached_l0_C_is_var = false;
    bool d_cached_l0_C_is_zero = true;
    double d_cached_l0_C_const = 0.0;
    bool d_cached_l0_D_is_var = false;
    double d_cached_l0_D_const = 0.0;
    int d_cached_l0_x_depth = -1;
    SAMRAI::hier::IntVector<NDIM> d_cached_l0_x_gcw = SAMRAI::hier::IntVector<NDIM>(-1);
    std::vector<double> d_cached_l0_C_vals;
    std::vector<double> d_cached_l0_D_vals;
    std::vector<double> d_cached_l0_bc_vals;
    long d_coarse_solver_reuse_count = 0;
    long d_coarse_solver_reinit_count = 0;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_SCPoissonPointRelaxationFACOperator
