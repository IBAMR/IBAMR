// Filename: CCHierarchyProjector.h
// Created on 18 Feb 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_CCHierarchyProjector
#define included_CCHierarchyProjector

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscksp.h>

// IBTK INCLUDES
#include <ibtk/CCDivGradHypreLevelSolver.h>
#include <ibtk/CCDivGradOperator.h>
#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/HierarchyMathOps.h>

// SAMRAI INCLUDES
#include <StandardTagAndInitStrategy.h>
#include <tbox/Serializable.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class CCHierarchyProjector provides so-called "exact" projection
 * functionality for cell-centered vector fields.
 *
 * \note The class implementation is intended to be used only for uniform grids
 * with periodic boundaries.  This class may also work in some cases for
 * locally-refined discretizations, but if it does, its success is purely
 * accidental.
 */
class CCHierarchyProjector
    : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public virtual SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for CCHierarchyProjector sets some default values, reads
     * in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     *
     * When assertion checking is active, passing in any null pointer or an
     * empty string will result in an unrecoverable exception.
     */
    CCHierarchyProjector(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        bool register_for_restart=true);

    /*!
     * The destructor for CCHierarchyProjector unregisters the integrator object
     * with the restart manager when so registered.
     */
    virtual
    ~CCHierarchyProjector();

    /*!
     * Return the name of the hierarchy projector object.
     */
    const std::string&
    getName() const;

    ///
    ///  The following routines:
    ///
    ///      getHierarchyMathOps(),
    ///      setHierarchyMathOps(),
    ///      isManagingHierarchyMathOps()
    ///
    ///  allow for the sharing of a single HierarchyMathOps object between
    ///  multiple HierarchyIntegrator objects.
    ///

    /*!
     * Return a pointer to the HierarchyMathOps object being used by this
     * integrator.
     *
     * The HierarchyMathOps object supplies discrete differential operations on
     * the patch hierarchy as well as cell weights used in computing discrete
     * norms of quantities defined on the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps>
    getHierarchyMathOps() const;

    /*!
     * Set the HierarchyMathOps object being used by this projector.
     *
     * When manage_ops is true, the HierarchyMathOps object is managed by the
     * integrator.  In particular, the integrator is responsible for invoking
     * HierarchyMathOps::setPatchHierarchy() and HierarchyMathOps::resetLevels()
     * following any changes to the configuration of the patch hierarchy.
     */
    void
    setHierarchyMathOps(
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        const bool manage_ops=false);

    /*!
     * Returns whether this projector is managing the state of its
     * HierarchyMathOps object.
     *
     * When the projector is managing the state of its HierarchyMathOps object,
     * the integrator is responsible for invoking
     * HierarchyMathOps::setPatchHierarchy() and HierarchyMathOps::resetLevels()
     * following any changes to the configuration of the patch hierarchy.
     */
    bool
    isManagingHierarchyMathOps() const;

    ///
    ///  The following routines:
    ///
    ///      getPoissonSolver()
    ///
    ///  allow other objects to access the Poisson solver and related data used
    ///  by this integrator.
    ///

    /*!
     * Returns a pointer to the concrete linear solver object employed by the
     * projector to solve the elliptic projection equation.
     */
    SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver>
    getPoissonSolver() const;

    ///
    ///  The following routines:
    ///
    ///      projectHierarchy()
    ///
    ///  provide the projection functionality.
    ///

    /*!
     * Project the cell centered co-located velocity W on the hierarchy.
     *
     * Computes U = W - Grad_Phi, where div U = Q.  If Q is not supplied, it is
     * assumed that div U = 0.
     */
    virtual void
    projectHierarchy(
        const double rho,
        const double dt,
        const double time,
        const int U_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& U_var,
        const int Phi_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Phi_var,
        const int Grad_Phi_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Grad_Phi_var,
        const int W_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& W_var,
        const int Q_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Q_var=NULL);

    ///
    ///  The following routines:
    ///
    ///      initializeLevelData(),
    ///      resetHierarchyConfiguration()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::mesh::StandardTagAndInitStrategy abstract base class.
    ///

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.  The level number indicates that of
     * the new level.  The old_level pointer corresponds to the level that
     * resided in the hierarchy before the level with the specified number was
     * introduced.  If the pointer is null, there was no level in the hierarchy
     * prior to the call and the level data is set based on the user routines
     * and the simulation time.  Otherwise, the specified level replaces the old
     * level and the new level receives data from the old level appropriately
     * before it is destroyed.
     *
     * The boolean argument initial_time indicates whether the level is being
     * introduced for the first time (i.e., at initialization time), or after
     * some regrid process during the calculation beyond the initial hierarchy
     * construction.  This information is provided since the initialization of
     * the data on a patch may be different in each of those circumstances.  The
     * can_be_refined boolean argument indicates whether the level is the finest
     * level allowed in the hierarchy.  This may or may not affect the data
     * initialization process depending on the problem.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, the level number does not match any
     * level in the hierarchy, or the old level number does not match the level
     * number (if the old level pointer is non-null).
     */
    virtual void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level=SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

    /*!
     * Reset cached communication schedules after the hierarchy has changed (for
     * example, due to regridding) and the data has been initialized on the new
     * levels.  The intent is that the cost of data movement on the hierarchy
     * will be amortized across multiple communication cycles, if possible.  The
     * level numbers indicate the range of levels in the hierarchy that have
     * changed.  However, this routine updates communication schedules every
     * level finer than and including that indexed by the coarsest level number
     * given.  When the integrator is managing the state of its HierarchyMathOps
     * object, the integrator also invokes HierarchyMathOps::setPatchHierarchy()
     * and HierarchyMathOps::resetLevels().
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, any pointer to a level in the hierarchy
     * that is coarser than the finest level is null, or the given level numbers
     * not specified properly; e.g., coarsest_level > finest_level.
     */
    virtual void
    resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    virtual void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CCHierarchyProjector();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCHierarchyProjector(
        const CCHierarchyProjector& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCHierarchyProjector&
    operator=(
        const CCHierarchyProjector& that);

    /*!
     * Read input values, indicated above, from given database.  The boolean
     * argument is_from_restart should be set to true if the simulation is
     * beginning from restart.  Otherwise it should be set to false.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *    -   The class version number and restart version number do not match.
     */
    void
    getFromRestart();

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * Indicates whether the projector should output logging messages.
     */
    bool d_do_log;

    /*
     * Pointers to the patch hierarchy object associated with this projector
     * object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    bool d_is_managing_hier_math_ops;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_var;
    int d_wgt_idx;
    double d_volume;

    /*
     * Cached communications algorithms and schedules.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_Phi_hier_bdry_fill_op, d_W_hier_bdry_fill_op, d_no_fill_op;

    /*
     * The Poisson solver and associated data including Poisson specifications,
     * boundary conditions, and the solver configuration database.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_F_var;
    int d_F_idx;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_null_space_var;
    std::vector<int> d_null_space_idxs;
    std::vector<Vec> d_petsc_null_space_vecs;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_W_var;
    int d_W_idx;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_sol_vec, d_rhs_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_null_space_vecs;

    int d_max_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    bool d_initial_guess_nonzero;

    SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> d_poisson_solver;
    SAMRAI::tbox::Pointer<IBTK::CCDivGradOperator> d_cc_div_grad_op;
    SAMRAI::tbox::Pointer<IBTK::CCDivGradHypreLevelSolver> d_cc_div_grad_hypre_solver;

    /*
     * Note that data is NOT allocated for these variables.  They are only used
     * to reset the Poisson solver.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_sol_var, d_rhs_var;
    int d_sol_idx, d_rhs_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/CCHierarchyProjector.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CCHierarchyProjector
