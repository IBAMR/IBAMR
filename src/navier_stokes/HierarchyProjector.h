#ifndef included_HierarchyProjector
#define included_HierarchyProjector

// Filename: HierarchyProjector.h
// Last modified: <15.Aug.2007 15:25:35 griffith@box221.cims.nyu.edu>
// Created on 30 Mar 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// STOOLS INCLUDES
#include <stools/LinearOperator.h>
#include <stools/CCLaplaceOperator.h>
#include <stools/CCPoissonFACOperator.h>
#include <stools/KrylovLinearSolver.h>
#include <stools/HierarchyMathOps.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <FACPreconditioner.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchyFaceDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <RobinBcCoefStrategy.h>
#include <SAMRAIVectorReal.h>
#include <StandardTagAndInitStrategy.h>
#include <VariableContext.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <map>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class HierarchyProjector provides MAC projection functionality for
 * face- and side-centered vector fields on AMR grids.
 */
class HierarchyProjector
    : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public virtual SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for HierarchyProjector sets some default values, reads in
     * configuration information from input and restart databases, and registers
     * the integrator object with the restart manager when requested.
     *
     * When assertion checking is active, passing in any null pointer or an
     * empty string will result in an unrecoverable exception.
     */
    HierarchyProjector(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        bool register_for_restart=true);

    /*!
     * The destructor for HierarchyProjector unregisters the integrator object
     * with the restart manager when so registered.
     */
    virtual
    ~HierarchyProjector();

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
    ///  mutiple HierarchyIntegrator objects.
    ///

    /*!
     * Return a pointer to the HierarchyMathOps object being used by this
     * integrator.
     *
     * The HierarchyMathOps object supplies discrete differential operations on
     * the patch hierarchy as well as cell weights used in computing discrete
     * norms of quantities defined on the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps>
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
        SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> hier_math_ops,
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
    ///      setPhysicalBcCoef(),
    ///      getPhysicalBcCoef(),
    ///      getPoissonSolver()
    ///
    ///  allow other objects to access the Poisson solver and related data used
    ///  by this integrator.
    ///

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note \a bc_coef may be NULL.  In this case, homogeneous Neumann boundary
     * conditions are employed.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition coefficients.
     */
    void
    setPhysicalBcCoef(
        const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef);

    /*!
     * Returns a pointer to the SAMRAI::solv::RobinBcCoefStrategy object
     * employed by the projector to solve the elliptic projection equation.
     */
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*
    getPhysicalBcCoef() const;

    /*!
     * Returns a pointer to the concrete linear solver object employed by the
     * projector to solve the elliptic projection equation.
     */
    SAMRAI::tbox::Pointer<STOOLS::KrylovLinearSolver>
    getPoissonSolver() const;

    ///
    ///  The following routines:
    ///
    ///      projectHierarchy()
    ///
    ///  provide the projection functionality.
    ///

    /*!
     * Project the face centered MAC velocity w on the hierarchy.
     *
     * Computes u = w - grad_Phi, where div u = Q.  If Q is not supplied, it is
     * assumed that div u = 0.
     *
     * Storage and communications schedules to fill ghost data for Phi are
     * required.
     *
     * Optional SAMRAI::xfer::RefineSchedule objects allow for the filling of
     * ghost cells in w at the specified time.  Additionally, w data along the
     * coarse-fine interface is synchronized when specified.
     */
    virtual void
    projectHierarchy(
        const double rho,
        const double dt,
        const int u_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& u_var,
        const int Phi_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Phi_var,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& Phi_bdry_fill,
        const double Phi_bdry_fill_time,
        const int grad_Phi_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& grad_Phi_var,
        const int w_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& w_var,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& w_bdry_fill,
        const double w_bdry_fill_time,
        const bool w_cf_bdry_synch=true,
        const int Q_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Q_var=NULL);

    /*!
     * Project the side centered MAC velocity w on the hierarchy.
     *
     * Computes u = w - grad_Phi, where div u = Q.  If Q is not supplied, it is
     * assumed that div u = 0.
     *
     * Storage and communications schedules to fill ghost data for Phi are
     * required.
     *
     * Optional SAMRAI::xfer::RefineSchedule objects allow for the filling of
     * ghost cells in w at the specified time.  Additionally, w data along the
     * coarse-fine interface is synchronized when specified.
     */
    virtual void
    projectHierarchy(
        const double rho,
        const double dt,
        const int u_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& u_var,
        const int Phi_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Phi_var,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& Phi_bdry_fill,
        const double Phi_bdry_fill_time,
        const int grad_Phi_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& grad_Phi_var,
        const int w_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& w_var,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& w_bdry_fill,
        const double w_bdry_fill_time,
        const bool w_cf_bdry_synch=true,
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
     * example, due to regidding) and the data has been initialized on the new
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

    ///
    ///  The following routines:
    ///
    ///      printClassData()
    ///
    ///  are provided for your viewing pleasure.
    ///

    /*!
     * Print out internal class data for debugging.
     */
    virtual void
    printClassData(
        std::ostream& os) const;

private:
    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > >               RefineAlgMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > >  RefineSchedMap;

    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > >              CoarsenAlgMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > > CoarsenSchedMap;

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyProjector();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyProjector(
        const HierarchyProjector& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyProjector&
    operator=(
        const HierarchyProjector& that);

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
     * Indicates whether the integrator should output logging messages.
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
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM,double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> d_hier_math_ops;
    bool d_is_managing_hier_math_ops;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_var;
    int d_wgt_idx;
    double d_volume;

    /*
     * The Poisson solver and associated data including Poisson specifications,
     * boundary conditions, and the solver configuation database.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_f_var;
    int d_f_idx;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_sol_vec, d_rhs_vec;

    int d_max_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;

    SAMRAI::solv::PoissonSpecifications d_poisson_spec;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_bc_coef;
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_bc_coef;

    SAMRAI::tbox::Pointer<STOOLS::KrylovLinearSolver> d_poisson_solver;
    SAMRAI::tbox::Pointer<STOOLS::CCLaplaceOperator> d_laplace_op;
    SAMRAI::tbox::Pointer<STOOLS::CCPoissonFACOperator> d_poisson_fac_op;
    SAMRAI::tbox::Pointer<SAMRAI::solv::FACPreconditioner<NDIM> > d_poisson_fac_pc;

    /*
     * Note that data is NOT allocated for these variables.  They are only used
     * to reset the Poisson solver.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_sol_var, d_rhs_var;
    int d_sol_idx, d_rhs_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/HierarchyProjector.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HierarchyProjector
