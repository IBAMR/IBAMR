#ifndef included_INSProjectionLaplaceOperator
#define included_INSProjectionLaplaceOperator

// Filename: INSProjectionLaplaceOperator.h
// Last modified: <28.Aug.2007 19:41:09 griffith@box221.cims.nyu.edu>
// Created on 19 Sep 2003 by Boyce Griffith (griffith@mstu1.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/INSProjectionBcCoef.h>

// STOOLS INCLUDES
#include <stools/CartRobinPhysBdryOp.h>
#include <stools/HierarchyMathOps.h>
#include <stools/LinearOperator.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>
#include <PoissonSpecifications.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SAMRAIVectorReal.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSProjectionLaplaceOperator is a reimplementation of class
 * STOOLS::CCLaplaceOperator that has been simplified and specialized for
 * solving the discrete Poisson problems that arrise in projection methods for
 * the incompressible Navier-Stokes equations.
 */
class INSProjectionLaplaceOperator
    : public STOOLS::LinearOperator
{
public:
    /*!
     * \brief Constructor for class INSProjectionLaplaceOperator initializes the
     * operator coefficients and boundary conditions.
     *
     * \param object_name     String used to register internal variables and for error reporting purposes.
     * \param bc_coef         Robin boundary conditions to use with this class.
     * \param homogeneous_bc  Whether to employ the homogeneous form of the boundary conditions.
     */
    INSProjectionLaplaceOperator(
        const std::string& object_name,
        INSProjectionBcCoef* const bc_coef,
        const bool homogeneous_bc=true);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSProjectionLaplaceOperator();

    /*!
     * \brief Set the INSProjectionBcCoef object used to specify physical
     * boundary conditions.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition coefficients
     */
    virtual void
    setPhysicalBcCoef(
        INSProjectionBcCoef* const bc_coef);

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    virtual void
    setHomogeneousBc(
        const bool homogeneous_bc);

    /*!
     * \brief Set the hierarchy time, for use with the refinement schedules and
     * boundary condition routines employed by the object.
     */
    void
    setTime(
        const double time);

    /*!
     * \brief Set the HierarchyMathOps object used by the operator.
     */
    void
    setHierarchyMathOps(
        SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> hier_math_ops);

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Modify y to account for inhomogeneous boundary conditions.
     *
     * Before calling this function, the form of the vector y should be set
     * properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in this vector should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * \note The operator MUST be initialized prior to calling
     * modifyRhsForInhomogeneousBc.
     *
     * \see initializeOperatorState
     *
     * \param y output: y=Ax
     *
     * \note A default implementation is provided which does nothing but warns
     * that inhomogeneous boundary conditions are not properly supported.
     */
    virtual void
    modifyRhsForInhomogeneousBc(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    virtual void
    apply(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute z=Ax+y.
     *
     * Before calling this function, the form of the vectors x, y, and z should
     * be set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * IMPORTANT \note In general, the vectors x and z \em cannot be the same.
     * Vectors y and z may be the same.
     *
     * Upon return from this function, the z vector will contain the result of
     * the application of A to x and adding y.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y input
     * \param z output: z=Ax+y
     */
    virtual void
    applyAdd(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& z);

    /*!
     * \brief Compute hierarchy-dependent data required for computing y=Ax (and
     * y=A'x).
     *
     * \param in input vector
     * \param out output vector
     *
     * \see KrylovLinearSolver::initializeSolverState
     */
    virtual void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out);

    /*!
     * \brief Remove all hierarchy-dependent data computed by
     * initializeOperatorState().
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() even state is already
     * deallocated.
     *
     * \see initializeOperatorState
     * \see KrylovLinearSolver::deallocateSolverState
     */
    virtual void
    deallocateOperatorState();

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     *
     * \param enabled logging state: true=on, false=off
     */
    virtual void
    enableLogging(
        bool enabled=true);

    /*!
     * \brief Print out internal class data for debugging.
     */
    virtual void
    printClassData(
        std::ostream& os) const;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSProjectionLaplaceOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSProjectionLaplaceOperator(
        const INSProjectionLaplaceOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSProjectionLaplaceOperator&
    operator=(
        const INSProjectionLaplaceOperator& that);

    // Housekeeping.
    std::string d_object_name;

    // Operator parameters.
    bool d_is_initialized;
    double d_apply_time;

    // Cached communications algorithms and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_homogeneous_bc_fill_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_homogeneous_bc_fill_scheds;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_inhomogeneous_bc_fill_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_inhomogeneous_bc_fill_scheds;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_no_fill;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_bc_refine_strategy;
    SAMRAI::tbox::Pointer<STOOLS::CartRobinPhysBdryOp> d_bc_op;

    // Scratch vectors.
    int d_scratch_idx;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_x, d_b;
    bool d_correcting_rhs;

    // Problem specification and mathematical operators.
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;
    INSProjectionBcCoef* d_bc_coef;
    bool d_homogeneous_bc;
    SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> d_hier_math_ops;
    bool d_hier_math_ops_external;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln, d_finest_ln;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <stools/INSProjectionLaplaceOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSProjectionLaplaceOperator
