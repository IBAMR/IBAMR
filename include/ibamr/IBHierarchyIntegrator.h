// ---------------------------------------------------------------------
//
// Copyright (c) 2006 - 2022 by the IBAMR developers
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

#ifndef included_IBAMR_IBHierarchyIntegrator
#define included_IBAMR_IBHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
template <int DIM>
class LoadBalancer;
} // namespace mesh
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHierarchyIntegrator provides an abstract interface for a time
 * integrator for various versions of the immersed boundary method on an AMR
 * grid hierarchy, along with basic data management for variables defined on
 * that hierarchy.
 *
 * <h2>Options Controlling Regridding</h2>
 *
 * Most IBAMR applications involve structure meshes defined on the finest level
 * of the patch hierarchy. These structures move: in particular, they can
 * potentially move off the finest grid level, causing the interaction routines
 * to no longer work.
 *
 * This class offers two different strategies for calculating how much the
 * structure has moved (which it then uses to specify that a regrid is
 * required):
 *
 * <ol>
 *   <li><em>Estimation based on the fluid</em>: The structure is assumed to
 *     move at the maximum velocity of the fluid at each time step. This class
 *     will regrid once a single fluid point has moved at least
 *     <code>regrid_fluid_cfl_interval</code> cell widths since the last
 *     regrid.</li>
 *
 *   <li><em>Estimation based on the structure</em>: The structure's
 *     displacement is calculated directly from its position vector and that
 *     value is used to determine the number of cell widths the structure has
 *     moved. In this case we will regrid once a point on the structure has
 *     moved at least <code>regrid_structure_cfl_interval</code> cell widths
 *     since the last regrid.</li>
 * </ol>
 *
 * Both <code>regrid_fluid_cfl_interval</code> and
 * <code>regrid_structure_cfl_interval</code> can be specified in the input
 * database. For backwards compatibility the value
 * <code>regrid_cfl_interval</code> is equivalent to
 * <code>regrid_fluid_cfl_interval</code>. <em>At the present time
 * <code>regrid_structure_cfl_interval</code> is not implemented for all
 * IBStrategy classes.</em>
 *
 * Alternatively, one can request that the solver (regardless of any computed
 * displacement or velocity) regrid every time a fixed number of timesteps have
 * occurred by specifying <code>regrid_interval</code> in the input database.
 * <em>If either <code>regrid_structure_cfl_interval</code> or
 * <code>regrid_fluid_cfl_interval</code> are provided in the input database
 * then <code>regrid_interval</code> is ignored.</em>
 */
class IBHierarchyIntegrator : public IBTK::HierarchyIntegrator
{
public:
    friend class IBStrategy;

    /*!
     * The destructor for class IBHierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~IBHierarchyIntegrator() = default;

    /*!
     * Return the time stepping scheme.
     */
    TimeSteppingType getTimeSteppingType() const;

    /*!
     * Return a pointer to the IBStrategy object registered with this
     * integrator.
     */
    SAMRAI::tbox::Pointer<IBStrategy> getIBStrategy() const;

    /*!
     * Supply a body force (optional).
     */
    void registerBodyForceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Register a load balancer for non-uniform load balancing.
     */
    virtual void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer) override;

    /*!
     * Return a pointer to the fluid velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getVelocityVariable() const;

    /*!
     * Return a pointer to the fluid pressure state variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getPressureVariable() const;

    /*!
     * Return a pointer to the body force variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getBodyForceVariable() const;

    /*!
     * Return a pointer to the source strength variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getFluidSourceVariable() const;

    /*!
     * Return a pointer to the velocity physical boundary conditions
     */
    IBTK::RobinPhysBdryPatchStrategy* getVelocityPhysBdryOp() const;

    /*!
     * Basic functions to prepare to advance data from current_time to new_time.
     *
     * A default implementation is provided that sets the current values of
     * num_cycles and the time step size and checks to see if the time step size
     * has changed.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Initialize the AMR patch hierarchy and data defined on the hierarchy at
     * the start of a computation.  If the computation is begun from a restart
     * file, the patch hierarchy and patch data are read from the hierarchy
     * database.  Otherwise, the patch hierarchy and patch data are initialized
     * by the gridding algorithm associated with the integrator object.
     *
     * The implementation of this function assumes that the hierarchy exists
     * upon entry to the function, but that it contains no patch levels.  On
     * return from this function, the state of the integrator object will be
     * such that it is possible to step through time via the advanceHierarchy()
     * function.
     */
    void initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

protected:
    /*!
     * Perform necessary data movement, workload estimation, and logging prior
     * to regridding.
     */
    void regridHierarchyBeginSpecialized() override;

    /*!
     * Perform necessary data movement and logging after regridding.
     */
    void regridHierarchyEndSpecialized() override;

    /*!
     * The constructor for class IBHierarchyIntegrator sets some default values,
     * reads in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     */
    IBHierarchyIntegrator(const std::string& object_name,
                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                          SAMRAI::tbox::Pointer<IBStrategy> ib_method_ops,
                          SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
                          bool register_for_restart = true);

    /*!
     * Function to determine whether regridding should occur at the current time
     * step.
     */
    bool atRegridPointSpecialized() const override;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                                        bool allocate_data) override;

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Add the work contributions (excluding the background grid) for the
     * current hierarchy into the variable with index
     * <code>workload_data_idx</code>. The only direct workload contribution
     * of this hierarchy manager is usually the work done by the IBStrategy
     * object.
     */
    virtual void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     const int workload_data_idx) override;

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized = false;

    /*!
     * Enum indicating the time integration employed for the IB equations.
     */
    TimeSteppingType d_time_stepping_type = MIDPOINT_RULE;

    /*!
     * Flag indicating whether to use an explicit predictor for the structure
     * configuration in the time stepping scheme.
     */
    bool d_use_structure_predictor;

    /*!
     * Flags to determine whether warnings or error messages should be emitted
     * when time step size changes are encountered.
     */
    bool d_error_on_dt_change = true, d_warn_on_dt_change = false;

    /*
     * The (optional) INSHierarchyIntegrator is used to provide time integration
     * capability for the incompressible Navier-Stokes equations.
     */
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> d_ins_hier_integrator;

    /**
     * The regrid CFL interval indicates the number of meshwidths a particle may
     * move in any coordinate direction between invocations of the regridding
     * process.
     *
     * @note Currently, when the CFL-based regrid interval is specified, it is
     * always used instead of the fixed-step regrid interval.
     */
    double d_regrid_fluid_cfl_interval = -1.0;

    /*
     * The regrid CFL interval indicates the number of meshwidths a particle may
     * move in any coordinate direction between invocations of the regridding
     * process.
     *
     * @note Currently, when the CFL-based regrid interval is specified, it is
     * always used instead of the fixed-step regrid interval.
     */
    double d_regrid_structure_cfl_interval = -1.0;

    /**
     * Estimation on the maximum fraction of fluid cells the structure has
     * moved based on the maximum fluid velocity.
     */
    double d_regrid_fluid_cfl_estimate = 0.0;

    /**
     * Estimation on the maximum fraction of fluid cells the structure has
     * moved based on the infinity norm of the structure's displacement.
     */
    double d_regrid_structure_cfl_estimate = 0.0;

    /*
     * IB method implementation object.
     */
    SAMRAI::tbox::Pointer<IBStrategy> d_ib_method_ops;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > d_hier_velocity_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > d_hier_pressure_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;

    /*
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_var, d_p_var, d_f_var, d_q_var;
    int d_u_idx = IBTK::invalid_index, d_p_idx = IBTK::invalid_index, d_f_idx = IBTK::invalid_index,
        d_f_current_idx = IBTK::invalid_index, d_q_idx = IBTK::invalid_index;

    /*!
     * Context containing all patch data indices relevant to IB operations.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ib_context;

    /*!
     * ComponentSelector corresponding to d_ib_context. Also contains patch data
     * indices for relevant cloned indices (which, as they are clones, cannot be
     * placed in the Context).
     */
    SAMRAI::hier::ComponentSelector d_ib_data;

    /*
     * Refine and coarsen algorithm data.
     * The base class, HierarchyIntegrator, is responsible for the d_u_phys_bdry_op and d_p_phys_bdry_op
     * objects as they are passed into d_ghostfill_strategies.
     */
    IBTK::RobinPhysBdryPatchStrategy *d_u_phys_bdry_op, *d_p_phys_bdry_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_u_ghostfill_alg, d_f_prolong_alg, d_p_ghostfill_alg,
        d_q_prolong_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_u_ghostfill_op, d_f_prolong_op, d_p_ghostfill_op,
        d_q_prolong_op;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_u_coarsen_alg, d_p_coarsen_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_u_coarsen_op, d_p_coarsen_op;

    /*
     * Body force functions.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_body_force_fcn;

    /*!
     * \brief A class to communicate the Eulerian body force computed by class
     * IBHierarchyIntegrator to the incompressible Navier-Stokes solver.
     */
    class IBEulerianForceFunction : public IBTK::CartGridFunction
    {
    public:
        /*!
         * \brief Constructor.
         */
        IBEulerianForceFunction(const IBHierarchyIntegrator* ib_solver);

        /*!
         * \brief Destructor.
         */
        ~IBEulerianForceFunction() = default;

        /*!
         * \name Methods to set the data.
         */
        //\{

        /*!
         * \note This concrete IBTK::CartGridFunction is time-dependent.
         */
        bool isTimeDependent() const override;

        /*!
         * \brief Set the data on the patch interiors on the specified levels of
         * the patch hierarchy.
         */
        void setDataOnPatchHierarchy(const int data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     const double data_time,
                                     const bool initial_time = false,
                                     const int coarsest_ln = IBTK::invalid_level_number,
                                     const int finest_ln = IBTK::invalid_level_number) override;

        /*!
         * Set the data on the patch interior.
         */
        void setDataOnPatch(int data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            double data_time,
                            bool initial_time = false,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        IBEulerianForceFunction() = delete;

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        IBEulerianForceFunction(const IBEulerianForceFunction& from) = delete;

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        IBEulerianForceFunction& operator=(const IBEulerianForceFunction& that) = delete;

        const IBHierarchyIntegrator* const d_ib_solver;
    };

    friend class IBEulerianForceFunction;

    /*!
     * \brief A class to communicate the Eulerian fluid source-sink distribution
     * computed by class IBHierarchyIntegrator to the incompressible
     * Navier-Stokes solver.
     */
    class IBEulerianSourceFunction : public IBTK::CartGridFunction
    {
    public:
        /*!
         * \brief Constructor.
         */
        IBEulerianSourceFunction(const IBHierarchyIntegrator* ib_solver);

        /*!
         * \brief Destructor.
         */
        ~IBEulerianSourceFunction() = default;

        /*!
         * \name Methods to set the data.
         */
        //\{

        /*!
         * \note This concrete IBTK::CartGridFunction is time-dependent.
         */
        bool isTimeDependent() const override;

        /*!
         * Set the data on the patch interior.
         */
        void setDataOnPatch(int data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            double data_time,
                            bool initial_time = false,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        IBEulerianSourceFunction() = delete;

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        IBEulerianSourceFunction(const IBEulerianSourceFunction& from) = delete;

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        IBEulerianSourceFunction& operator=(const IBEulerianSourceFunction& that) = delete;

        const IBHierarchyIntegrator* const d_ib_solver;
    };

    friend class IBEulerianSourceFunction;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHierarchyIntegrator(const IBHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHierarchyIntegrator& operator=(const IBHierarchyIntegrator& that) = delete;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_IBHierarchyIntegrator
