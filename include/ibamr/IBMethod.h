// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_IBMethod
#define included_IBAMR_IBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBImplicitStrategy.h"
#include "ibamr/IBInstrumentPanel.h"
#include "ibamr/IBLagrangianForceStrategy.h"
#include "ibamr/IBLagrangianSourceStrategy.h"
#include "ibamr/IBMethodPostProcessStrategy.h"

#include "ibtk/LInitStrategy.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ibtk_utilities.h"

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "tbox/Pointer.h"

#include "petscmat.h"
#include "petscsys.h"
#include "petscvec.h"

#include <limits>
#include <set>
#include <string>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace IBTK
{
class LData;
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace tbox
{
class Database;
template <class TYPE>
class Array;
} // namespace tbox
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
 * \brief Class IBMethod is an implementation of the abstract base class
 * IBImplicitStrategy that provides functionality required by the standard IB
 * method.
 */
class IBMethod : public IBImplicitStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBMethod(std::string object_name,
             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
             bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~IBMethod();

    /*!
     * Supply a Lagrangian force object.
     */
    void registerIBLagrangianForceFunction(SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> ib_force_fcn);

    /*!
     * Supply a Lagrangian source object.
     */
    void registerIBLagrangianSourceFunction(SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> ib_source_fcn);

    /*!
     * Supply a Lagrangian initialization object.
     */
    void registerLInitStrategy(SAMRAI::tbox::Pointer<IBTK::LInitStrategy> l_initializer);

    /*!
     * Free references to Lagrangian initialization objects.
     */
    void freeLInitStrategy();

    /*!
     * Supply a post processor object.
     */
    void registerIBMethodPostProcessor(SAMRAI::tbox::Pointer<IBMethodPostProcessStrategy> post_processor);

    /*!
     * Return a pointer to the Lagrangian data manager object.
     */
    IBTK::LDataManager* getLDataManager() const;

    /*!
     * Return a pointer to the instrumentation manager object.
     */
    SAMRAI::tbox::Pointer<IBInstrumentPanel> getIBInstrumentPanel() const;

    /*!
     * Register a Lagrangian Silo data writer so this class will write plot
     * files that may be postprocessed with the VisIt visualization tool.
     */
    void registerLSiloDataWriter(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer);

    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>& getMinimumGhostCellWidth() const override;

    /*!
     * Setup the tag buffer.
     */
    void setupTagBuffer(SAMRAI::tbox::Array<int>& tag_buffer,
                        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) const override;

    /*!
     * Inactivate a structure/part. See IBAMR::IBStrategy::inactivateLagrangianStructure().
     */
    virtual void inactivateLagrangianStructure(int structure_number = 0,
                                               int level_number = std::numeric_limits<int>::max()) override;

    /*!
     * Activate a previously inactivated structure/part to be used again in
     * FSI calculations. See IBAMR::IBStrategy::activateLagrangianStructure().
     */
    virtual void activateLagrangianStructure(int structure_number = 0,
                                             int level_number = std::numeric_limits<int>::max()) override;

    /*!
     * Determine whether or not the given structure or part is currently
     * activated. See IBAMR::IBStrategy::getLagrangianStructureIsActivated().
     */
    virtual bool getLagrangianStructureIsActivated(int structure_number = 0,
                                                   int level_number = std::numeric_limits<int>::max()) const override;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Create solution and rhs data on the specified level of the patch
     * hierarchy.
     */
    void createSolverVecs(Vec* X_vec, Vec* F_vec) override;

    /*!
     * Setup solution and rhs data on the specified level of the patch
     * hierarchy.
     */
    void setupSolverVecs(Vec* X_vec, Vec* F_vec) override;

    /*!
     * Set the value of the updated position vector.
     */
    void setUpdatedPosition(Vec& X_new_vec) override;

    /*!
     * Get the value of the updated position vector.
     */
    void getUpdatedPosition(Vec& X_new_vec) override;

    /*!
     * Compute the nonlinear residual for backward Euler time stepping.
     */
    void computeResidualBackwardEuler(Vec& R_vec) override;

    /*!
     * Compute the nonlinear residual for midpoint rule time stepping.
     */
    void computeResidualMidpointRule(Vec& R_vec) override;

    /*!
     * Compute the nonlinear residual for trapezoidal rule time stepping.
     */
    void computeResidualTrapezoidalRule(Vec& R_vec) override;

    /*!
     * Update the positions used for the "fixed" interpolation and spreading
     * operators.
     */
    void updateFixedLEOperators() override;

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the backward Euler
     * method.
     */
    void backwardEulerStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the midpoint rule.
     */
    void midpointStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the trapezoidal
     * rule.
     */
    void trapezoidalStep(double current_time, double new_time) override;

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time) override;

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time) override;

    /*!
     * Construct the IB interpolation operator.
     */
    void constructInterpOp(Mat& J,
                           void (*spread_fnc)(const double, double*),
                           int stencil_width,
                           const std::vector<int>& num_dofs_per_proc,
                           int dof_index_idx,
                           double data_time);

    /*!
     * Indicate whether there are any internal fluid sources/sinks.
     */
    bool hasFluidSources() const override;

    /*!
     * Compute the Lagrangian source/sink density at the specified time within
     * the current time interval.
     */
    void computeLagrangianFluidSource(double data_time) override;

    /*!
     * Spread the Lagrangian source/sink density to the Cartesian grid at the
     * specified time within the current time interval.
     */
    void spreadFluidSource(
        int q_data_idx,
        IBTK::RobinPhysBdryPatchStrategy* q_phys_bdry_op,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& q_prolongation_scheds,
        double data_time) override;

    /*!
     * Compute the pressures at the positions of any distributed internal fluid
     * sources or sinks.
     */
    void interpolatePressure(
        int p_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& p_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& p_ghost_fill_scheds,
        double data_time) override;

    /*!
     * Execute user-defined post-processing operations.
     */
    void postprocessData() override;

    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time) override;

    /*!
     * Register a load balancer and work load patch data index with the IB
     * strategy object.
     *
     * @deprecated This method is no longer necessary with the current
     * workload estimation scheme.
     */
    void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
                              int workload_data_idx) override;

    /*!
     * Add the estimated computational work from the current object per cell
     * into the specified <code>workload_data_idx</code>.
     */
    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const int workload_data_idx) override;

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    void beginDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    void endDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                             bool allocate_data) override;

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_level,
                                     int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too) override;

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * Get the current structure position data.
     */
    void getPositionData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >** X_data,
                         bool** X_needs_ghost_fill,
                         double data_time);

    /*!
     * Get the current interpolation/spreading position data.
     */
    void getLECouplingPositionData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >** X_LE_data,
                                   bool** X_LE_needs_ghost_fill,
                                   double data_time);

    /*!
     * Get the current structure velocity data.
     */
    void getVelocityData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >** U_data, double data_time);

    /*!
     * Get the current structure force data.
     */
    void getForceData(std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >** F_data,
                      bool** F_needs_ghost_fill,
                      double data_time);

    /*!
     * Interpolate the current and new data to obtain values at the midpoint of
     * the time interval.
     */
    void reinitMidpointData(const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& current_data,
                            const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& new_data,
                            const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& half_data);

    /*!
     * Set the elements of the Lagrangian vector to zero at anchored nodes of
     * the curvilinear mesh.
     */
    void
    resetAnchorPointValues(std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > U_data, int coarsest_ln, int finest_ln);

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log = false;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * The current time step interval.
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_half_time = std::numeric_limits<double>::quiet_NaN();

    /*
     * Boolean values tracking whether certain quantities need to be
     * reinitialized.
     */
    bool d_X_current_needs_ghost_fill = true, d_X_new_needs_ghost_fill = true, d_X_half_needs_ghost_fill = true,
         d_X_LE_new_needs_ghost_fill = true, d_X_LE_half_needs_ghost_fill = true;
    bool d_F_current_needs_ghost_fill = true, d_F_new_needs_ghost_fill = true, d_F_half_needs_ghost_fill = true;

    /*
     * The LDataManager is used to coordinate the distribution of Lagrangian
     * data on the patch hierarchy.
     */
    IBTK::LDataManager* d_l_data_manager;
    std::string d_interp_kernel_fcn = "IB_4", d_spread_kernel_fcn = "IB_4";
    bool d_error_if_points_leave_domain = false;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;

    /*
     * Lagrangian variables.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_X_current_data, d_X_new_data, d_X_half_data, d_X_LE_new_data,
        d_X_LE_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_U_current_data, d_U_new_data, d_U_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_F_current_data, d_F_new_data, d_F_half_data;

    /*
     * List of local indices of local anchor points.
     *
     * NOTE: IB points are automatically considered to be anchored if they are
     * within 2.0*sqrt(epsilon_mach) of the physical boundary.
     */
    std::vector<std::set<int> > d_anchor_point_local_idxs;

    /*
     * Instrumentation (flow meter and pressure gauge) algorithms and data
     * structures.
     */
    SAMRAI::tbox::Pointer<IBInstrumentPanel> d_instrument_panel;
    std::vector<double> d_total_flow_volume;

    /*
     * The specification and initialization information for the Lagrangian data
     * used by the integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::LInitStrategy> d_l_initializer;

    /*
     * The force generators.
     */
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> d_ib_force_fcn;
    bool d_ib_force_fcn_needs_init = true;

    /*
     * The source/sink generators.
     */
    SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> d_ib_source_fcn;
    bool d_ib_source_fcn_needs_init = true;
    std::vector<std::vector<IBTK::Point> > d_X_src;
    std::vector<std::vector<double> > d_r_src, d_P_src, d_Q_src;
    std::vector<int> d_n_src;
    bool d_normalize_source_strength = false;

    /*
     * Post-processor object.
     */
    SAMRAI::tbox::Pointer<IBMethodPostProcessStrategy> d_post_processor;

    /*
     * Visualization data writers.
     */
    SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> d_silo_writer;

    /*
     * Nonuniform load balancing data structures.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;
    int d_workload_idx = IBTK::invalid_index;

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBMethod() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMethod(const IBMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMethod& operator=(const IBMethod& that) = delete;

    /*!
     * Reset the Lagrangian force function object.
     */
    void resetLagrangianForceFunction(double init_data_time, bool initial_time);

    /*!
     * Reset the Lagrangian source function object.
     */
    void resetLagrangianSourceFunction(double init_data_time, bool initial_time);

    /*!
     * Compute the flow rates and pressures in the internal flow meters and
     * pressure gauges.
     */
    void updateIBInstrumentationData(int timestep_num, double data_time);

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

#endif // #ifndef included_IBAMR_IBMethod
