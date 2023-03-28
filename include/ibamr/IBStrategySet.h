// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_IBStrategySet
#define included_IBAMR_IBStrategySet

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBStrategy.h"

#include "IntVector.h"
#include "tbox/Pointer.h"

#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace IBAMR
{
class IBHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
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
 * \brief Class IBStrategySet is a utility class that allows multiple IBStrategy
 * objects to be employed by a single IBHierarchyIntegrator.
 */
class IBStrategySet : public IBStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    template <typename InputIterator>
    IBStrategySet(InputIterator first, InputIterator last) : d_strategy_set(first, last)
    {
        // intentionally blank
        return;
    } // IBStrategySet

    /*!
     * Register the IBHierarchyIntegrator object that is using this strategy
     * class.
     */
    void registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver) override;

    /*!
     * Register Eulerian variables with the parent IBHierarchyIntegrator with
     * the VariableDatabase, or via the various versions of the protected method
     * IBStrategy::registerVariable().
     */
    void registerEulerianVariables() override;

    /*!
     * Register Eulerian refinement or coarsening algorithms with the parent
     * IBHierarchyIntegrator using the two versions of the protected methods
     * IBStrategy::registerGhostfillRefineAlgorithm(),
     * IBStrategy::registerProlongRefineAlgorithm(), and
     * IBStrategy::registerCoarsenAlgorithm().
     */
    void registerEulerianCommunicationAlgorithms() override;

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
     * Same as the base class: considers all owned IBStrategy objects.
     */
    virtual double getMaxPointDisplacement() const override;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

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
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
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
     * Execute user-defined routines just before solving the fluid equations.
     */
    void preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num) override;

    /*!
     * Execute user-defined routines just after solving the fluid equations.
     */
    void postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num) override;

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

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBStrategySet(const IBStrategySet& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBStrategySet& operator=(const IBStrategySet& that) = delete;

    /*!
     * \brief The set of IBStrategy objects.
     */
    std::vector<SAMRAI::tbox::Pointer<IBStrategy> > d_strategy_set;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBStrategySet
