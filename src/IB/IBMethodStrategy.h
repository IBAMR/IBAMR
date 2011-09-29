// Filename: IBMethodStrategy.h
// Created on 21 Sep 2011 by Boyce Griffith
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

#ifndef included_IBMethodStrategy
#define included_IBMethodStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <CoarsenSchedule.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitStrategy.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// FORWARD DECLARATION //////////////////////////

namespace IBAMR
{
class IBHierarchyIntegrator;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMethodStrategy provides a generic interface for specifying the
 * implementation details of a particular version of the IB method.
 */
class IBMethodStrategy
    : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Constructor.
     */
    IBMethodStrategy();

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBMethodStrategy();

    /*!
     * Register the IBHierarchyIntegrator object that is using this strategy
     * class.
     */
    void
    registerIBHierarchyIntegrator(
        IBHierarchyIntegrator* ib_solver);

    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    virtual const SAMRAI::hier::IntVector<NDIM>&
    getMinimumGhostCellWidth() const = 0;

    /*!
     * Virtual method to prepare to advance data from current_time to new_time.
     *
     * An empty default implementation is provided.
     */
    virtual void
    preprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);

    /*!
     * Virtual method to clean up data following call(s) to
     * integrateHierarchy().
     *
     * An empty default implementation is provided.
     */
    virtual void
    postprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    virtual void
    interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) = 0;

    /*!
     * Advance the positions of the Lagrangian structure using forward Euler.
     */
    virtual void
    eulerStep(
        double current_time,
        double new_time) = 0;

    /*!
     * Advance the positions of the Lagrangian structure using the midpoint
     * rule.
     */
    virtual void
    midpointStep(
        double current_time,
        double new_time) = 0;

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    virtual void
    computeLagrangianForce(
        double data_time) = 0;

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    virtual void
    spreadForce(
        int f_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
        double data_time) = 0;

    /*!
     * Indicate whether there are any internal fluid sources/sinks.
     *
     * A default implementation is provided that returns false.
     */
    virtual bool
    hasFluidSources() const;

    /*!
     * Compute the Lagrangian source/sink density at the specified time within
     * the current time interval.
     *
     * An empty default implementation is provided.
     */
    virtual void
    computeLagrangianFluidSource(
        double data_time);

    /*!
     * Spread the Lagrangian source/sink density to the Cartesian grid at the
     * specified time within the current time interval.
     *
     * An empty default implementation is provided.
     */
    virtual void
    spreadFluidSource(
        int q_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& q_prolongation_scheds,
        double data_time);

    /*!
     * Compute the pressures at the positions of any distributed internal fluid
     * sources or sinks.
     *
     * An empty default implementation is provided.
     */
    virtual void
    interpolatePressure(
        int p_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& p_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& p_ghost_fill_scheds,
        double data_time);

    /*!
     * Execute user-defined post-processing operations.
     *
     * An empty default implementation is provided.
     */
    virtual void
    postprocessData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     *
     * An empty default implementation is provided.
     */
    virtual void
    initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time);

    /*!
     * Register a load balancer and work load patch data index with the IB
     * strategy object.
     *
     * An empty default implementation is provided.
     */
    virtual void
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
        int workload_data_idx);

    /*!
     * Update work load estimates on each level of the patch hierarchy.
     *
     * An empty default implementation is provided.
     */
    virtual void
    updateWorkloadEstimates(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        int workload_data_idx);

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     *
     * An empty default implementation is provided.
     */
    virtual void
    beginDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     *
     * An empty default implementation is provided.
     */
    virtual void
    endDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     *
     * An empty default implementation is provided.
     */
    void
    initializeLevelData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double init_data_time,
        bool can_be_refined,
        bool initial_time,
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     *
     * An empty default implementation is provided.
     */
    void
    resetHierarchyConfiguration(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     *
     * An empty default implementation is provided.
     */
    void
    applyGradientDetector(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double error_data_time,
        int tag_index,
        bool initial_time,
        bool uses_richardson_extrapolation_too);

    /*!
     * Write out object state to the given database.
     *
     * An empty default implementation is provided.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * The IBHierarchyIntegrator object that is using this strategy class.
     */
    IBHierarchyIntegrator* d_ib_solver;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMethodStrategy(
        const IBMethodStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMethodStrategy&
    operator=(
        const IBMethodStrategy& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBMethodStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMethodStrategy
