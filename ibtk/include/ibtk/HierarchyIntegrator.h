// Filename: HierarchyIntegrator.h
// Created on 10 Aug 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#ifndef included_HierarchyIntegrator
#define included_HierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <deque>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "StandardTagAndInitStrategy.h"
#include "VariableContext.h"
#include "VisItDataWriter.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/ibtk_enums.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenPatchStrategy;
template <int DIM>
class RefinePatchStrategy;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class HierarchyIntegrator provides an abstract interface for a time
 * integrator for a system of equations defined on an AMR grid hierarchy, along
 * with basic data management for variables defined on that hierarchy.
 */
class HierarchyIntegrator : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>, public SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for class HierarchyIntegrator sets some default values,
     * reads in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     */
    HierarchyIntegrator(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        bool register_for_restart);

    /*!
     * The destructor for class HierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~HierarchyIntegrator();

    /*!
     * Return the name of the hierarchy integrator object.
     */
    const std::string& getName() const;

    /*!
     * Virtual method to initialize the variables, basic communications
     * algorithms, solvers, and other data structures used by a concrete time
     * integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.
     *
     * \note This method should be implemented so that it is safe to make
     * multiple calls to this method.  Implementations should indicate whether
     * it is possible for users to be able to make an explicit call to
     * initializeHierarchyIntegrator() prior to calling
     * initializePatchHierarchy().
     *
     * \note This method is called \em prior to the initial construction of the
     * patch hierarchy.  Consequently, implementations of this method are unable
     * to initialize patch data associated with the variables managed by the
     * integrator object, nor can they initialize hierarchy-dependent
     * communications schedules associated with the integrator.
     */
    virtual void
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) = 0;

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
    virtual void initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                          SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Integrate data on all patches on all levels of the patch hierarchy over
     * the specified time increment.
     */
    virtual void advanceHierarchy(double dt);

    /*!
     * Return the current value of the minimum time step size for the integrator
     * object.
     *
     * Subclasses can control the method used to determined the time step size
     * by overriding the protected virtual member function
     * getMinimumTimeStepSizeSpecialized().
     */
    double getMinimumTimeStepSize();

    /*!
     * Return the current value of the maximum time step size for the integrator
     * object.
     *
     * Subclasses can control the method used to determined the time step size
     * by overriding the protected virtual member function
     * getMaximumTimeStepSizeSpecialized().
     */
    double getMaximumTimeStepSize();

    /*!
     * Synchronize data defined on the grid hierarchy.
     *
     * Subclasses can control the method used to synchronize data on the grid
     * hierarchy by overriding the protected virtual member function
     * synchronizeHierarchyDataSpecialized().
     */
    void synchronizeHierarchyData(VariableContextType ctx_type);

    /*!
     * Reset the current data to equal the new data, update the time level of
     * the current data, and deallocate the scratch and new data.
     *
     * Subclasses can control the method used to reset data on the grid
     * hierarchy by overriding the protected virtual member function
     * resetTimeDependentHierarchyDataSpecialized().
     */
    void resetTimeDependentHierarchyData(double new_time);

    /*!
     * Reset the hierarchy integrator to the state at the beginning of the
     * current time step.
     *
     * Subclasses can control the method used to reset data on the grid
     * hierarchy by overriding the protected virtual member function
     * resetTimeDependentHierarchyDataSpecialized().
     */
    void resetIntegratorToPreadvanceState();

    /*!
     * Virtual method to regrid the patch hierarchy.
     *
     * A default implementation is provided that calls
     * GriddingAlgorithm::regidAllFinerLevels() to regrid the patch hierarchy.
     * Subclasses can control the method used to regrid the patch hierarchy by
     * overriding this public virtual member function.
     */
    virtual void regridHierarchy();

    /*!
     * Return a boolean value that indicates whether regridding should occur at
     * the current time step.
     *
     * Subclasses can control the method used to trigger adaptive regridding by
     * overriding the protected virtual member function
     * atRegridPointSpecialized().
     */
    bool atRegridPoint() const;

    /*!
     * Return the current integration time.
     */
    double getIntegratorTime() const;

    /*!
     * Return the initial integration time.
     */
    double getStartTime() const;

    /*!
     * Return the final integration time.
     */
    double getEndTime() const;

    /*!
     * Return the number of time steps taken by the integrator.
     */
    int getIntegratorStep() const;

    /*!
     * Return the maximum number of time steps allowed by the integrator.
     */
    int getMaxIntegratorSteps() const;

    /*!
     * Return true if any time steps remain, false otherwise.
     */
    bool stepsRemaining() const;

    /*!
     * Return a pointer to the patch hierarchy managed by the integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > getPatchHierarchy() const;

    /*!
     * Return a pointer to the gridding algorithm object managed by the integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > getGriddingAlgorithm() const;

    /*!
     * Register a VisIt data writer so the integrator can output data files that
     * may be postprocessed with the VisIt visualization tool.
     */
    void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * Get a pointer to the VisIt data writer registered with the solver.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > getVisItDataWriter() const;

    /*!
     * Prepare variables for plotting.
     *
     * Subclasses can control the method used to setup plot data by overriding
     * the protected virtual member function setupPlotData().
     *
     * \note Subclasses are allowed to require that this function be called
     * immediately before writing visualization data.
     */
    void setupPlotData();

    ///
    ///  Routines to implement the time integration scheme.
    ///

    /*!
     * Virtual method to return the number of cycles to perform for the present
     * time step.
     */
    virtual int getNumberOfCycles() const;

    /*!
     * Virtual method to return the current cycle number within the present time
     * step.
     *
     * The default implementation returns a value of -1 when it is not advancing
     * the hierarchy.
     */
    virtual int getCurrentCycleNumber() const;

    /*!
     * Virtual method to return the current time step size.
     *
     * The default implementation returns the value
     * numeric_limits<>::quiet_NaN() when it is not advancing the hierarchy.
     */
    virtual double getCurrentTimeStepSize() const;

    /*!
     * Virtual method to prepare to advance data from current_time to new_time.
     *
     * A default implementation is provided that sets the current values of
     * num_cycles and the time step size.
     */
    virtual void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Pure virtual method to advance data from current_time to new_time.
     *
     * Implementations of this virtual function are not required to synchronize
     * data on the patch hierarchy.  Data synchronization may be done
     * (optionally) in a specialization of the public virtual member function
     * postprocessIntegrateHierarchy().
     */
    virtual void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) = 0;

    /*!
     * Method to skip a cycle of the time integration scheme (e.g. for cases in
     * which time integration is handled by another class).
     */
    void skipCycle(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Virtual method to clean up data following call(s) to
     * integrateHierarchy().
     *
     * A default implementation is provided that resets the current values of
     * num_cycles and the time step size.
     */
    virtual void postprocessIntegrateHierarchy(double current_time,
                                               double new_time,
                                               bool skip_synchronize_new_state_data,
                                               int num_cycles = 1);

    /*!
     * Callback function specification to enable further specialization of
     * preprocessIntegrateHierarchy().
     */
    typedef void (*PreprocessIntegrateHierarchyCallbackFcnPtr)(double current_time,
                                                               double new_time,
                                                               int num_cycles,
                                                               void* ctx);

    /*!
     * Register a callback function to enable further specialization of
     * preprocessIntegrateHierarchy().
     */
    void registerPreprocessIntegrateHierarchyCallback(PreprocessIntegrateHierarchyCallbackFcnPtr callback,
                                                      void* ctx = NULL);

    /*!
     * Callback function specification to enable further specialization of
     * integrateHierarchy().
     */
    typedef void (*IntegrateHierarchyCallbackFcnPtr)(double current_time, double new_time, int cycle_num, void* ctx);

    /*!
     * Register a callback function to enable further specialization of
     * integrateHierarchy().
     */
    void registerIntegrateHierarchyCallback(IntegrateHierarchyCallbackFcnPtr callback, void* ctx = NULL);

    /*!
     * Callback function specification to enable further specialization of
     * postprocessIntegrateHierarchy().
     */
    typedef void (*PostprocessIntegrateHierarchyCallbackFcnPtr)(double current_time,
                                                                double new_time,
                                                                bool skip_synchronize_new_state_data,
                                                                int num_cycles,
                                                                void* ctx);

    /*!
     * Register a callback function to enable further specialization of
     * postprocessIntegrateHierarchy().
     */
    void registerPostprocessIntegrateHierarchyCallback(PostprocessIntegrateHierarchyCallbackFcnPtr callback,
                                                       void* ctx = NULL);

    /*!
     * Callback function specification to enable further specialization of
     * applyGradientDetector().
     */
    typedef void (*ApplyGradientDetectorCallbackFcnPtr)(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double error_data_time,
        int tag_index,
        bool initial_time,
        bool uses_richardson_extrapolation_too,
        void* ctx);

    /*!
     * Register a callback function to enable further specialization of
     * applyGradientDetector().
     */
    void registerApplyGradientDetectorCallback(ApplyGradientDetectorCallbackFcnPtr callback, void* ctx = NULL);

    ///
    ///  Implementations of functions declared in the
    ///  SAMRAI::mesh::StandardTagAndInitStrategy abstract base class.
    ///

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \note Subclasses should not override the implementation of this function
     * provided by class HierarchyIntegrator.  Instead, they should override the
     * protected virtual member function initializeLevelDataSpecialized().
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level =
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
                             bool allocate_data = true);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \note Subclasses should not override the implementation of this function
     * provided by class HierarchyIntegrator.  Instead, they should override the
     * protected virtual member function
     * resetHierarchyConfigurationSpecialized().
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_level,
                                     int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \note Subclasses should not override the implementation of this function
     * provided by class HierarchyIntegrator.  Instead, they should override the
     * protected virtual member function applyGradientDetectorSpecialized().
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too);

    ///
    ///  Routines to access to the variable contexts maintained by the
    ///  integrator.
    ///

    /*!
     * Return a pointer the variable context corresponding to the specified
     * variable context type.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> getContext(VariableContextType ctx_type) const;

    /*!
     * Return a pointer to the "current" variable context used by integrator.
     * Current data corresponds to state data at the beginning of a time step,
     * or when a new level is initialized.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> getCurrentContext() const;

    /*!
     * Return a pointer to the "new" variable context used by integrator.  New
     * data corresponds to advanced state data at the end of a time step.  The
     * data is one time step later than the "current" data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> getNewContext() const;

    /*!
     * Return a pointer to the "scratch" variable context used by integrator.
     * Scratch data typically corresponds to storage that user-routines in the
     * concrete GodunovAdvector object manipulate; in particular, scratch data
     * contains ghost cells.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> getScratchContext() const;

    /*!
     * Check whether a patch data index corresponds to allocated data over the
     * specified range of patch level numbers.
     *
     * NOTE: This method will return "false" without error for invalid (i.e.,
     * negative) patch data indices.
     */
    bool isAllocatedPatchData(int data_idx, int coarsest_ln = -1, int finest_ln = -1) const;

    /*!
     * Allocate a patch data index over the specified range of patch level
     * numbers.
     */
    void allocatePatchData(int data_idx, double data_time, int coarsest_ln = -1, int finest_ln = -1) const;

    /*!
     * Deallocate a patch data index over the specified range of patch level
     * numbers.
     */
    void deallocatePatchData(int data_idx, int coarsest_ln = -1, int finest_ln = -1) const;

    ///
    ///  Routines to access utility classeses managed by the integrator.
    ///

    SAMRAI::tbox::Pointer<HierarchyMathOps> getHierarchyMathOps() const;

    ///
    ///  Implementations of functions declared in the SAMRAI::tbox::Serializable
    ///  abstract base class.
    ///

    /*!
     * Write out object state to the given database.
     *
     * \note Subclasses should not override the implementation of this function
     * provided by class HierarchyIntegrator.  Instead, they should override the
     * protected virtual member function putToDatabaseSpecialized().
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * Virtual method to compute an implementation-specific minimum stable time
     * step size.
     *
     * A default implementation is provided that returns
     * min(dt_max,dt_growth_factor*dt_current).  The growth condition prevents
     * excessive changes in the time step size as the computation progresses.
     */
    virtual double getMinimumTimeStepSizeSpecialized();

    /*!
     * Virtual method to compute an implementation-specific maximum stable time
     * step size.
     *
     * A default implementation is provided that returns
     * min(dt_max,dt_growth_factor*dt_current).  The growth condition prevents
     * excessive changes in the time step size as the computation progresses.
     */
    virtual double getMaximumTimeStepSizeSpecialized();

    /*!
     * Virtual method to perform implementation-specific data synchronization.
     *
     * A default implementation is provided that synchronizes state data
     * registered with the HierarchyIntegrator object using the coarsen
     * operations specified by calls to registerVariable().
     */
    virtual void synchronizeHierarchyDataSpecialized(VariableContextType ctx_type);

    /*!
     * Virtual method to perform implementation-specific data reset operations.
     *
     * A default implementation is provided that first swaps the current and new
     * PatchData pointers, and then deallocates the new and scratch data
     * contexts.
     */
    virtual void resetTimeDependentHierarchyDataSpecialized(double new_time);

    /*!
     * Virtual method to perform implementation-specific data reset operations.
     *
     * A default implementation is provided that deallocates the new and scratch
     * data contexts when data associated with those contexts have been
     * allocated.
     */
    virtual void resetIntegratorToPreadvanceStateSpecialized();

    /*!
     * Virtual method to provide implementation-specific function to determine
     * whether regridding should occur at the current time step.
     *
     * A default implementation is provided that indicates that the hierarchy
     * should be regridded at a fixed integer interval of time steps unless a
     * parent integrator has been registered with this integrator.  If a parent
     * integrator has been registered with this integrator,
     * atRegridPointSpecialized() returns false, in order to allow the parent
     * integrator to control the timing of regridding.
     */
    virtual bool atRegridPointSpecialized() const;

    /*!
     * Virtual method to perform implementation-specific visualization setup.
     *
     * An empty default implementation is provided.
     */
    virtual void setupPlotDataSpecialized();

    /*!
     * Virtual method to perform implementation-specific data initialization on
     * a new level after it is inserted into an AMR patch hierarchy by the
     * gridding algorithm.
     *
     * An empty default implementation is provided.
     */
    virtual void
    initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                   int level_number,
                                   double init_data_time,
                                   bool can_be_refined,
                                   bool initial_time,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                                   bool allocate_data);

    /*!
     * Virtual method to perform implementation-specific data reset operations.
     *
     * An empty default implementation is provided.
     */
    virtual void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level);

    /*!
     * Virtual method to perform implementation-specific cell tagging
     * operations.
     *
     * An empty default implementation is provided.
     */
    virtual void
    applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int level_number,
                                     double error_data_time,
                                     int tag_index,
                                     bool initial_time,
                                     bool uses_richardson_extrapolation_too);

    /*!
     * Protecethod to write implementation-specific object state to a database.
     *
     * An empty default implementation is provided.
     */
    virtual void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Execute any user-specified preprocessIntegrateHierarchy callback
     * functions.
     */
    virtual void executePreprocessIntegrateHierarchyCallbackFcns(double current_time, double new_time, int num_cycles);

    /*!
     * Execute any user-specified integrateHierarchy callback functions.
     */
    virtual void executeIntegrateHierarchyCallbackFcns(double current_time, double new_time, int cycle_num);

    /*!
     * Execute any user-specified postprocessIntegrateHierarchy callback
     * functions.
     */
    virtual void executePostprocessIntegrateHierarchyCallbackFcns(double current_time,
                                                                  double new_time,
                                                                  bool skip_synchronize_new_state_data,
                                                                  int num_cycles);

    /*!
     * Execute any user-specified applyGradientDetector callback functions.
     */
    virtual void
    executeApplyGradientDetectorCallbackFcns(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                             int level_number,
                                             double error_data_time,
                                             int tag_index,
                                             bool initial_time,
                                             bool uses_richardson_extrapolation_too);

    /*!
     * Register a state variable with the integrator.  When a refine operator is
     * specified, the data for the variable are automatically maintained as the
     * patch hierarchy evolves.
     *
     * All state variables are registered with three contexts: current, new, and
     * scratch.  The current context of a state variable is maintained from time
     * step to time step and, if the necessary coarsen and refine operators are
     * specified, as the patch hierarchy evolves.
     */
    void
    registerVariable(int& current_idx,
                     int& new_idx,
                     int& scratch_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
                     const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts = SAMRAI::hier::IntVector<NDIM>(0),
                     const std::string& coarsen_name = "NO_COARSEN",
                     const std::string& refine_name = "NO_REFINE",
                     SAMRAI::tbox::Pointer<CartGridFunction> init_fcn = SAMRAI::tbox::Pointer<CartGridFunction>(NULL));

    /*!
     * Register a variable with the integrator that may not be maintained from
     * time step to time step.
     *
     * By default, variables are registered with the scratch context, which is
     * deallocated after each time step.
     */
    void registerVariable(int& idx,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
                          const SAMRAI::hier::IntVector<NDIM>& ghosts = SAMRAI::hier::IntVector<NDIM>(0),
                          SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx =
                              SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>(NULL));

    /*!
     * Register a ghost cell-filling refine algorithm.
     */
    void registerGhostfillRefineAlgorithm(const std::string& name,
                                          SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > ghostfill_alg,
                                          SAMRAI::xfer::RefinePatchStrategy<NDIM>* ghostfill_patch_strategy = NULL);

    /*!
     * Register a data-prolonging refine algorithm.
     */
    void registerProlongRefineAlgorithm(const std::string& name,
                                        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > prolong_alg,
                                        SAMRAI::xfer::RefinePatchStrategy<NDIM>* prolong_patch_strategy = NULL);

    /*!
     * Register a coarsen algorithm.
     */
    void registerCoarsenAlgorithm(const std::string& name,
                                  SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > coarsen_alg,
                                  SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* coarsen_patch_strategy = NULL);

    /*!
     * Get ghost cell-filling refine algorithm.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> >
    getGhostfillRefineAlgorithm(const std::string& name) const;

    /*!
     * Get data-prolonging refine algorithm.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> >
    getProlongRefineAlgorithm(const std::string& name) const;

    /*!
     * Get coarsen algorithm.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > getCoarsenAlgorithm(const std::string& name) const;

    /*!
     * Get ghost cell-filling refine schedules.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
    getGhostfillRefineSchedules(const std::string& name) const;

    /*!
     * Get data-prolonging refine schedules.
     *
     * \note These schedules are allocated only for level numbers >= 1.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >&
    getProlongRefineSchedules(const std::string& name) const;

    /*!
     * Get coarsen schedules.
     *
     * \note These schedules are allocated only for level numbers >= 1.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >&
    getCoarsenSchedules(const std::string& name) const;

    /*!
     * Register a "child" integrator object with this integrator object.
     *
     * \note Multiple child integrator objects may be registered with a single
     * parent integrator object.
     */
    void registerChildHierarchyIntegrator(HierarchyIntegrator* child_integrator);

    /*!
     * Register a "parent" integrator object with this integrator object.
     *
     * \note Only a single parent integrator object may be registered with a
     * particular child integrator object.
     */
    void registerParentHierarchyIntegrator(HierarchyIntegrator* parent_integrator);

    /*!
     * Build the HierarchyMathOps object.
     */
    SAMRAI::tbox::Pointer<HierarchyMathOps>
    buildHierarchyMathOps(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Setup the tag buffer.
     */
    void setupTagBuffer(SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Returns true when we are regridding the patch hierarchy.
     */
    bool regriddingHierarchy() const
    {
        return d_regridding_hierarchy;
    }

    /*!
     * Returns true when we are executing a time step in which a regridding
     * operation was performed.
     */
    bool atRegridTimeStep() const
    {
        return d_at_regrid_time_step;
    }

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

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this time integration object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * Indicates whether the hierarchy has been initialized.
     */
    bool d_hierarchy_is_initialized;

    /*
     * Collection of child integrator objects.
     */
    HierarchyIntegrator* d_parent_integrator;
    std::set<HierarchyIntegrator*> d_child_integrators;

    /*
     * The object used to write out data for postprocessing by the VisIt
     * visualization tool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;

    /*
     * Time and time step size data read from input or set at initialization.
     */
    double d_integrator_time, d_start_time, d_end_time;
    double d_dt_init, d_dt_min, d_dt_max, d_dt_growth_factor;
    int d_integrator_step, d_max_integrator_steps;
    std::deque<double> d_dt_previous;

    /*
     * The number of cycles of fixed-point iteration to use per timestep.
     */
    int d_num_cycles;

    /*
     * The number of cycles for the current time step, the current cycle number,
     * and the current time step size.
     */
    int d_current_num_cycles, d_current_cycle_num;
    double d_current_dt;

    /*
     * The number of integration steps taken between invocations of the
     * regridding process.
     */
    int d_regrid_interval;

    /*
     * The regrid mode.  "Standard" regridding involves only one call to
     * SAMRAI::mesh::GriddingAlgorithm::regridAllFinerLevels().  This limits the
     * amount that the AMR grid hierarchy can change within a single call to
     * regridHierarchy().  "Agressive" regridding involes multiple calls to
     * SAMRAI::mesh::GriddingAlgorithm::regridAllFinerLevels(), effectively
     * allowing arbitrary changes to the grid hierarchy configuration within a
     * single call to regridHierarchy().
     */
    RegridMode d_regrid_mode;

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_enable_logging;

    /*
     * The type of extrapolation to use at physical boundaries when prolonging
     * data during regridding.
     */
    std::string d_bdry_extrap_type;

    /*
     * The number of cells on each level by which tagged cells will be buffered
     * after they have selected for refinement.  The tag buffer helps to
     * guarantee that refined cells near important features in the solution will
     * remain refined until the level is regridded next.
     */
    SAMRAI::tbox::Array<int> d_tag_buffer;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<HierarchyMathOps> d_hier_math_ops;
    bool d_manage_hier_math_ops;

    /*
     * SAMRAI::hier::Variable lists and SAMRAI::hier::ComponentSelector objects
     * are used for data management.
     */
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_state_variables;
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_scratch_variables;

    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_copy_scratch_to_current_fast;
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_copy_scratch_to_current_slow;

    SAMRAI::hier::ComponentSelector d_current_data, d_new_data, d_scratch_data;

    std::map<SAMRAI::hier::Variable<NDIM>*, SAMRAI::tbox::Pointer<CartGridFunction> > d_state_var_init_fcns;

    /*!
     * Variable contexts.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_current_context, d_new_context, d_scratch_context;

    /*!
     * Names of special coarsen algorithms/schedules.
     */
    static const std::string SYNCH_CURRENT_DATA_ALG, SYNCH_NEW_DATA_ALG;

    /*!
     * Regridding-related communications algorithms and other data structures.
     */
    SAMRAI::hier::ComponentSelector d_fill_after_regrid_bc_idxs;
    SAMRAI::xfer::RefineAlgorithm<NDIM> d_fill_after_regrid_prolong_alg;
    SAMRAI::xfer::RefinePatchStrategy<NDIM>* d_fill_after_regrid_phys_bdry_bc_op;

    /*!
     * Callback functions and callback function contexts.
     */
    std::vector<PreprocessIntegrateHierarchyCallbackFcnPtr> d_preprocess_integrate_hierarchy_callbacks;
    std::vector<void*> d_preprocess_integrate_hierarchy_callback_ctxs;
    std::vector<IntegrateHierarchyCallbackFcnPtr> d_integrate_hierarchy_callbacks;
    std::vector<void*> d_integrate_hierarchy_callback_ctxs;
    std::vector<PostprocessIntegrateHierarchyCallbackFcnPtr> d_postprocess_integrate_hierarchy_callbacks;
    std::vector<void*> d_postprocess_integrate_hierarchy_callback_ctxs;
    std::vector<ApplyGradientDetectorCallbackFcnPtr> d_apply_gradient_detector_callbacks;
    std::vector<void*> d_apply_gradient_detector_callback_ctxs;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyIntegrator(const HierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyIntegrator& operator=(const HierarchyIntegrator& that);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*
     * Indicates whether we are currently regridding the hierarchy, or whether
     * the time step began by regridding the hierarchy.
     */
    bool d_regridding_hierarchy; // true only when we are regridding
    bool d_at_regrid_time_step;  // true for the duration of a time step that included a regrid
                                 // operation

    /*
     * Cached communications algorithms, strategies, and schedules.
     */
    typedef std::map<std::string, SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > > RefineAlgorithmMap;
    typedef std::map<std::string, SAMRAI::xfer::RefinePatchStrategy<NDIM>*> RefinePatchStrategyMap;
    typedef std::map<std::string, std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > >
        RefineScheduleMap;

    RefineAlgorithmMap d_ghostfill_algs;
    RefinePatchStrategyMap d_ghostfill_strategies;
    RefineScheduleMap d_ghostfill_scheds;

    RefineAlgorithmMap d_prolong_algs;
    RefinePatchStrategyMap d_prolong_strategies;
    RefineScheduleMap d_prolong_scheds;

    typedef std::map<std::string, SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > > CoarsenAlgorithmMap;
    typedef std::map<std::string, SAMRAI::xfer::CoarsenPatchStrategy<NDIM>*> CoarsenPatchStrategyMap;
    typedef std::map<std::string, std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > >
        CoarsenScheduleMap;

    CoarsenAlgorithmMap d_coarsen_algs;
    CoarsenPatchStrategyMap d_coarsen_strategies;
    CoarsenScheduleMap d_coarsen_scheds;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HierarchyIntegrator
