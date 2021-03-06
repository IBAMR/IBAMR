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

#ifndef included_IBAMR_AdvDiffPredictorCorrectorHierarchyIntegrator
#define included_IBAMR_AdvDiffPredictorCorrectorHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"

#include "HyperbolicLevelIntegrator.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class BasePatchLevel;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffPredictorCorrectorHierarchyIntegrator manages the spatial
 * discretization and time integration of scalar- and vector-valued quantities
 * whose dynamics are governed by the advection-diffusion equation.
 *
 * Each quantity \f$ Q \f$ managed by the integrator may have a unique diffusion
 * coefficient \f$ \kappa \f$ and damping coefficient \f$ \lambda \f$, and may
 * optionally have a forcing term \f$ F \f$.  Additionally, a different
 * advection velocity may be used with each quantity registered with the
 * integrator.
 *
 * This hierarchy integrator advances all levels of the patch hierarchy
 * synchronously in time.  In particular, subcycling in time is \em not
 * performed.
 *
 * Either Crank-Nicolson (i.e., the trapezoidal rule) or backward Euler is used
 * for the linearly implicit treatment of the diffusive terms.  The advective
 * terms are discretized by the AdvectorExplicitPredictorPatchOps object supplied to the class
 * constructor.
 *
 * \see AdvDiffPredictorCorrectorHyperbolicPatchOps
 * \see HierarchyIntegrator
 * \see AdvectorExplicitPredictorPatchOps
 * \see SAMRAI::algs::HyperbolicLevelIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class AdvDiffPredictorCorrectorHierarchyIntegrator : public AdvDiffHierarchyIntegrator
{
public:
    /*!
     * The constructor for class AdvDiffPredictorCorrectorHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<AdvectorExplicitPredictorPatchOps> explicit_predictor,
        bool register_for_restart = true);

    /*!
     * The destructor for class AdvDiffPredictorCorrectorHierarchyIntegrator unregisters
     * the integrator object with the restart manager when the object is so
     * registered.
     */
    ~AdvDiffPredictorCorrectorHierarchyIntegrator() = default;

    /*!
     * Return a pointer to the level integrator object used to integrate the
     * advective terms.
     */
    SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> > getHyperbolicLevelIntegrator() const;

    /*!
     * Return a pointer to the patch strategy object used to specify the
     * numerical routines used to integrate the advective terms.
     */
    SAMRAI::tbox::Pointer<AdvDiffPredictorCorrectorHyperbolicPatchOps> getHyperbolicPatchStrategy() const;

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
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

protected:
    /*!
     * Return the maximum stable time step size.
     */
    double getMaximumTimeStepSizeSpecialized() override;

    /*!
     * Reset the current data to equal the new data, update the time level of
     * the current data, and deallocate the scratch and new data.
     */
    void resetTimeDependentHierarchyDataSpecialized(double new_time) override;

    /*!
     * Reset the hierarchy integrator to the state at the beginning of the
     * current time step.
     */
    void resetIntegratorToPreadvanceStateSpecialized() override;

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
     * should occur according to gradient criteria specified by the
     * AdvectorExplicitPredictorPatchOps object.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator(const AdvDiffPredictorCorrectorHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator&
    operator=(const AdvDiffPredictorCorrectorHierarchyIntegrator& that) = delete;

    /*
     * The SAMRAI::algs::HyperbolicLevelIntegrator supplies generic operations
     * use to handle the explicit integration of advection terms.
     *
     * The advection patch strategy supplies the advection-specific operations
     * needed to treat data on patches in the AMR grid hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> > d_hyp_level_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_hyp_level_integrator_db;
    SAMRAI::tbox::Pointer<AdvDiffPredictorCorrectorHyperbolicPatchOps> d_hyp_patch_ops;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_hyp_patch_ops_db;
    SAMRAI::tbox::Pointer<AdvectorExplicitPredictorPatchOps> d_explicit_predictor;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_AdvDiffPredictorCorrectorHierarchyIntegrator
