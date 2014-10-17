// Filename: AdvDiffPredictorCorrectorHierarchyIntegrator.h
// Created on 16 Mar 2004 by Boyce Griffith
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

#ifndef included_AdvDiffPredictorCorrectorHierarchyIntegrator
#define included_AdvDiffPredictorCorrectorHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "HyperbolicLevelIntegrator.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/AdvDiffPredictorCorrectorHyperbolicPatchOps.h"
#include "ibamr/AdvectorExplicitPredictorPatchOps.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

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
    ~AdvDiffPredictorCorrectorHierarchyIntegrator();

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
    void initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1);

protected:
    /*!
     * Return the maximum stable time step size.
     */
    double getMaximumTimeStepSizeSpecialized();

    /*!
     * Reset the current data to equal the new data, update the time level of
     * the current data, and deallocate the scratch and new data.
     */
    void resetTimeDependentHierarchyDataSpecialized(double new_time);

    /*!
     * Reset the hierarchy integrator to the state at the beginning of the
     * current time step.
     */
    void resetIntegratorToPreadvanceStateSpecialized();

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
                                        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level);

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
                                          bool uses_richardson_extrapolation_too);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator(const AdvDiffPredictorCorrectorHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffPredictorCorrectorHierarchyIntegrator& operator=(const AdvDiffPredictorCorrectorHierarchyIntegrator& that);

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

#endif //#ifndef included_AdvDiffPredictorCorrectorHierarchyIntegrator
