// Filename: INSCollocatedHierarchyIntegrator.h
// Created on 24 Aug 2011 by Boyce Griffith
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

#ifndef included_INSCollocatedHierarchyIntegrator
#define included_INSCollocatedHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "CellVariable.h"
#include "FaceVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
class ConvectiveOperator;
} // namespace IBAMR
namespace IBTK
{
class PoissonSolver;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class BasePatchLevel;
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
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
 * \brief Class INSCollocatedHierarchyIntegrator provides a collocated
 * projection method-based solver for the incompressible Navier-Stokes equations
 * on an AMR grid hierarchy.
 */
class INSCollocatedHierarchyIntegrator : public INSHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSCollocatedHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSCollocatedHierarchyIntegrator(const std::string& object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                     bool register_for_restart = true);

    /*!
     * The destructor for class INSCollocatedHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSCollocatedHierarchyIntegrator();

    /*!
     * Get the convective operator being used by this solver class.
     *
     * If the time integrator is configured to solve the time-dependent
     * (creeping) Stokes equations, then the returned pointer will be NULL.
     *
     * If the convective operator has not already been constructed, and if the
     * time integrator is not configured to solve the time-dependent (creeping)
     * Stokes equations, then this function will initialize the default type of
     * convective operator, which may be set in the class input database.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator();

    /*!
     * Get the subdomain solver for the velocity subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getVelocitySubdomainSolver();

    /*!
     * Get the subdomain solver for the pressure subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getPressureSubdomainSolver();

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

    /*!
     * Regrid the patch hierarchy.
     */
    void regridHierarchy();

protected:
    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

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
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too);

    /*!
     * Prepare variables for plotting.
     */
    void setupPlotDataSpecialized();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSCollocatedHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSCollocatedHierarchyIntegrator(const INSCollocatedHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSCollocatedHierarchyIntegrator& operator=(const INSCollocatedHierarchyIntegrator& that);

    /*!
     * Compute the appropriate source term that must be added to the momentum
     * equation when the fluid contains internal sources and sinks.
     */
    void computeDivSourceTerm(int F_idx, int Q_idx, int U_idx);

    /*!
     * Reinitialize the operators and solvers used by the hierarchy integrator.
     */
    void reinitializeOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Project the velocity field following a regridding operation.
     */
    void regridProjection();

    /*!
     * Value determining the type of projection method to use.
     */
    ProjectionMethodType d_projection_method_type;

    /*!
     * Boolean indicating whether to use the Brown-Cortez-Minion 2nd-order
     * pressure update.
     */
    bool d_using_2nd_order_pressure_update;

    /*!
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;

    /*
     * Hierarchy operators and solvers.
     */
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_adv_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_N_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_Phi_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_Phi_rhs_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_U_nul_vecs;
    bool d_vectors_need_init;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_Phi_bdry_bc_fill_op;

    /*!
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_U_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_u_ADV_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_N_old_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_U_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_u_ADV_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_Norm_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Grad_P_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Grad_Phi_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_Grad_Phi_fc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_F_div_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_U_current_idx, d_U_new_idx, d_U_scratch_idx;
    int d_u_ADV_current_idx, d_u_ADV_new_idx, d_u_ADV_scratch_idx;
    int d_P_current_idx, d_P_new_idx, d_P_scratch_idx;
    int d_F_current_idx, d_F_new_idx, d_F_scratch_idx;
    int d_Q_current_idx, d_Q_new_idx, d_Q_scratch_idx;
    int d_N_old_current_idx, d_N_old_new_idx, d_N_old_scratch_idx;

    /*
     * Patch data descriptor indices for all "plot" variables managed by the
     * integrator.
     *
     * Plot variables have one context: current.
     */
    int d_Omega_idx, d_Div_U_idx, d_Div_u_ADV_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_Omega_Norm_idx, d_Grad_P_idx, d_Phi_idx, d_Grad_Phi_cc_idx, d_Grad_Phi_fc_idx, d_F_div_idx;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSCollocatedHierarchyIntegrator
