// Filename: IBHierarchyIntegrator.h
// Created on 12 Jul 2004 by Boyce Griffith
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

#ifndef included_IBHierarchyIntegrator
#define included_IBHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/IBEulerianForceFunction.h>
#include <ibamr/IBEulerianSourceFunction.h>
#include <ibamr/IBMethodStrategy.h>
#include <ibamr/INSHierarchyIntegrator.h>

// IBTK INCLUDES
#include <ibtk/LDataManager.h>

/////////////////////////////// FORWARD DECLARATION //////////////////////////

namespace IBAMR
{
class IBMethod;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHierarchyIntegrator is an implementation of a formally
 * second-order accurate, semi-implicit version of the immersed boundary method.
 */
class IBHierarchyIntegrator
    : public HierarchyIntegrator
{
public:
    friend class IBMethod;

    /*!
     * The constructor for class IBHierarchyIntegrator sets some default values,
     * reads in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     */
    IBHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<IBMethodStrategy> ib_method_ops,
        SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
        bool register_for_restart=true);

    /*!
     * The destructor for class IBHierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~IBHierarchyIntegrator();

    /*!
     * Supply a body force (optional).
     */
    void
    registerBodyForceFunction(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Register a load balancer for non-uniform load balancing.
     */
    void
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer);

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
    initializeHierarchyIntegrator(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
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
    void
    initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int
    getNumberOfCycles();

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void
    preprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        int num_cycles=1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        double current_time,
        double new_time,
        int cycle_num=0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        bool skip_synchronize_new_state_data,
        int num_cycles=1);

    /*!
     * Regrid the hierarchy.
     */
    void
    regridHierarchy();

    /*!
     * Execute user-defined post-processing operations.
     */
    void
    postProcessData();

protected:
    /*!
     * Function to determine whether to regridding should occur at the current
     * time step.
     */
    bool
    atRegridPointSpecialized() const;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void
    initializeLevelDataSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
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
    resetHierarchyConfigurationSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void
    applyGradientDetectorSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double error_data_time,
        int tag_index,
        bool initial_time,
        bool uses_richardson_extrapolation_too);

    /*!
     * Write out specialized object state to the given database.
     */
    void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHierarchyIntegrator(
        const IBHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHierarchyIntegrator&
    operator=(
        const IBHierarchyIntegrator& that);

    /*!
     * Read input values from a given database.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void
    getFromRestart();

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*
     * The INSHierarchyIntegrator is used to provide time integration
     * capability for the incompressible Navier-Stokes equations.
     */
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> d_ins_hier_integrator;

    /*
     * The regrid CFL interval indicates the number of meshwidths a particle may
     * move in any coordinate direction between invocations of the regridding
     * process.
     *
     * NOTE: Currently, when the CFL-based regrid interval is specified, it is
     * used instead of the fixed regrid interval.
     */
    double d_regrid_cfl_interval, d_regrid_cfl_estimate;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM,double> > d_hier_velocity_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_pressure_cc_data_ops;

    /*
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_var, d_f_var, d_q_var;
    int d_u_idx, d_f_idx, d_q_idx;

    /*
     * Body force functions.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_body_force_fcn;
    SAMRAI::tbox::Pointer<IBEulerianForceFunction> d_eulerian_force_fcn;

    /*
     * The source/sink distribution functions.
     */
    SAMRAI::tbox::Pointer<IBEulerianSourceFunction> d_eulerian_source_fcn;

    /*
     * IB method strategy object.
     */
    SAMRAI::tbox::Pointer<IBMethodStrategy> d_ib_method_ops;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHierarchyIntegrator
