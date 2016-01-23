// Filename: IBHierarchyIntegrator.h
// Created on 12 Jul 2004 by Boyce Griffith
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

#ifndef included_IBHierarchyIntegrator
#define included_IBHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "MultiblockDataTranslator.h"
#include "PatchLevel.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "Variable.h"
#include "VariableContext.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/LMarkerSetVariable.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"

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
 */
class IBHierarchyIntegrator : public IBTK::HierarchyIntegrator
{
public:
    friend class IBStrategy;

    /*!
     * The destructor for class IBHierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~IBHierarchyIntegrator();

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
    void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer);

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
     * Basic functions to prepare to advance data from current_time to new_time.
     *
     * A default implementation is provided that sets the current values of
     * num_cycles and the time step size and checks to see if the time step size
     * has changed.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

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
     * Regrid the hierarchy.
     */
    void regridHierarchy();

protected:
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
    bool atRegridPointSpecialized() const;

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
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*!
     * Enum indicating the time integration employed for the IB equations.
     */
    TimeSteppingType d_time_stepping_type;

    /*!
     * Flags to determine whether warnings or error messages should be emitted
     * when time step size changes are encountered.
     */
    bool d_error_on_dt_change, d_warn_on_dt_change;

    /*
     * The (optional) INSHierarchyIntegrator is used to provide time integration
     * capability for the incompressible Navier-Stokes equations.
     */
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> d_ins_hier_integrator;

    /*
     * The regrid CFL interval indicates the number of meshwidths a particle may
     * move in any coordinate direction between invocations of the regridding
     * process.
     *
     * NOTE: Currently, when the CFL-based regrid interval is specified, it is
     * always used instead of the fixed-step regrid interval.
     */
    double d_regrid_cfl_interval, d_regrid_cfl_estimate;

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
    int d_u_idx, d_p_idx, d_f_idx, d_f_current_idx, d_q_idx;
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_ib_context;

    /*
     * Refine and coarsen algorithm data.
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

    /*
     * Nonuniform load balancing data structures.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_workload_var;
    int d_workload_idx;

    /*
     * Lagrangian marker data structures.
     */
    SAMRAI::tbox::Pointer<IBTK::LMarkerSetVariable> d_mark_var;
    int d_mark_current_idx, d_mark_new_idx, d_mark_scratch_idx;
    std::vector<IBTK::Point> d_mark_init_posns;
    std::string d_mark_file_name;

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
        ~IBEulerianForceFunction();

        /*!
         * \name Methods to set the data.
         */
        //\{

        /*!
         * \note This concrete IBTK::CartGridFunction is time-dependent.
         */
        bool isTimeDependent() const;

        /*!
         * \brief Set the data on the patch interiors on the specified levels of
         * the patch hierarchy.
         */
        void setDataOnPatchHierarchy(const int data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     const double data_time,
                                     const bool initial_time = false,
                                     const int coarsest_ln = -1,
                                     const int finest_ln = -1);

        /*!
         * Set the data on the patch interior.
         */
        void setDataOnPatch(int data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            double data_time,
                            bool initial_time = false,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        IBEulerianForceFunction();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        IBEulerianForceFunction(const IBEulerianForceFunction& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        IBEulerianForceFunction& operator=(const IBEulerianForceFunction& that);

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
        ~IBEulerianSourceFunction();

        /*!
         * \name Methods to set the data.
         */
        //\{

        /*!
         * \note This concrete IBTK::CartGridFunction is time-dependent.
         */
        bool isTimeDependent() const;

        /*!
         * Set the data on the patch interior.
         */
        void setDataOnPatch(int data_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            double data_time,
                            bool initial_time = false,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        IBEulerianSourceFunction();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        IBEulerianSourceFunction(const IBEulerianSourceFunction& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        IBEulerianSourceFunction& operator=(const IBEulerianSourceFunction& that);

        const IBHierarchyIntegrator* const d_ib_solver;
    };

    friend class IBEulerianSourceFunction;

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
    IBHierarchyIntegrator(const IBHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHierarchyIntegrator& operator=(const IBHierarchyIntegrator& that);

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

#endif //#ifndef included_IBHierarchyIntegrator
