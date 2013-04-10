// Filename: IBFEMethod.h
// Created on 5 Oct 2011 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#ifndef included_IBFEMethod
#define included_IBFEMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/IBStrategy.h>

// IBTK INCLUDES
#include <ibtk/FEDataManager.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEMethod is an implementation of the abstract base class
 * IBStrategy that provides functionality required by the IB method with finite
 * element elasticity.
 */
class IBFEMethod
    : public IBStrategy
{
public:
    static const std::string        COORDS_SYSTEM_NAME;
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string         FORCE_SYSTEM_NAME;
    static const std::string      VELOCITY_SYSTEM_NAME;
    static const std::string BODY_VELOCITY_SYSTEM_NAME;

    /*!
     * \brief Constructor.
     */
    IBFEMethod(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        libMesh::Mesh* mesh,
        int max_level_number,
        bool register_for_restart=true);

    /*!
     * \brief Constructor.
     */
    IBFEMethod(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::vector<libMesh::Mesh*>& meshes,
        int max_level_number,
        bool register_for_restart=true);

    /*!
     * \brief Destructor.
     */
    ~IBFEMethod();

    /*!
     * Return a pointer to the finite element data manager object for the
     * specified part.
     */
    IBTK::FEDataManager*
    getFEDataManager(
        unsigned int part=0) const;

    /*!
     * Indicate that a part is constrained.
     */
    void
    registerConstrainedPart(
        unsigned int part=0);

    /*!
     * Typedef specifying interface for specifying constrained body velocities.
     */
    typedef
    void
    (*ConstrainedPartVelocityFcnPtr)(
        libMesh::NumericVector<double>& U_b,
        libMesh::NumericVector<double>& U,
        libMesh::NumericVector<double>& X,
        libMesh::MeshBase& mesh,
        double time,
        void* ctx);

    /*!
     * Register a constrained body velocity function.
     */
    void
    registerConstrainedPartVelocityFunction(
        ConstrainedPartVelocityFcnPtr constrained_part_velocity_fcn,
        void* constrained_part_velocity_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for coordinate mapping function.
     */
    typedef
    void
    (*CoordinateMappingFcnPtr)(
        libMesh::Point& X,
        const libMesh::Point& s,
        void* ctx);

    /*!
     * Register the (optional) function used to initialize the physical
     * coordinates from the Lagrangian coordinates.
     *
     * \note If no function is provided, the initial physical coordinates are
     * taken to be the same as the Lagrangian coordinate system, i.e., the
     * initial coordinate mapping is assumed to be the identity mapping.
     */
    void
    registerInitialCoordinateMappingFunction(
        CoordinateMappingFcnPtr coordinate_mapping_fcn,
        void* coordinate_mapping_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for PK1 stress tensor function.
     */
    typedef
    void
    (*PK1StressFcnPtr)(
        libMesh::TensorValue<double>& PP,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute the first Piola-Kirchhoff
     * stress tensor, used to compute the forces on the Lagrangian finite
     * element mesh.
     */
    void
    registerPK1StressTensorFunction(
        PK1StressFcnPtr PK1_stress_fcn,
        std::vector<unsigned int> PK1_stress_fcn_systems=std::vector<unsigned int>(),
        void* PK1_stress_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for Lagrangian body force distribution
     * function.
     */
    typedef
    void
    (*LagBodyForceFcnPtr)(
        libMesh::VectorValue<double>& F,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     */
    void
    registerLagBodyForceFunction(
        LagBodyForceFcnPtr lag_body_force_fcn,
        std::vector<unsigned int> lag_body_force_fcn_systems=std::vector<unsigned int>(),
        void* lag_body_force_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for Lagrangian pressure force distribution
     * function.
     */
    typedef
    void
    (*LagPressureFcnPtr)(
        double& P,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        unsigned short int side,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute surface pressure
     * distributions on the Lagrangian finite element mesh.
     */
    void
    registerLagPressureFunction(
        LagPressureFcnPtr lag_pressure_fcn,
        std::vector<unsigned int> lag_pressure_fcn_systems=std::vector<unsigned int>(),
        void* lag_pressure_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Typedef specifying interface for Lagrangian surface force distribution
     * function.
     */
    typedef
    void
    (*LagSurfaceForceFcnPtr)(
        libMesh::VectorValue<double>& F,
        const libMesh::TensorValue<double>& FF,
        const libMesh::Point& X,
        const libMesh::Point& s,
        libMesh::Elem* elem,
        unsigned short int side,
        libMesh::NumericVector<double>& X_vec,
        const std::vector<libMesh::NumericVector<double>*>& system_data,
        double time,
        void* ctx);

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     */
    void
    registerLagSurfaceForceFunction(
        LagSurfaceForceFcnPtr lag_surface_force_fcn,
        std::vector<unsigned int> lag_surface_force_fcn_systems=std::vector<unsigned int>(),
        void* lag_surface_force_fcn_ctx=NULL,
        unsigned int part=0);

    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>&
    getMinimumGhostCellWidth() const;

    /*!
     * Setup the tag buffer.
     */
    void
    setupTagBuffer(
        SAMRAI::tbox::Array<int>& tag_buffer,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) const;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void
    preprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void
    interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time);

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void
    eulerStep(
        double current_time,
        double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void
    midpointStep(
        double current_time,
        double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    void
    trapezoidalStep(
        double current_time,
        double new_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void
    computeLagrangianForce(
        double data_time);

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(
        int f_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
        double data_time);

    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    void
    initializeFEData();

    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void
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
     */
    void
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
        int workload_data_idx);

    /*!
     * Update work load estimates on each level of the patch hierarchy.
     */
    void
    updateWorkloadEstimates(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        int workload_data_idx);

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    void
    beginDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    void
    endDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
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
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*
     * \brief Compute the constraint force density.
     */
    void
    computeConstraintForceDensity(
        libMesh::PetscVector<double>& F_vec,
        libMesh::PetscVector<double>& U_b_vec,
        libMesh::PetscVector<double>& U_vec,
        double time,
        unsigned int part);

    /*
     * \brief Compute the interior elastic density, possibly splitting off the
     * normal component of the transmission force along the physical boundary of
     * the Lagrangian structure.
     */
    void
    computeInteriorForceDensity(
        libMesh::PetscVector<double>& G_vec,
        libMesh::PetscVector<double>& X_vec,
        double time,
        unsigned int part);

    /*!
     * \brief Spread the transmission force density along the physical boundary
     * of the Lagrangian structure.
     */
    void
    spreadTransmissionForceDensity(
        int f_data_idx,
        libMesh::PetscVector<double>& X_ghost_vec,
        double time,
        unsigned int part);

    /*!
     * \brief Impose jump conditions determined from the interior and
     * transmission force densities along the physical boundary of the
     * Lagrangian structure.
     */
    void
    imposeJumpConditions(
        int f_data_idx,
        libMesh::PetscVector<double>& F_ghost_vec,
        libMesh::PetscVector<double>& X_ghost_vec,
        double time,
        unsigned int part);

    /*!
     * \brief Initialize the physical coordinates using the supplied coordinate
     * mapping function.  If no function is provided, the initial coordinates
     * are taken to be the Lagrangian coordinates.
     */
    void
    initializeCoordinates(
        unsigned int part);

    /*!
     * \brief Compute dX = X - s, useful mainly for visualization purposes.
     */
    void
    updateCoordinateMapping(
        unsigned int part);

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;
    bool d_is_initialized;

    /*
     * The current time step interval.
     */
    double d_current_time, d_new_time, d_half_time;

    /*
     * FE data associated with this object.
     */
    std::vector<libMesh::Mesh*> d_meshes;
    std::vector<libMesh::EquationSystems*> d_equation_systems;

    const unsigned int d_num_parts;
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;
    std::vector<libMesh::System*> d_X_systems, d_U_systems, d_U_b_systems, d_F_systems;
    std::vector<libMesh::PetscVector<double>*> d_X_current_vecs, d_X_new_vecs, d_X_half_vecs, d_X_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_current_vecs, d_U_new_vecs, d_U_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_b_current_vecs;
    std::vector<libMesh::PetscVector<double>*> d_F_half_vecs, d_F_IB_ghost_vecs;

    /*
     * Method paramters.
     */
    bool d_use_IB_spread_operator;
    std::string d_spread_delta_fcn;
    bool d_use_IB_interp_operator;
    std::string d_interp_delta_fcn;
    bool d_split_forces;
    bool d_use_jump_conditions;
    bool d_use_consistent_mass_matrix;
    libMeshEnums::FEFamily d_fe_family;
    libMeshEnums::Order d_fe_order;
    libMeshEnums::QuadratureType d_quad_type;
    libMeshEnums::Order d_quad_order;

    /*
     * Data related to handling constrained body constraints.
     */
    bool d_has_constrained_parts;
    std::vector<bool> d_constrained_part;
    std::vector<ConstrainedPartVelocityFcnPtr> d_constrained_part_velocity_fcns;
    std::vector<void*> d_constrained_part_velocity_fcn_ctxs;
    double d_constraint_omega;

    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnPtr> d_coordinate_mapping_fcns;
    std::vector<void*> d_coordinate_mapping_fcn_ctxs;

    /*
     * Functions used to compute the first Piola-Kirchhoff stress tensor.
     */
    std::vector<PK1StressFcnPtr> d_PK1_stress_fcns;
    std::vector<std::vector<unsigned int> > d_PK1_stress_fcn_systems;
    std::vector<void*> d_PK1_stress_fcn_ctxs;

    /*
     * Functions used to compute additional body and surface forces on the
     * Lagrangian mesh.
     */
    std::vector<LagBodyForceFcnPtr> d_lag_body_force_fcns;
    std::vector<std::vector<unsigned int> > d_lag_body_force_fcn_systems;
    std::vector<void*> d_lag_body_force_fcn_ctxs;

    std::vector<LagPressureFcnPtr> d_lag_pressure_fcns;
    std::vector<std::vector<unsigned int> > d_lag_pressure_fcn_systems;
    std::vector<void*> d_lag_pressure_fcn_ctxs;

    std::vector<LagSurfaceForceFcnPtr> d_lag_surface_force_fcns;
    std::vector<std::vector<unsigned int> > d_lag_surface_force_fcn_systems;
    std::vector<void*> d_lag_surface_force_fcn_ctxs;

    /*
     * Nonuniform load balancing data structures.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;
    int d_workload_idx;

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
    IBFEMethod();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEMethod(
        const IBFEMethod& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEMethod&
    operator=(
        const IBFEMethod& that);

    /*!
     * Implementation of class constructor.
     */
    void
    commonConstructor(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::vector<libMesh::Mesh*>& meshes,
        int max_level_number,
        bool register_for_restart);

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
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBFEMethod.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBFEMethod
