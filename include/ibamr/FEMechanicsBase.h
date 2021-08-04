// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
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

#ifndef included_IBAMR_FEMechanicsBase
#define included_IBAMR_FEMechanicsBase

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibamr/ibamr_enums.h"

#include "ibtk/FEDataManager.h"
#include "ibtk/LibMeshSystemVectors.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/Serializable.h"

#include "libmesh/coupling_matrix.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"

#include <string>
#include <utility>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * @brief      Class FEMechanicsBase provides core finite element mechanics
 *             functionality and data management.
 *
 * <h2>Parameters read from the input database</h2>
 * <ol>
 *   <li>FEProjector: Input database passed along to the object responsible for
 *     computing projections onto the finite element space. See the
 *     documentation of IBTK::FEProjector for more information.</li>
 * </ol>
 */
class FEMechanicsBase : public SAMRAI::tbox::Serializable
{
public:
    // TODO: Update the names of these systems:
    // COORDS_SYSTEM_NAME --> CURRENT_COORDINATES_SYSTEM_NAME
    // COORD_MAPPING_SYSTEM_NAME --> DISPLACEMENT_SYSTEM_NAME
    static const std::string COORDS_SYSTEM_NAME;
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string FORCE_SYSTEM_NAME;
    static const std::string PRESSURE_SYSTEM_NAME;
    static const std::string VELOCITY_SYSTEM_NAME;

    /*!
     * Constructor for a single-part model.
     */
    FEMechanicsBase(const std::string& object_name,
                    const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                    libMesh::MeshBase* mesh,
                    bool register_for_restart = true,
                    const std::string& restart_read_dirname = "",
                    unsigned int restart_restore_number = 0);

    /*!
     * Constructor for a multi-part model.
     */
    FEMechanicsBase(const std::string& object_name,
                    const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                    const std::vector<libMesh::MeshBase*>& meshes,
                    bool register_for_restart = true,
                    const std::string& restart_read_dirname = "",
                    unsigned int restart_restore_number = 0);

    /*!
     * Deleted default constructor.
     */
    FEMechanicsBase() = delete;

    /*!
     * Deleted copy constructor.
     */
    FEMechanicsBase(const FEMechanicsBase& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    FEMechanicsBase& operator=(const FEMechanicsBase& that) = delete;

    /*!
     * Destructor.
     */
    virtual ~FEMechanicsBase() override;

    /*!
     * Return a pointer to the equations system object for the specified part.
     */
    libMesh::EquationSystems* getEquationSystems(unsigned int part = 0) const;

    /*!
     * Return a pointer to the FEData object for the specified part.
     */
    std::shared_ptr<IBTK::FEData> getFEData(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for coordinate mapping function.
     */
    using CoordinateMappingFcnPtr = void (*)(libMesh::Point& x, const libMesh::Point& X, void* ctx);

    /*!
     * Struct encapsulating coordinate mapping function data.
     */
    struct CoordinateMappingFcnData
    {
        CoordinateMappingFcnData(CoordinateMappingFcnPtr fcn = nullptr, void* ctx = nullptr) : fcn(fcn), ctx(ctx)
        {
        }

        CoordinateMappingFcnPtr fcn;
        void* ctx;
    };

    /*!
     * Register the (optional) function used to initialize the physical
     * coordinates from the Lagrangian coordinates.
     *
     * @note       If no function is provided, the initial physical coordinates
     *             are taken to be the same as the Lagrangian coordinate system,
     *             i.e., the initial coordinate mapping is assumed to be the
     *             identity mapping.
     */
    virtual void registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, unsigned int part = 0);

    /*!
     * Get the initial coordinate mapping function data.
     */
    CoordinateMappingFcnData getInitialCoordinateMappingFunction(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for initial velocity specification function.
     */
    using InitialVelocityFcnPtr = void (*)(libMesh::VectorValue<double>& U0, const libMesh::Point& X0, void* ctx);

    /*!
     * Struct encapsulating initial velocity specification function data.
     */
    struct InitialVelocityFcnData
    {
        InitialVelocityFcnData(InitialVelocityFcnPtr fcn = nullptr, void* ctx = nullptr) : fcn(fcn), ctx(ctx)
        {
        }

        InitialVelocityFcnPtr fcn;
        void* ctx;
    };

    /*!
     * Register the (optional) function used to initialize the velocity of the
     * solid mesh.
     *
     * @note       If no function is provided, the initial velocity is taken to
     *             be zero.
     */
    virtual void registerInitialVelocityFunction(const InitialVelocityFcnData& data, unsigned int part = 0);

    /*!
     * Get the initial velocity function data.
     */
    InitialVelocityFcnData getInitialVelocityFunction(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for PK1 stress tensor function.
     */
    using PK1StressFcnPtr = IBTK::TensorMeshFcnPtr;

    /*!
     * Struct encapsulating PK1 stress tensor function data.
     */
    struct PK1StressFcnData
    {
        PK1StressFcnData(PK1StressFcnPtr fcn = nullptr,
                         std::vector<IBTK::SystemData> system_data = {},
                         void* const ctx = nullptr,
                         const libMesh::QuadratureType& quad_type = libMesh::INVALID_Q_RULE,
                         const libMesh::Order& quad_order = libMesh::INVALID_ORDER)
            : fcn(fcn), system_data(std::move(system_data)), ctx(ctx), quad_type(quad_type), quad_order(quad_order)
        {
        }

        PK1StressFcnPtr fcn;
        std::vector<IBTK::SystemData> system_data;
        void* ctx;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
    };

    /*!
     * Register the (optional) function to compute the first Piola-Kirchhoff
     * stress tensor, used to compute the forces on the Lagrangian finite
     * element mesh.
     *
     * @note       It is possible to register multiple PK1 stress functions with
     *             this class.  This is intended to be used to implement
     *             selective reduced integration.
     */
    virtual void registerPK1StressFunction(const PK1StressFcnData& data, unsigned int part = 0);

    /*!
     * Get the PK1 stress function data.
     */
    std::vector<PK1StressFcnData> getPK1StressFunction(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for Lagrangian body force distribution
     * function.
     */
    using LagBodyForceFcnPtr = IBTK::VectorMeshFcnPtr;

    /*!
     * Struct encapsulating Lagrangian body force distribution data.
     */
    struct LagBodyForceFcnData
    {
        LagBodyForceFcnData(LagBodyForceFcnPtr fcn = nullptr,
                            std::vector<IBTK::SystemData> system_data = {},
                            void* const ctx = nullptr)
            : fcn(fcn), system_data(std::move(system_data)), ctx(ctx)
        {
        }

        LagBodyForceFcnPtr fcn;
        std::vector<IBTK::SystemData> system_data;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     *
     * @note       It is @em NOT possible to register multiple body force
     *             functions with this class.
     */
    virtual void registerLagBodyForceFunction(const LagBodyForceFcnData& data, unsigned int part = 0);

    /*!
     * Get the Lagrangian body force function data.
     */
    LagBodyForceFcnData getLagBodyForceFunction(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for Lagrangian pressure force distribution
     * function.
     */
    using LagSurfacePressureFcnPtr = IBTK::ScalarSurfaceFcnPtr;

    /*!
     * Struct encapsulating Lagrangian surface pressure distribution data.
     */
    struct LagSurfacePressureFcnData
    {
        LagSurfacePressureFcnData(LagSurfacePressureFcnPtr fcn = nullptr,
                                  std::vector<IBTK::SystemData> system_data = {},
                                  void* const ctx = nullptr)
            : fcn(fcn), system_data(std::move(system_data)), ctx(ctx)
        {
        }

        LagSurfacePressureFcnPtr fcn;
        std::vector<IBTK::SystemData> system_data;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute surface pressure
     * distributions on the Lagrangian finite element mesh.
     *
     * @note       It is @em NOT possible to register multiple pressure
     *             functions with this class.
     */
    virtual void registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, unsigned int part = 0);

    /*!
     * Get the Lagrangian surface pressure function data.
     */
    LagSurfacePressureFcnData getLagSurfacePressureFunction(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for Lagrangian surface force distribution
     * function.
     */
    using LagSurfaceForceFcnPtr = IBTK::VectorSurfaceFcnPtr;

    /*!
     * Struct encapsulating Lagrangian surface force distribution data.
     */
    struct LagSurfaceForceFcnData
    {
        LagSurfaceForceFcnData(LagSurfaceForceFcnPtr fcn = nullptr,
                               std::vector<IBTK::SystemData> system_data = {},
                               void* const ctx = nullptr)
            : fcn(fcn), system_data(std::move(system_data)), ctx(ctx)
        {
        }

        LagSurfaceForceFcnPtr fcn;
        std::vector<IBTK::SystemData> system_data;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     *
     * @note       It is @em NOT possible to register multiple surface force
     *             functions with this class.
     */
    virtual void registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, unsigned int part = 0);

    /*!
     * Get the Lagrangian surface force function data.
     */
    LagSurfaceForceFcnData getLagSurfaceForceFunction(unsigned int part = 0) const;

    /*!
     * Function signature for specifying the energy functional that determines
     * the pressure.
     */
    using VolumetricEnergyDerivativeFcn = double (*)(double);

    /*!
     * Indicate that a part should include a static pressure.
     *
     * The pressure is determined via (P, Q) = (U'(J), Q), using either a
     * consistent or lumped mass matrix, or via a locally stabilized projection
     * of the form (P, Q) + epsilon (P - Pi P, Q - Pi Q) = (U'(J), Q), in which
     * P is the pressure and Q is an arbitrary test function.
     *
     * Users can provide a function to evaluate U'(J).  If no function is
     * provided, we default to using U(J) = -kappa (J ln(J) − J + 1), so that
     * U'(J) = -kappa ln J. (Ref: C.H. Liu, G. Hofstetter, H.A. Mang, 3D finite
     * element analysis of rubber-like materials at finite strains, Eng. Comput.
     * 11 (2) (1994) 111–128.)
     *
     * The sign convention used in the implementation generates a PK1 stress of
     * the form PP = -J P FF^{-T}.
     *
     * @note       The same part cannot have both static and dynamic pressures.
     *
     * @see        registerDynamicPressurePart
     */
    virtual void registerStaticPressurePart(PressureProjectionType projection_type = CONSISTENT_PROJECTION,
                                            VolumetricEnergyDerivativeFcn dU_dJ_fcn = nullptr,
                                            unsigned int part = 0);

    /*!
     * Indicate that a part should include a dynamic pressure.
     *
     * The pressure is determined via (P-dot, Q) = (J U''(J) (FF : Grad U), Q),
     * using either a consistent or lumped mass matrix, or via a locally
     * stabilized projection of the form (P, Q) + epsilon (P - Pi P, Q - Pi Q) =
     * (U'(J), Q), in which P is the pressure and Q is an arbitrary test
     * function.
     *
     * Users can provide a function to evaluate U''(J).  If no function is
     * provided, we default to using U(J) = -kappa (J ln(J) − J + 1), so that
     * U'(J) = -kappa ln J and U''(J) = -kappa J^{-1}. (Ref: C.H. Liu, G.
     * Hofstetter, H.A. Mang, 3D finite element analysis of rubber-like
     * materials at finite strains, Eng. Comput. 11 (2) (1994) 111–128.)
     *
     * The sign convention used in the implementation generates a PK1 stress of
     * the form PP = -J P FF^{-T}.
     *
     * @note       The same part cannot have both static and dynamic pressures.
     *
     * @see        registerStaticPressurePart
     */
    virtual void registerDynamicPressurePart(PressureProjectionType projection_type = CONSISTENT_PROJECTION,
                                             VolumetricEnergyDerivativeFcn d2U_dJ2_fcn = nullptr,
                                             unsigned int part = 0);

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    virtual void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    virtual void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Initialize the FE equation systems objects.  This method must be called
     * prior to calling initializeFEData().
     */
    virtual void initializeFEEquationSystems();

    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    virtual void initializeFEData();

    /*!
     * Reinitialize FE data by calling `reinit` on each part's EquationSystem,
     * reassembling the system matrices, and setting boundary conditions.
     */
    virtual void reinitializeFEData();

    /*!
     * Write out object state to the given database.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * For technical reasons this class does not use SAMRAI's RestartManager, so
     * restart files must be separately written for the FE objects. This
     * function saves the solutions to the defined EquationSystems in an xdr
     * file in restart_dump_dirname for each FE part. An example snippet is
     * included below to show the distinct FE restart data saving step. The data
     * will then be automatically read back into the system along with the
     * RestartManager data during restart.
     *
     * @code
     * if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
     * {
     *     RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
     *     fe_mechanics_base->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
     * }
     * @endcode
     */
    virtual void writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number);

protected:
    /*!
     * Do the actual work in initializeFEEquationSystems.
     */
    virtual void doInitializeFEEquationSystems();

    /*!
     * Do the actual work of setting up libMesh system vectors.
     */
    virtual void doInitializeFESystemVectors();

    /*!
     * Do the actual work in reinitializeFEData and initializeFEData. if @p
     * use_present_data is `true` then the current content of the solution
     * vectors is used: more exactly, the coordinates and velocities (computed
     * by initializeCoordinates and initializeVelocity) are considered as being
     * up to date.
     */
    virtual void doInitializeFEData(bool use_present_data);

    /*!
     * \brief Compute the static pressure field.
     */
    void computeStaticPressure(libMesh::PetscVector<double>& P_vec,
                               libMesh::PetscVector<double>& X_vec,
                               double data_time,
                               unsigned int part);

    /*!
     * \brief Compute the dynamic pressure time rate of change.
     */
    void computeDynamicPressureRateOfChange(libMesh::PetscVector<double>& dP_dt_vec,
                                            libMesh::PetscVector<double>& X_vec,
                                            libMesh::PetscVector<double>& U_vec,
                                            double data_time,
                                            unsigned int part);

    /*!
     * Assemble the RHS for the interior elastic density, possibly splitting off
     * the normal component of the transmission force along the physical
     * boundary of the Lagrangian structure.
     */
    virtual void assembleInteriorForceDensityRHS(libMesh::PetscVector<double>& F_rhs_vec,
                                                 libMesh::PetscVector<double>& X_vec,
                                                 libMesh::PetscVector<double>* P_vec,
                                                 double data_time,
                                                 unsigned int part);

    /*!
     * Initialize the physical coordinates using the supplied coordinate mapping
     * function.  If no function is provided, the initial coordinates are taken
     * to be the Lagrangian coordinates.
     */
    virtual void initializeCoordinates(unsigned int part);

    /*!
     * Compute dX = x - X, useful mainly for visualization purposes.
     */
    virtual void updateCoordinateMapping(unsigned int part);

    /*!
     * Initialize the velocity field using the supplied initial velocity
     * specification function.  If no function is provided, the initial velocity
     * is taken to be zero.
     */
    virtual void initializeVelocity(unsigned int part);

    /*!
     * Convenience function to setup system vectors and, if necessary, convert
     * PARALLEL vectors into GHOSTED vectors for a collection of Systems.
     *
     * @deprecated use IBTK::setup_system_vectors instead.
     */
    static void setup_system_vectors(libMesh::EquationSystems* equation_systems,
                                     const std::vector<std::string>& system_names,
                                     const std::vector<std::string>& vector_names);

    /*!
     * Convenience function to setup a system vector and, if necessary, convert
     * a PARALLEL vector into a GHOSTED vector.
     *
     * @deprecated use IBTK::setup_system_vector instead.
     */
    static void setup_system_vector(libMesh::System& system, const std::string& vector_name);

    /*!
     * Get the libMesh restart file name.
     */
    static std::string libmesh_restart_file_name(const std::string& restart_dump_dirname,
                                                 unsigned int time_step_number,
                                                 unsigned int part,
                                                 const std::string& extension);

    /*!
     * Cached input databases.
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_fe_projector_db;

    /*!
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log = false;

    /*
     * The current time step interval.
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_half_time = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Meshes provided to this object. These are set up and managed outside this
     * class. These meshes are modified by FEMechanicsBase since this class
     * creates several libMesh Systems (and hence stores DoF information in
     * these meshes).
     */
    std::vector<libMesh::MeshBase*> d_meshes;

    /*!
     * The libMesh Systems set up by this system (for example, for velocity
     * projection) consist of one variable per spatial component. By default,
     * libMesh assumes that all variables in a given System couple to each other
     * which, since we only ever solve projection problems in this class, is not
     * the case. Hence we can save some memory by explicitly informing libMesh
     * that the variables in a system only couple to themselves by providing a
     * diagonal coupling matrix to each System.
     */
    libMesh::CouplingMatrix d_diagonal_system_coupling;

    /*!
     * EquationSystems objects, one per part. These contain the actual matrices
     * and solution vectors for each relevant libMesh system.
     */
    std::vector<std::unique_ptr<libMesh::EquationSystems> > d_equation_systems;

    /// FEData objects provide key FE data management.
    std::vector<std::shared_ptr<IBTK::FEData> > d_fe_data;

    /// FEProjector objects provide L2 projection functionality.
    std::vector<std::shared_ptr<IBTK::FEProjector> > d_fe_projectors;

    /// Vectors of pointers to the systems for each part (for position,
    /// velocity, force, and pressure).
    std::vector<libMesh::ExplicitSystem*> d_X_systems, d_U_systems, d_F_systems, d_P_systems;

    /*!
     * Object managing access to libMesh system vectors for the position.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_X_vecs;

    /*!
     * Object managing access to libMesh system vectors for the velocity.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_U_vecs;

    /*!
     * Object managing access to libMesh system vectors for the force.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_F_vecs;

    /*!
     * Object managing access to libMesh system vectors for the pressure.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_P_vecs;

    /*!
     * Whether or not the libMesh equation systems objects have been initialized
     * (i.e., whether or not initializeFEEquationSystems has been called).
     */
    bool d_fe_equation_systems_initialized = false;

    /*!
     * Whether or not all finite element data (including that initialized by
     * initializeFEEquationSystems), such system matrices, is available.
     */
    bool d_fe_data_initialized = false;

    /*!
     * Type of partitioner to use. See the main documentation of this class for
     * more information.
     */
    LibmeshPartitionerType d_libmesh_partitioner_type = LIBMESH_DEFAULT;

    /*!
     * Whether or not to use AMR in the finite element discretization. This
     * feature is not yet implemented and currently defaults to false.
     */
    bool d_libmesh_use_amr = false;

    /*!
     * Method parameters.
     */
    std::vector<libMesh::Order> d_fe_order_position, d_fe_order_force, d_fe_order_pressure;
    std::vector<libMesh::FEFamily> d_fe_family_position, d_fe_family_force, d_fe_family_pressure;
    std::vector<libMesh::QuadratureType> d_default_quad_type_stress, d_default_quad_type_force,
        d_default_quad_type_pressure;
    std::vector<libMesh::Order> d_default_quad_order_stress, d_default_quad_order_force, d_default_quad_order_pressure;
    bool d_use_consistent_mass_matrix = true;
    bool d_allow_rules_with_negative_weights = true;
    bool d_include_normal_stress_in_weak_form = false;
    bool d_include_tangential_stress_in_weak_form = false;
    bool d_include_normal_surface_forces_in_weak_form = true;
    bool d_include_tangential_surface_forces_in_weak_form = true;

    /*!
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnData> d_coordinate_mapping_fcn_data;

    /*!
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<InitialVelocityFcnData> d_initial_velocity_fcn_data;

    /*!
     * Functions used to compute the first Piola-Kirchhoff stress tensor.
     */
    std::vector<std::vector<PK1StressFcnData> > d_PK1_stress_fcn_data;

    /*!
     * Functions used to compute additional body and surface forces on the
     * Lagrangian mesh.
     */
    std::vector<LagBodyForceFcnData> d_lag_body_force_fcn_data;
    std::vector<LagSurfacePressureFcnData> d_lag_surface_pressure_fcn_data;
    std::vector<LagSurfaceForceFcnData> d_lag_surface_force_fcn_data;

    /*!
     * Data related to handling static pressures.
     */
    double d_static_pressure_kappa = 0.0, d_static_pressure_stab_param = 0.0;
    bool d_has_static_pressure_parts = false;
    std::vector<bool> d_static_pressure_part;
    std::vector<PressureProjectionType> d_static_pressure_proj_type;
    std::vector<VolumetricEnergyDerivativeFcn> d_static_pressure_dU_dJ_fcn;

    /*!
     * Data related to handling dynamic pressures.
     */
    double d_dynamic_pressure_kappa = 0.0, d_dynamic_pressure_stab_param = 0.0;
    bool d_has_dynamic_pressure_parts = false;
    std::vector<bool> d_dynamic_pressure_part;
    std::vector<PressureProjectionType> d_dynamic_pressure_proj_type;
    std::vector<VolumetricEnergyDerivativeFcn> d_dynamic_pressure_d2U_dJ2_fcn;

    /*!
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*!
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;

    /*!
     * Directory and time step number to use when restarting.
     */
    std::string d_libmesh_restart_read_dir;
    unsigned int d_libmesh_restart_restore_number;

    /*!
     * Restart file type for libMesh equation systems (e.g. xda or xdr).
     */
    std::string d_libmesh_restart_file_extension;

private:
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const std::string& object_name,
                           const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                           const std::vector<libMesh::MeshBase*>& meshes,
                           bool register_for_restart,
                           const std::string& restart_read_dirname,
                           unsigned int restart_restore_number);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_FEMechanicsBase
