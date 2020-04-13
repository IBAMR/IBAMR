// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_FEMechanicsExplicitIntegrator
#define included_IBAMR_FEMechanicsExplicitIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/ibamr_enums.h"

#include "ibtk/FEDataManager.h"
#include "ibtk/FEProjector.h"
#include "ibtk/LibMeshSystemVectors.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "tbox/Serializable.h"

#include "libmesh/coupling_matrix.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/explicit_system.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class FEMechanicsExplicitIntegrator is an explicit FE elastodynamics
 * solver.
 */
class FEMechanicsExplicitIntegrator : SAMRAI::tbox::Serializable
{
public:
    static const std::string COORDS_SYSTEM_NAME;
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string FORCE_SYSTEM_NAME;
    static const std::string PRESSURE_SYSTEM_NAME;
    static const std::string VELOCITY_SYSTEM_NAME;

    /// Constructor.
    FEMechanicsExplicitIntegrator(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                  libMesh::MeshBase* mesh,
                                  bool register_for_restart = true,
                                  const std::string& restart_read_dirname = "",
                                  unsigned int restart_restore_number = 0);

        /// Constructor.
    FEMechanicsExplicitIntegrator(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                  const std::vector<libMesh::MeshBase*>& meshes,
                                  bool register_for_restart = true,
                                  const std::string& restart_read_dirname = "",
                                  unsigned int restart_restore_number = 0);

    /// Deleted default constructor.
    FEMechanicsExplicitIntegrator() = delete;

    /// Deleted copy constructor.
    FEMechanicsExplicitIntegrator(const FEMechanicsExplicitIntegrator& from) = delete;

    /// Deleted assignment operator.
    FEMechanicsExplicitIntegrator& operator=(const FEMechanicsExplicitIntegrator& that) = delete;

     /// Destructor.
    ~FEMechanicsExplicitIntegrator() override;

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
     * \note If no function is provided, the initial physical coordinates are
     * taken to be the same as the Lagrangian coordinate system, i.e., the
     * initial coordinate mapping is assumed to be the identity mapping.
     */
    void registerInitialCoordinateMappingFunction(const CoordinateMappingFcnData& data, unsigned int part = 0);

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
     * \note If no function is provided, the initial velocity is taken to be
     * zero.
     */
    void registerInitialVelocityFunction(const InitialVelocityFcnData& data, unsigned int part = 0);

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
                         libMesh::QuadratureType quad_type = libMesh::INVALID_Q_RULE,
                         libMesh::Order quad_order = libMesh::INVALID_ORDER)
            : fcn(fcn), system_data(std::move(system_data)), ctx(ctx), quad_type(std::move(quad_type)), quad_order(std::move(quad_order))
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
     * \note It is possible to register multiple PK1 stress functions with this
     * class.  This is intended to be used to implement selective reduced
     * integration.
     */
    void registerPK1StressFunction(const PK1StressFcnData& data, unsigned int part = 0);

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
     * \note It is \em NOT possible to register multiple body force functions
     * with this class.
     */
    void registerLagBodyForceFunction(const LagBodyForceFcnData& data, unsigned int part = 0);

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
     * \note It is \em NOT possible to register multiple pressure functions with
     * this class.
     */
    void registerLagSurfacePressureFunction(const LagSurfacePressureFcnData& data, unsigned int part = 0);

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
     * \note It is \em NOT possible to register multiple surface force functions
     * with this class.
     */
    void registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, unsigned int part = 0);

    /*!
     * Get the Lagrangian surface force function data.
     */
    LagSurfaceForceFcnData getLagSurfaceForceFunction(unsigned int part = 0) const;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time);

    /*!
     * Initialize the FE equation systems objects.  This method must be called
     * prior to calling initializeFEData().
     */
    void initializeFEEquationSystems();

    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    void initializeFEData();

    /*!
     * Reinitialize FE data by calling `reinit` on each part's EquationSystem,
     * reassembling the system matrices, and setting boundary conditions.
     */
    void reinitializeFEData();

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * For technical reasons this class does not use SAMRAI's RestartManager, so
     * restart files must be separately written for the IBFE objects. This function
     * saves the solutions to the defined EquationSystems in an xdr file in
     * restart_dump_dirname for each FE part. An example snippet is included below to show
     * the distinct IBFE restart data saving step. The data will then be automatically
     * read back into the system along with the RestartManager data during restart.
     *
     * @code
     * if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
     * {
     *     RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
     *     ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
     * }
     * @endcode
     */
    void writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number);

protected:
    /*!
     * \brief Assemble the RHS for the structural force density.
     */
    void assembleInteriorForceDensityRHS(libMesh::PetscVector<double>& F_rhs_vec,
                                         libMesh::PetscVector<double>& X_vec,
                                         libMesh::PetscVector<double>& U_vec,
                                         libMesh::PetscVector<double>& P_vec,
                                         double data_time,
                                         unsigned int part);

    /*!
     * \brief Initialize the physical coordinates using the supplied coordinate
     * mapping function.  If no function is provided, the initial coordinates
     * are taken to be the Lagrangian coordinates.
     */
    void initializeCoordinates(unsigned int part);

    /*!
     * \brief Compute dX = x - X, useful mainly for visualization purposes.
     */
    void updateCoordinateMapping(unsigned int part);

    /*!
     * \brief Initialize the velocity field using the supplied initial velocity
     * specification function.  If no function is provided, the initial
     * velocity is taken to be zero.
     */
    void initializeVelocity(unsigned int part);

    /*!
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log = false;

    /*!
     * Whether or not we have started time integration. This is only used to
     * determine whether or not we print some initial logging output: see
     * d_skip_initial_workload_log for more information.
     */
    bool d_started_time_integration = false;

    /*
     * The current time step interval.
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_half_time = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Meshes provided to this object. These are set up and managed outside
     * this class. These meshes are modified by FEMechanicsExplicitIntegrator since this class
     * creates several libMesh Systems (and hence stores DoF information in
     * these meshes).
     */
    std::vector<libMesh::MeshBase*> d_meshes;

    /*!
     * The libMesh Systems set up by this system (for example, for velocity
     * projection) consist of one variable per spatial component. By default,
     * libMesh assumes that all variables in a given System couple to
     * eachother which, since we only ever solve projection problems in this
     * class, is not the case. Hence we can save some memory by explicitly
     * informing libMesh that the variables in a system only couple to
     * themselves by providing a diagonal coupling matrix to each System.
     */
    libMesh::CouplingMatrix d_diagonal_system_coupling;

    /*!
     * EquationSystems objects, one per part. These contain the actual
     * matrices and solution vectors for each relevant libMesh system.
     */
    std::vector<std::unique_ptr<libMesh::EquationSystems> > d_equation_systems;

    /// Number of parts owned by the present object.
    const unsigned int d_num_parts = 1;

    /// (Uniform) mass density in the reference configuration.
    std::vector<double> d_rho0;

    /// FE data.
    std::vector<std::shared_ptr<IBTK::FEData> > d_fe_data;

    /// FE projection functionality.
    std::vector<std::shared_ptr<IBTK::FEProjector> > d_fe_projector;

    /// Vectors of pointers to the systems for each part (for position, velocity, force
    /// density, sources, and body stress normalization).
    std::vector<libMesh::ExplicitSystem*> d_X_systems, d_U_systems, d_F_systems, d_P_systems;

    /*!
     * Object managing access to libMesh system vectors for the structure position.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_X_vecs;

    /*!
     * Object managing access to libMesh system vectors for the structure velocity.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_U_vecs;

    /*!
     * Object managing access to libMesh system vectors for the structure force.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_F_vecs;

    /*!
     * Object managing access to libMesh system vectors for the structure pressure.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_P_vecs;

    /*!
     * Whether or not the libMesh equation systems objects have been
     * initialized (i.e., whether or not initializeFEEquationSystems has been
     * called).
     */
    bool d_fe_equation_systems_initialized = false;

    /*!
     * Whether or not all finite element data (including that initialized by
     * initializeFEEquationSystems), such system matrices, is available.
     */
    bool d_fe_data_initialized = false;

    /*!
     * Type of partitioner to use. See the main documentation of this class
     * for more information.
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
    std::vector<libMesh::FEFamily> d_fe_family;
    std::vector<libMesh::Order> d_fe_order;
    std::vector<libMesh::QuadratureType> d_default_quad_type;
    std::vector<libMesh::Order> d_default_quad_order;
    bool d_use_consistent_mass_matrix = true;

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
    int d_libmesh_restart_restore_number;

    /*!
     * Restart file type for libMesh equation systems (e.g. xda or xdr).
     */
    std::string d_libmesh_restart_file_extension;

private:
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const std::string& object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           const std::vector<libMesh::MeshBase*>& meshes,
                           bool register_for_restart,
                           const std::string& restart_read_dirname,
                           unsigned int restart_restore_number);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Do the actual work in reinitializeFEData and initializeFEData.
     */
    void doInitializeFEData(const bool use_present_data);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_FEMechanicsExplicitIntegrator
