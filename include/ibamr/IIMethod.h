// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBAMR_IIMethod
#define included_IBAMR_IIMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBStrategy.h"

#include "ibtk/FEDataManager.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "SideIndex.h"
#include "tbox/Pointer.h"

#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/vector_value.h"

#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
class SAMRAIDataCache;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
template <int DIM>
class BasePatchLevel;
} // namespace hier
namespace tbox
{
class Database;
template <class TYPE>
class Array;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI
namespace libMesh
{
class EquationSystems;
class Mesh;
class Point;
class System;
template <typename T>
class NumericVector;
template <typename T>
class PetscVector;
class MeshBase;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IIMethod is an implementation of the abstract base
 * class IBStrategy that provides functionality required by the IB method with
 * a finite element representation of a surface mesh.
 *
 * Coupling schemes include both IB formulations (integral operations with
 * regularized delta function kernels) and an immersed interface method (IIM)
 * scheme (E. M. Kolahdouz, A. P. S. Bhalla, B. A. Craven, and B. E. Griffith.
 * An immersed interface method for discrete surfaces. J Comput Phys,
 * 400:108854 (37 pages), 2020).
 *
 * \note When using the IIM implementation, it is recommended that users set
 * all linear solvers to use tight relative tolerances (1e-10).
 */
class IIMethod : public IBStrategy
{
public:
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string COORDS_SYSTEM_NAME;
    static const std::string FORCE_SYSTEM_NAME;
    static const std::string NORMAL_VELOCITY_SYSTEM_NAME;
    static const std::string PRESSURE_IN_SYSTEM_NAME;
    static const std::string PRESSURE_JUMP_SYSTEM_NAME;
    static const std::string PRESSURE_OUT_SYSTEM_NAME;
    static const std::string TANGENTIAL_VELOCITY_SYSTEM_NAME;
    static const std::string TAU_IN_SYSTEM_NAME;
    static const std::string TAU_OUT_SYSTEM_NAME;
    static const std::string VELOCITY_SYSTEM_NAME;
    static const std::string WSS_IN_SYSTEM_NAME;
    static const std::string WSS_OUT_SYSTEM_NAME;
    static const std::array<std::string, NDIM> VELOCITY_JUMP_SYSTEM_NAME;

    /*!
     * \brief Constructor.
     */
    IIMethod(const std::string& object_name,
             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
             libMesh::MeshBase* mesh,
             int max_level_number,
             bool register_for_restart = true,
             const std::string& restart_read_dirname = "",
             unsigned int restart_restore_number = 0);

    /*!
     * \brief Constructor.
     */
    IIMethod(const std::string& object_name,
             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
             const std::vector<libMesh::MeshBase*>& meshes,
             int max_level_number,
             bool register_for_restart = true,
             const std::string& restart_read_dirname = "",
             unsigned int restart_restore_number = 0);

    /*!
     * \brief Destructor.
     */
    ~IIMethod();

    /*!
     * Return a pointer to the finite element data manager object for the
     * specified part.
     */
    IBTK::FEDataManager* getFEDataManager(unsigned int part = 0) const;

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
     * Register relevant part to use discontinuous element type family (L2_LAGRANGE or MONOMIAL)
     * for the calculation of jump/traction quantities. This option should be used for geometries with
     * sharp corners.
     *
     *
     * \note The relveant FE family is provided via input file. The default value is set to regular LAGRANGE).
     */
    void registerDisconElemFamilyForJumps(unsigned int part = 0);

    /*!
     * Register relevant part to only use the tangential component of the Lagrangian velocity
     * for displacement. The lack of normal velocity needs to be compensated by a large
     * frictional penalty force (direct forcing for the normal velocity)
     *
     * \note If no function is provided, all parts will be moving according to the full Lagrangian velocity
     */
    void registerTangentialVelocityMotion(unsigned int part = 0);

    /*!
     * Register relevant part for which we would like to normalize the pressure jump.
     *
     * \note This is only recommended for enclosed bodies.
     */
    void registerPressureJumpNormalization(unsigned int part = 0);

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
                                  const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                                  void* const ctx = nullptr)
            : fcn(fcn), system_data(system_data), ctx(ctx)
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
                               const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                               void* const ctx = nullptr)
            : fcn(fcn), system_data(system_data), ctx(ctx)
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
     * The current value of the integrated surface force.
     */
    const libMesh::VectorValue<double>& getSurfaceForceIntegral(unsigned int part = 0) const;

    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>& getMinimumGhostCellWidth() const override;

    /*!
     * Setup the tag buffer.
     */
    void setupTagBuffer(SAMRAI::tbox::Array<int>& tag_buffer,
                        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) const override;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time) override;

    /*!
     * Compute the fluid traction if all jump conditions are applied.
     */
    void computeFluidTraction(double current_time, unsigned int part = 0);

    /*!
     * Compute the interior/exterior pressure used in the calculation of the fluid traction (pressure jump condition
     * needs to be applied).
     */
    void extrapolatePressureForTraction(int p_data_idx, double data_time, unsigned int part = 0);

    /*!
     * A wrapper to compute the interfacial pressure and fluid traction from the jumps.
     */
    void calculateInterfacialFluidForces(int p_data_idx, double data_time);

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void forwardEulerStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time) override;

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time) override;

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time) override;

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time) override;

    /*!
     * Get the default interpolation spec object used by the class.
     */
    IBTK::FEDataManager::InterpSpec getDefaultInterpSpec() const;

    /*!
     * Get the default spread spec object used by the class.
     */
    IBTK::FEDataManager::SpreadSpec getDefaultSpreadSpec() const;

    /*!
     * Set the interpolation spec object used with a particular mesh part.
     */
    void setInterpSpec(const IBTK::FEDataManager::InterpSpec& interp_spec, unsigned int part = 0);

    /*!
     * Set the spread spec object used with a particular mesh part.
     */
    void setSpreadSpec(const IBTK::FEDataManager::SpreadSpec& spread_spec, unsigned int part = 0);

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
     * \brief Register Eulerian variables with the parent IBHierarchyIntegrator.
     */
    void registerEulerianVariables() override;

    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time) override;

    /*!
     * Register a load balancer and work load patch data index with the IB
     * strategy object.
     *
     * @deprecated This method is no longer necessary with the current
     * workload estimation scheme.
     */
    void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
                              int workload_data_idx) override;

    /*!
     * Add the estimated computational work from the current object (i.e., the
     * work required by the owned Lagrangian objects) per cell into the
     * specified <code>workload_data_idx</code>.
     */
    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const int workload_data_idx) override;

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    void beginDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    void endDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                             bool allocate_data) override;

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_level,
                                     int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too) override;

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Write the equation_systems data to a restart file in the specified directory.
     */
    void writeFEDataToRestartFile(const std::string& restart_dump_dirname, unsigned int time_step_number);

protected:
    /*!
     * Impose the jump conditions.
     */
    void imposeJumpConditions(const int f_data_idx,
                              libMesh::PetscVector<double>& P_jump_ghost_vec,
                              std::array<libMesh::PetscVector<double>*, NDIM>& DU_jump_ghost_vec,
                              libMesh::PetscVector<double>& X_ghost_vec,
                              const double data_time,
                              const unsigned int part);

    /*!
     * \brief Helper function for checking possible double-counting
     *  intesection points
     */
    bool checkDoubleCountingIntersection(int axis,
                                         const double* dx,
                                         const libMesh::VectorValue<double>& n,
                                         const libMesh::Point& x,
                                         const libMesh::Point& xi,
                                         const SAMRAI::pdat::SideIndex<NDIM>& i_s,
                                         const SAMRAI::pdat::SideIndex<NDIM>& i_s_prime,
                                         const std::vector<libMesh::Point>& candidate_coords,
                                         const std::vector<libMesh::Point>& candidate_ref_coords,
                                         const std::vector<libMesh::VectorValue<double> >& candidate_normals);

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

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log = false;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;
    bool d_is_initialized = false;

    /*
     * Scratch data caching objects.
     */
    std::shared_ptr<IBTK::SAMRAIDataCache> d_eulerian_data_cache;

    /*
     * The current time step interval.
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_half_time = std::numeric_limits<double>::quiet_NaN();

    /*
     * FE data associated with this object.
     * d_X: coordinates system
     * d_F: IB force system
     * d_U: velocity system
     * d_U_n: normal velocity system
     * d_U_t: tangential velocity system
     * d_P_jump: pressure jump system
     * d_DU_jump: velocity gradient jump system
     * d_WSS_in: one sided interior shear stress system
     * d_WSS_out: one sided exterior shear stress system
     * d_P_in: one sided interior pressure system
     * d_P_out: one sided exterior pressure system
     * d_TAU_in: one sided interior fluid traction system
     * d_TAU_out: one sided exterior fluid traction system
     */
    std::vector<libMesh::MeshBase*> d_meshes;
    int d_max_level_number;
    std::vector<libMesh::EquationSystems*> d_equation_systems;
    std::vector<bool> d_use_discon_elem_for_jumps = { false };
    std::vector<bool> d_use_tangential_velocity = { false };
    std::vector<bool> d_normalize_pressure_jump = { false };
    const unsigned int d_num_parts = 1;
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;
    SAMRAI::hier::IntVector<NDIM> d_ghosts = 0;
    std::vector<libMesh::System*> d_X_systems, d_U_systems, d_U_n_systems, d_U_t_systems, d_F_systems, d_P_jump_systems,
        d_WSS_in_systems, d_WSS_out_systems, d_P_in_systems, d_P_out_systems, d_TAU_in_systems, d_TAU_out_systems;
    std::vector<std::array<libMesh::System*, NDIM> > d_DU_jump_systems;
    std::vector<libMesh::PetscVector<double>*> d_F_half_vecs, d_F_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_X_current_vecs, d_X_new_vecs, d_X_half_vecs, d_X0_vecs,
        d_X_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_current_vecs, d_U_new_vecs, d_U_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_n_current_vecs, d_U_n_new_vecs, d_U_n_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_t_current_vecs, d_U_t_new_vecs, d_U_t_half_vecs;
    std::vector<std::array<libMesh::PetscVector<double>*, NDIM> > d_DU_jump_half_vecs, d_DU_jump_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_P_jump_half_vecs, d_P_jump_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_P_in_half_vecs, d_P_in_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_P_out_half_vecs, d_P_out_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_WSS_in_half_vecs, d_WSS_in_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_WSS_out_half_vecs, d_WSS_out_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_TAU_in_half_vecs, d_TAU_in_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_TAU_out_half_vecs, d_TAU_out_IB_ghost_vecs;

    bool d_fe_equation_systems_initialized = false, d_fe_data_initialized = false;

    /*
     * Method parameters.
     */
    IBTK::FEDataManager::InterpSpec d_default_interp_spec;
    IBTK::FEDataManager::SpreadSpec d_default_spread_spec;
    IBTK::FEDataManager::WorkloadSpec d_default_workload_spec;
    std::vector<IBTK::FEDataManager::InterpSpec> d_interp_spec;
    std::vector<IBTK::FEDataManager::SpreadSpec> d_spread_spec;
    bool d_use_pressure_jump_conditions = false;
    libMesh::FEFamily d_pressure_jump_fe_family = libMesh::LAGRANGE;
    bool d_use_velocity_jump_conditions = false;
    libMesh::FEFamily d_velocity_jump_fe_family = libMesh::LAGRANGE;
    bool d_compute_fluid_traction = false;
    libMesh::FEFamily d_wss_fe_family = libMesh::LAGRANGE;
    libMesh::FEFamily d_tau_fe_family = libMesh::LAGRANGE;
    bool d_perturb_fe_mesh_nodes = true;
    std::vector<libMesh::FEFamily> d_fe_family;
    std::vector<libMesh::Order> d_fe_order;
    std::vector<libMesh::QuadratureType> d_default_quad_type;
    std::vector<libMesh::Order> d_default_quad_order;
    bool d_use_consistent_mass_matrix = true;
    bool d_use_direct_forcing = false;
    double d_exterior_calc_coef = 1.0;
    double d_wss_calc_width = 1.05;
    double d_p_calc_width = 1.3;

    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnData> d_coordinate_mapping_fcn_data;

    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<InitialVelocityFcnData> d_initial_velocity_fcn_data;

    /*
     * Functions used to compute surface forces on the Lagrangian mesh.
     */
    std::vector<LagSurfacePressureFcnData> d_lag_surface_pressure_fcn_data;
    std::vector<LagSurfaceForceFcnData> d_lag_surface_force_fcn_data;
    std::vector<libMesh::VectorValue<double> > d_lag_surface_force_integral;

    /*
     * Eulerian data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_p_var;
    int d_p_scratch_idx = IBTK::invalid_index;

    /*
     * Nonuniform load balancing data structures.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;
    int d_workload_idx = IBTK::invalid_index;

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
     * Directory and time step number to use when restarting.
     */
    std::string d_libmesh_restart_read_dir;
    int d_libmesh_restart_restore_number;

    /*
     * Restart file type for libMesh equation systems (e.g. xda or xdr).
     */
    std::string d_libmesh_restart_file_extension = "xdr";

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IIMethod() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IIMethod(const IIMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IIMethod& operator=(const IIMethod& that) = delete;

    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const std::string& object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           const std::vector<libMesh::MeshBase*>& meshes,
                           int max_level_number,
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
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IIMethod
