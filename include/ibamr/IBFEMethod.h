// Filename: IBFEMethod.h
// Created on 5 Oct 2011 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBAMR_IBFEMethod
#define included_IBAMR_IBFEMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <array>
#include <set>
#include <string>
#include <vector>

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBStrategy.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/explicit_system.h"
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
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBFEMethod is an implementation of the abstract base class
 * IBStrategy that provides functionality required by the IB method with finite
 * element elasticity.
 */
class IBFEMethod : public IBStrategy
{
public:
    static const std::string COORDS_SYSTEM_NAME;
    static const std::string COORD_MAPPING_SYSTEM_NAME;
    static const std::string FORCE_SYSTEM_NAME;
    static const std::string PHI_SYSTEM_NAME;
    static const std::string SOURCE_SYSTEM_NAME;
    static const std::string VELOCITY_SYSTEM_NAME;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > mask_var;
    int mask_current_idx, mask_new_idx, mask_scratch_idx;

    /*!
     * \brief Constructor.
     */
    IBFEMethod(const std::string& object_name,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
               libMesh::MeshBase* mesh,
               int max_level_number,
               bool register_for_restart = true,
               const std::string& restart_read_dirname = "",
               unsigned int restart_restore_number = 0);

    /*!
     * \brief Constructor.
     */
    IBFEMethod(const std::string& object_name,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
               const std::vector<libMesh::MeshBase*>& meshes,
               int max_level_number,
               bool register_for_restart = true,
               const std::string& restart_read_dirname = "",
               unsigned int restart_restore_number = 0);

    /*!
     * \brief Destructor.
     */
    ~IBFEMethod();

    /*!
     * Return a pointer to the finite element data manager object for the
     * specified part.
     */
    IBTK::FEDataManager* getFEDataManager(unsigned int part = 0) const;

    /*!
     * Indicate that a part should use stress normalization.
     */
    void registerStressNormalizationPart(unsigned int part = 0);

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
                         const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                         void* const ctx = nullptr,
                         const libMesh::QuadratureType& quad_type = libMesh::INVALID_Q_RULE,
                         const libMesh::Order& quad_order = libMesh::INVALID_ORDER)
            : fcn(fcn), system_data(system_data), ctx(ctx), quad_type(quad_type), quad_order(quad_order)
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
                            const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                            void* const ctx = nullptr)
            : fcn(fcn), system_data(system_data), ctx(ctx)
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
     * Get the Lagrangian surface force function data.
     */
    LagSurfaceForceFcnData getLagSurfaceForceFunction(unsigned int part = 0) const;

    /*!
     * Typedef specifying interface for Lagrangian mass source/sink distribution
     * function.
     */
    using LagBodySourceFcnPtr = IBTK::ScalarMeshFcnPtr;

    /*!
     * Struct encapsulating Lagrangian mass source/sink distribution data.
     */
    struct LagBodySourceFcnData
    {
        LagBodySourceFcnData(LagBodySourceFcnPtr fcn = nullptr,
                             const std::vector<IBTK::SystemData>& system_data = std::vector<IBTK::SystemData>(),
                             void* const ctx = nullptr)
            : fcn(fcn), system_data(system_data), ctx(ctx)
        {
        }

        LagBodySourceFcnPtr fcn;
        std::vector<IBTK::SystemData> system_data;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute a mass source/sink
     * distribution on the Lagrangian finite element mesh.
     */
    void registerLagBodySourceFunction(const LagBodySourceFcnData& data, unsigned int part = 0);

    /*!
     * Get the Lagrangian body source function data.
     */
    LagBodySourceFcnData getLagBodySourceFunction(unsigned int part = 0) const;

    /*!
     * Use tether forces to constrain the motion of a pair of parts.
     */
    void constrainPartOverlap(unsigned int part1,
                              unsigned int part2,
                              double kappa,
                              libMesh::QBase* qrule1 = nullptr,
                              libMesh::QBase* qrule2 = nullptr);

    /*!
     * Always reset the velocity of the nodes of part1 that overlap part2 to
     * equal the velocity of part2.
     */
    void registerOverlappingVelocityReset(unsigned int part1, unsigned int part2);

    /*!
     * Set up forces to penalize relative motion between part1 and part2.
     */
    void registerOverlappingForceConstraint(unsigned int part1,
                                            unsigned int part2,
                                            double kappa,
                                            libMesh::QBase* qrule1 = nullptr,
                                            libMesh::QBase* qrule2 = nullptr);

    /*!
     * Register the (optional) direct forcing kinematics object with the finite
     * element mesh.
     */
    void registerDirectForcingKinematics(SAMRAI::tbox::Pointer<IBAMR::IBFEDirectForcingKinematics> data,
                                         unsigned int part = 0);

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
     * Indicate whether there are any internal fluid sources/sinks.
     */
    bool hasFluidSources() const override;

    /*!
     * Compute the Lagrangian source/sink density at the specified time within
     * the current time interval.
     */
    void computeLagrangianFluidSource(double data_time) override;

    /*!
     * Spread the Lagrangian source/sink density to the Cartesian grid at the
     * specified time within the current time interval.
     */
    void spreadFluidSource(
        int q_data_idx,
        IBTK::RobinPhysBdryPatchStrategy* q_phys_bdry_op,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& q_prolongation_scheds,
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
     * Set the workload spec object used with a particular mesh part.
     */
    void setWorkloadSpec(const IBTK::FEDataManager::WorkloadSpec& workload_spec, unsigned int part = 0);

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
     * Reinitialize FE data by calling `reinit` on each part's EquationSystem,
     * reassembling the system matrices, and setting boundary conditions.
     */
    void reinitializeFEData();

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
     * This function should be called before IBFEMethod::initializeFEData().
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
     * \brief Compute the stress normalization field Phi.
     */
    void computeStressNormalization(libMesh::PetscVector<double>& Phi_vec,
                                    libMesh::PetscVector<double>& X_vec,
                                    double data_time,
                                    unsigned int part);

    /*!
     * \brief Assemble the RHS for the interior elastic density, possibly
     * splitting off the normal component of the transmission force along the
     * physical boundary of the Lagrangian structure.
     */
    void assembleInteriorForceDensityRHS(libMesh::PetscVector<double>& G_rhs_vec,
                                         libMesh::PetscVector<double>& X_vec,
                                         libMesh::PetscVector<double>* Phi_vec,
                                         double data_time,
                                         unsigned int part);

    /*!
     * \brief Reset positions in overlap regions.
     */
    void resetOverlapNodalValues(const std::string& system_name,
                                 const std::vector<libMesh::NumericVector<double>*>& F_vecs);

    /*!
     * \brief Reset positions in overlap regions.
     */
    void resetOverlapNodalValues(const std::string& system_name,
                                 const std::vector<libMesh::PetscVector<double>*>& F_vecs);

    /*!
     * \brief Reset values in the specified part to equal the
     * corresponding velocities of the overlapping part.
     */
    void resetOverlapNodalValues(unsigned int part,
                                 const std::string& system_name,
                                 libMesh::NumericVector<double>* F_vec,
                                 libMesh::NumericVector<double>* F_master_vec);

    /*!
     * \brief Compute constraint forces between pairs of overlapping bodies.
     */
    void computeOverlapConstraintForceDensity(std::vector<libMesh::PetscVector<double>*>& G_vec,
                                              std::vector<libMesh::PetscVector<double>*>& X_vec);

    /*!
     * \brief Spread the transmission force density along the physical boundary
     * of the Lagrangian structure.
     */
    void spreadTransmissionForceDensity(int f_data_idx,
                                        libMesh::PetscVector<double>& X_ghost_vec,
                                        IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                                        double data_time,
                                        unsigned int part);

    /*!
     * \brief Impose jump conditions determined from the interior and
     * transmission force densities along the physical boundary of the
     * Lagrangian structure.
     */
    void imposeJumpConditions(int f_data_idx,
                              libMesh::PetscVector<double>& F_ghost_vec,
                              libMesh::PetscVector<double>& X_ghost_vec,
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
     * The current time step interval.
     */
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_half_time = std::numeric_limits<double>::quiet_NaN();

    /*
     * FE data associated with this object.
     */
    std::vector<libMesh::MeshBase*> d_meshes;
    int d_max_level_number;
    std::vector<std::unique_ptr<libMesh::EquationSystems>> d_equation_systems;

    const unsigned int d_num_parts = 1;
    std::vector<IBTK::FEDataManager*> d_fe_data_managers;
    SAMRAI::hier::IntVector<NDIM> d_ghosts = 0;
    std::vector<libMesh::ExplicitSystem*> d_X_systems, d_U_systems, d_F_systems, d_Q_systems, d_Phi_systems;
    std::vector<libMesh::PetscVector<double>*> d_X_current_vecs, d_X_rhs_vecs, d_X_new_vecs, d_X_half_vecs,
        d_X_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_current_vecs, d_U_rhs_vecs, d_U_new_vecs, d_U_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_F_half_vecs, d_F_rhs_vecs, d_F_tmp_vecs, d_F_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_Q_half_vecs, d_Q_rhs_vecs, d_Q_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_Phi_half_vecs, d_Phi_rhs_vecs;

    bool d_fe_equation_systems_initialized = false, d_fe_data_initialized = false;

    /*
     * Method paramters.
     */
    IBTK::FEDataManager::InterpSpec d_default_interp_spec;
    IBTK::FEDataManager::SpreadSpec d_default_spread_spec;
    std::vector<IBTK::FEDataManager::WorkloadSpec> d_workload_spec;
    std::vector<IBTK::FEDataManager::InterpSpec> d_interp_spec;
    std::vector<IBTK::FEDataManager::SpreadSpec> d_spread_spec;
    bool d_split_normal_force = false, d_split_tangential_force = false;
    bool d_use_jump_conditions = false;
    std::vector<libMesh::FEFamily> d_fe_family;
    std::vector<libMesh::Order> d_fe_order;
    std::vector<libMesh::QuadratureType> d_default_quad_type;
    std::vector<libMesh::Order> d_default_quad_order;
    bool d_use_consistent_mass_matrix = true;

    /*
     * Data related to handling stress normalization.
     */
    double d_epsilon = 0.0;
    bool d_has_stress_normalization_parts = false;
    std::vector<bool> d_is_stress_normalization_part;

    /*
     * Data related to constraining overlaps between pairs of parts.
     */
    double d_overlap_tolerance = 0.0;

    bool d_has_overlap_velocity_parts = false;
    std::vector<bool> d_is_overlap_velocity_part, d_is_overlap_velocity_master_part;
    std::vector<int> d_overlap_velocity_master_part;
    std::vector<std::map<libMesh::dof_id_type, libMesh::dof_id_type> > d_overlap_velocity_part_node_to_elem_map;
    std::vector<std::set<libMesh::dof_id_type> > d_overlap_velocity_part_ghost_idxs;

    bool d_has_overlap_force_parts = false;
    std::vector<bool> d_is_overlap_force_part;
    std::vector<std::set<libMesh::dof_id_type> > d_overlap_force_part_ghost_idxs;
    std::vector<std::array<unsigned int, 2> > d_overlap_force_part_idxs;
    std::vector<std::array<std::map<libMesh::dof_id_type, std::map<unsigned int, libMesh::dof_id_type> >, 2> >
        d_overlapping_elem_map;
    std::vector<double> d_overlap_force_part_kappa;
    std::vector<std::array<libMesh::QBase*, 2> > d_overlap_force_part_qrule; // \todo let's try to fix this when we switch to C++11!
    std::vector<std::vector<double> > d_overlap_force_part_max_displacement;

    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnData> d_coordinate_mapping_fcn_data;

    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<InitialVelocityFcnData> d_initial_velocity_fcn_data;

    /*
     * Functions used to compute the first Piola-Kirchhoff stress tensor.
     */
    std::vector<std::vector<PK1StressFcnData> > d_PK1_stress_fcn_data;

    /*
     * Objects used to impose direct forcing kinematics.
     */
    std::vector<SAMRAI::tbox::Pointer<IBAMR::IBFEDirectForcingKinematics> > d_direct_forcing_kinematics_data;

    /*
     * Functions used to compute additional body and surface forces on the
     * Lagrangian mesh.
     */
    std::vector<LagBodyForceFcnData> d_lag_body_force_fcn_data;
    std::vector<LagSurfacePressureFcnData> d_lag_surface_pressure_fcn_data;
    std::vector<LagSurfaceForceFcnData> d_lag_surface_force_fcn_data;

    /*
     * Functions used to compute source/sink strength on the Lagrangian mesh.
     */
    bool d_has_lag_body_source_parts = false;
    std::vector<bool> d_lag_body_source_part;
    std::vector<LagBodySourceFcnData> d_lag_body_source_fcn_data;

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
    std::string d_libmesh_restart_file_extension;

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
    IBFEMethod(const IBFEMethod& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEMethod& operator=(const IBFEMethod& that) = delete;

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

    /*!
     * Do the actual work in reinitializeFEData and initializeFEData. if @p
     * use_present_data is `true` then the current content of the solution
     * vectors is used: more exactly, the coordinates and velocities (computed
     * by initializeCoordinates and initializeVelocity) are considered as
     * being up to date, as is the direct forcing kinematic data.
     */
    void doInitializeFEData(const bool use_present_data);
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBFEMethod
