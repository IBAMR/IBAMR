// Filename: IBFEMethod.h
// Created on 5 Oct 2011 by Boyce Griffith
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

#ifndef included_IBFEMethod
#define included_IBFEMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stdbool.h>
#include <stddef.h>
#include <set>
#include <string>
#include <vector>

#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "ibamr/IBStrategy.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/libmesh_utilities.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
class CCLaplaceOperator;
class PETScKrylovPoissonSolver;
class CCPoissonPointRelaxationFACOperator;
class FACPreconditioner;
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
namespace solv
{
class PoissonSpecifications;
template <int DIM>
class LocationIndexRobinBcCoefs;
} // namespace solv
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
    static const std::string VELOCITY_SYSTEM_NAME;
    static const std::string BODY_VELOCITY_SYSTEM_NAME;

    /*!
     * \brief Constructor.
     */
    IBFEMethod(const std::string& object_name,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
               libMesh::Mesh* mesh,
               int max_level_number,
               bool register_for_restart = true);

    /*!
     * \brief Constructor.
     */
    IBFEMethod(const std::string& object_name,
               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
               const std::vector<libMesh::Mesh*>& meshes,
               int max_level_number,
               bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    ~IBFEMethod();
	
	/*!
	 * \brief Register Eulerian variables for velocity correction step.
	 */
	void registerEulerianVariables();

    /*!
     * Return a pointer to the finite element data manager object for the
     * specified part.
     */
    IBTK::FEDataManager* getFEDataManager(unsigned int part = 0) const;

    /*!
     * Indicate that a part is constrained.
     */
    void registerConstrainedPart(unsigned int part = 0);

    /*!
     * Typedef specifying interface for specifying constrained body velocities.
     */
    typedef void (*ConstrainedVelocityFcnPtr)(libMesh::NumericVector<double>& U_b,
											  libMesh::NumericVector<double>& U,
											  libMesh::NumericVector<double>& X,
											  const Eigen::Vector3d& X_com,
											  Eigen::Vector3d& U_com,
											  Eigen::Vector3d& W_com,
											  libMesh::EquationSystems* equation_systems,
											  double data_time,
											  void* ctx);
    /*!
     * Struct encapsulating constrained velocity function data.
     */
    struct ConstrainedVelocityFcnData
	{
        ConstrainedVelocityFcnData(ConstrainedVelocityFcnPtr fcn = NULL,
								   void* ctx = NULL) : fcn(fcn), ctx(ctx)
        {
        }

		ConstrainedVelocityFcnPtr fcn;
        void* ctx;
    };

    /*!
     * Register a constrained body velocity function.
     */
	void registerConstrainedVelocityFunction(ConstrainedVelocityFcnPtr fcn,
											 void* ctx = NULL,
											 unsigned int part = 0,
											 const Eigen::Vector3d& u_com_initial = Eigen::Vector3d::Zero(),
											 const Eigen::Vector3d& w_com_initial = Eigen::Vector3d::Zero());

    /*!
     * Register a constrained body velocity function data.
     */
    void registerConstrainedVelocityFunction(const ConstrainedVelocityFcnData& data,
											 unsigned int part = 0,
											 const Eigen::Vector3d& u_com_initial = Eigen::Vector3d::Zero(),
											 const Eigen::Vector3d& w_com_initial = Eigen::Vector3d::Zero());
	
	/*!
	 * Typedef specifying interface for specifying constrained position update.
	 */
	typedef void (*ConstrainedPositionFcnPtr)(libMesh::NumericVector<double>& X_new,
											  libMesh::NumericVector<double>& X_current,
											  const Eigen::Vector3d& X_com,
											  const Eigen::Vector3d& u_com,
											  const Eigen::Vector3d& w_com,
											  libMesh::EquationSystems* equation_systems,
											  double data_time,
											  double dt,
	                                          void* ctx);
	/*!
	 * Struct encapsulating constrained position function data.
	 */
	struct ConstrainedPositionFcnData
	{
		ConstrainedPositionFcnData(ConstrainedPositionFcnPtr fcn = NULL,
								   void* ctx = NULL) : fcn(fcn), ctx(ctx)
		{
		}
		
		ConstrainedPositionFcnPtr fcn;
		void* ctx;
	};
	
	/*!
	 * Register a constrained position update function.
	 */
	void registerConstrainedPositionFunction(ConstrainedPositionFcnPtr fcn,
											 void* ctx = NULL,
											 unsigned int part = 0);
	
	/*!
	 * Register a constrained position update function data.
	 */
	void registerConstrainedPositionFunction(const ConstrainedPositionFcnData& data,
											 unsigned int part = 0);
	
	/*!
     * Typedef specifying interface for coordinate mapping function.
     */
    typedef void (*CoordinateMappingFcnPtr)(libMesh::Point& X, const libMesh::Point& s, void* ctx);

    /*!
     * Struct encapsulating coordinate mapping function data.
     */
    struct CoordinateMappingFcnData
    {
        CoordinateMappingFcnData(CoordinateMappingFcnPtr fcn = NULL, void* ctx = NULL) : fcn(fcn), ctx(ctx)
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
    void registerInitialCoordinateMappingFunction(CoordinateMappingFcnPtr fcn, void* ctx = NULL, unsigned int part = 0);

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
     * Typedef specifying interface for PK1 stress tensor function.
     */
    typedef IBTK::TensorMeshFcnPtr PK1StressFcnPtr;

    /*!
     * Struct encapsulating PK1 stress tensor function data.
     */
    struct PK1StressFcnData
    {
        PK1StressFcnData(PK1StressFcnPtr fcn = NULL,
                         const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                         void* ctx = NULL,
                         libMesh::QuadratureType quad_type = libMesh::INVALID_Q_RULE,
                         libMesh::Order quad_order = libMesh::INVALID_ORDER)
            : fcn(fcn), systems(systems), ctx(ctx), quad_type(quad_type), quad_order(quad_order)
        {
        }

        PK1StressFcnPtr fcn;
        std::vector<unsigned int> systems;
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
    void registerPK1StressFunction(PK1StressFcnPtr fcn,
                                   const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                                   void* ctx = NULL,
                                   libMesh::QuadratureType quad_type = libMesh::INVALID_Q_RULE,
                                   libMesh::Order quad_order = libMesh::INVALID_ORDER,
                                   unsigned int part = 0);

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
     * Typedef specifying interface for Lagrangian body force distribution
     * function.
     */
    typedef IBTK::VectorMeshFcnPtr LagBodyForceFcnPtr;

    /*!
     * Struct encapsulating Lagrangian body force distribution data.
     */
    struct LagBodyForceFcnData
    {
        LagBodyForceFcnData(LagBodyForceFcnPtr fcn = NULL,
                            const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                            void* ctx = NULL)
            : fcn(fcn), systems(systems), ctx(ctx)
        {
        }

        LagBodyForceFcnPtr fcn;
        std::vector<unsigned int> systems;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple body force functions
     * with this class.
     */
    void registerLagBodyForceFunction(LagBodyForceFcnPtr fcn,
                                      const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                                      void* ctx = NULL,
                                      unsigned int part = 0);

    /*!
     * Register the (optional) function to compute body force distributions on
     * the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple body force functions
     * with this class.
     */
    void registerLagBodyForceFunction(const LagBodyForceFcnData& data, unsigned int part = 0);

    /*!
     * Typedef specifying interface for Lagrangian pressure force distribution
     * function.
     */
    typedef IBTK::ScalarSurfaceFcnPtr LagSurfacePressureFcnPtr;

    /*!
     * Struct encapsulating Lagrangian surface pressure distribution data.
     */
    struct LagSurfacePressureFcnData
    {
        LagSurfacePressureFcnData(LagSurfacePressureFcnPtr fcn = NULL,
                                  const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                                  void* ctx = NULL)
            : fcn(fcn), systems(systems), ctx(ctx)
        {
        }

        LagSurfacePressureFcnPtr fcn;
        std::vector<unsigned int> systems;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute surface pressure
     * distributions on the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple pressure functions with
     * this class.
     */
    void registerLagSurfacePressureFunction(LagSurfacePressureFcnPtr fcn,
                                            const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                                            void* ctx = NULL,
                                            unsigned int part = 0);

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
    typedef IBTK::VectorSurfaceFcnPtr LagSurfaceForceFcnPtr;

    /*!
     * Struct encapsulating Lagrangian surface force distribution data.
     */
    struct LagSurfaceForceFcnData
    {
        LagSurfaceForceFcnData(LagSurfaceForceFcnPtr fcn = NULL,
                               const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                               void* ctx = NULL)
            : fcn(fcn), systems(systems), ctx(ctx)
        {
        }

        LagSurfaceForceFcnPtr fcn;
        std::vector<unsigned int> systems;
        void* ctx;
    };

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple surface force functions
     * with this class.
     */
    void registerLagSurfaceForceFunction(LagSurfaceForceFcnPtr fcn,
                                         const std::vector<unsigned int>& systems = std::vector<unsigned int>(),
                                         void* ctx = NULL,
                                         unsigned int part = 0);

    /*!
     * Register the (optional) function to compute surface force distributions
     * on the Lagrangian finite element mesh.
     *
     * \note It is \em NOT possible to register multiple surface force functions
     * with this class.
     */
    void registerLagSurfaceForceFunction(const LagSurfaceForceFcnData& data, unsigned int part = 0);

    /*!
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>& getMinimumGhostCellWidth() const;

    /*!
     * Setup the tag buffer.
     */
    void setupTagBuffer(SAMRAI::tbox::Array<int>& tag_buffer,
                        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) const;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void preprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateData(double current_time, double new_time, int num_cycles);

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time);

    /*!
     * Advance the positions of the Lagrangian structure using the forward Euler
     * method.
     */
    void eulerStep(double current_time, double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * midpoint rule.
     */
    void midpointStep(double current_time, double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the (explicit)
     * trapezoidal rule.
     */
    void trapezoidalStep(double current_time, double new_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void computeLagrangianForce(double data_time);

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(int f_data_idx,
                IBTK::RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
                double data_time);
	
	void
	postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num);

    /*!
     * Initialize FE data.  This method must be called prior to calling
     * IBHierarchyIntegrator::initializePatchHierarchy().
     */
    void initializeFEData();

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
        bool initial_time);

    /*!
     * Register a load balancer and work load patch data index with the IB
     * strategy object.
     */
    void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
                              int workload_data_idx);

    /*!
     * Update work load estimates on each level of the patch hierarchy.
     */
    void updateWorkloadEstimates(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 int workload_data_idx);

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    void beginDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    void endDataRedistribution(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                               SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

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
                             bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_level,
                                     int finest_level);

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
                               bool uses_richardson_extrapolation_too);

    /*!
     * Write out object state to the given database.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
	
	/*!
	 * \brief Fill the rotation matrix.
	 * \param rot_mat Matrix to set.
	 * \param dt Time interval of rotation.
	 */
	static void setRotationMatrix(const Eigen::Vector3d& rot_vel,
						          Eigen::Matrix3d& R,
						          const double dt);

protected:
    /*
     * \brief Compute the constraint force density.
     */
    void computeConstraintForceDensity(libMesh::PetscVector<double>& F_vec,
                                       libMesh::PetscVector<double>& X_vec,
                                       libMesh::PetscVector<double>& U_vec,
                                       libMesh::PetscVector<double>& U_b_vec,
                                       double data_time,
                                       unsigned int part);

    /*
     * \brief Compute the interior elastic density, possibly splitting off the
     * normal component of the transmission force along the physical boundary of
     * the Lagrangian structure.
     */
    void computeInteriorForceDensity(libMesh::PetscVector<double>& G_vec,
                                     libMesh::PetscVector<double>& X_vec,
                                     double data_time,
                                     unsigned int part);

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
     * \brief Compute dX = X - s, useful mainly for visualization purposes.
     */
    void updateCoordinateMapping(unsigned int part);
	
	/*!
	 * \brief Compute center of mass and moment of inertia of the structure.
	 */
	void computeCOMandMOI(const unsigned part,
						  Eigen::Vector3d& center_of_mass,
						  Eigen::Matrix3d& moment_of_inertia,
						  libMesh::PetscVector<double>* X);
	
	/*!
	 * \brief Copy fluid variable from solver to a widened variable
	 * and fill the ghost cells.
	 */
	void copyFluidVariable(const int from_idx,
						   const int to_idx);
	
	/*!
	 * \brief Build projection solver for imposing divergence free constraint on
	 * velocity. This is needed for certain cases of imposing rigidity constraint.
	 */
	void buildProjectionSolver();
	
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
    std::vector<libMesh::System*> d_X_systems, d_U_systems, d_F_systems, d_U_b_systems;
    std::vector<libMesh::PetscVector<double>*> d_X_current_vecs, d_X_new_vecs, d_X_half_vecs, d_X_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_current_vecs, d_U_new_vecs, d_U_half_vecs;
    std::vector<libMesh::PetscVector<double>*> d_F_half_vecs, d_F_IB_ghost_vecs;
    std::vector<libMesh::PetscVector<double>*> d_U_b_current_vecs, d_U_b_new_vecs, d_U_b_half_vecs;

    bool d_fe_data_initialized;

    /*
     * Method paramters.
     */
    bool d_use_IB_interp_operator;
    IBTK::FEDataManager::InterpSpec d_interp_spec;
    bool d_use_IB_spread_operator;
    IBTK::FEDataManager::SpreadSpec d_spread_spec;
    bool d_split_forces;
    bool d_use_jump_conditions;
    libMesh::FEFamily d_fe_family;
    libMesh::Order d_fe_order;
    libMesh::QuadratureType d_quad_type;
    libMesh::Order d_quad_order;
    bool d_use_consistent_mass_matrix;

    /*
     * Data related to handling constrained body constraints.
     */
    bool d_has_constrained_parts;
    std::vector<bool> d_constrained_part;
    std::vector<ConstrainedVelocityFcnData> d_constrained_velocity_fcn_data;
	std::vector<ConstrainedPositionFcnData> d_constrained_position_fcn_data;
	
	
	/*!
	 * Rigid body velocity of the structures.
	 */
	std::vector<Eigen::Vector3d> d_com_u_current, d_com_u_half, d_com_u_new;
	std::vector<Eigen::Vector3d> d_com_w_current, d_com_w_half, d_com_w_new;
	
	/*!
	 * Center of mass and moment of inertia.
	 */
	std::vector<Eigen::Vector3d> d_com_current, d_com_half;
	std::vector<Eigen::Matrix3d> d_moi_current, d_moi_half;
	
	/*!
     * Eulerian data for fluid velocity and momentum correction.
	 */
	SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_ins_var;
	int d_u_ins_idx, d_u_ins_cib_idx;
	
    /*
     * Functions used to compute the initial coordinates of the Lagrangian mesh.
     */
    std::vector<CoordinateMappingFcnData> d_coordinate_mapping_fcn_data;

    /*
     * Functions used to compute the first Piola-Kirchhoff stress tensor.
     */
    std::vector<std::vector<PK1StressFcnData> > d_PK1_stress_fcn_data;

    /*
     * Functions used to compute additional body and surface forces on the
     * Lagrangian mesh.
     */
    std::vector<LagBodyForceFcnData> d_lag_body_force_fcn_data;
    std::vector<LagSurfacePressureFcnData> d_lag_surface_pressure_fcn_data;
    std::vector<LagSurfaceForceFcnData> d_lag_surface_force_fcn_data;

    /*
     * Collection of all systems required to evaluate various quantities.
     */
    std::vector<std::set<unsigned int> > d_fcn_systems, d_body_fcn_systems, d_surface_fcn_systems;

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
	
	/*!
	 * Data structures for applying divergence free projection.
	 */
	bool d_impose_div_free_constraint;
	SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* d_proj_bc_coef;
	SAMRAI::solv::PoissonSpecifications* d_proj_spec;
	SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> d_proj_solver_op;
	SAMRAI::tbox::Pointer<IBTK::PETScKrylovPoissonSolver> d_proj_solver;
	SAMRAI::tbox::Pointer<IBTK::CCPoissonPointRelaxationFACOperator> d_proj_pc_op;
	SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_proj_solver_db, d_proj_pc_db;
	SAMRAI::tbox::Pointer<IBTK::FACPreconditioner> d_proj_pc;
	SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_proj_var;
	int d_u_proj_idx, d_phi_proj_idx, d_div_u_proj_idx;
	
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBFEMethod(const IBFEMethod& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBFEMethod& operator=(const IBFEMethod& that);

    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const std::string& object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           const std::vector<libMesh::Mesh*>& meshes,
                           int max_level_number,
                           bool register_for_restart);
	
	/*!
	 * Apply divergence free projection.
	 */
	void applyProjection();
	
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

#endif //#ifndef included_IBFEMethod
