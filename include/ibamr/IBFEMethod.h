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

#ifndef included_IBAMR_IBFEMethod
#define included_IBAMR_IBFEMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibamr/FEMechanicsBase.h"
#include "ibamr/IBFEDirectForcingKinematics.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/FEDataManager.h"
#include "ibtk/LibMeshSystemIBVectors.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/SAMRAIGhostDataAccumulator.h"
#include "ibtk/SecondaryHierarchy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "BoxGeneratorStrategy.h"
#include "GriddingAlgorithm.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "TagAndInitializeStrategy.h"
#include "Variable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "libmesh/coupling_matrix.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/explicit_system.h"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
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
template <int DIM>
class RefinePatchStrategy;
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
class ExplicitSystem;
class MeshBase;
template <typename T>
class VectorValue;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
class IBFEDirectForcingKinematics;

/*!
 * \brief Class IBFEMethod is an implementation of the abstract base class
 * IBStrategy that provides functionality required by the IB method with finite
 * element elasticity. Much of the setup for finite element computations is done
 * by the IBAMR::FEMechanicsBase: see the documentation of that class for
 * additional information on input parameters.
 *
 * By default, the libMesh data is partitioned once at the beginning of the
 * computation by libMesh's default partitioner.
 *
 * <h2>Options Controlling Interpolation and Spreading</h2>
 * Like other classes inheriting from IBStrategy, most options regarding the
 * actual IB method implementation can be specified with the provided input
 * database. Parameters starting with <code>IB_</code> set and override those
 * with the same name starting with <code>interp_</code> or
 * <code>spread_</code>: e.g., <code>IB_delta_fcn</code> overrides both
 * <code>interp_delta_fcn</code> and <code>spread_delta_fcn</code>.
 * <ul>
 *   <li><code>interp_quad_type</code>: Quadrature type for interpolation,
 *   provided as a string. Can be any quadrature type known to libMesh.
 *   Defaults to <code>"QGAUSS"</code>.</li>
 *   <li><code>spread_quad_type</code>: Quadrature type for spreading,
 *   provided as a string. Parsed in the same was as <code>interp_quad_type</code>.</li>
 *   <li><code>IB_quad_type</code>: overriding alias for the two previous
 *   entries - has the same default.</li>
 *   <li><code>interp_use_adaptive_quadrature</code>: Whether or not the current
 *   deformation of each element should be considered when determining which
 *   quadrature rule to use. Defaults to <code>TRUE</code>.</li>
 *   <li><code>interp_point_density</code>: Multiplier on the number of points
 *   to use for doing IB calculations. Defaults to <code>2</code>.
 *   <li><code>spread_point_density</code>: Same as above, but for spreading.
 *   <li><code>IB_point_density</code>: overriding alias for the previous
 *   entry - has the same default.</li>
 *   <li><code>interp_point_density</code>: Parameter for adaptively computing the
 *   number of quadrature points in a quadrature rule. Defaults to
 *   <code>2.0</code>. See IBTK::getQuadratureKey() for a detailed
 *   description.</li>
 *   <li><code>spread_point_density</code>: Same as above, but for spreading.
 *   <li><code>IB_point_density</code>: overriding alias for the two previous
 *   entries - has the same default.</li>
 *   <li><code>interp_use_consistent_mass_matrix</code>: Whether or not mass
 *   lumping should be applied when solving the L2 projection for computing
 *   the velocity of the structure. Defaults to FALSE. Note that no linear
 *   system is solved when computing forces so this parameter does not have a
 *   spreading equivalent.</li>
 *   <li><code>use_consistent_mass_matrix</code>: Overriding alias of the
 *   previous entry.
 *   <li><code>IB_use_consistent_mass_matrix</code>: Overriding alias of
 *   the previous entry.</li>
 *   <li><code>interp_use_nodal_quadrature</code>: Whether or not nodal
 *   quadrature should be used, which is essentially interpolation instead of
 *   projection. This is an experimental feature. Defaults to
 *   <code>FALSE</code>.</li>
 *   <li><code>spread_use_nodal_quadrature</code>: Same as above, but for spreading.
 *   <li><code>IB_use_nodal_quadrature</code>: overriding alias for the two previous
 *   entries - has the same default.</li>
 * </ul>
 *
 * <h2>Options Controlling libMesh Partitioning</h2>
 * <em>This feature is experimental: at the present time the default settings
 * have the best performance and are the correct choice.</em>
 *
 * This class can repartition libMesh data in a way that matches SAMRAI's
 * distribution of patches; put another way, if a certain region of space on
 * the finest level is assigned to processor N, then all libMesh nodes and
 * elements within that region will also be assigned to processor N. The
 * actual partitioning here is done by the IBTK::BoxPartitioner class. See the
 * discussion in IBTK::HierarchyIntegrator and IBTK::FEDataManager for
 * descriptions on how this partitioning is performed.
 *
 * The choice of libMesh partitioner depends on the libmesh_partitioner_type
 * parameter in the input database:
 * <ul>
 *  <li>If <code>libmesh_partitioner_type</code> is
 *      <code>LIBMESH_DEFAULT</code> then this class will never repartition
 *      libMesh data, since the default libMesh partitioner is already used at
 *      the beginning of the computation and, since no degrees of freedom are
 *      added or removed, any subsequent partitioning would have no
 *      effect.</li>
 *
 *  <li>If <code>libmesh_partitioner_type</code> is <code>SAMRAI_BOX</code>
 *      then this class will always repartition the libMesh data with
 *      IBTK::BoxPartitioner every time the Eulerian data is regridded.</li>
 * </ul>
 * The default value for <code>libmesh_partitioner_type</code> is
 * <code>LIBMESH_DEFAULT</code>. The intent of these choices is to
 * automatically use the fairest (that is, partitioning based on equal work
 * when computing force densities and L2 projections) partitioner.
 *
 * <h2>Options Controlling IB Data Partitioning</h2>
 *
 * The main computational expenses of this class are
 * IBFEMethod::interpolateVelocity() and IBFEMethod::spreadForce(). These two
 * methods compute at IB points placed inside the patches owned on the current
 * processor: i.e., they use the Eulerian partitioning of the domain. This
 * partitioning scales very poorly at higher processor counts with some
 * Lagrangian geometries since the Eulerian partitioning places equal number
 * of cells, which do not necessarily coincide with IB points, on different
 * processors: i.e., some processors will have a large number of IB points and
 * some may have zero.
 *
 * To get around this, this class can optionally work with a different
 * partitioning of the Eulerian data that is partitioned so that each
 * processor has roughly the same number of IB points, or some more elaborate
 * partitioning scheme that takes into account the number of mesh nodes as
 * well. This class will set up this scratch hierarchy and manage its state
 * (see IBFEMethod::d_scratch_hierarchy). The scratch hierarchy can be set up
 * by adding the following parameters to the input database:
 *
 * @code
 * IBFEMethod {
 *
 *    // Place whatever database entries you typically use with IBFEMethod here,
 *    // e.g., define IB_delta_fcn, split_forces, use_consistent_mass_matrix etc.
 *    // as usual. The parameters listed below solely pertain to the scratch
 *    // hierarchy.
 *
 *    use_scratch_hierarchy = TRUE
 *    workload_quad_point_weight = 1.0
 *    workload_duplicated_node_weight = 0.8
 *
 *    // The values supplied here should usually be the same as those provided to
 *    // the top-level GriddingAlgorithm.
 *    GriddingAlgorithm
 *    {
 *        max_levels = MAX_LEVELS
 *        ratio_to_coarser
 *        {
 *            level_1 = REF_RATIO,REF_RATIO
 *        }
 *
 *        largest_patch_size
 *        {
 *            // We recommend using very large values here: large patches
 *            // are more efficient, especially with the merging load balancer.
 *            level_0 = 512,512
 *        }
 *
 *        smallest_patch_size
 *        {
 *            // on the other hand, smaller patch sizes here typically enable
 *            // better load balancing at the cost of creating more total work
 *            // due to an increased number of ghost cells (and, therefore,
 *            // an increased number of elements in more than one patch).
 *            // We recommend adjusting this number to be as large as possible
 *            // by examining log output - if higher-rank processors do not have
 *            // enough work then it should be slightly decreased. 8 x 8 is usually
 *            // too small, but 16 x 16 is reasonable for many setups.
 *            level_0 = 16, 16
 *        }
 *
 *        efficiency_tolerance = 0.80e0
 *        combine_efficiency   = 0.80e0
 *        coalesce_boxes = TRUE
 *        allow_patches_smaller_than_minimum_size_to_prevent_overlaps = TRUE
 *    }
 *
 *    // Smaller workload factors improve load balancing but increase the total
 *    // amount of work since more elements will end up on multiple patches.
 *    // This value is a good compromise.
 *    // Similarly, since intraprocessor patch communication is less of a concern
 *    // here than in the fluid solver, we recommend using the greedy load
 *    // balancer bin packing method.
 *    LoadBalancer
 *    {
 *       type                = "MERGING"
 *       bin_pack_method     = "GREEDY"
 *       max_workload_factor = 0.5
 *    }
 * }
 * @endcode
 *
 * i.e., providing <code>use_scratch_hierarchy = TRUE</code> (the default is
 * <code>FALSE</code>) turns on the scratch hierarchy and the remaining
 * parameters determine how patches are generated and load balanced. The extra
 * argument <code>type</code> to <code>LoadBalancer</code> specifies whether
 * an IBTK::MergingLoadBalancer (chosen by <code>"MERGING"</code>) or the
 * default SAMRAI LoadBalancer (chosen by <code>"DEFAULT"</code>) is
 * used. Since IBTK::MergingLoadBalancer is usually what one wants
 * <code>"MERGING"</code> is the default. The merging option is better since
 * it reduces the total number of elements which end up in patch ghost
 * regions since some patches will be merged together.
 *
 * The parameter <code>workload_quad_point_weight</code> is the multiplier
 * assigned to an IB point when calculating the work per processor. Similarly,
 * the parameter <code>workload_duplicated_node_weight</code> is the multiplier
 * assigned to each node of every element (i.e., each node is counted more than
 * once): see IBTK::FEDataManager::WorkloadSpec for more information.
 *
 * For efficiency reasons this class only associates elements with patches if
 * they can interact with the patches (e.g., the points used for velocity
 * interpolation can lie inside the patch). Hence, over time, this association
 * needs to change as the structure moves. The frequency of reassociation is
 * controlled by two things:
 * <ol>
 *   <li>The rate at which the hierarchy is regridded, as if the hierarchy itself
 *   changes then we must recompute the association.</li>
 *   <li>The value d_patch_association_cfl, which determines how frequently we regrid
 *   purely based on the displacement of the structure. This parameter can be set
 *   by providing <code>patch_association_cfl</code> to the input database.
 * </ol>
 *
 * If you set IBHierarchyIntegrator::d_regrid_cfl_interval to a value larger
 * than 1 and use AMR this class may not work correctly since, in that case,
 * elements may move outside the patch level they are associated with before new
 * cells are tagged to be on that level.
 *
 * <h2>Options Controlling Logging</h2>
 * The logging options set by this class are propagated to the owned
 * IBTK::FEDataManager objects.
 * <ol>
 *   <li><code>enable_logging</code>: set to <code>TRUE</code> to enable logging.</code>.
 *   Defaults to <code>false</code>.</li>
 *   <li><code>skip_initial_workload_log</code>: For testing purposes (see
 *   d_skip_initial_workload_log) it is necessary to disable some output: this
 *   option disables logging of workload data (quadrature point counts, etc.)
 *   before the first time step if set to <code>TRUE</code>. Defaults to
 *   <code>false</code>.</li>
 * </ol>
 *
 * <h2>Using IBFEMethod with your own libMesh System objects</h2>
 *
 * It is often useful to add your own libMesh data to the EquationSystems
 * object used by IBAMR objects. One such example would be defining fields of
 * fibers to give direction to add anisotropy to solid models. Since a libMesh
 * Mesh object also stores all the degrees of freedom of all the systems, there
 * can only be one EquationSystems object for each Mesh object. If you set up
 * your own systems in this way then they will be automatically saved and loaded
 * from restart data.
 *
 * Since initialization of this class occurs in multiple stages and IBFEMethod
 * assumes that it is ultimately responsible for setting up the EquationSystems
 * object (by calling <code>EquationSystems::init()</code>) other libMesh Systems
 * must be added in a specific order. Additionally, some parts must only be set
 * up if restart information is not available.
 *
 * This code is based on an IBFE example and assumes that the IBAMR objects
 * objects are already set up in the usual way.
 *
 * @code
 * // Actually create the EquationSystems objects.
 * ib_method_ops->initializeFEEquationSystems();
 *
 * // This code assumes we only have one part, so there is only one
 * // EquationSystems object.
 * libMesh::EquationSystems *equation_systems =
 *     ib_method_ops->getFEDataManager()->getEquationSystems();
 *
 * const bool from_restart = RestartManager::getManager()->isFromRestart();
 * // Set up the external System. It is not necessary to create the System
 * // if restart data is available, since in that case libMesh will create
 * // the system itself from loaded data.
 * if (!from_restart)
 * {
 *     auto& fiber =
 *         equation_systems->add_system<ExplicitSystem>("fiber");
 *     fiber.add_variable("f_0");
 *     fiber.add_variable("f_1");
 *
 *     // do not call init here - the IBFEMethod object will do that itself
 *     // later in the call to initializeFEData().
 * }
 *
 * // do more initialization of IBAMR objects here...
 *
 * // IBFEMethod::initializeFEData() will initialize all System objects in each
 * // EquationSystems object: i.e., libMesh will fail with an assertion check if
 * // you added your own System to that EquationSystems object and initialized it
 * // it before this point. Since this includes all Systems, i.e., also
 * // external systems, any data vectors you added to the System will be
 * // loaded in from the same restart file as the rest of IBAMR's data.
 * ib_method_ops->initializeFEData();
 *
 * // initialize the patch hierarchy, postprocessor, etc. here too, as usual.
 *
 * // Continue configuring the external system. Like before, its not necessary
 * // to add vectors if restart data is available, since in that case libMesh
 * // will already know the vectors exist.
 * if (!from_restart)
 * {
 *     auto& fiber =
 *         equation_systems->get_system<ExplicitSystem>("fiber");
 *     fiber.add_vector("NF");
 *     NumericVector<double> &vec = fiber.get_vector("NF");
 *
 *     vec = 42.0;
 *     vec.close();
 * }
 * // For testing purposes, lets make sure that on a restarted run we have the
 * // correct system and vector available:
 * else
 * {
 *     auto& fiber =
 *         equation_systems->get_system<ExplicitSystem>("fiber");
 *     plog << "number of variables: " << fiber.n_vars() << '\n';
 *     NumericVector<double> &vec = fiber.get_vector("NF");
 *     // This will print the vector to plog: since each entry is 42 we know
 *     // that the vector was loaded correctly.
 *     vec.print(plog);
 * }
 * @endcode
 *
 * <h2>Passing Data to the FEDataManager class</h2>
 * IBFEMethod uses IBTK::FEDataManager to actually perform IB calculations with
 * the finite element mesh. IBTK::FEDataManager objects are configured by
 * setting spreading and interpolation parameters in the usual way (i.e., by
 * providing the parameters described above in the input database). Options
 * specific to the behavior of FEDataManager can be set by defining a database
 * named <code>FEDataManager</code> inside the database provided to this class -
 * see the documentation of FEDataManager for more information.
 *
 * <h2>Handling Restart Data</h2>
 * The caching of the IBFE restart data is not managed by SAMRAI's
 * SAMRAI::tbox::RestartManager. It is instead handled by
 * FEMechanicsBase::writeFEDataToRestartFile() given a restart_dump_dirname and
 * time_step_number. Each instance of IBFEMethod is registered for restart by
 * default, but the this option can be turned off. During a restart, the data
 * is handled by the SAMRAI::tbox::RestartManager automatically to reinitiate
 * the IBFEMethod.
 */
class IBFEMethod : public FEMechanicsBase, public IBStrategy
{
public:
    static const std::string SOURCE_SYSTEM_NAME;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > mask_var;
    int mask_current_idx, mask_new_idx, mask_scratch_idx;

    /*!
     * \brief Constructor for a single-part model.
     */
    IBFEMethod(const std::string& object_name,
               const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
               libMesh::MeshBase* mesh,
               int max_levels,
               bool register_for_restart = true,
               const std::string& restart_read_dirname = "",
               unsigned int restart_restore_number = 0);

    /*!
     * \brief Constructor for a multi-part model.
     */
    IBFEMethod(const std::string& object_name,
               const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
               const std::vector<libMesh::MeshBase*>& meshes,
               int max_levels,
               bool register_for_restart = true,
               const std::string& restart_read_dirname = "",
               unsigned int restart_restore_number = 0);

    /*!
     * \brief Deleted default constructor.
     */
    IBFEMethod() = delete;

    /*!
     * \brief Deleted copy constructor.
     */
    IBFEMethod(const IBFEMethod& from) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    IBFEMethod& operator=(const IBFEMethod& that) = delete;

    /*!
     * \brief Defaulted destructor.
     */
    ~IBFEMethod() override = default;

    /*!
     * Return a pointer to the finite element data manager object for the
     * specified part.
     */
    IBTK::FEDataManager* getFEDataManager(unsigned int part = 0) const;

    /*!
     * Indicate that a part should include a static pressure.
     *
     * @see FEMechanicsBase::registerStaticPressurePart for more details.
     *
     * @note A given part cannot be registered both for stress normalization and
     * also to have a static pressure. Attempting to do so will generate a fatal
     * error.
     */
    virtual void registerStaticPressurePart(PressureProjectionType projection_type = CONSISTENT_PROJECTION,
                                            FEMechanicsBase::VolumetricEnergyDerivativeFcn U_prime_fcn = nullptr,
                                            unsigned int part = 0) override;

    /*!
     * Indicate that a part should use stress normalization.
     *
     * @note A given part cannot be registered both for stress normalization and
     * also to have a static pressure. Attempting to do so will generate a fatal
     * error.
     */
    virtual void registerStressNormalizationPart(unsigned int part = 0);

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
                             std::vector<IBTK::SystemData> system_data = std::vector<IBTK::SystemData>(),
                             void* const ctx = nullptr)
            : fcn(fcn), system_data(std::move(system_data)), ctx(ctx)
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
    virtual void registerLagBodySourceFunction(const LagBodySourceFcnData& data, unsigned int part = 0);

    /*!
     * Get the Lagrangian body source function data.
     */
    LagBodySourceFcnData getLagBodySourceFunction(unsigned int part = 0) const;

    /*!
     * Register the (optional) direct forcing kinematics object with the finite
     * element mesh.
     */
    virtual void registerDirectForcingKinematics(const SAMRAI::tbox::Pointer<IBAMR::IBFEDirectForcingKinematics>& data,
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
     * Same as the base class.
     */
    virtual double getMaxPointDisplacement() const override;

    /*!
     * Inactivate a structure/part. See IBAMR::IBStrategy::inactivateLagrangianStructure().
     *
     * @note Since this class assumes that structures live on the finest grid
     * level the second argument is ignored.
     */
    virtual void inactivateLagrangianStructure(int structure_number = 0,
                                               int level_number = std::numeric_limits<int>::max()) override;

    /*!
     * Activate a previously inactivated structure/part to be used again in
     * FSI calculations. See IBAMR::IBStrategy::activateLagrangianStructure().
     *
     * @note Since this class assumes that structures live on the finest grid
     * level the second argument is ignored.
     */
    virtual void activateLagrangianStructure(int structure_number = 0,
                                             int level_number = std::numeric_limits<int>::max()) override;

    /*!
     * Determine whether or not the given structure or part is currently
     * activated. See IBAMR::IBStrategy::getLagrangianStructureIsActivated().
     *
     * @note Since this class assumes that structures live on the finest grid
     * level the second argument is ignored.
     */
    virtual bool getLagrangianStructureIsActivated(int structure_number = 0,
                                                   int level_number = std::numeric_limits<int>::max()) const override;

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
     * backward Euler method.
     */
    void backwardEulerStep(double current_time, double new_time) override;

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
                             int workload_data_idx) override;

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
     * this function only exists for compatibility with the base class and
     * does nothing: data reinitialization is handled by
     * endDataRedistribution() instead.
     *
     * The reasoning is this: since this class stores data only on particular
     * levels (at the present time, the structure is always on the finest
     * level) setting up level data is nontrivial when generating the initial
     * grid (i.e., when tagging cells that contain interaction points for
     * refinement). In a sense there is no level data to compute until we are
     * done regridding.
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
     * The name here is misleading, but SAMRAI expects us to use one of two
     * tagging methods to refine the grid, and IBAMR consistently uses
     * gradient detection: hence this function has the same name but tags
     * cells in a different way.
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
     * Return a pointer to the scratch hierarchy used by this object. See the
     * main documentation of this class for more information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > getScratchHierarchy();

protected:
    /*!
     * Do the actual work in initializeFEEquationSystems.
     */
    virtual void doInitializeFEEquationSystems() override;

    /*!
     * Do the actual work of setting up libMesh system vectors.
     */
    void doInitializeFESystemVectors() override;

    /*!
     * Do the actual work in reinitializeFEData and initializeFEData. if @p
     * use_present_data is `true` then the current content of the solution
     * vectors is used: more exactly, the coordinates and velocities (computed
     * by initializeCoordinates and initializeVelocity) are considered as
     * being up to date, as is the direct forcing kinematic data.
     */
    virtual void doInitializeFEData(bool use_present_data) override;

    /*!
     * \brief Compute the stress normalization field.
     */
    void computeStressNormalization(libMesh::PetscVector<double>& P_vec,
                                    libMesh::PetscVector<double>& X_vec,
                                    double data_time,
                                    unsigned int part);

    /*!
     * \brief Spread the transmission force density along the physical boundary
     * of the Lagrangian structure.
     */
    void spreadTransmissionForceDensity(int f_data_idx,
                                        libMesh::PetscVector<double>& X_ghost_vec,
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
     * Get the coarsest patch level number on which elements (including all
     * parts) are assigned.
     */
    int getCoarsestPatchLevelNumber() const;

    /*!
     * Get the finest patch level number on which elements (including all parts) are assigned.
     */
    int getFinestPatchLevelNumber() const;

    /*!
     * Get the schedule used to prolong force values. Data is read from @p
     * coarse_data_idx on level @p level_number and written into fine_data_idx
     * on level level_number + 1.
     */
    SAMRAI::xfer::RefineSchedule<NDIM>&
    getProlongationSchedule(int level_number, int coarse_data_idx, int fine_data_idx);

    /*!
     * Cached input databases.
     */
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_fe_data_manager_db;

    /*!
     * Whether or not the initial (i.e., before the regrid prior to
     * timestepping) workload calculations should be logged. This output is
     * generally not stable between machines and so this is usually disabled
     * in tests.
     */
    bool d_skip_initial_workload_log = false;

    /*!
     * Whether or not we have started time integration. This is only used to
     * determine whether or not we print some initial logging output: see
     * d_skip_initial_workload_log for more information.
     */
    bool d_started_time_integration = false;

    /*!
     * Boolean controlling whether or not the scratch hierarchy should be
     * used.
     */
    bool d_use_scratch_hierarchy = false;

    /*!
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;
    bool d_is_initialized = false;

    /*!
     * Scratch data caching object. Used both by this class and the
     * FEDataManager objects responsible for working on the primary
     * partitioning.
     */
    std::shared_ptr<IBTK::SAMRAIDataCache> d_eulerian_data_cache;

    int d_lagrangian_workload_current_idx = IBTK::invalid_index;
    int d_lagrangian_workload_new_idx = IBTK::invalid_index;
    int d_lagrangian_workload_scratch_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_lagrangian_workload_var;
    const std::string d_lagrangian_workload_coarsen_type = "CONSERVATIVE_COARSEN";
    const std::string d_lagrangian_workload_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    /*!
     * Maximum level number in the patch hierarchy.
     */
    int d_max_level_number = IBTK::invalid_level_number;

    /// CFL-like number used to determine when we should call
    /// reinitElementMappings() based on maximum structure point displacement.
    /// More exactly: this class will call that function if the maximum
    /// displacement of the structure (calculated by comparing the position
    /// vector as of the last reassociation to the current position vector)
    /// exceeds dx * d_patch_association_cfl, where dx is the smallest Eulerian
    /// cell width.
    ///
    /// Note that this is not a regridding, in the sense that the grid changes:
    /// instead only the association between patches and elements changes.
    ///
    /// @seealso IBHierarchyIntegrator::d_regrid_cfl_interval
    ///
    /// @note Most applications use a fluid solver regrid value of 0.5 - i.e.,
    /// the default value given here is a conservative choice.
    double d_patch_association_cfl = 0.75;

    /// Indexing information determining whether a given part is active or not.
    /// The default state for each part is to be active. Parts are active
    /// unless inactivated via inactivateLagrangianStructure().
    std::vector<bool> d_part_is_active{ true };

    /// FEDataManager objects associated with the primary hierarchy (i.e.,
    /// d_hierarchy). These are used by some other objects (such as
    /// IBFEPostProcessor); IBFEMethod keeps them up to date (i.e.,
    /// reinitializing data after regrids).
    std::vector<IBTK::FEDataManager*> d_primary_fe_data_managers;

    /// FEDataManager objects that use the scratch hierarchy instead of
    /// d_hierarchy. These are only used internally by IBFEMethod and are not
    /// intended to be accessed by any other object.
    std::vector<IBTK::FEDataManager*> d_scratch_fe_data_managers;

    /// The FEDataManager objects that are actually used in computations. This
    /// vector will be equal to either d_primary_fe_data_managers or
    /// d_scratch_fe_data_managers, dependent on which is actually used in IB
    /// calculations.
    std::vector<IBTK::FEDataManager*> d_active_fe_data_managers;

    /// Pointer to object used to accumulate forces during spreading.
    std::unique_ptr<IBTK::SAMRAIGhostDataAccumulator> d_ghost_data_accumulator;

    /*!
     * Schedules for prolonging data during spreading. The keys are the level
     * number, the patch data index for the coarse level data, and the patch
     * data index which will be filled with fine level data.
     */
    std::map<std::pair<int, std::pair<int, int> >, SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >
        d_prolongation_schedules;

    /// Minimum ghost cell width.
    SAMRAI::hier::IntVector<NDIM> d_ghosts = 0;

    /// Vectors of pointers to the systems for each part (for sources; other
    /// systems are handled by the base class).
    std::vector<libMesh::ExplicitSystem*> d_Q_systems;

    /*!
     * Object managing access to libMesh system vectors for the source/sink
     * strength.
     */
    std::unique_ptr<IBTK::LibMeshSystemVectors> d_Q_vecs;

    /*!
     * Object managing access to libMesh system vectors for the position that
     * include IB ghosting information.
     */
    std::unique_ptr<IBTK::LibMeshSystemIBVectors> d_X_IB_vecs;

    /*!
     * Object managing access to libMesh system vectors for the velocity that
     * include IB ghosting information.
     */
    std::unique_ptr<IBTK::LibMeshSystemIBVectors> d_U_IB_vecs;

    /*!
     * Object managing access to libMesh system vectors for the force that
     * include IB ghosting information.
     */
    std::unique_ptr<IBTK::LibMeshSystemIBVectors> d_F_IB_vecs;

    /*!
     * Object managing access to libMesh system vectors for the source/sink
     * strength that include IB ghosting information.
     */
    std::unique_ptr<IBTK::LibMeshSystemIBVectors> d_Q_IB_vecs;

    /*!
     * IBFE method parameters.
     */
    IBTK::FEDataManager::InterpSpec d_default_interp_spec;
    IBTK::FEDataManager::SpreadSpec d_default_spread_spec;
    IBTK::FEDataManager::WorkloadSpec d_default_workload_spec;
    std::vector<IBTK::FEDataManager::WorkloadSpec> d_workload_spec;
    std::vector<IBTK::FEDataManager::InterpSpec> d_interp_spec;
    std::vector<IBTK::FEDataManager::SpreadSpec> d_spread_spec;
    bool d_split_normal_force = false, d_split_tangential_force = false;
    bool d_use_jump_conditions = false;

    /*!
     * Data related to handling stress normalization.
     */
    double d_epsilon = 0.0;
    bool d_has_stress_normalization_parts = false;
    std::vector<bool> d_stress_normalization_part;

    /*!
     * Objects used to impose direct forcing kinematics.
     */
    std::vector<SAMRAI::tbox::Pointer<IBAMR::IBFEDirectForcingKinematics> > d_direct_forcing_kinematics_data;

    /*!
     * Functions used to compute source/sink strength on the Lagrangian mesh.
     */
    bool d_has_lag_body_source_parts = false;
    std::vector<bool> d_lag_body_source_part;
    std::vector<LagBodySourceFcnData> d_lag_body_source_fcn_data;

    /*!
     * Nonuniform load balancing data structures.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;
    int d_workload_idx = IBTK::invalid_index;

    /**
     * The optional second (i.e., scratch) hierarchy. This hierarchy is only
     * used for evaluation of IB terms, i.e., in
     * IBFEMethod::interpolateVelocity(), IBFEMethod::spreadForce(), and
     * IBFEMethod::spreadFluidSource().
     */
    std::unique_ptr<IBTK::SecondaryHierarchy> d_secondary_hierarchy;

private:
    /*!
     * Implementation of class constructor.
     */
    void commonConstructor(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db, int max_levels);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * At the present time this class and FEDataManager assume that the finite
     * element mesh is always on the finest grid level. This function
     * explicitly asserts that this condition is met.
     */
    void assertStructureOnFinestLevel() const;

    /*!
     * Convenience function that reinitializes the patch-to-element mappings on
     * all relevant FEDataManagers (i.e., for all parts and, if enabled, on the
     * scratch hierarchy).
     */
    void reinitElementMappings();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif //#ifndef included_IBAMR_IBFEMethod
