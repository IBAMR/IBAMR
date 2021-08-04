// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_FEDataManager
#define included_IBTK_FEDataManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#ifdef IBTK_HAVE_LIBMESH

#include "ibtk/QuadratureCache.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/ibtk_utilities.h"

#include "CellVariable.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "RefineSchedule.h"
#include "VariableContext.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/system.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/multi_array.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace IBTK
{
class FEDataManager;
class FEProjector;
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace SAMRAI
{
namespace geom
{
template <int DIM>
class CartesianPatchGeometry;
} // namespace geom
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class containing all of the finite element data structures that will be
 * used by FEDataManager and IBFEMethod.
 */
class FEData : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * Class which enables fast lookup of dofs on a given libMesh <code>elem</code>.
     *
     * @note The contents of this cache are invalidated when we regrid and the
     * caches should be reset at that point. Copies of this class should
     * always be retrieved via FEData::getDofCache() to avoid this
     * problem.
     */
    class SystemDofMapCache
    {
    public:
        /*!
         * Constructor.
         */
        inline SystemDofMapCache(libMesh::System& system) : d_dof_map(system.get_dof_map())
        {
        }

        /*!
         * Populate the vector @p dof_indices with the dofs corresponding to
         * variable @var on element @elem. The dof indices on each cell are
         * cached: i.e., the second call to this function with the same @p
         * elem and @p var is much faster than the first.
         */
        inline void
        dof_indices(const libMesh::Elem* const elem, std::vector<unsigned int>& dof_indices, const unsigned int var = 0)
        {
            const boost::multi_array<libMesh::dof_id_type, 2>& elem_dof_indices = d_dof_cache[elem->id()];
            copy_dof_ids_to_vector(var, elem_dof_indices, dof_indices);
            return;
        }

        /*!
         * Alternative indexing operation: retrieve all dof indices of all
         * variables in the given system at once by reference.
         */
        inline const boost::multi_array<libMesh::dof_id_type, 2>& dof_indices(const libMesh::Elem* const elem)
        {
            boost::multi_array<libMesh::dof_id_type, 2>& elem_dof_indices = d_dof_cache[elem->id()];
            if (elem_dof_indices.shape()[0] == 0)
            {
                d_scratch_dofs.resize(d_dof_map.n_variables());
                // we have to be careful - different elements may have different
                // numbers of dofs, so grab all of them first:
                for (unsigned int var_n = 0; var_n < d_dof_map.n_variables(); ++var_n)
                {
                    d_scratch_dofs[var_n].clear();
                    d_dof_map.dof_indices(elem, d_scratch_dofs[var_n], var_n);
                }
                // We assume in a lot of other places that all variables in a
                // system have the same finite element (i.e., dofs per element)
                // so assert that here too
                for (unsigned int var_n = 0; var_n < d_dof_map.n_variables(); ++var_n)
                {
                    TBOX_ASSERT(d_scratch_dofs[var_n].size() == d_scratch_dofs[0].size());
                }

                boost::multi_array<libMesh::dof_id_type, 2>::extent_gen extents;
                elem_dof_indices.resize(extents[d_scratch_dofs.size()][d_scratch_dofs[0].size()]);
                for (unsigned int var_n = 0; var_n < d_dof_map.n_variables(); ++var_n)
                {
                    std::copy(d_scratch_dofs[var_n].begin(), d_scratch_dofs[var_n].end(), &elem_dof_indices[var_n][0]);
                }
            }
            return elem_dof_indices;
        }

    private:
        libMesh::DofMap& d_dof_map;
        std::unordered_map<libMesh::dof_id_type, boost::multi_array<libMesh::dof_id_type, 2> > d_dof_cache;
        std::vector<std::vector<libMesh::dof_id_type> > d_scratch_dofs;
    };

    /*!
     * Constructor. Registers the object with the restart database: i.e.,
     * inheriting classes should not also register themselves.
     */
    FEData(std::string object_name, libMesh::EquationSystems& equation_systems, const bool register_for_restart);

    /*!
     * Destructor.
     */
    ~FEData();

    /*!
     * Set up the object with data stored in the restart database.
     *
     * @note this is called inside the constructor.
     */
    virtual void getFromRestart();

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \brief Set the equations systems object that is associated with the
     * FEData object.
     *
     * @deprecated The equation systems object should be set by the constructor.
     */
    void setEquationSystems(libMesh::EquationSystems* equation_systems, int level_number);

    /*!
     * \return A pointer to the equations systems object that is associated with
     * the FEData object.
     */
    libMesh::EquationSystems* getEquationSystems() const;

    /*!
     * \return The DofMapCache for a specified system.
     */
    SystemDofMapCache* getDofMapCache(const std::string& system_name);

    /*!
     * \return The DofMapCache for a specified system.
     */
    SystemDofMapCache* getDofMapCache(unsigned int system_num);

    /*!
     * Clear all cached (i.e., computed at first request and then stored for
     * future calls) data that depends on the Eulerian data partitioning.
     */
    void clearPatchHierarchyDependentData();

protected:
    /*!
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*!
     * Name of the coordinates system.
     *
     * @note For backwards compatibility reasons this string may be reassigned
     * to another value by assignment to
     * FEDataManager::COORDINATES_SYSTEM_NAME.
     */
    std::string d_coordinates_system_name = "coordinates system";

    /*!
     * FE equation system associated with this data manager object.
     */
    libMesh::EquationSystems* d_es = nullptr;

    /*!
     * Cache of libMesh quadrature objects. Defaults to being an NDIM cache
     * but is overwritten once a libMesh::EquationSystems object is attached.
     */
    QuadratureCache d_quadrature_cache;

    /*!
     * Mapping between system numbers and SystemDofMapCache objects.
     */
    std::map<std::pair<unsigned int, libMesh::FEType>, std::unique_ptr<SystemDofMapCache> > d_system_dof_map_cache;

    /**
     * Permit FEDataManager to directly examine the internals of this class.
     */
    friend class FEDataManager;
};

/*!
 * Class that can translate libMesh subdomain IDs into patch level numbers.
 *
 * The primary use of this class is multilevel IBFE - i.e., enabling a finite
 * element mesh to interact with multiple patch levels.
 */
class SubdomainToPatchLevelTranslation
{
public:
    /*!
     * Constructor. Takes as argument the maximum level number in the Cartesian
     * grid patch hierarchy, a set of all subdomain ids (including those not on
     * the current processor), and an input database enumerating the
     * level-to-subdomain mapping, e.g.,
     * @code
     * {
     *   level_1 = 1, 2
     *   level_2 = 3
     * }
     * @endcode
     * Any unspecified subdomain ids will be assigned to the finest patch level.
     * Duplicated assignments are not permitted.
     */
    SubdomainToPatchLevelTranslation(const int max_level_number,
                                     const std::set<libMesh::subdomain_id_type>& subdomain_ids,
                                     const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db);

    /*!
     * Given a libMesh subdomain id, return the patch level of the Cartesian
     * grid hierarchy on which that subdomain id interacts.
     */
    const int& operator[](const libMesh::subdomain_id_type id) const
    {
        if (id < fixed_array_size) return d_fixed[id];
        const auto it = d_map.find(id);
        if (it == d_map.end()) return d_max_level_number;
        return it->second;
    }

    /*!
     * Return whether or not there are any elements on a given patch level.
     */
    bool levelHasElements(const int level_number) const
    {
        const auto array_it = std::find(d_fixed.begin(), d_fixed.end(), level_number);
        if (array_it != d_fixed.end())
        {
            return true;
        }
        else
        {
            return std::find_if(d_map.begin(),
                                d_map.end(),
                                [&](const std::pair<libMesh::subdomain_id_type, int>& pair)
                                { return pair.second == level_number; }) != d_map.end();
        }
    }

private:
    /*!
     * like operator[], but returns a mutable reference. Used to set up the object.
     */
    int& get(const libMesh::subdomain_id_type id)
    {
        if (id < d_fixed.size()) return d_fixed[id];
        return d_map[id];
    }

    /*!
     * Maximum level number.
     */
    int d_max_level_number = IBTK::invalid_level_number;

    /*!
     * Size of the fixed-size array.
     */
    static constexpr int fixed_array_size = 1024;

    /*!
     * The overwhelming majority of subdomain IDs used with IBAMR come from
     * block IDs set by ExodusII - these are numbered sequentially from zero.
     * However, in principle, a subdomain ID could be any signed 64-bit integer
     * so we cannot use a small fixed length array.
     *
     * Since we look up element levels a lot optimize for the common case by
     * using a fixed-length array and a map for everything else.
     */
    std::array<int, fixed_array_size> d_fixed;

    /*!
     * The map used for everything else.
     */
    std::map<libMesh::subdomain_id_type, int> d_map;
};

/*!
 * \brief Class FEDataManager coordinates data required for
 * Lagrangian-Eulerian interaction between a Lagrangian finite element (FE)
 * mesh. In particular, the FEData member object stores the necessary finite
 * element data while this class stores additional data dependent on the
 * Eulerian grid.
 *
 * <h2>Parameters read from the input database</h2>
 *
 * <code>node_outside_patch_check</code>: parameter controlling how this class
 * responds to mesh nodes outside the finest patch level. In all cases, for
 * backwards compatibility, nodes that are outside the computational domain are
 * permitted and are ignored by this check. Possible values are:
 * <ol>
 *   <li><code>node_outside_permit</code>: Permit nodes to be outside the finest
 *   patch level.</li>
 *   <li><code>node_outside_warn</code>: Permit nodes to be outside the finest
 *   patch level, but log a warning whenever this is detected.
 *   <li><code>node_outside_error</code>: Abort the program when nodes are detected
 *   outside the finest patch level.
 * </ol>
 * The default value is <code>node_outside_error</code>.
 *
 * <code>subdomain_ids_on_levels</code>: a database correlating libMesh subdomain
 * IDs to patch levels. A possible value for this is
 * @code
 * subdomain_ids_on_levels
 * {
 *   level_1 = 4
 *   level_3 = 10, 12, 14
 *   level_5 = 1000, 1003, 1006, 1009
 * }
 * @endcode
 * This particular input will associate all elements with subdomain id 4 with
 * patch level 1, all elements with subdomain ids 10, 12, or 14 with patch level
 * 3, etc. All unspecified subdomain ids will be associated with the finest
 * patch level. All inputs in this database for levels finer than the finest
 * level are ignored (e.g., if the maximum patch level number is 4, then the
 * values given in the example for level 5 ultimately end up on level 4).
 * <em>This feature is experimental</em>: at the current time it is known that
 * it produces some artifacts at the coarse-fine interface, but that these
 * generally don't effect the overall solution quality.
 *
 * <h2>Parameters effecting workload estimate calculations</h2>
 * FEDataManager can estimate the amount of work done in IBFE calculations
 * (such as FEDataManager::spread). Since most calculations use a variable
 * number of quadrature points on each libMesh element this estimate can vary
 * quite a bit over different Eulerian cells corresponding to a single
 * mesh. The current implementation estimates the workload on each cell of the
 * background Eulerian grid by applying a background value representing the
 * work on the Eulerian cell itself and a weight times the number of
 * quadrature points on that cell. These values are set at the time of object
 * construction through the FEDataManager::WorkloadSpec object, which contains
 * reasonable defaults.
 *
 * \note Multiple FEDataManager objects may be instantiated simultaneously.
 */
class FEDataManager : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * Alias FEData::SystemDofMapCache for backwards compatibility.
     */
    using SystemDofMapCache = FEData::SystemDofMapCache;

    /*!
     * \brief Struct InterpSpec encapsulates data needed to specify the manner
     * in which Eulerian-to-Lagrangian interpolation is performed when using an
     * FE structural discretization.
     */
    struct InterpSpec
    {
        InterpSpec() = default;

        InterpSpec(const std::string& kernel_fcn,
                   const libMesh::QuadratureType& quad_type,
                   const libMesh::Order& quad_order,
                   bool use_adaptive_quadrature,
                   double point_density,
                   bool use_consistent_mass_matrix,
                   bool use_nodal_quadrature,
                   bool allow_rules_with_negative_weights = true)
            : kernel_fcn(kernel_fcn),
              quad_type(quad_type),
              quad_order(quad_order),
              use_adaptive_quadrature(use_adaptive_quadrature),
              point_density(point_density),
              use_consistent_mass_matrix(use_consistent_mass_matrix),
              use_nodal_quadrature(use_nodal_quadrature),
              allow_rules_with_negative_weights(allow_rules_with_negative_weights)
        {
        }

        std::string kernel_fcn;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
        bool use_adaptive_quadrature;
        double point_density;
        bool use_consistent_mass_matrix;
        bool use_nodal_quadrature;
        bool allow_rules_with_negative_weights;
    };

    /*!
     * \brief Struct SpreadSpec encapsulates data needed to specify the manner
     * in which Lagrangian-to-Eulerian spreading is performed when using an FE
     * structural discretization.
     */
    struct SpreadSpec
    {
        SpreadSpec() = default;

        SpreadSpec(const std::string& kernel_fcn,
                   const libMesh::QuadratureType& quad_type,
                   const libMesh::Order& quad_order,
                   bool use_adaptive_quadrature,
                   double point_density,
                   bool use_nodal_quadrature,
                   bool allow_rules_with_negative_weights = true)
            : kernel_fcn(kernel_fcn),
              quad_type(quad_type),
              quad_order(quad_order),
              use_adaptive_quadrature(use_adaptive_quadrature),
              point_density(point_density),
              use_nodal_quadrature(use_nodal_quadrature),
              allow_rules_with_negative_weights(allow_rules_with_negative_weights)
        {
        }

        std::string kernel_fcn;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
        bool use_adaptive_quadrature;
        double point_density;
        bool use_nodal_quadrature;
        bool allow_rules_with_negative_weights;
    };

    /*!
     * \brief Struct WorkloadSpec encapsulates the parameters used to
     * calculate workload estimates (i.e., the input to the load balancing
     * algorithm).
     */
    struct WorkloadSpec
    {
        /// The multiplier applied to each quadrature point. This value accounts
        /// for work done at each IB point (e.g., the work done inside the
        /// Fortran spreading and interpolation kernels).
        ///
        /// If nodal quadrature is used then this value simply corresponds to
        /// counting the nodes since those are the quadrature points.
        double q_point_weight = 1.0;

        /// The multiplier applied to the nodes of elements. This value accounts
        /// for the work done for each element regardless of the number of
        /// quadrature points (e.g. calculating the size of the element
        /// in the deformed configuration).
        ///
        /// These work values are calculated in an unusual way in the sense that
        /// if a node is attached to N elements, then we count that node N
        /// times. This accounts for the fact that work is done on an element
        /// level but elements may exist on multiple patches at once: i.e., it
        /// makes more sense to compute the work associated with an element by
        /// assigning work to its nodes rather than the element centroid.
        ///
        /// This value should not be set to anything but zero when doing nodal
        /// quadrature.
        ///
        /// A good value for this is 0.8.
        double duplicated_node_weight = 0.0;
    };

protected:
    /*!
     * FEData object that contains the libMesh data structures.
     *
     * @note multiple FEDataManager objects may use the same FEData object,
     * usually combined with different hierarchies.
     */
    std::shared_ptr<FEData> d_fe_data;

    /*!
     * FEProjector object that handles L2 projection functionality.
     */
    std::shared_ptr<FEProjector> d_fe_projector;

    /*!
     * IB ghosted diagonal mass matrix representations.
     */
    std::map<std::string, std::unique_ptr<libMesh::PetscVector<double> > > d_L2_proj_matrix_diag_ghost;

public:
    /*!
     * \brief The name of the equation system which stores the spatial position
     * data. The actual string is stored by FEData.
     *
     * \note The default value for this string is "coordinates system".
     */
    std::string& COORDINATES_SYSTEM_NAME;

    /*!
     * \brief The libMesh boundary IDs to use for specifying essential boundary
     * conditions.
     */
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_X_BDRY_ID;
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_Y_BDRY_ID;
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_Z_BDRY_ID;
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_XY_BDRY_ID;
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_XZ_BDRY_ID;
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_YZ_BDRY_ID;
    static const libMesh::boundary_id_type ZERO_DISPLACEMENT_XYZ_BDRY_ID;

    /*!
     * Return a pointer to the instance of the Lagrangian data manager
     * corresponding to the specified name.  Access to FEDataManager objects is
     * mediated by the getManager() function.
     *
     * @note When a manager is accessed for the first time, the
     * FEDataManager::freeAllManagers() static method is registered with the
     * SAMRAI::tbox::ShutdownRegistry class. Consequently, all allocated
     * managers are freed at program completion. Thus, users of this class do
     * not explicitly allocate or deallocate the FEDataManager instances.
     *
     * \return A pointer to the data manager instance.
     */
    static FEDataManager*
    getManager(std::shared_ptr<FEData> fe_data,
               const std::string& name,
               const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
               const int max_levels,
               const InterpSpec& default_interp_spec,
               const SpreadSpec& default_spread_spec,
               const WorkloadSpec& default_workload_spec,
               const SAMRAI::hier::IntVector<NDIM>& min_ghost_width = SAMRAI::hier::IntVector<NDIM>(0),
               std::shared_ptr<SAMRAIDataCache> eulerian_data_cache = nullptr,
               bool register_for_restart = true);

    /*!
     * Deallocate all of the FEDataManager instances.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the SAMRAI::tbox::ShutdownRegistry class.
     */
    static void freeAllManagers();

    /*!
     * \brief Set the equations systems object that is associated with the
     * FEData object. Currently, each set of equation systems must be assigned
     * to a particular level of the AMR grid.
     *
     * @deprecated This function is deprecated since the FEData constructor now
     * requires an EquationSystems argument.
     */
    void setEquationSystems(libMesh::EquationSystems* equation_systems, int level_number);

    /*!
     * \return A pointer to the equations systems object that is associated with
     * the FEData object.
     */
    libMesh::EquationSystems* getEquationSystems() const;

    /*!
     * \return The DofMapCache for a specified system.
     */
    SystemDofMapCache* getDofMapCache(const std::string& system_name);

    /*!
     * \return The DofMapCache for a specified system.
     */
    SystemDofMapCache* getDofMapCache(unsigned int system_num);

    /*!
     * \brief Enable or disable logging.
     *
     * @note This is usually set by the IBFEMethod which owns the current
     * FEDataManager, which reads the relevant boolean from the database.
     */
    void setLoggingEnabled(bool enable_logging = true);

    /*!
     * \brief Determine whether logging is enabled or disabled.
     */
    bool getLoggingEnabled() const;

    /*!
     * \name Methods to set and get the patch hierarchy and range of patch
     * levels associated with this manager class.
     */
    //\{

    /*!
     * \brief Reset patch hierarchy over which operations occur.
     *
     * The patch hierarchy must be fully set up (i.e., contain all the levels it
     * is expected to have) at the point this function is called. If you need to
     * tag cells for refinement to create the initial hierarchy then use
     * applyGradientDetector, which does not use the stored patch hierarchy.
     */
    void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Get the patch hierarchy used by this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > getPatchHierarchy() const;

    /*!
     * Get the coarsest patch level number on which elements are assigned.
     */
    int getCoarsestPatchLevelNumber() const;

    /*!
     * Get the finest patch level number on which elements are assigned.
     */
    int getFinestPatchLevelNumber() const;

    /*!
     * \return The ghost cell width used for quantities that are to be
     * interpolated from the Cartesian grid to the FE mesh.
     */
    const SAMRAI::hier::IntVector<NDIM>& getGhostCellWidth() const;

    /*!
     * \return The specifications of the scheme used for interpolating from the
     * Cartesian grid to the FE mesh.
     */
    const InterpSpec& getDefaultInterpSpec() const;

    /*!
     * \return The specifications of the scheme used for spreading densities
     * from the FE mesh to the Cartesian grid
     */
    const SpreadSpec& getDefaultSpreadSpec() const;

    /*!
     * \return A const reference to the map from local patch number to local
     * active elements.
     */
    const std::vector<std::vector<libMesh::Elem*> >& getActivePatchElementMap() const;

    /*!
     * \return A const reference to the map from local patch number to local
     * active nodes.
     *
     * \note The local active nodes are the nodes of the local active elements.
     */
    const std::vector<std::vector<libMesh::Node*> >& getActivePatchNodeMap() const;

    /*!
     * \brief Reinitialize the mappings from elements to Cartesian grid patches.
     */
    void reinitElementMappings();

    /*!
     * \return A pointer to the unghosted solution vector associated with the
     * specified system.
     */
    libMesh::NumericVector<double>* getSolutionVector(const std::string& system_name) const;

    /*!
     * \return A pointer to the ghosted solution vector associated with the
     * specified system. The vector contains positions for values in the
     * relevant IB ghost region which are populated if @p localize_data is
     * <code>true</code>.
     *
     * @note The vector returned by pointer is owned by this class (i.e., no
     * copying is done).
     *
     * @deprecated Use buildIBGhostedVector() instead which clones a vector
     * with the same ghost region.
     */
    libMesh::NumericVector<double>* buildGhostedSolutionVector(const std::string& system_name,
                                                               bool localize_data = true);

    /*!
     * \return A pointer to a vector, with ghost entries corresponding to
     * relevant IB data, associated with the specified system.
     */
    std::unique_ptr<libMesh::PetscVector<double> > buildIBGhostedVector(const std::string& system_name);

    /*!
     * \return A pointer to the unghosted coordinates (nodal position) vector.
     */
    libMesh::NumericVector<double>* getCoordsVector() const;

    /*!
     * \return A pointer to the ghosted coordinates (nodal position) vector.
     *
     * @deprecated Use buildIBGhostedVector() instead.
     */
    libMesh::NumericVector<double>* buildGhostedCoordsVector(bool localize_data = true);

    /*!
     * \return The shared pointer to the object managing the Lagrangian data.
     */
    std::shared_ptr<FEData>& getFEData();

    /*!
     * \return The shared pointer to the object managing the Lagrangian data.
     */
    const std::shared_ptr<FEData>& getFEData() const;

    /*!
     * \brief Spread a density from the FE mesh to the Cartesian grid using the
     * default spreading spec.
     *
     * @note This function may spread forces from points near the physical
     * boundary into ghost cells outside the physical domain. To account for
     * these forces one should usually call
     * IBTK::RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData()
     * after this function to move said force values into the physical domain.
     */
    void spread(int f_data_idx,
                libMesh::NumericVector<double>& F,
                libMesh::NumericVector<double>& X,
                const std::string& system_name);

    /*!
     * \brief Spread a density from the FE mesh to the Cartesian grid.
     *
     * @param[in] f_data_idx SAMRAI index into the
     * SAMRAI::hier::PatchHierarchy object owned by this class. The spread
     * force values will be added into this index.
     * @param[in] F Finite element solution vector containing the field which
     * will be spread onto the Eulerian grid.
     * @param[in] X Finite element solution vector containing the current
     * position of the mesh.
     * @param[in] system_name name of the system corresponding to @p F.
     *
     * Both @p X and @p F should contain ghost values corresponding to the IB
     * partitioning of the Lagrangian data, i.e., vectors returned from
     * buildIBGhostedVector.
     *
     * @note This function may spread forces from points near the physical
     * boundary into ghost cells outside the physical domain. To account for
     * these forces one should usually call
     * IBTK::RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData()
     * after this function to move said force values into the physical domain.
     */
    void spread(int f_data_idx,
                libMesh::NumericVector<double>& F,
                libMesh::NumericVector<double>& X,
                const std::string& system_name,
                const SpreadSpec& spread_spec);

    /*!
     * \brief Spread a density from the FE mesh to the Cartesian grid using
     * the default spreading spec. This is a convenience overload of spread
     * where @p f_phys_bdry_op is used to correctly accumulate force values
     * spread outside the physical domain into it.
     */
    void spread(int f_data_idx,
                libMesh::NumericVector<double>& F,
                libMesh::NumericVector<double>& X,
                const std::string& system_name,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                double fill_data_time,
                bool close_F = true,
                bool close_X = true);

    /*!
     * \brief Spread a density from the FE mesh to the Cartesian grid using a
     * specified spreading spec. This is a convenience overload of spread
     * where @p f_phys_bdry_op is used to correctly accumulate force values
     * spread outside the physical domain into it.
     */
    void spread(int f_data_idx,
                libMesh::NumericVector<double>& F,
                libMesh::NumericVector<double>& X,
                const std::string& system_name,
                const SpreadSpec& spread_spec,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                double fill_data_time,
                bool close_F = true,
                bool close_X = true);

    /*!
     * \brief Prolong a value or a density from the FE mesh to the Cartesian
     * grid.
     */
    void prolongData(int f_data_idx,
                     libMesh::NumericVector<double>& F,
                     libMesh::NumericVector<double>& X,
                     const std::string& system_name,
                     bool is_density = true,
                     bool accumulate_on_grid = true,
                     bool close_F = true,
                     bool close_X = true);

    /*!
     * \brief Set up the right-hand side @p F of an L2 projection problem
     * where Eulerian data given by @p f_data_idx is projected onto the finite
     * element space given by @p system_name.
     *
     * @param[in] f_data_idx Index of the variable being projected.
     *
     * @param[in] IB-ghosted position vector containing the current position
     * of all dofs whose basis functions have support (in the deformed
     * configuration) on the locally owned set of patches.
     *
     * @param[out] F Vector into which the L2-projection RHS vector will be
     * assembled. If this vector is ghosted then off-processor entries are
     * assembled into the ghost value region of the vector and callers should
     * complete the assembly process with, e.g.,
     * @code
     * ierr = VecGhostUpdateBegin(F.vec(), ADD_VALUES, SCATTER_REVERSE);
     * IBTK_CHKERRQ(ierr);
     * ierr = VecGhostUpdateEnd(F.vec(), INSERT_VALUES, SCATTER_FORWARD);
     * IBTK_CHKERRQ(ierr);
     * @endcode
     * otherwise values are added into the vector directly (which, for a
     * libMesh::PetscVector, will involve the use of the PETSc VecCache object
     * to track off-processor entries). It is strongly recommended, for
     * performance reasons, that one assemble into a ghosted vector instead of
     * the approach with VecCache.
     *
     * @param[in] system_name Name of the libMesh system corresponding to the
     * vector @p F.
     *
     * @param[in] f_refine_scheds Refinement schedules to process before
     * actually running this function.
     *
     * @param[in] fill_data_time Time at which the data in @p f_data_idx was filled.
     *
     * @param[in] Whether or not to close @p F after assembly.
     *
     * @param[in] Whether or not to close @p X before assembly.
     *
     * @note We recommend against using the last two booleans: it is usually
     * better to do any vector communication before calling this function.
     *
     * @note This function is poorly named: it actually sets up a dual-space
     * vector @p F which is the right-hand side of an L2 projection
     * problem. Callers will still need to solve the resulting linear
     * system. The result, therefore, is projected, not interpolated.
     *
     * @note The vector set up by this function is set up without applying any
     * constraints (e.g., hanging node or periodicity constraints) during
     * assembly. This is because we do not store the constraints on the IB
     * data partitioning, only on libMesh's own data partitioning. For more
     * information on how to assemble RHS vectors in this case, see the
     * documentation of the function apply_transposed_constraint_matrix().
     */
    void
    interpWeighted(int f_data_idx,
                   libMesh::NumericVector<double>& F,
                   libMesh::NumericVector<double>& X,
                   const std::string& system_name,
                   const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_refine_scheds =
                       std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                   double fill_data_time = 0.0,
                   bool close_F = true,
                   bool close_X = true);

    /*!
     * \brief Interpolate a value from the Cartesian grid to the FE mesh using a
     * specified interpolation spec. This interpolation function does NOT do
     * an L2-projection of the interpolated quantity. It does however weighs/filters
     * the interpolated quantity at the quadrature points to the nodes. Here, the
     * basis functions of the deformational field is used as the filter.
     */
    void
    interpWeighted(int f_data_idx,
                   libMesh::NumericVector<double>& F,
                   libMesh::NumericVector<double>& X,
                   const std::string& system_name,
                   const InterpSpec& interp_spec,
                   const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_refine_scheds =
                       std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                   double fill_data_time = 0.0,
                   bool close_F = true,
                   bool close_X = true);

    /*!
     * \brief Interpolate a value from the Cartesian grid to the FE mesh using
     * the default interpolation spec.
     */
    void interp(int f_data_idx,
                libMesh::NumericVector<double>& F,
                libMesh::NumericVector<double>& X,
                const std::string& system_name,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_refine_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool close_X = true);

    /*!
     * \brief Interpolate a value from the Cartesian grid to the FE mesh using a
     * specified interpolation spec.
     */
    void interp(int f_data_idx,
                libMesh::NumericVector<double>& F,
                libMesh::NumericVector<double>& X,
                const std::string& system_name,
                const InterpSpec& interp_spec,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_refine_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool close_X = true);

    /*!
     * \brief Restrict a value from the Cartesian grid to the FE mesh.
     */
    void restrictData(int f_data_idx,
                      libMesh::NumericVector<double>& F,
                      libMesh::NumericVector<double>& X,
                      const std::string& system_name,
                      bool use_consistent_mass_matrix = true,
                      bool close_X = true);

    /*!
     * \return Pointer to vector representation of diagonal L2 mass matrix.
     */
    libMesh::NumericVector<double>* buildDiagonalL2MassMatrix(const std::string& system_name);

    /*!
     * \return Pointer to IB ghosted vector representation of diagonal L2 mass matrix.
     */
    libMesh::PetscVector<double>* buildIBGhostedDiagonalL2MassMatrix(const std::string& system_name);

    /*!
     * \brief Set U to be the L2 projection of F.
     */
    bool computeL2Projection(libMesh::NumericVector<double>& U,
                             libMesh::NumericVector<double>& F,
                             const std::string& system_name,
                             bool consistent_mass_matrix = true,
                             bool close_U = true,
                             bool close_F = true,
                             double tol = 1.0e-6,
                             unsigned int max_its = 100);

    /*!
     * Update the quadrature rule for the current element.  If the provided
     * qrule is already configured appropriately, it is not modified.
     *
     * \return true if the quadrature rule is updated or otherwise requires
     * reinitialization (e.g. because the element type or p_level changed);
     * false otherwise.
     */
    static bool updateQuadratureRule(std::unique_ptr<libMesh::QBase>& qrule,
                                     libMesh::QuadratureType quad_type,
                                     libMesh::Order quad_order,
                                     bool use_adaptive_quadrature,
                                     double point_density,
                                     bool allow_rules_with_negative_weights,
                                     const libMesh::Elem* elem,
                                     const boost::multi_array<double, 2>& X_node,
                                     double dx_min);

    /*!
     * Update the quadrature rule for the current element used by the
     * Lagrangian-Eulerian interaction scheme.  If the provided qrule is already
     * configured appropriately, it is not modified.
     *
     * \return true if the quadrature rule is updated or otherwise requires
     * reinitialization (e.g. because the element type or p_level changed);
     * false otherwise.
     */
    static bool updateInterpQuadratureRule(std::unique_ptr<libMesh::QBase>& qrule,
                                           const InterpSpec& spec,
                                           const libMesh::Elem* elem,
                                           const boost::multi_array<double, 2>& X_node,
                                           double dx_min);

    /*!
     * Update the quadrature rule for the current element used by the
     * Lagrangian-Eulerian interaction scheme.  If the provided qrule is already
     * configured appropriately, it is not modified.
     *
     * \return true if the quadrature rule is updated or otherwise requires
     * reinitialization (e.g. because the element type or p_level changed);
     * false otherwise.
     */
    static bool updateSpreadQuadratureRule(std::unique_ptr<libMesh::QBase>& qrule,
                                           const SpreadSpec& spec,
                                           const libMesh::Elem* elem,
                                           const boost::multi_array<double, 2>& X_node,
                                           double dx_min);

    /*!
     * \brief Update the cell workload estimate by adding a value (see the
     * main documentation of this class for information on how this is
     * computed) to the <code>d_workload_idx</code> cell variable.
     */
    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const int workload_data_idx,
                             const int coarsest_ln = -1,
                             const int finest_ln = -1);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur due to the presence of Lagrangian data.  The double time
     * argument is the regrid time.  The integer "tag_index" argument is the
     * patch descriptor index of the cell centered integer tag array on each
     * patch in the hierarchy.  The boolean argument initial_time indicates
     * whether the level is being subject to refinement at the initial
     * simulation time.  If it is false, then the error estimation process is
     * being invoked at some later time after the AMR hierarchy was initially
     * constructed.  The boolean argument uses_richardson_extrapolation_too is
     * true when Richardson extrapolation error estimation is used in addition
     * to the gradient detector, and false otherwise.  This argument helps the
     * user to manage multiple regridding criteria.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null or the level number does not match any
     * existing level in the hierarchy.
     *
     * @note This function is analogous to
     * SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector() and is
     * only meant to be called from IBAMR::IBFEMethod::applyGradientDetector().
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too);

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \brief Zero the values corresponding to points on the given patch
     * (described by @p patch_geometry) that might be duplicated on another
     * patch.
     *
     * In the process of computing either forces on the Lagrangian mesh to
     * spread to the Eulerian grid or interpolating values from the Eulerian
     * grid to use in the Lagrangian structure we may, under extraordinary
     * circumstances, double-count a point that lies on the boundary between
     * two patches. Note that due to the CFL condition assumed by this class
     * and the width of the tag buffer used by IBAMR::IBFEMethod we can assume
     * that no IB point will ever lie on the boundary of the finest grid
     * level: hence, if a point lies on the boundary of a patch, it must also
     * lie on the boundary of a neighboring patch.
     *
     * Hence, to establish uniqueness, we zero data that lies on the 'upper'
     * face (corresponding to the unit normal vector pointing towards
     * infinity) since we are guaranteed, by the previous assumptions, that
     * that data must also lie on the lower face of a neighboring patch.
     */
    static void zeroExteriorValues(const SAMRAI::geom::CartesianPatchGeometry<NDIM>& patch_geom,
                                   const std::vector<double>& X_qp,
                                   std::vector<double>& F_qp,
                                   int n_vars);

protected:
    /*!
     * \brief Constructor, where the FEData object owned by this class may be
     * co-owned by other objects.
     */
    FEDataManager(std::string object_name,
                  const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db,
                  const int max_levels,
                  InterpSpec default_interp_spec,
                  SpreadSpec default_spread_spec,
                  WorkloadSpec default_workload_spec,
                  SAMRAI::hier::IntVector<NDIM> ghost_width,
                  std::shared_ptr<SAMRAIDataCache> eulerian_data_cache,
                  std::shared_ptr<FEData> fe_data,
                  bool register_for_restart = true);

    /*!
     * \brief The FEDataManager destructor cleans up any allocated data objects.
     */
    ~FEDataManager();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FEDataManager() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FEDataManager(const FEDataManager& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FEDataManager& operator=(const FEDataManager& that) = delete;

    /*!
     * Compute the quadrature point counts in each cell of the level in which
     * the FE mesh is embedded.  Also zeros out node count data for other levels
     * within the specified range of level numbers.
     */
    void updateQuadPointCountData(int coarsest_ln, int finest_ln);

    /*!
     * Collect all of the active elements which are located within a local
     * Cartesian grid patch grown by a ghost width of 1 (like
     * IBTK::LEInteractor::getMinimumGhostWidth(), we assume that IB points
     * are allowed to move no more than one cell width between regridding
     * operations).
     *
     * The parameters refer to the levels of different objects:
     * <ol>
     *   <li>@p level_number - the level number in the patch hierarchy on which
     *     we are identifying intersections.</li>
     *   <li>@p coarsest_elem_ln - The minimum level number of elements we should
     *     consider (see the main documentation of this class for an explanation
     *     on how elements are assigned to particular levels)</li>
     *   <li>@p finest_elem_ln - The maximum level number of elements we should
     *     consider.</li>
     * </ol>
     *
     * All three parameters are necessary because we use this function both to
     * tag cells for refinement (i.e., we want to refine cells containing
     * elements on levels higher than the present level) and to do IB
     * calculations (where all three numbers will be the same).
     *
     * In this method, the determination as to whether an element is local or
     * not is based on the position of the bounding box of the element.
     */
    void collectActivePatchElements(std::vector<std::vector<libMesh::Elem*> >& active_patch_elems,
                                    int level_number,
                                    int coarsest_elem_ln,
                                    int finest_elem_ln);

    /*!
     * Collect all of the nodes of the active elements that are located within a
     * local Cartesian grid patch grown by the specified ghost cell width.
     */
    void collectActivePatchNodes(std::vector<std::vector<libMesh::Node*> >& active_patch_nodes,
                                 const std::vector<std::vector<libMesh::Elem*> >& active_patch_elems);

    /*!
     * Store the association between subdomain ids and patch levels.
     */
    SubdomainToPatchLevelTranslation d_level_lookup;

    /*!
     * Get the patch level on which an element lives.
     */
    int getPatchLevel(const libMesh::Elem* elem) const;

    /*!
     * Collect all ghost DOF indices for the specified collection of elements.
     */
    void collectGhostDOFIndices(std::vector<unsigned int>& ghost_dofs,
                                const std::vector<libMesh::Elem*>& active_elems,
                                const std::string& system_name);

    /*!
     * Reinitialize IB ghosted DOF data structures for the specified system.
     */
    void reinitializeIBGhostedDOFs(const std::string& system_name);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *
     *    -   The class version number and restart version number do not match.
     *
     */
    virtual void getFromRestart();

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static std::map<std::string, FEDataManager*> s_data_manager_instances;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*!
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*!
     * Whether or not to log data to the screen: see
     * FEDataManager::setLoggingEnabled() and
     * FEDataManager::getLoggingEnabled().
     *
     * @note This is usually set by IBFEMethod, which reads the relevant
     * boolean from the database.
     */
    bool d_enable_logging = false;

    /*!
     * Grid hierarchy information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    /*!
     * Maximum possible level number in the patch hierarchy.
     */
    int d_max_level_number = IBTK::invalid_level_number;

    /*!
     * Cached Eulerian data to reduce the number of allocations/deallocations.
     */
    std::shared_ptr<SAMRAIDataCache> d_eulerian_data_cache;

    /*!
     * SAMRAI::hier::VariableContext object used for data management.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    /*!
     * SAMRAI::hier::Variable pointer and patch data descriptor indices for the
     * cell variable used to keep track of the count of the quadrature points in
     * each cell.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_qp_count_var;
    int d_qp_count_idx;

    /*!
     * The default parameters used during workload calculations.
     */
    const WorkloadSpec d_default_workload_spec;

    /*!
     * The default kernel functions and quadrature rule used to mediate
     * Lagrangian-Eulerian interaction.
     */
    const InterpSpec d_default_interp_spec;
    const SpreadSpec d_default_spread_spec;

    /*!
     * after reassociating patches with elements a node may still lie
     * outside all patches on the finest level in unusual circumstances
     * (like when the parent integrator class does not regrid sufficiently
     * frequently and has more than one patch level). This enum controls
     * what we do when this problem is detected.
     */
    NodeOutsidePatchCheckType d_node_patch_check = NODE_OUTSIDE_ERROR;

    /*!
     * SAMRAI::hier::IntVector object which determines the required ghost cell
     * width of this class.
     */
    const SAMRAI::hier::IntVector<NDIM> d_ghost_width;

    /*!
     * SAMRAI::hier::IntVector object which determines how many ghost cells we
     * should enlarge a patch by when associating an element with a patch. An
     * element is associated with a patch when its bounding box (defined as
     * the bounding box of both its nodes and quadrature points) intersects
     * the bounding box (including ghost cells) of that patch.
     *
     * @note At the present time this is always 1, which matches the
     * assumption made by IBTK::LEInteractor::getMinimumGhostWidth().
     */
    const SAMRAI::hier::IntVector<NDIM> d_associated_elem_ghost_width = SAMRAI::hier::IntVector<NDIM>(1);

    /*!
     * Data to manage mappings between mesh elements and grid patches.
     */
    std::vector<std::vector<std::vector<libMesh::Elem*> > > d_active_patch_elem_map;
    std::vector<std::vector<std::vector<libMesh::Node*> > > d_active_patch_node_map;
    std::map<std::string, std::vector<unsigned int> > d_active_patch_ghost_dofs;
    std::vector<libMesh::Elem*> d_active_elems;

    /*!
     * Ghost vectors for the various equation systems.
     */
    std::map<std::string, std::unique_ptr<libMesh::NumericVector<double> > > d_system_ghost_vec;

    /*!
     * Exemplar relevant IB-ghosted vectors for the various equation
     * systems. These vectors are cloned for fast initialization in
     * buildIBGhostedVector.
     */
    std::map<std::string, std::unique_ptr<libMesh::PetscVector<double> > > d_system_ib_ghost_vec;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBTK_HAVE_LIBMESH
#endif //#ifndef included_IBTK_FEDataManager
