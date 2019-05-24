// Filename: FEDataManager.h
// Created on 19 Apr 2010 by Boyce Griffith
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

#ifndef included_IBTK_FEDataManager
#define included_IBTK_FEDataManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <ostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BasePatchLevel.h"
#include "CellVariable.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "RefineSchedule.h"
#include "SideVariable.h"
#include "StandardTagAndInitStrategy.h"
#include "VariableContext.h"
#include "boost/multi_array.hpp"
#include "ibtk/QuadratureCache.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/system.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

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
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
namespace libMesh
{
class Elem;
class EquationSystems;
class QBase;
template <typename T>
class LinearSolver;
template <typename T>
class NumericVector;
template <typename T>
class SparseMatrix;
} // namespace libMesh

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FEDataManager coordinates data required for Lagrangian-Eulerian
 * interaction between a Lagrangian finite element (FE) mesh.
 *
 * <h3>Parameters effecting workload estimate calculations</h3>
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
class FEDataManager : public SAMRAI::tbox::Serializable, public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    /*!
     * Class which enables fast lookup of dofs on a given libMesh <code>elem</code>.
     *
     * @note The contents of this cache are invalidated when we regrid and the
     * caches should be reset at that point. Copies of this class should
     * always be retrieved via FEDataManager::getDofCache() to avoid this
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
            dof_indices = this->dof_indices(elem)[var];
            return;
        }

        /*!
         * Alternative indexing operation: retrieve all dof indices of all
         * variables in the given system at once by reference.
         */
        inline const std::vector<std::vector<libMesh::dof_id_type> >& dof_indices(const libMesh::Elem* const elem)
        {
            std::vector<std::vector<libMesh::dof_id_type> >& elem_dof_indices = d_dof_cache[elem->id()];
            if (elem_dof_indices.empty())
            {
                elem_dof_indices.resize(d_dof_map.n_variables());
                for (unsigned int var_n = 0; var_n < d_dof_map.n_variables(); ++var_n)
                {
                    d_dof_map.dof_indices(elem, elem_dof_indices[var_n], var_n);
                }
            }
            return elem_dof_indices;
        }

    private:
        libMesh::DofMap& d_dof_map;
        std::unordered_map<libMesh::dof_id_type, std::vector<std::vector<unsigned int> > > d_dof_cache;
    };

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
                   bool use_nodal_quadrature)
            : kernel_fcn(kernel_fcn),
              quad_type(quad_type),
              quad_order(quad_order),
              use_adaptive_quadrature(use_adaptive_quadrature),
              point_density(point_density),
              use_consistent_mass_matrix(use_consistent_mass_matrix),
              use_nodal_quadrature(use_nodal_quadrature)
        {
        }

        std::string kernel_fcn;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
        bool use_adaptive_quadrature;
        double point_density;
        bool use_consistent_mass_matrix;
        bool use_nodal_quadrature;
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
                   bool use_nodal_quadrature)
            : kernel_fcn(kernel_fcn),
              quad_type(quad_type),
              quad_order(quad_order),
              use_adaptive_quadrature(use_adaptive_quadrature),
              point_density(point_density),
              use_nodal_quadrature(use_nodal_quadrature)
        {
        }

        std::string kernel_fcn;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
        bool use_adaptive_quadrature;
        double point_density;
        bool use_nodal_quadrature;
    };

    /*!
     * \brief Struct WorkloadSpec encapsulates the parameters used to
     * calculate workload estimates (i.e., the input to the load balancing
     * algorithm).
     */
    struct WorkloadSpec
    {
        /// The multiplier applied to each quadrature point.
        double q_point_weight = 2.0;
    };

    /*!
     * \brief The name of the equation system which stores the spatial position
     * data.
     *
     * \note The default value for this string is "coordinates system".
     */
    std::string COORDINATES_SYSTEM_NAME = "coordinates system";

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
     * Note that when a manager is accessed for the first time, the
     * freeAllManagers static method is registered with the ShutdownRegistry
     * class.  Consequently, all allocated managers are freed at program
     * completion.  Thus, users of this class do not explicitly allocate or
     * deallocate the FEDataManager instances.
     *
     * \return A pointer to the data manager instance.
     */
    static FEDataManager*
    getManager(const std::string& name,
               const InterpSpec& default_interp_spec,
               const SpreadSpec& default_spread_spec,
               const WorkloadSpec& default_workload_spec,
               const SAMRAI::hier::IntVector<NDIM>& min_ghost_width = SAMRAI::hier::IntVector<NDIM>(0),
               bool register_for_restart = true);

    /*!
     * Same as the last function, but uses a default workload specification
     * for compatibility.
     *
     * \return A pointer to the data manager instance.
     */
    static FEDataManager*
    getManager(const std::string& name,
               const InterpSpec& default_interp_spec,
               const SpreadSpec& default_spread_spec,
               const SAMRAI::hier::IntVector<NDIM>& min_ghost_width = SAMRAI::hier::IntVector<NDIM>(0),
               bool register_for_restart = true);

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
     * Deallocate all of the FEDataManager instances.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeAllManagers();

    /*!
     * \brief Register a load balancer for non-uniform load balancing.
     *
     * @deprecated Since the correct workload index is passed in via
     * addWorkloadEstimate this function is no longer necessary.
     */
    void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
                              int workload_data_idx);

    /*!
     * \name Methods to set and get the patch hierarchy and range of patch
     * levels associated with this manager class.
     */
    //\{

    /*!
     * \brief Reset patch hierarchy over which operations occur.
     */
    void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Get the patch hierarchy used by this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > getPatchHierarchy() const;

    /*!
     * \brief Reset range of patch levels over which operations occur.
     *
     * The levels must exist in the hierarchy or an assertion failure will
     * result.
     */
    void setPatchLevels(int coarsest_ln, int finest_ln);

    /*!
     * \brief Get the range of patch levels used by this object.
     *
     * \note Returns [coarsest_ln,finest_ln+1).
     */
    std::pair<int, int> getPatchLevels() const;

    //\}

    /*!
     * \brief Set the equations systems object that is associated with the
     * FEDataManager.  Currently, each set of equation systems must be assigned
     * to a particular level of the AMR grid.
     */
    void setEquationSystems(libMesh::EquationSystems* equation_systems, int level_number);

    /*!
     * \return A pointer to the equations systems object that is associated with
     * the FEDataManager.
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
     * \return The level number to which the equations system object managed by
     * the FEDataManager is assigned.
     */
    int getLevelNumber() const;

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
     * @deprecated Use buildIBGhostedVector instead which clones a vector with
     * the same ghost region.
     */
    libMesh::NumericVector<double>* buildGhostedSolutionVector(const std::string& system_name,
                                                               bool localize_data = true);

    /*!
     * \return A pointer to a vector, with ghost entries corresponding to
     * relevant IB data, associated with the specified system.
     */
    std::unique_ptr<libMesh::PetscVector<double> > buildIBGhostedVector(const std::string& system_name) const;

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
     * \brief Spread a density from the FE mesh to the Cartesian grid using the
     * default spreading spec.
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
     * specified spreading spec.
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
     * \return Pointers to a linear solver and sparse matrix corresponding to a
     * L2 projection operator.
     */
    std::pair<libMesh::LinearSolver<double>*, libMesh::SparseMatrix<double>*>
    buildL2ProjectionSolver(const std::string& system_name);

    /*!
     * \return Pointer to vector representation of diagonal L2 mass matrix.
     */
    libMesh::NumericVector<double>* buildDiagonalL2MassMatrix(const std::string& system_name);

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
     * Update the quarature rule for the current element.  If the provided
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
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.  The level number indicates that of
     * the new level.  The old_level pointer corresponds to the level that
     * resided in the hierarchy before the level with the specified number was
     * introduced.  If the pointer is null, there was no level in the hierarchy
     * prior to the call and the level data is set based on the user routines
     * and the simulation time.  Otherwise, the specified level replaces the old
     * level and the new level receives data from the old level appropriately
     * before it is destroyed.
     *
     * The boolean argument initial_time indicates whether the level is being
     * introduced for the first time (i.e., at initialization time) or after
     * some regrid process during the calculation beyond the initial hierarchy
     * construction.  This information is provided since the initialization of
     * the data on a patch may be different in each of those circumstances.  The
     * can_be_refined boolean argument indicates whether the level is the finest
     * level allowed in the hierarchy.  This may or may not affect the data
     * initialization process depending on the problem.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, the level number does not match any
     * level in the hierarchy, or the old level number does not match the level
     * number (if the old level pointer is non-null).
     */
    void initializeLevelData(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                             int level_number,
                             double init_data_time,
                             bool can_be_refined,
                             bool initial_time,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level =
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
                             bool allocate_data = true) override;

    /*!
     * Reset cached communication schedules after the hierarchy has changed (for
     * example, due to regridding) and the data has been initialized on the new
     * levels.  The intent is that the cost of data movement on the hierarchy
     * will be amortized across multiple communication cycles, if possible.  The
     * level numbers indicate the range of levels in the hierarchy that have
     * changed.  However, this routine updates communication schedules every
     * level finer than and including that indexed by the coarsest level number
     * given.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, any pointer to a level in the hierarchy
     * that is coarser than the finest level is null, or the given level numbers
     * not specified properly; e.g., coarsest_ln > finest_ln.
     */
    void resetHierarchyConfiguration(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                     int coarsest_ln,
                                     int finest_ln) override;

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
     */
    void applyGradientDetector(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                               int level_number,
                               double error_data_time,
                               int tag_index,
                               bool initial_time,
                               bool uses_richardson_extrapolation_too) override;

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * \brief Constructor.
     */
    FEDataManager(std::string object_name,
                  InterpSpec default_interp_spec,
                  SpreadSpec default_spread_spec,
                  WorkloadSpec default_workload_spec,
                  SAMRAI::hier::IntVector<NDIM> ghost_width,
                  bool register_for_restart = true);

    /*!
     * \brief The FEDataManager destructor cleans up any allocated data objects.
     */
    ~FEDataManager() = default;

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
     * Compute the bounding boxes of all active elements.
     *
     * \note For inactive elements, the lower and upper bound values will be
     * identically zero.
     */
    std::vector<std::pair<Point, Point> >* computeActiveElementBoundingBoxes();

    /*!
     * Collect all of the active elements which are located within a local
     * Cartesian grid patch grown by the specified ghost cell width.
     *
     * In this method, the determination as to whether an element is local or
     * not is based on the position of the bounding box of the element.
     */
    void collectActivePatchElements(std::vector<std::vector<libMesh::Elem*> >& active_patch_elems,
                                    int level_number,
                                    const SAMRAI::hier::IntVector<NDIM>& ghost_width);

    /*!
     * Collect all of the nodes of the active elements that are located within a
     * local Cartesian grid patch grown by the specified ghost cell width.
     */
    void collectActivePatchNodes(std::vector<std::vector<libMesh::Node*> >& active_patch_nodes,
                                 const std::vector<std::vector<libMesh::Elem*> >& active_patch_elems);

    /*!
     * Collect all ghost DOF indices for the specified collection of elements.
     */
    void collectGhostDOFIndices(std::vector<unsigned int>& ghost_dofs,
                                const std::vector<libMesh::Elem*>& active_elems,
                                const std::string& system_name);

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
    void getFromRestart();

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static std::map<std::string, FEDataManager*> s_data_manager_instances;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * Whether or not to log data to the screen: see
     * FEDataManager::setLoggingEnabled() and
     * FEDataManager::getLoggingEnabled().
     *
     * @note This is usually set by IBFEMethod, which reads the relevant
     * boolean from the database.
     */
    bool d_enable_logging = false;

    /*
     * We cache a pointer to the load balancer.
     *
     * @deprecated This pointer is never used and will be removed in the
     * future.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /*
     * Grid hierarchy information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    /*
     * Cached Eulerian data to reduce the number of allocations/deallocations.
     */
    SAMRAIDataCache d_cached_eulerian_data;

    /*
     * SAMRAI::hier::VariableContext object used for data management.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    /*
     * SAMRAI::hier::Variable pointer and patch data descriptor indices for the
     * cell variable used to keep track of the count of the quadrature points in
     * each cell.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_qp_count_var;
    int d_qp_count_idx;

    /*
     * SAMRAI::xfer::RefineAlgorithm pointer to fill the ghost cell region of
     * SAMRAI variables.
     */
    SAMRAI::xfer::RefineAlgorithm<NDIM> d_ghost_fill_alg;

    /*
     * SAMRAI::hier::Variable pointer and patch data descriptor indices for the
     * cell variable used to determine the workload for nonuniform load
     * balancing.
     *
     * @deprecated d_workload_var and d_workload_idx will not be stored in a
     * future release since the correct workload index will be provided as an
     * argument to addWorkloadEstimate.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_workload_var;
    int d_workload_idx = IBTK::invalid_index;

    /*
     * The default parameters used during workload calculations.
     */
    const WorkloadSpec d_default_workload_spec;

    /*
     * The default kernel functions and quadrature rule used to mediate
     * Lagrangian-Eulerian interaction.
     */
    const InterpSpec d_default_interp_spec;
    const SpreadSpec d_default_spread_spec;

    /*
     * SAMRAI::hier::IntVector object which determines the ghost cell width used
     * to determine elements that are associated with each Cartesian grid patch.
     */
    const SAMRAI::hier::IntVector<NDIM> d_ghost_width;

    /*
     * FE equation system associated with this data manager object.
     */
    libMesh::EquationSystems* d_es = nullptr;
    int d_level_number = IBTK::invalid_level_number;
    std::map<unsigned int, std::unique_ptr<SystemDofMapCache> > d_system_dof_map_cache;

    /*
     * Data to manage mappings between mesh elements and grid patches.
     */
    std::vector<std::vector<libMesh::Elem*> > d_active_patch_elem_map;
    std::vector<std::vector<libMesh::Node*> > d_active_patch_node_map;
    std::map<std::string, std::vector<unsigned int> > d_active_patch_ghost_dofs;
    std::vector<std::pair<Point, Point> > d_active_elem_bboxes;

    /*
     * Cache of libMesh quadrature objects. Defaults to being an NDIM cache
     * but is overwritten once a libMesh::EquationSystems object is attached.
     */
    QuadratureCache d_quadrature_cache;

    /*
     * Ghost vectors for the various equation systems.
     */
    std::map<std::string, std::unique_ptr<libMesh::NumericVector<double> > > d_system_ghost_vec;

    /*
     * Exemplar relevant IB-ghosted vectors for the various equation
     * systems. These vectors are cloned for fast initialization in
     * buildIBGhostedVector.
     */
    std::map<std::string, std::unique_ptr<libMesh::PetscVector<double> > > d_system_ib_ghost_vec;

    /*
     * Linear solvers and related data for performing interpolation in the IB-FE
     * framework.
     */
    std::map<std::string, std::unique_ptr<libMesh::LinearSolver<double> > > d_L2_proj_solver;
    std::map<std::string, std::unique_ptr<libMesh::SparseMatrix<double> > > d_L2_proj_matrix;
    std::map<std::string, std::unique_ptr<libMesh::NumericVector<double> > > d_L2_proj_matrix_diag;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FEDataManager
