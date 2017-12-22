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

#include <stdbool.h>
#include <stddef.h>
#include <map>
#include <ostream>
#include <string>
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
#include "boost/unordered_map.hpp"
#include "ibtk/ibtk_utilities.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
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
 * \note Multiple FEDataManager objects may be instantiated simultaneously.
 */
class FEDataManager : public SAMRAI::tbox::Serializable, public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    class SystemDofMapCache
    {
    public:
        inline SystemDofMapCache(libMesh::System& system) : d_dof_map(system.get_dof_map())
        {
        }

        inline ~SystemDofMapCache()
        {
        }

        inline void
        dof_indices(const libMesh::Elem* const elem, std::vector<unsigned int>& dof_indices, const unsigned int var = 0)
        {
            const libMesh::dof_id_type elem_id = elem->id();
            std::vector<std::vector<unsigned int> >& elem_dof_indices = d_dof_cache[elem_id];
            if (elem_dof_indices.size() <= var)
            {
                elem_dof_indices.resize(var + 1);
            }
            if (elem_dof_indices[var].empty())
            {
                d_dof_map.dof_indices(elem, elem_dof_indices[var], var);
            }
            dof_indices = elem_dof_indices[var];
            return;
        }

    private:
        libMesh::DofMap& d_dof_map;
        boost::unordered_map<libMesh::dof_id_type, std::vector<std::vector<unsigned int> > > d_dof_cache;
    };

    /*!
     * \brief Struct InterpSpec encapsulates data needed to specify the manner
     * in which Eulerian-to-Lagrangian interpolation is performed when using an
     * FE structural discretization.
     */
    struct InterpSpec
    {
        InterpSpec()
        {
        }

        InterpSpec(const std::string& kernel_fcn,
                   const libMesh::QuadratureType& quad_type,
                   const libMesh::Order& quad_order,
                   bool use_adaptive_quadrature,
                   double point_density,
                   bool use_consistent_mass_matrix)
            : kernel_fcn(kernel_fcn),
              quad_type(quad_type),
              quad_order(quad_order),
              use_adaptive_quadrature(use_adaptive_quadrature),
              point_density(point_density),
              use_consistent_mass_matrix(use_consistent_mass_matrix)
        {
        }

        std::string kernel_fcn;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
        bool use_adaptive_quadrature;
        double point_density;
        bool use_consistent_mass_matrix;
    };

    /*!
     * \brief Struct SpreadSpec encapsulates data needed to specify the manner
     * in which Lagrangian-to-Eulerian spreading is performed when using an FE
     * structural discretization.
     */
    struct SpreadSpec
    {
        SpreadSpec()
        {
        }

        SpreadSpec(const std::string& kernel_fcn,
                   const libMesh::QuadratureType& quad_type,
                   const libMesh::Order& quad_order,
                   bool use_adaptive_quadrature,
                   double point_density)
            : kernel_fcn(kernel_fcn),
              quad_type(quad_type),
              quad_order(quad_order),
              use_adaptive_quadrature(use_adaptive_quadrature),
              point_density(point_density)
        {
        }

        std::string kernel_fcn;
        libMesh::QuadratureType quad_type;
        libMesh::Order quad_order;
        bool use_adaptive_quadrature;
        double point_density;
    };

    /*!
     * \brief The name of the equation system which stores the spatial position
     * data.
     *
     * \note The default value for this string is "coordinates system".
     */
    std::string COORDINATES_SYSTEM_NAME;

    /*!
     * \brief The libMesh boundary IDs to use for specifying essential boundary
     * conditions.
     */
    static const short int ZERO_DISPLACEMENT_X_BDRY_ID;
    static const short int ZERO_DISPLACEMENT_Y_BDRY_ID;
    static const short int ZERO_DISPLACEMENT_Z_BDRY_ID;
    static const short int ZERO_DISPLACEMENT_XY_BDRY_ID;
    static const short int ZERO_DISPLACEMENT_XZ_BDRY_ID;
    static const short int ZERO_DISPLACEMENT_YZ_BDRY_ID;
    static const short int ZERO_DISPLACEMENT_XYZ_BDRY_ID;

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
               const SAMRAI::hier::IntVector<NDIM>& min_ghost_width = SAMRAI::hier::IntVector<NDIM>(0),
               bool register_for_restart = true);

    /*!
     * Deallocate all of the FEDataManager instances.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeAllManagers();

    /*!
     * \brief Register a load balancer for non-uniform load balancing.
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
     * specified system.
     */
    libMesh::NumericVector<double>* buildGhostedSolutionVector(const std::string& system_name,
                                                               bool localize_data = true);

    /*!
     * \return A pointer to the unghosted coordinates (nodal position) vector.
     */
    libMesh::NumericVector<double>* getCoordsVector() const;

    /*!
     * \return A pointer to the ghosted coordinates (nodal position) vector.
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
                double fill_data_time);

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
                double fill_data_time);

    /*!
     * \brief Prolong a value or a density from the FE mesh to the Cartesian
     * grid.
     */
    void prolongData(int f_data_idx,
                     libMesh::NumericVector<double>& F,
                     libMesh::NumericVector<double>& X,
                     const std::string& system_name,
                     bool is_density = true,
                     bool accumulate_on_grid = true);

    /*!
     * \brief Interpolate a value from the Cartesian grid to the FE mesh using
     * the default interpolation spec. This interpolation function does NOT do
     * an L2-projection of the interpolated quantity. It does however weighs/filters
     * the interpolated quantity at the quadrature points to the nodes. Here, the
     * basis functions of the deformational field is used as the filter.
     */
    void
    interpWeighted(int f_data_idx,
                   libMesh::NumericVector<double>& F,
                   libMesh::NumericVector<double>& X,
                   const std::string& system_name,
                   const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_refine_scheds =
                       std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                   double fill_data_time = 0.0);

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
                   double fill_data_time = 0.0);

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
                double fill_data_time = 0.0);

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
                double fill_data_time = 0.0);

    /*!
     * \brief Restrict a value from the Cartesian grid to the FE mesh.
     */
    void restrictData(int f_data_idx,
                      libMesh::NumericVector<double>& F,
                      libMesh::NumericVector<double>& X,
                      const std::string& system_name,
                      bool use_consistent_mass_matrix = true);

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
    static bool updateQuadratureRule(libMesh::UniquePtr<libMesh::QBase>& qrule,
                                     libMesh::QuadratureType quad_type,
                                     libMesh::Order quad_order,
                                     bool use_adaptive_quadrature,
                                     double point_density,
                                     libMesh::Elem* elem,
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
    static bool updateInterpQuadratureRule(libMesh::UniquePtr<libMesh::QBase>& qrule,
                                           const InterpSpec& spec,
                                           libMesh::Elem* elem,
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
    static bool updateSpreadQuadratureRule(libMesh::UniquePtr<libMesh::QBase>& qrule,
                                           const SpreadSpec& spec,
                                           libMesh::Elem* elem,
                                           const boost::multi_array<double, 2>& X_node,
                                           double dx_min);

    /*!
     * \brief Update the cell workload estimate.
     */
    void updateWorkloadEstimates(int coarsest_ln = -1, int finest_ln = -1);

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
                             bool allocate_data = true);

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
                                     int finest_ln);

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
                               bool uses_richardson_extrapolation_too);

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * \brief Constructor.
     */
    FEDataManager(const std::string& object_name,
                  const InterpSpec& default_interp_spec,
                  const SpreadSpec& default_spread_spec,
                  const SAMRAI::hier::IntVector<NDIM>& ghost_width,
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
    FEDataManager();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FEDataManager(const FEDataManager& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FEDataManager& operator=(const FEDataManager& that);

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
     * We cache a pointer to the load balancer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /*
     * Grid hierarchy information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

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
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_workload_var;
    int d_workload_idx;

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
    libMesh::EquationSystems* d_es;
    int d_level_number;
    std::map<unsigned int, SAMRAI::tbox::Pointer<SystemDofMapCache> > d_system_dof_map_cache;

    /*
     * Data to manage mappings between mesh elements and grid patches.
     */
    std::vector<std::vector<libMesh::Elem*> > d_active_patch_elem_map;
    std::map<std::string, std::vector<unsigned int> > d_active_patch_ghost_dofs;
    std::vector<std::pair<Point, Point> > d_active_elem_bboxes;

    /*
     * Ghost vectors for the various equation systems.
     */
    std::map<std::string, libMesh::NumericVector<double>*> d_system_ghost_vec;

    /*
     * Linear solvers and related data for performing interpolation in the IB-FE
     * framework.
     */
    std::map<std::string, libMesh::LinearSolver<double>*> d_L2_proj_solver;
    std::map<std::string, libMesh::SparseMatrix<double>*> d_L2_proj_matrix;
    std::map<std::string, libMesh::NumericVector<double>*> d_L2_proj_matrix_diag;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FEDataManager
