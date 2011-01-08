// Filename: FEDataManager.h
// Created on 19 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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

#ifndef included_FEDataManager
#define included_FEDataManager

/////////////////////////////// INCLUDES /////////////////////////////////////

// LIBMESH INCLUDES
#define LIBMESH_REQUIRE_SEPARATE_NAMESPACE
#include <../base/variable.h>
#include <enum_order.h>
#include <enum_quadrature_type.h>
#include <equation_systems.h>
#include <linear_solver.h>
#include <sparse_matrix.h>

// PETSC INCLUDES
#include <petsc.h>

// IBTK INCLUDES
#include <ibtk/libmesh_utilities.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitStrategy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FEDataManager coordinates data required for Lagrangian-Eulerian
 * interaction between a Lagrangian finite element (FE) mesh.
 *
 * \note Multiple FEDataManager objects may be instantiated simultaneously.
 */
class FEDataManager
    : public SAMRAI::tbox::Serializable,
      public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
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
//  static const short int     NORMAL_DIRICHLET_BDRY_ID = 256;
//  static const short int TANGENTIAL_DIRICHLET_BDRY_ID = 512;
    static const short int            DIRICHLET_BDRY_ID = 256 | 512;

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
    getManager(
        const std::string& name,
        const std::string& interp_weighting_fcn,
        const std::string& spread_weighting_fcn,
        libMesh::QBase* const qrule,
        const bool interp_uses_consistent_mass_matrix,
        bool register_for_restart=true);

    /*!
     * Deallocate all of the FEDataManager instances.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void
    freeAllManagers();

    /*!
     * \name Methods to set the hierarchy and range of levels.
     */
    //\{

    /*!
     * \brief Reset patch hierarchy over which operations occur.
     */
    void
    setPatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Reset range of patch levels over which operations occur.
     *
     * The levels must exist in the hierarchy or an assertion failure will
     * result.
     */
    void
    resetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //\}

    /*!
     * \brief Set the equations systems object that is associated with the
     * FEDataManager.  Currently, each set of equation systems must be assigned
     * to a particular level of the AMR grid.
     */
    void
    setEquationSystems(
        libMesh::EquationSystems* const equation_systems,
        const int level_number);

    /*!
     * \return A pointer to the equations systems object that is associated with
     * the FEDataManager.
     */
    libMesh::EquationSystems*
    getEquationSystems() const;

    /*!
     * \return The level number to which the equations system object managed by
     * the FEDataManager is assigned.
     */
    int
    getLevelNumber() const;

    /*!
     * \return The ghost cell width used for quantities that are to be
     * interpolated from the Cartesian grid to the FE mesh.
     */
    const SAMRAI::hier::IntVector<NDIM>&
    getGhostCellWidth() const;

    /*!
     * \return The name of the weighting function used for interpolating from
     * the Cartesian grid to the FE mesh.
     */
    const std::string&
    getInterpWeightingFunction() const;

    /*!
     * \return The name of the weighting function used for spreading densities
     * from the FE mesh to the Cartesian grid.
     */
    const std::string&
    getSpreadWeightingFunction() const;

    /*!
     * \return A pointer to the quadrature rule used to construct the discrete
     * Lagrangian-Eulerian interation operators.
     */
    libMesh::QBase*
    getQuadratureRule() const;

    /*!
     * \return A boolean value indicating whether the interpolation operator is
     * defined in terms of a consistent mass matrix.
     */
    bool
    getInterpUsesConsistentMassMatrix() const;

    /*!
     * \return A const reference to the map from patch number to local elements.
     */
    const std::map<int,std::set<libMesh::Elem*> >&
    getActivePatchElements();

    /*!
     * \brief Reinitialize the mappings from elements to Cartesian grid patches.
     */
    void
    reinitElementMappings();

    /*!
     * \return A pointer to the unghosted solution vector associated with the
     * specified system.
     */
    libMesh::NumericVector<double>*
    getSolutionVector(
        const std::string& system_name);

    /*!
     * \return A pointer to the ghosted solution vector associated with the
     * specified system.
     */
    libMesh::NumericVector<double>*
    getGhostedSolutionVector(
        const std::string& system_name);

    /*!
     * \return A pointer to the unghosted coordinates (nodal position) vector.
     */
    libMesh::NumericVector<double>*
    getCoordsVector();

    /*!
     * \return A pointer to the ghosted coordinates (nodal position) vector.
     */
    libMesh::NumericVector<double>*
    getGhostedCoordsVector();

    /*!
     * \brief Spread a density from the FE mesh to the Cartesian grid.
     */
    void
    spread(
        const int f_data_idx,
        libMesh::NumericVector<double>& F,
        libMesh::NumericVector<double>& X,
        const std::string& system_name,
        const bool close_F=true,
        const bool close_X=true);

    /*!
     * \brief Interpolate a value from the Cartesian grid to the FE mesh.
     */
    void
    interp(
        const int f_data_idx,
        libMesh::NumericVector<double>& F,
        libMesh::NumericVector<double>& X,
        const std::string& system_name,
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > f_refine_scheds=std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
        const double fill_data_time=0.0,
        const bool close_X=true);

    /*!
     * \return Pointers to a linear solver and sparse matrix corresponding to a
     * L2 projection.
     */
    std::pair<libMesh::LinearSolver<double>*,libMesh::SparseMatrix<double>*>
    getL2ProjectionSolver(
        const std::string& system_name,
        const bool consistent_mass_matrix=true,
        const libMeshEnums::QuadratureType quad_type=QGAUSS,
        const libMeshEnums::Order quad_order=FIFTH);

    ///
    ///  The following routines:
    ///
    ///      initializeLevelData(),
    ///      resetHierarchyConfiguration(),
    ///      applyGradientDetector()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::mesh::StandardTagAndInitStrategy abstract base class.
    ///

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
    void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level=SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

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
    void
    resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_ln,
        const int finest_ln);

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
    virtual void
    applyGradientDetector(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);

    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * \brief Constructor.
     */
    FEDataManager(
        const std::string& object_name,
        const std::string& interp_weighting_fcn,
        const std::string& spread_weighting_fcn,
        libMesh::QBase* const qrule,
        const bool interp_uses_consistent_mass_matrix,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width,
        bool register_for_restart=true);

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
    FEDataManager(
        const FEDataManager& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FEDataManager&
    operator=(
        const FEDataManager& that);

    /*!
     * Compute the quadrature point counts in each cell of the level in which
     * the FE mesh is embedded.  Also zeros out node count data for other levels
     * within the specified range of level numbers.
     */
    void
    updateQuadPointCountData(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * Compute the bounding boxes of all active elements.
     *
     * \note For inactive elements, the lower and upper bound values will be
     * identically zero.
     */
    void
    computeActiveElementBoundingBoxes(
        std::vector<double>& elem_bounds);

    /*!
     * Collect all of the active elements which are located within a local
     * Cartesian grid patch grown by the specified ghost cell width.
     *
     * In this method, the determination as to whether an element is local or
     * not is based on the position of the bounding box of the element.
     */
    void
    collectActivePatchElements(
        std::map<int,std::set<libMesh::Elem*> >& active_patch_elems,
        const int level_number,
        const SAMRAI::hier::IntVector<NDIM>& ghost_width);

    /*!
     * Collect all ghost DOF indices for the specified collection of elements.
     */
    void
    collectGhostDOFIndices(
        std::vector<unsigned int>& ghost_dofs,
        const std::map<int,std::set<libMesh::Elem*> >& active_patch_elems,
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
    void
    getFromRestart();

    /*!
     * Static data members used to control access to and destruction of
     * singleton data manager instance.
     */
    static std::map<std::string,FEDataManager*> s_data_manager_instances;
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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_qp_count_var;
    int d_qp_count_idx;

    /*
     * The weighting functions and quadrature rule used to mediate
     * Lagrangian-Eulerian interaction.
     */
    const std::string d_interp_weighting_fcn;
    const std::string d_spread_weighting_fcn;
    libMesh::QBase* const d_qrule;
    const bool d_interp_uses_consistent_mass_matrix;

    /*
     * SAMRAI::hier::IntVector object which determines the ghost cell width of
     * the LNodeIndexData SAMRAI::hier::PatchData objects.
     */
    const SAMRAI::hier::IntVector<NDIM> d_ghost_width;

    /*
     * FE equation system associated with this data manager object.
     */
    libMesh::EquationSystems* d_es;
    int d_level_number;

    /*
     * Data to manage mappings between mesh elements and grid patches.
     */
    std::map<int,std::set<libMesh::Elem*> > d_active_patch_elems;
    std::map<std::string,std::vector<unsigned int> > d_active_patch_ghost_dofs;

    /*
     * Ghost vectors for the various equation systems.
     */
    std::map<std::string,libMesh::NumericVector<double>*> d_system_ghost_vec;

    /*
     * Linear solvers and related data for performing interpolation in the IB-FE
     * framework.
     */
    std::map<std::string,libMesh::LinearSolver<double>*> d_L2_projection_solvers;
    std::map<std::string,libMesh::SparseMatrix<double>*> d_L2_mass_matrices;
    std::map<std::string,bool> d_L2_consistent_mass_matrix;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/FEDataManager.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FEDataManager
