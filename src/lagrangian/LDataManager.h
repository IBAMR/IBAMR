#ifndef included_LDataManager
#define included_LDataManager

// Filename: LDataManager.h
// Created on 01 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <07.Nov.2006 00:58:00 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeIndexVariable.h>
#include <ibamr/LNodeJacobianInitStrategy.h>
#include <ibamr/LNodePosnInitStrategy.h>

// SAMRAI INCLUDES
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CellIndex.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <ComponentSelector.h>
#include <Index.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitStrategy.h>
#include <VariableContext.h>
#include <VisItDataWriter.h>
#include <tbox/Arena.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// PETSC INCLUDES
#include <petscvec.h>
#include <petscao.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBAMR
{
class LNodeIndexSet;
class LNodeLevelData;
class LagSiloDataWriter;
}// namespace IBAMR

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * Class LDataManager coordinates the irregular distribution of
 * LNodeIndexData and LNodeLevelData on the patch hierarchy.  The
 * manager class is responsible for maintaining this data distribution
 * and for all interprocessor communications.
 */
class LDataManager
    : public SAMRAI::tbox::Serializable,
      public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    /*!
     * The name of the LNodeLevelData that specifies the location of
     * the curvilinear mesh nodes.
     */
    static const string COORDS_DATA_NAME;

    /*!
     * The name of the LNodeLevelData that specifies the Jacobian
     * determinant of the curvilinear mesh nodes.
     */
    static const string JACOBIAN_DATA_NAME;

    /*!
     * Return a pointer to the instance of the Lagrangian data manager
     * corresponding to the specified name.  Access to LDataManager
     * objects is mediated by the getManager() function.
     *
     * Note that when a manager is accessed for the first time, the
     * freeAllManagers static method is registered with the
     * ShutdownRegistry class.  Consequently, all allocated managers
     * are freed at program completion.  Thus, users of this class do
     * not explicitly allocate or deallocate the LDataManager
     * instances.
     *
     * @return A pointer to the data manager instance.
     *
     * @todo fix this interface
     */
    static LDataManager* getManager(
        const string& name,
        const SAMRAI::hier::IntVector<NDIM>& ghosts=3,
        bool register_for_restart=true);

    /*!
     * Deallocate all of the LDataManager instances.
     *
     * It is not necessary to call this function at program
     * termination since it is automatically called by the
     * ShutdownRegistry class.
     */
    static void freeAllManagers();

    //@{ @name Methods to set the hierarchy and range of levels.

    /*!
     * @brief Reset patch hierarchy over which operations occur.
     */
    void setPatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * @brief Reset range of patch levels over which operations occur.
     *
     * The levels must exist in the hierarchy or an assertion failure
     * will result.
     */
    void resetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //@}

    /*!
     * Register a concrete strategy object with the integrator that
     * specifies the Jacobian determinant at the curvilinear mesh
     * nodes.
     */
    void registerLNodeJacobianInitStrategy(
        SAMRAI::tbox::Pointer<LNodeJacobianInitStrategy> lag_jac_init);

    /*!
     * Free the concrete Jacobian strategy object.
     *
     * NOTE: Be sure to call this method only once the initialization
     * object is no longer needed.
     */
    void freeLNodeJacobianInitStrategy();

    /*!
     * Register a concrete strategy object with the integrator that
     * specifies the initial configuration of the curvilinear mesh
     * nodes.
     */
    void registerLNodePosnInitStrategy(
        SAMRAI::tbox::Pointer<LNodePosnInitStrategy> lag_posn_init);

    /*!
     * Free the concrete initialization strategy object.
     *
     * NOTE: Be sure to call this method only once the initialization
     * object is no longer needed.
     */
    void freeLNodePosnInitStrategy();

    /*!
     * @brief Register a VisIt data writer with the manager.
     */
    void registerVisItDataWriter(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * @brief Register a Silo data writer with the manager.
     */
    void registerLagSiloDataWriter(
        SAMRAI::tbox::Pointer<LagSiloDataWriter> silo_writer);

    /*!
     * @brief Register a load balancer for non-uniform load balancing.
     */
    void registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer);

    /*!
     * @brief Indicates whether there is Lagrangian data on the given
     * patch hierarchy level.
     */
    bool levelContainsLagrangianData(
        const int level_number) const;

    /*!
     * @return The number of total nodes of the Lagrangian data for
     * the specified level of the patch hierarchy.
     */
    int getNumberOfNodes(
        const int level_number) const;

    /*!
     * @return The number of local (i.e., on processor) nodes of the
     * Lagrangian data for the specified level of the patch hierarchy.
     *
     * NOTE: This count does not include nodes which only lie in ghost
     * cells for the current process.
     */
    int getNumberOfLocalNodes(
        const int level_number) const;

    /*!
     * @return The number of nodes on all processors with MPI rank
     * less than the current process on the specified level of the
     * patch hierarchy.
     *
     * NOTE: This count does not include nodes which only lie in ghost
     * cells for the current process.
     */
    int getGlobalNodeOffset(
        const int level_number) const;

    /*!
     * @brief Get the specified Lagrangian quantity data on the given
     * patch hierarchy level.
     */
    SAMRAI::tbox::Pointer<LNodeLevelData> getLNodeLevelData(
        const string& quantity_name,
        const int level_number);

    /*!
     * @brief Allocate new Lagrangian level data with the specified
     * name and depth.  If specified, the quantity is maintained as
     * the patch hierarchy evolves.
     *
     * NOTE: Quantities maintained by the LDataManager must have
     * unique names.  The name "X" is reserved for the nodal
     * coordinates.
     */
    SAMRAI::tbox::Pointer<LNodeLevelData> createLNodeLevelData(
        const string& quantity_name,
        const int level_number,
        const int depth=1,
        const bool maintain_data=false);

    /*!
     * @brief Get the patch data descriptor index for the Lagrangian
     * index data.
     */
    int getLNodeIndexPatchDescriptorIndex() const;

    /*!
     * @brief Get the patch data descriptor index for the workload
     * cell data.
     */
    int getWorkloadPatchDescriptorIndex() const;

    /*!
     * @brief Get the patch data descriptor index for the Lagrangian
     * node count cell data.
     */
    int getNodeCountPatchDescriptorIndex() const;

    /*!
     * @brief Get the patch data descriptor index for the MPI process
     * mapping cell data.
     */
    int getProcMappingPatchDescriptorIndex() const;

    /*!
     * @brief Map the collection of Lagrangian indices to the
     * corresponding global PETSc indices.
     */
    void mapLagrangianToPETSc(
        std::vector<int>& inds,
        const int level_number) const;

    /*!
     * @brief Map the collection of global PETSc indices to the
     * corresponding Lagrangian indices.
     */
    void mapPETScToLagrangian(
        std::vector<int>& inds,
        const int level_number) const;

    /*!
     * @brief Start the process of redistributing the Lagrangian data.
     *
     * This method uses the present location of each Lagrangian mesh
     * node to redistribue the LNodeIndexData managed by this object.
     *
     * IMPORTANT NOTE: This routine assumes that the time interval
     * between node redistribution satisfies a timestep restriction of
     * the form dt <= C*dx*|U| with C <= 1.  This restriction prevents
     * nodes from moving more than one cell width per timestep.
     *
     * @see endDataRedistribution
     */
    void beginDataRedistribution(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * @brief Finish the process of redistributing the Lagrangian
     * data.
     *
     * This method redistribues the quantites associated with each
     * node in the Lagrangian mesh according to the data distribution
     * defined by the LNodeIndexData managed by this object.  This
     * routine potentially involves SUBSTANTIAL interprocessor
     * communication.
     *
     * IMPORTANT NOTE: Since this routine potentially results in a
     * large amount of interprocessor communication, it may be worth
     * putting it off for as long as possible.  If the timestep dt
     * satisfies a condition of the form dt <= C*dx*|U| with C << 1,
     * it may be possible to redistribute the Lagrangian data less
     * frequently than every timestep.
     *
     * @see beginDataRedistribution
     */
    void endDataRedistribution(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * @brief Update the workload and count of nodes per cell.
     *
     * This routine updates cell data which is maintained on the patch
     * hierarchy in order to keep track of the number of nodes in each
     * cell of the AMR index space.  The node count data is used to
     * tag cells for refinement, and to specify non-uniform load
     * balancing.  The workload per cell is defined by
     *
     *    workload(i) = alpha_work + beta_work*node_count(i)
     *
     * where alpha and beta are parameters that each default to the
     * value 1.
     */
    void updateWorkloadAndNodeCount(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * @brief Each LNodeIndex object owns a pointer to its nodal
     * location.  This routine updates these pointers based on the
     * current state of the Lagrangian nodal position data.
     *
     * NOTE: It is important to note that any operation on the
     * LNodeLevelData which results in the restoration of the local
     * form of the underlying PETSc Vec object potentially invalidates
     * these pointers.
     */
    void restoreLocationPointers(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * @brief Each LNodeIndex object owns a pointer to its nodal
     * location.  This routine invalidates these pointers (an action
     * which is mainly useful for debugging purposes).
     */
    void invalidateLocationPointers(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    ///
    ///  The following routines:
    ///
    ///      initializeLevelData(),
    ///      resetHierarchyConfiguration(),
    ///      applyGradientDetector()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
    ///

    /*!
     * Initialize data on a new level after it is inserted into an AMR
     * patch hierarchy by the gridding algorithm.  The level number
     * indicates that of the new level.  The old_level pointer
     * corresponds to the level that resided in the hierarchy before
     * the level with the specified number was introduced.  If the
     * pointer is null, there was no level in the hierarchy prior to
     * the call and the level data is set based on the user routines
     * and the simulation time.  Otherwise, the specified level
     * replaces the old level and the new level receives data from the
     * old level appropriately before it is destroyed.
     *
     * The boolean argument initial_time indicates whether the level
     * is being introduced for the first time (i.e., at initialization
     * time) or after some regrid process during the calculation
     * beyond the initial hierarchy construction.  This information is
     * provided since the initialization of the data on a patch may be
     * different in each of those circumstances.  The can_be_refined
     * boolean argument indicates whether the level is the finest
     * level allowed in the hierarchy.  This may or may not affect the
     * data initialization process depending on the problem.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null, the level number
     * does not match any level in the hierarchy, or the old level
     * number does not match the level number (if the old level
     * pointer is non-null).
     */
    void initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level=SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

    /*!
     * Reset cached communication schedules after the hierarchy has
     * changed (for example, due to regidding) and the data has been
     * initialized on the new levels.  The intent is that the cost of
     * data movement on the hierarchy will be amortized across
     * multiple communication cycles, if possible.  The level numbers
     * indicate the range of levels in the hierarchy that have
     * changed.  However, this routine updates communication schedules
     * every level finer than and including that indexed by the
     * coarsest level number given.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null, any pointer to a
     * level in the hierarchy that is coarser than the finest level is
     * null, or the given level numbers not specified properly; e.g.,
     * coarsest_level > finest_level.
     */
    void resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the
     * given level should occur due to the presence of Lagrangian
     * data.  The double time argument is the regrid time.  The
     * integer "tag_index" argument is the patch descriptor index of
     * the cell centered integer tag array on each patch in the
     * hierarchy.  The boolean argument initial_time indicates whether
     * the level is being subject to refinement at the initial
     * simulation time.  If it is false, then the error estimation
     * process is being invoked at some later time after the AMR
     * hierarchy was initially constructed.  The boolean argument
     * uses_richardson_extrapolation_too is true when Richardson
     * extrapolation error estimation is used in addition to the
     * gradient detector, and false otherwise.  This argument helps
     * the user to manage multiple regridding criteria.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null or the level
     * number does not match any existing level in the hierarchy.
     */
    virtual void applyGradientDetector(
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
     * When assertion checking is active, database pointer must be
     * non-null.
     */
    void putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * @brief Constructor.
     *
     * XXXX: Need mechanism for specifiying GCW!
     */
    LDataManager(
        const string& object_name,
        const SAMRAI::hier::IntVector<NDIM>& ghosts=3,
        bool register_for_restart=true);

    /*!
     * @brief The LDataManager destructor cleans up any remaining
     * PETSc AO objects.
     */
    ~LDataManager();

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    LDataManager();

    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * @param from The value to copy to this object.
     */
    LDataManager(
        const LDataManager& from);

    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * @param that The value to assign to this object.
     *
     * @return A reference to this object.
     */
    LDataManager& operator=(
        const LDataManager& that);

    /*!
     * @brief Begin the process of refilling nonlocal Lagrangian
     * quantites over the specified range of levels in the patch
     * hierarchy.
     *
     * The operation is essentially equivalent to refilling ghost
     * cells for structured (SAMRAI native) data.
     */
    void beginNonlocalDataFill(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * @brief End the process of refilling nonlocal Lagrangian
     * quantites over the specified range of levels in the patch
     * hierarchy.
     *
     * The operation is essentially equivalent to refilling ghost
     * cells for structured (SAMRAI native) data.
     */
    void endNonlocalDataFill(
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * Determines the global Lagrangian and PETSc indices of the local
     * and nonlocal nodes associated with the processor as well as the
     * local PETSc indices of the interior and ghost nodes in each
     * patch of the specified level.
     *
     * NOTE: The set of local Lagrangian indices lists all the nodes
     * which are owned by this processor.  The set of nonlocal
     * Lagrangian indices lists all of the nodes which are not owned
     * by this processor but which appear in the ghost cell region of
     * some patch on this processor.  Both of these sets of node
     * indices use the fixed, global Lagrangian indexing scheme.
     *
     * NOTE: The set of interior local indices lists the nodes which
     * live on the interior on each patch.  The set of ghost local
     * indices lists the nodes which live in the ghost cell region of
     * each patch.  Both of these sets of node indices use the local
     * PETSc indexing scheme, determined by the present distribution
     * of data across the processors.
     *
     * Since each processor may own multiple patches in a given level,
     * nodes appearing in the ghost cell region of a patch may or may
     * not be owned by this processor.
     */
    int computeNodeDistribution(
        std::vector<int>& local_lag_indices,
        std::vector<int>& nonlocal_lag_indices,
        AO& ao,
        std::vector<int>& local_petsc_indices,
        std::vector<int>& nonlocal_petsc_indices,
        int& num_nodes,
        int& node_offset,
        std::map<int,std::vector<int>*>& patch_interior_local_indices,
        std::map<int,std::vector<int>*>& patch_ghost_local_indices,
        const int ln);

    /*!
     * Determine the number of local Lagrangian nodes on all MPI
     * processes with rank less than the rank of the current MPI
     * process.
     */
    static void computeNodeOffsets(
        int& num_nodes,
        int& node_offset,
        const int& num_local_nodes);

    /*!
     * Read object state from the restart file and initialize class
     * data members.  The database from which the restart data is read
     * is determined by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found
     *        in the restart file.
     *
     *    -   The class version number and restart version number do not
     *        match.
     *
     */
    void getFromRestart();

    /*!
     * Static data members used to control access to and destruction
     * of singleton data manager instance.
     */
    static std::map<string,LDataManager*> s_data_manager_instances;
    static bool s_registered_callback;
    static unsigned char s_shutdown_priority;

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.  The boolean is
     * used to control restart file writing operations.
     */
    string d_object_name;
    bool d_registered_for_restart;

    /*
     * Grid hierarchy information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln, d_finest_ln;

    /*
     * We cache a pointer to the VisIt and Silo data writers to
     * register plot variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;
    SAMRAI::tbox::Pointer<LagSiloDataWriter> d_silo_writer;

    /*
     * We cache a pointer to the load balancer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /*
     * Objects used to specify the Jacobian determinant of the
     * curvilinear coordinate system on the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<LNodeJacobianInitStrategy> d_lag_jac_init;

    /*
     * Objects used to specify and initialize the Lagrangian data on
     * the patch hierarchy.
     */
    SAMRAI::tbox::Pointer<LNodePosnInitStrategy> d_lag_posn_init;
    std::vector<bool> d_level_contains_lag_data;

    /*
     * SAMRAI::hier::Variable<NDIM> pointer and patch data descriptor indices
     * for the LNodeIndexData used to define the data distribution.
     */
    SAMRAI::tbox::Pointer<LNodeIndexVariable> d_lag_node_index_var;
    int d_lag_node_index_current_idx, d_lag_node_index_scratch_idx;

    /*
     * SAMRAI::hier::Variable<NDIM> pointer and patch data descriptor indices
     * for the cell variable used to determine the workload for
     * nonuniform load balancing.
     */
    double d_alpha_work, d_beta_work;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_workload_var;
    int d_workload_idx;

    /*
     * SAMRAI::hier::Variable<NDIM> pointer and patch data descriptor indices
     * for the cell variable used to keep track of the count of the
     * nodes in each cell for visualization and tagging purposes.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_node_count_var;
    int d_node_count_idx;
    bool d_output_node_count;

    /*
     * SAMRAI::hier::Variable<NDIM> pointer and patch data descriptor indices
     * for the cell variable used to keep track of the MPI process
     * assigned to each patch.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,int> > d_mpi_proc_var;
    int d_mpi_proc_idx;
    bool d_output_mpi_proc;

    /*
     * SAMRAI::hier::IntVector<NDIM> object which determines the ghost cell
     * width of the LNodeIndexData SAMRAI::hier::PatchData<NDIM> objects.
     */
    SAMRAI::hier::IntVector<NDIM> d_ghosts;

    /*
     * Communications algorithms and schedules.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_lag_node_index_bdry_fill_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_lag_node_index_bdry_fill_scheds;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_node_count_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_node_count_coarsen_scheds;

    /*
     * SAMRAI::hier::VariableContext objects are used for data management.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_current_context, d_scratch_context;

    /*
     * ComponenetSelector object allow for the collective allocation
     * and deallocation of SAMRAI::hier::PatchData<NDIM>.
     */
    SAMRAI::hier::ComponentSelector d_current_data, d_scratch_data;

    //@{ @Data that is separately maintained for each level of the
    //patch hierarchy.

    /*!
     * The Lagrangian quantity data owned by the manager object.
     */
    std::vector<std::map<string,SAMRAI::tbox::Pointer<LNodeLevelData> > > d_lag_quantity_data;

    /*!
     * Indicates whether the LNodeLevelData is in synch with the
     * LNodeIndexData.
     */
    std::vector<bool> d_needs_synch;

    /*!
     * PETSc AO objects provide mappings between the fixed global
     * Lagrangian node IDs and the ever-changing global PETSc
     * ordering.
     */
    std::vector<AO> d_ao;
    static std::vector<int> s_ao_dummy;

    /*!
     * The total number of nodes for all processors.
     */
    std::vector<int> d_num_nodes;

    /*!
     * The total number of local nodes for all processors with rank
     * less than the rank of the current processor.
     */
    std::vector<int> d_node_offset;

    /*!
     * The Lagrangian node indices of all local and nonlocal nodes on
     * each level of the patch hierarchy.
     *
     * A local node is one that is owned by a patch on this processor,
     * while a nonlocal node is one that is owned by a patch on
     * another processor, but found in the ghost region of some patch
     * owned by this processor.
     *
     * Note that these sets of indices provide the information
     * necessary to determine the local PETSc index for all nodes.
     * Local node d_local_lag_indices[ln][j] has local PETSc index j,
     * while nonlocal node d_nonlocal_lag_indices[ln][k] has local
     * PETSc index d_local_lag_indices.size()+j.
     *
     * It is possible to determine the global PETSc index of a local
     * node by making use of d_node_offset.  Local node
     * d_local_lag_indices[ln][j] has global PETSc index
     * j+d_node_offset[ln].  A similar mapping for nonlocal nodes is
     * not well defined.
     */
    std::vector<std::vector<int> > d_local_lag_indices;
    std::vector<std::vector<int> > d_nonlocal_lag_indices;

    /*!
     * The node indices of all local nodes (i.e. the nodes owned by
     * this processor) on each level of the hierarchy.  The indices
     * are in the global PETSc ordering corresponding to a depth of 1.
     */
    std::vector<std::vector<int> > d_local_petsc_indices;

    /*!
     * The node indices of all nonlocal nodes (i.e. the nodes owned by
     * another processor which appear in the ghost region of some
     * patch owned by this processor) on each level of the hierarchy.
     * The indices are in the global PETSc ordering corresponding to a
     * depth of 1.
     *
     * NOTE: These sets are used to create the VecScatter objects used
     * to transfer data from the old PETSc ordering to the new PETSc
     * ordering.  Since the ordering is different for different depths
     * of LNodeLevelData, we compute one set of indices for each depth
     * which is being reordered.
     *
     * TODO: Is this data really necessary for depth != 1?  It doesn't
     * seem to be used in any essential way by the class.
     */
    std::vector<std::map<int,std::vector<int> > > d_nonlocal_petsc_indices;

    //@}
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LDataManager.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LDataManager
