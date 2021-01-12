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

#ifndef included_IBTK_LDataManager
#define included_IBTK_LDataManager

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LInitStrategy.h"
#include "ibtk/LNodeSet.h"
#include "ibtk/LNodeSetVariable.h"
#include "ibtk/LSiloDataWriter.h"
#include "ibtk/ParallelSet.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_utilities.h"

#include "BasePatchLevel.h"
#include "CartesianGridGeometry.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenSchedule.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "LoadBalancer.h"
#include "PatchHierarchy.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "StandardTagAndInitStrategy.h"
#include "VariableContext.h"
#include "VisItDataWriter.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include "petscao.h"
#include "petscvec.h"

#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
class LData;
class LMesh;
class LNode;
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

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LDataManager coordinates the irregular distribution of LNode and
 * LData on the patch hierarchy.
 *
 * The manager class is responsible for maintaining this data distribution and
 * for all inter-processor communications.  All access to instantiated
 * LDataManager objects is via the static method getManager().
 *
 * Each Lagrangian point is associated with an LNode object which provides an interface between the Lagrangian and PETSc
 * ordering as well as the data storage for force specification objects and other objects associated with the point.
 * From each LNode object, we can retrieve the Lagrangian index as well as any associated node data. The LNode objects
 * for each processor are contained in an LMesh object. The LMesh object contains two vectors of LNodes: one containing
 * just the local LNode objects, and another containing 'ghost' LNodes. The 'ghost' LNodes are Lagrangian points
 * assigned to other processors but have data needed for calculations. Interacting with Lagrangian data is mediated
 * through the static LDataManager object. From this class, we can get the local LMesh object for each level of the
 * Cartesian grid, and then loop through all associated Lagrangian nodes.
 *
 * As an example, suppose we have a circle of Lagrangian points that are tethered to target points, and we wish to
 * specify the motion of the target points, and hence the circle, to move in a straight line. We can do this by looping
 * through all the nodes, pull out the IBAMR::IBTargetPointForceSpec, and then update the target point location. In
 * parallel, we also need to update the target point locations of the ghost data.
 *
 * \code{.cpp}
 * void
 * update_target_points(Pointer<PatchHierarchy<NDIM>> hierarchy,
 *                      LDataManager* const l_data_manager,
 *                      const double current_time,
 *                      const double dt)
 * {
 *     const int finest_ln = hierarchy->getFinestLevelNumber();
 *     // Note we assume the circle is the 0th structure and on the finest level.
 *     const std::pair<int, int>& circle_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);
 *     Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
 *     const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
 *     const std::vector<LNode*>& ghost_nodes = mesh->getGhostNodes();
 *     std::vector<LNode*> nodes = local_nodes;
 *     nodes.insert(nodes.end(), ghost_nodes.begin(), ghost_nodes.end());
 *     for (auto& node : nodes)
 *     {
 *         auto force_spec = node->getNodeDataItem<IBTargetPointForceSpec>();
 *         if (force_spec == nullptr) continue;
 *         const int lag_idx = node->getLagrangianIndex();
 *         Point& X_target = force_spec->getTargetPointPosition();
 *         if (circle_lag_idxs.first <= lag_idx && lag_idx < circle_lag_idxs.second)
 *         {
 *             X_target[0] += 0.1 * dt; // Move point to the right with speed 0.1
 *         }
 *     }
 *     return;
 * }
 * \endcode
 *
 * For some applications, it might be useful to get the Eulerian data of the Lagrangian points. To do this, it is
 * important to understand the Lagrangian and PETSc index ordering. The Lagrangian ordering is fixed at the beginning of
 * the simulation and is set by the reading of the vertices, whether through an input file or programmatically. The
 * PETSc ordering, however, will change over the course of the simulation. Whenever a regridding operation is triggered,
 * the PETSc indexing is changed to ensure efficient memory usage. The Eulerian data is stored in PETSc ordering, so in
 * order to access that data, we need to map Lagrangian indices to their corresponding PETSc indices. We can then access
 * the corresponding components of the PETSc data vector which is wrapped in an LData object. As an example, here we
 * loop over and print out the physical coordinates of all the Lagrangian points on the finest level.
 *
 * \code{.cpp}
 * void
 * print_eul_data(Pointer<PatchHierarchy<NDIM>> hierarchy,
 *                LDataManager* const l_data_manager)
 * {
 *     const int finest_ln = hierarchy->getFinestLevelNumber();
 *     Pointer<LData> X_data = l_data_manager->getLData("X", finest_ln);
 *     Vec X_vec = X_data->getVec();
 *     double* x_vals;
 *     int ierr = VecGetArray(X_vec, &x_vals);
 *     IBTK_CHKERRQ(ierr);
 *
 *     Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
 *     const std::vector<LNode*>& nodes = mesh->getLocalNodes();
 *     for (const auto& node : nodes)
 *     {
 *         const int lag_idx = node->getLagrangianIndex();
 *         // Note we use the local index instead of the global index.
 *         // The global index is the local index plust the total indices on all processors of lower rank.
 *         const int petsc_idx = node->getLocalPETScIndex();
 *         Eigen::Map<VectorNd> X(&x_vals[petsc_idx * NDIM]);
 *         pout << "Euerian location of node " << lag_idx << ":\n"
 *              << X << "\n";
 *     }
 *     return;
 * }
 * \endcode
 *
 * \note Multiple LDataManager objects may be instantiated simultaneously.
 *
 * \see LMesh, LNode, LData
 */
class LDataManager : public SAMRAI::tbox::Serializable, public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>
{
public:
    /*!
     * The name of the LData that specifies the current positions of
     * the curvilinear mesh nodes.
     */
    static const std::string POSN_DATA_NAME;

    /*!
     * The name of the LData that specifies the initial positions of
     * the curvilinear mesh nodes.
     */
    static const std::string INIT_POSN_DATA_NAME;

    /*!
     * The name of the LData that specifies the velocities of the
     * curvilinear mesh nodes.
     */
    static const std::string VEL_DATA_NAME;

    /*!
     * Return a pointer to the instance of the Lagrangian data manager
     * corresponding to the specified name.  Access to LDataManager objects is
     * mediated by the getManager() function.
     *
     * Note that when a manager is accessed for the first time, the
     * freeAllManagers static method is registered with the ShutdownRegistry
     * class.  Consequently, all allocated managers are freed at program
     * completion.  Thus, users of this class do not explicitly allocate or
     * deallocate the LDataManager instances.
     *
     * \return A pointer to the data manager instance.
     *
     * \note By default, the ghost cell width is set according to the
     * interpolation and spreading kernel functions.
     */
    static LDataManager*
    getManager(const std::string& name,
               const std::string& default_interp_kernel_fcn,
               const std::string& default_spread_kernel_fcn,
               bool error_if_points_leave_domain = false,
               const SAMRAI::hier::IntVector<NDIM>& min_ghost_width = SAMRAI::hier::IntVector<NDIM>(0),
               bool register_for_restart = true);

    /*!
     * Deallocate all of the LDataManager instances.
     *
     * It is not necessary to call this function at program termination since it
     * is automatically called by the ShutdownRegistry class.
     */
    static void freeAllManagers();

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
     * \brief Return the ghost cell width associated with the interaction
     * scheme.
     */
    const SAMRAI::hier::IntVector<NDIM>& getGhostCellWidth() const;

    /*!
     * \brief Return the default kernel function associated with the
     * Eulerian-to-Lagrangian interpolation scheme.
     */
    const std::string& getDefaultInterpKernelFunction() const;

    /*!
     * \brief Return the default kernel function associated with the
     * Lagrangian-to-Eulerian spreading scheme.
     */
    const std::string& getDefaultSpreadKernelFunction() const;

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using the default spreading kernel function.
     *
     * \note This spreading operation does include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s)) ds(q,r,s)
     *
     * This is the standard regularized delta function spreading operation,
     * which spreads densities, \em NOT values.
     */
    void spread(int f_data_idx,
                SAMRAI::tbox::Pointer<LData> F_data,
                SAMRAI::tbox::Pointer<LData> X_data,
                SAMRAI::tbox::Pointer<LData> ds_data,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                int level_num,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true,
                bool ds_data_ghost_node_update = true);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using a specified spreading kernel function.
     *
     * \note This spreading operation does include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s)) ds(q,r,s)
     *
     * This is the standard regularized delta function spreading operation,
     * which spreads densities, \em NOT values.
     */
    void spread(int f_data_idx,
                SAMRAI::tbox::Pointer<LData> F_data,
                SAMRAI::tbox::Pointer<LData> X_data,
                SAMRAI::tbox::Pointer<LData> ds_data,
                const std::string& spread_kernel_fcn,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                int level_num,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true,
                bool ds_data_ghost_node_update = true);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using the default spreading kernel function.
     *
     * \note This spreading operation does include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s)) ds(q,r,s)
     *
     * This is the standard regularized delta function spreading operation,
     * which spreads densities, \em NOT values.
     */
    void spread(int f_data_idx,
                std::vector<SAMRAI::tbox::Pointer<LData> >& F_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& X_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& ds_data,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true,
                bool ds_data_ghost_node_update = true,
                int coarsest_ln = -1,
                int finest_ln = -1);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using a specified spreading kernel function.
     *
     * \note This spreading operation does include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s)) ds(q,r,s)
     *
     * This is the standard regularized delta function spreading operation,
     * which spreads densities, \em NOT values.
     */
    void spread(int f_data_idx,
                std::vector<SAMRAI::tbox::Pointer<LData> >& F_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& X_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& ds_data,
                const std::string& spread_kernel_fcn,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true,
                bool ds_data_ghost_node_update = true,
                int coarsest_ln = -1,
                int finest_ln = -1);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using the default spreading kernel function.
     *
     * \note This spreading operation does NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operation spreads values, \em NOT densities.
     */
    void spread(int f_data_idx,
                SAMRAI::tbox::Pointer<LData> F_data,
                SAMRAI::tbox::Pointer<LData> X_data,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                int level_num,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using the specified spreading kernel function.
     *
     * \note This spreading operation does NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operation spreads values, \em NOT densities.
     */
    void spread(int f_data_idx,
                SAMRAI::tbox::Pointer<LData> F_data,
                SAMRAI::tbox::Pointer<LData> X_data,
                const std::string& spread_kernel_fcn,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                int level_num,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using the default spreading kernel function.
     *
     * \note This spreading operation does NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operation spreads values, \em NOT densities.
     */
    void spread(int f_data_idx,
                std::vector<SAMRAI::tbox::Pointer<LData> >& F_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& X_data,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true,
                int coarsest_ln = -1,
                int finest_ln = -1);

    /*!
     * \brief Spread a quantity from the Lagrangian mesh to the Eulerian grid
     * using the specified spreading kernel function.
     *
     * \note This spreading operation does NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     f(i,j,k) = f(i,j,k) + Sum_{q,r,s} F(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operation spreads values, \em NOT densities.
     */
    void spread(int f_data_idx,
                std::vector<SAMRAI::tbox::Pointer<LData> >& F_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& X_data,
                const std::string& spread_kernel_fcn,
                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                bool F_data_ghost_node_update = true,
                bool X_data_ghost_node_update = true,
                int coarsest_ln = -1,
                int finest_ln = -1);

    /*!
     * \brief Interpolate a quantity from the Eulerian grid to the Lagrangian
     * mesh using the default interpolation kernel function.
     */
    void interp(int f_data_idx,
                SAMRAI::tbox::Pointer<LData> F_data,
                SAMRAI::tbox::Pointer<LData> X_data,
                int level_num,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& f_synch_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >(),
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_ghost_fill_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0);

    /*!
     * \brief Interpolate a quantity from the Eulerian grid to the Lagrangian
     * mesh using the specified interpolation kernel function.
     */
    void interp(int f_data_idx,
                SAMRAI::tbox::Pointer<LData> F_data,
                SAMRAI::tbox::Pointer<LData> X_data,
                const std::string& interp_kernel_fcn,
                int level_num,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& f_synch_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >(),
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_ghost_fill_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0);

    /*!
     * \brief Interpolate a quantity from the Eulerian grid to the Lagrangian
     * mesh using the default interpolation kernel function.
     */
    void interp(int f_data_idx,
                std::vector<SAMRAI::tbox::Pointer<LData> >& F_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& X_data,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& f_synch_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >(),
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_ghost_fill_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                int coarsest_ln = -1,
                int finest_ln = -1);

    /*!
     * \brief Interpolate a quantity from the Eulerian grid to the Lagrangian
     * mesh using the specified interpolation kernel function.
     */
    void interp(int f_data_idx,
                std::vector<SAMRAI::tbox::Pointer<LData> >& F_data,
                std::vector<SAMRAI::tbox::Pointer<LData> >& X_data,
                const std::string& interp_kernel_fcn,
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& f_synch_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >(),
                const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_ghost_fill_scheds =
                    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(),
                double fill_data_time = 0.0,
                int coarsest_ln = -1,
                int finest_ln = -1);

    /*!
     * Register a concrete strategy object with the integrator that specifies
     * the initial configuration of the curvilinear mesh nodes.
     *
     * \note This function calls LInitStrategy::init(). All preprocessing should be completed before
     * registering a LInitStrategy object.
     */
    void registerLInitStrategy(SAMRAI::tbox::Pointer<LInitStrategy> lag_init);

    /*!
     * Free the concrete initialization strategy object.
     *
     * \note Be sure to call this method only once the initialization object is
     * no longer needed.
     */
    void freeLInitStrategy();

    /*!
     * \brief Register a VisIt data writer with the manager.
     */
    void registerVisItDataWriter(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * \brief Register a Silo data writer with the manager.
     */
    void registerLSiloDataWriter(SAMRAI::tbox::Pointer<LSiloDataWriter> silo_writer);

    /*!
     * \brief Register a load balancer for non-uniform load balancing.
     *
     * @deprecated This method is deprecated since the current strategy for
     * handling non-uniform load balancing does not require that this object
     * store a pointer to the load balancer.
     */
    void registerLoadBalancer(SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
                              int workload_data_idx);

    /*!
     * \brief Indicates whether there is Lagrangian data on the given patch
     * hierarchy level.
     */
    bool levelContainsLagrangianData(int level_number) const;

    /*!
     * \return The number of total nodes of the Lagrangian data for the
     * specified level of the patch hierarchy.
     */
    unsigned int getNumberOfNodes(int level_number) const;

    /*!
     * \return The number of local (i.e., on-processor) nodes of the Lagrangian
     * data for the specified level of the patch hierarchy.
     *
     * \note This count does not include nodes that only lie in ghost cells for
     * the current process.
     *
     * \see getNumberOfNodes
     * \see getNumberOfGhostNodes
     */
    unsigned int getNumberOfLocalNodes(int level_number) const;

    /*!
     * \return The number of ghost (i.e., off-processor) nodes of the Lagrangian
     * data for the specified level of the patch hierarchy.
     *
     * \see getNumberOfNodes
     * \see getNumberOfLocalNodes
     */
    unsigned int getNumberOfGhostNodes(int level_number) const;

    /*!
     * \return The number of nodes on all processors with MPI rank less than the
     * current process on the specified level of the patch hierarchy.
     *
     * \note This count does not include nodes that only lie in ghost cells for
     * the current process.
     */
    unsigned int getGlobalNodeOffset(int level_number) const;

    /*!
     * \brief Get the Lagrangian mesh associated with the given patch hierarchy
     * level.
     */
    SAMRAI::tbox::Pointer<LMesh> getLMesh(int level_number) const;

    /*!
     * \brief Get the specified Lagrangian quantity data on the given patch
     * hierarchy level.
     */
    SAMRAI::tbox::Pointer<LData> getLData(const std::string& quantity_name, int level_number) const;

    /*!
     * \brief Allocate new Lagrangian level data with the specified name and
     * depth.  If specified, the quantity is maintained as the patch hierarchy
     * evolves.
     *
     * \note Quantities maintained by the LDataManager must have unique names.
     * The name "X" is reserved for the nodal coordinates.
     */
    SAMRAI::tbox::Pointer<LData>
    createLData(const std::string& quantity_name, int level_number, unsigned int depth = 1, bool maintain_data = false);

    /*!
     * \brief Get the patch data descriptor index for the Lagrangian index data.
     */
    int getLNodePatchDescriptorIndex() const;

    /*!
     * \brief Get the patch data descriptor index for the workload cell data.
     *
     * @deprecated This method is deprecated since, in future versions of
     * IBAMR, this value will no longer be stored and will only be available
     * via the parent hierarchy integrator.
     */
    int getWorkloadPatchDescriptorIndex() const;

    /*!
     * \brief Get the patch data descriptor index for the Lagrangian node count
     * cell data.
     */
    int getNodeCountPatchDescriptorIndex() const;

    /*!
     * \brief Get a list of Lagrangian structure names for the specified level
     * of the patch hierarchy.
     */
    std::vector<std::string> getLagrangianStructureNames(int level_number) const;

    /*!
     * \brief Get a list of Lagrangian structure IDs for the specified level of
     * the patch hierarchy.
     */
    std::vector<int> getLagrangianStructureIDs(int level_number) const;

    /*!
     * \brief Get the ID of the Lagrangian structure associated with the
     * specified Lagrangian index.
     *
     * \note Returns -1 in the case that the Lagrangian index is not associated
     * with any Lagrangian structure.
     */
    int getLagrangianStructureID(int lagrangian_index, int level_number) const;

    /*!
     * \brief Get the ID of the Lagrangian structure with the specified name.
     *
     * \note Returns -1 in the case that the Lagrangian structure name is not
     * associated with any Lagrangian structure.
     */
    int getLagrangianStructureID(const std::string& structure_name, int level_number) const;

    /*!
     * \brief Get the name of the Lagrangian structure with the specified ID.
     *
     * \note Returns "UNKNOWN" in the case that the Lagrangian structure ID is
     * not associated with any Lagrangian structure.
     */
    std::string getLagrangianStructureName(int structure_id, int level_number) const;

    /*!
     * \brief Get the range of Lagrangian indices for the Lagrangian structure
     * with the specified ID.
     *
     * \return A pair of indices such that if pair.first <= lag_idx <
     * pair.second, then lag_idx is associated with the specified structure;
     * otherwise, lag_idx is not associated with the specified structure.
     *
     * \note Returns std::make_pair(-1,-1) in the case that the Lagrangian
     * structure ID is not associated with any Lagrangian structure.
     */
    std::pair<int, int> getLagrangianStructureIndexRange(int structure_id, int level_number) const;

    /*!
     * \brief Get the center of mass of the Lagrangian structure with the
     * specified ID.
     *
     * \note The center of mass X of a particular structure is computed as
     *
     *    X = (1/N) Sum_{k in structure} X_k
     *
     * in which N is the number of nodes associated with that structure.
     *
     * \note Returns Point::Zero() in the case that the Lagrangian structure
     * ID is not associated with any Lagrangian structure.
     */
    Point computeLagrangianStructureCenterOfMass(int structure_id, int level_number);

    /*!
     * \brief Get the bounding box of the Lagrangian structure with the
     * specified ID.
     *
     * \note Returns the entire range of double precision values in the case
     * that the Lagrangian structure ID is not associated with any Lagrangian
     * structure.
     */
    std::pair<Point, Point> computeLagrangianStructureBoundingBox(int structure_id, int level_number);

    /*!
     * \brief Reset the positions of the nodes of the Lagrangian structure with
     * the specified ID to be equal to the initial positions but shifted so that
     * the bounding box of the structure is centered about X_center.
     *
     * \note This operation must be performed immediately before a regridding
     * operation, otherwise the results are undefined.
     */
    void reinitLagrangianStructure(const Point& X_center, int structure_id, int level_number);

    /*!
     * \brief Shift the positions of the nodes of the Lagrangian structure with
     * the specified ID by a displacement dX.
     *
     * \note This operation must be performed immediately before a regridding
     * operation, otherwise the results are undefined.
     *
     * \warning All displacements must involve shifts that do \em not cross
     * periodic boundaries.
     */
    void displaceLagrangianStructure(const Vector& dX, int structure_id, int level_number);

    /*!
     * \brief Activate the Lagrangian structures with the specified ID numbers.
     *
     * \note This method is collective (i.e., must be called by all MPI
     * processes); however, each MPI process may provide a different collection
     * of structures to activate.
     */
    void activateLagrangianStructures(const std::vector<int>& structure_ids, int level_number);

    /*!
     * \brief Inactivate the Lagrangian structures with the specified ID
     * numbers.
     *
     * \note This method is collective (i.e., must be called by all MPI
     * processes); however, each MPI process may provide a different collection
     * of structures to inactivate.
     */
    void inactivateLagrangianStructures(const std::vector<int>& structure_ids, int level_number);

    /*!
     * \brief Determine whether the Lagrangian structure with the specified ID
     * number is activated.
     */
    bool getLagrangianStructureIsActivated(int structure_id, int level_number) const;

    /*!
     * \brief Set the components of the supplied LData object to zero
     * for those entries that correspond to inactivated structures.
     */
    void zeroInactivatedComponents(SAMRAI::tbox::Pointer<LData> lag_data, int level_number) const;

    /*!
     * \brief Map the collection of Lagrangian indices to the corresponding
     * global PETSc indices.
     */
    void mapLagrangianToPETSc(std::vector<int>& inds, int level_number) const;

    /*!
     * \brief Map the collection of global PETSc indices to the corresponding
     * Lagrangian indices.
     */
    void mapPETScToLagrangian(std::vector<int>& inds, int level_number) const;

    /*!
     * \brief Scatter data from the Lagrangian ordering to the global PETSc
     * ordering.
     *
     * \todo Optimize the implementation of this method.
     */
    void scatterLagrangianToPETSc(Vec& lagrangian_vec, Vec& petsc_vec, int level_number) const;

    /*!
     * \brief Scatter data from the global PETSc ordering to the Lagrangian
     * ordering.
     *
     * \todo Optimize the implementation of this method.
     */
    void scatterPETScToLagrangian(Vec& petsc_vec, Vec& lagrangian_vec, int level_number) const;

    /*!
     * \brief Scatter data from a distributed PETSc vector to all processors.
     *
     * \todo Optimize the implementation of this method.
     */
    void scatterToAll(Vec& parallel_vec, Vec& sequential_vec) const;

    /*!
     * \brief Scatter data from a distributed PETSc vector to processor zero.
     *
     * \todo Optimize the implementation of this method.
     */
    void scatterToZero(Vec& parallel_vec, Vec& sequential_vec) const;

    /*!
     * \brief Start the process of redistributing the Lagrangian data.
     *
     * This method uses the present location of each Lagrangian mesh node to
     * redistribute the LNodeData managed by this object.
     *
     * \note This routine assumes that the time interval between node
     * redistribution satisfies a timestep restriction of the form dt <=
     * C*dx*|U| with C <= 1.  This restriction prevents nodes from moving more
     * than one cell width per timestep.
     *
     * \see endDataRedistribution
     */
    void beginDataRedistribution(int coarsest_ln = -1, int finest_ln = -1);

    /*!
     * \brief Finish the process of redistributing the Lagrangian data.
     *
     * This method redistributes the quantities associated with each node in the
     * Lagrangian mesh according to the data distribution defined by the
     * LNodeData managed by this object.  This routine potentially
     * involves SUBSTANTIAL inter-processor communication.
     *
     * \note Since this routine potentially results in a large amount of
     * inter-processor communication, it may be worth putting it off for as long
     * as possible.  If the timestep dt satisfies a condition of the form dt <=
     * C*dx*|U| with C << 1, it may be possible to redistribute the Lagrangian
     * data less frequently than every timestep.
     *
     * \see beginDataRedistribution
     */
    void endDataRedistribution(int coarsest_ln = -1, int finest_ln = -1);

    /*!
     * \brief Update the workload and count of nodes per cell.
     *
     * This routine updates cell data that is maintained on the patch hierarchy
     * to track the number of nodes in each cell of the AMR index space.  The
     * node count data is used to tag cells for refinement, and to specify
     * non-uniform load balancing.  The workload per cell is defined by
     *
     *    workload(i) = 1 + beta_work*node_count(i)
     *
     * in which alpha and beta are parameters that each default to the value 1.
     */
    void addWorkloadEstimate(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const int workload_data_idx,
                             const int coarsest_ln = -1,
                             const int finest_ln = -1);

    /*!
     * \brief Update the count of nodes per cell.
     *
     * This routine updates cell data that is maintained on the patch hierarchy
     * to track the number of nodes in each cell of the AMR index space.  The
     * node count data is used to tag cells for refinement, and to specify
     * non-uniform load balancing.
     */
    void updateNodeCountData(int coarsest_ln = -1, int finest_ln = -1);

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

    /*!
     * Register user defined Lagrangian data to be maintained
     *
     */
    void registerUserDefinedLData(const std::string& data_name, int depth);

protected:
    /*!
     * \brief Constructor.
     */
    LDataManager(std::string object_name,
                 std::string default_interp_kernel_fcn,
                 std::string default_spread_kernel_fcn,
                 bool error_if_points_leave_domain,
                 SAMRAI::hier::IntVector<NDIM> ghost_width,
                 bool register_for_restart = true);

    /*!
     * \brief The LDataManager destructor cleans up any remaining PETSc AO
     * objects.
     */
    ~LDataManager();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LDataManager();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LDataManager(const LDataManager& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LDataManager& operator=(const LDataManager& that) = delete;

    /*!
     * \brief Common implementation of scatterPETScToLagrangian() and
     * scatterLagrangianToPETSc().
     */
    void scatterData(Vec& lagrangian_vec, Vec& petsc_vec, int level_number, ScatterMode mode) const;

    /*!
     * \brief Begin the process of refilling nonlocal Lagrangian quantities over
     * the specified range of levels in the patch hierarchy.
     *
     * The operation is essentially equivalent to refilling ghost cells for
     * structured (SAMRAI native) data.
     */
    void beginNonlocalDataFill(int coarsest_ln = -1, int finest_ln = -1);

    /*!
     * \brief End the process of refilling nonlocal Lagrangian quantities over
     * the specified range of levels in the patch hierarchy.
     *
     * The operation is essentially equivalent to refilling ghost cells for
     * structured (SAMRAI native) data.
     */
    void endNonlocalDataFill(int coarsest_ln = -1, int finest_ln = -1);

    /*!
     * Determines the global Lagrangian and PETSc indices of the local and
     * nonlocal nodes associated with the processor as well as the local PETSc
     * indices of the interior and ghost nodes in each patch of the specified
     * level.
     *
     * \note The set of local Lagrangian indices lists all the nodes that are
     * owned by this processor.  The set of nonlocal Lagrangian indices lists
     * all of the nodes that are not owned by this processor but that appear in
     * the ghost cell region of some patch on this processor.  Both of these
     * sets of node indices use the fixed, global Lagrangian indexing scheme.
     *
     * \note The set of interior local indices lists the nodes that live on the
     * interior on each patch.  The set of ghost local indices lists the nodes
     * that live in the ghost cell region of each patch.  Both of these sets of
     * node indices use the local PETSc indexing scheme, determined by the
     * present distribution of data across the processors.
     *
     * Since each processor may own multiple patches in a given level, nodes
     * appearing in the ghost cell region of a patch may or may not be owned by
     * this processor.
     */
    void computeNodeDistribution(AO& ao,
                                 std::vector<int>& local_lag_indices,
                                 std::vector<int>& nonlocal_lag_indices,
                                 std::vector<int>& local_petsc_indices,
                                 std::vector<int>& nonlocal_petsc_indices,
                                 unsigned int& num_nodes,
                                 unsigned int& node_offset,
                                 int level_number);

    /*!
     * Determine the number of local Lagrangian nodes on all MPI processes with
     * rank less than the rank of the current MPI process.
     */
    static void computeNodeOffsets(unsigned int& num_nodes, unsigned int& node_offset, unsigned int num_local_nodes);

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
    static std::map<std::string, LDataManager*> s_data_manager_instances;
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
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    /*
     * Cached Eulerian data to reduce the number of allocations/deallocations.
     */
    SAMRAIDataCache d_cached_eulerian_data;

    /*
     * We cache a pointer to the visualization data writers to register plot
     * variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;
    SAMRAI::tbox::Pointer<LSiloDataWriter> d_silo_writer;

    /*
     * We cache a pointer to the load balancer.
     *
     * @deprecated This will be removed in a future release since it is no longer necessary.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /*
     * Objects used to specify and initialize the Lagrangian data on the patch
     * hierarchy.
     */
    SAMRAI::tbox::Pointer<LInitStrategy> d_lag_init;
    std::vector<bool> d_level_contains_lag_data;

    /*
     * SAMRAI::hier::Variable pointer and patch data descriptor indices for the
     * LNodeData used to define the data distribution.
     */
    SAMRAI::tbox::Pointer<LNodeSetVariable> d_lag_node_index_var;
    int d_lag_node_index_current_idx = IBTK::invalid_index, d_lag_node_index_scratch_idx = IBTK::invalid_index;
    std::vector<SAMRAI::tbox::Pointer<std::vector<LNode> > > d_local_and_ghost_nodes;

    /*
     * SAMRAI::hier::Variable pointer and patch data descriptor indices for the
     * cell variable used to determine the workload for nonuniform load
     * balancing.
     *
     * @deprecated d_workload_var and d_workload_idx will not be stored in a
     * future release since the correct workload index will be provided as an
     * argument to addWorkloadEstimate.
     */
    double d_beta_work = 1.0;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_workload_var;
    int d_workload_idx = IBTK::invalid_index;
    bool d_output_workload = false;

    /*
     * SAMRAI::hier::Variable pointer and patch data descriptor indices for the
     * cell variable used to keep track of the count of the nodes in each cell
     * for visualization and tagging purposes.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_node_count_var;
    int d_node_count_idx = IBTK::invalid_index;
    bool d_output_node_count = false;

    /*
     * The kernel functions used to mediate Lagrangian-Eulerian interaction.
     */
    const std::string d_default_interp_kernel_fcn;
    const std::string d_default_spread_kernel_fcn;

    /*
     * Whether to emit an error message if IB points "escape" from the computational
     * domain.
     */
    bool d_error_if_points_leave_domain;

    /*
     * SAMRAI::hier::IntVector object that determines the ghost cell width of
     * the LNodeData SAMRAI::hier::PatchData objects.
     */
    const SAMRAI::hier::IntVector<NDIM> d_ghost_width;

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
     * ComponenetSelector object allow for the collective allocation and
     * deallocation of SAMRAI::hier::PatchData.
     */
    SAMRAI::hier::ComponentSelector d_current_data, d_scratch_data;

    /*!
     * \name Data that is separately maintained for each level of the patch
     * hierarchy.
     */
    //\{

    /*!
     * Information about the names and IDs of the various Lagrangian structures.
     */
    std::vector<std::map<std::string, int> > d_strct_name_to_strct_id_map;
    std::vector<std::map<int, std::string> > d_strct_id_to_strct_name_map;
    std::vector<std::map<int, std::pair<int, int> > > d_strct_id_to_lag_idx_range_map;
    std::vector<std::map<int, int> > d_last_lag_idx_to_strct_id_map;
    std::vector<ParallelSet> d_inactive_strcts;
    std::vector<std::vector<int> > d_displaced_strct_ids;
    std::vector<std::vector<std::pair<Point, Point> > > d_displaced_strct_bounding_boxes;
    std::vector<std::vector<LNodeSet::value_type> > d_displaced_strct_lnode_idxs;
    std::vector<std::vector<Point> > d_displaced_strct_lnode_posns;

    /*!
     * Lagrangian mesh data.
     */
    std::vector<SAMRAI::tbox::Pointer<LMesh> > d_lag_mesh;

    /*!
     * The Lagrangian mesh data owned by the manager object.
     */
    std::vector<std::map<std::string, SAMRAI::tbox::Pointer<LData> > > d_lag_mesh_data;

    /*!
     * Indicates whether the LData is in synch with the LNodeData.
     */
    std::vector<bool> d_needs_synch;

    /*!
     * PETSc AO objects provide mappings between the fixed global Lagrangian
     * node IDs and the ever-changing global PETSc ordering.
     */
    std::vector<AO> d_ao;
    static std::vector<int> s_ao_dummy;

    /*!
     * The total number of nodes for all processors.
     */
    std::vector<unsigned int> d_num_nodes;

    /*!
     * The total number of local nodes for all processors with rank less than
     * the rank of the current processor.
     */
    std::vector<unsigned int> d_node_offset;

    /*!
     * The Lagrangian node indices of all local and nonlocal nodes on each level
     * of the patch hierarchy.
     *
     * A local node is one that is owned by a patch on this processor, while a
     * nonlocal node is one that is owned by a patch on another processor, but
     * found in the ghost region of some patch owned by this processor.
     *
     * Note that these sets of indices provide the information necessary to
     * determine the local PETSc index for all nodes.  Local node
     * d_local_lag_indices[ln][j] has local PETSc index j, while nonlocal node
     * d_nonlocal_lag_indices[ln][k] has local PETSc index
     * d_local_lag_indices.size()+j.
     *
     * It is possible to determine the global PETSc index of a local node by
     * making use of d_node_offset.  Local node d_local_lag_indices[ln][j] has
     * global PETSc index j+d_node_offset[ln].  A similar mapping for nonlocal
     * nodes is not well defined.
     */
    std::vector<std::vector<int> > d_local_lag_indices;
    std::vector<std::vector<int> > d_nonlocal_lag_indices;

    /*!
     * The node indices of all local nodes (i.e. the nodes owned by this
     * processor) on each level of the hierarchy.  The indices are in the global
     * PETSc ordering corresponding to a depth of 1.
     */
    std::vector<std::vector<int> > d_local_petsc_indices;

    /*!
     * The node indices of all nonlocal nodes (i.e. the nodes owned by another
     * processor that appear in the ghost region of some patch owned by this
     * processor) on each level of the hierarchy.  The indices are in the global
     * PETSc ordering corresponding to a depth of 1.
     */
    std::vector<std::vector<int> > d_nonlocal_petsc_indices;

    /*!
     * Container for additional user defined Lagrangian data
     */
    std::map<std::string, int> d_user_defined_ldata;

    //\}
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/LDataManager-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LDataManager
