// ---------------------------------------------------------------------
//
// Copyright (c) 2025 by Vinod Thale
// Part of IBAMR Vinod Extensions
//
// ---------------------------------------------------------------------

#ifndef included_MultiStructureForceTracker
#define included_MultiStructureForceTracker

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/IBHydrodynamicForceEvaluator.h>

#include <ibtk/ibtk_utilities.h>

#include <PatchHierarchy.h>
#include <tbox/Pointer.h>
#include <tbox/Database.h>
#include <appu/VisItDataWriter.h>
#include <solv/RobinBcCoefStrategy.h>

#include <string>
#include <vector>

namespace VINOD
{
/*!
 * \brief MultiStructureForceTracker manages hydrodynamic force evaluation
 * for multiple immersed structures.
 *
 * This utility class simplifies the setup and management of control volumes
 * for computing hydrodynamic forces and torques on multiple immersed bodies.
 * It eliminates code duplication when working with multi-body simulations.
 *
 * Key Features:
 * - Automatic control volume registration for multiple structures
 * - Input validation with clear error messages
 * - Simplified force and torque computation
 * - Center of mass tracking
 *
 * Example Usage:
 * \code
 *   Pointer<MultiStructureForceTracker> force_tracker =
 *       new MultiStructureForceTracker("force_tracker", rho, mu, start_time);
 *
 *   // Register all structures from input database
 *   force_tracker->registerStructuresFromDatabase(
 *       input_db, "InitHydroForceBox", num_structures, patch_hierarchy);
 *
 *   // Set torque origins from center of mass
 *   force_tracker->setTorqueOriginsFromCOM(structure_COM);
 *
 *   // Register for visualization
 *   force_tracker->registerVisualization(visit_data_writer, patch_hierarchy);
 *
 *   // During time stepping:
 *   force_tracker->updateAllControlVolumes(COM_velocities, dt, patch_hierarchy);
 *   force_tracker->computeAllForces(u_idx, p_idx, patch_hierarchy, dt, bc_u, bc_p);
 * \endcode
 */
class MultiStructureForceTracker
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name Name for this force tracker
     * \param rho_fluid Fluid density
     * \param mu_fluid Fluid dynamic viscosity
     * \param start_time Initial simulation time
     * \param use_adaptive_control_volumes Whether to use adaptive CVs (default: true)
     */
    MultiStructureForceTracker(const std::string& object_name,
                              double rho_fluid,
                              double mu_fluid,
                              double start_time,
                              bool use_adaptive_control_volumes = true);

    /*!
     * \brief Destructor.
     */
    ~MultiStructureForceTracker();

    /*!
     * \brief Register control volumes for all structures from input database.
     *
     * Reads control volume configuration from input database with naming pattern:
     * db_name_prefix + "_" + structure_id
     * Example: "InitHydroForceBox_0", "InitHydroForceBox_1", etc.
     *
     * \param input_db Input database containing CV configurations
     * \param db_name_prefix Prefix for CV database names
     * \param num_structures Number of structures
     * \param patch_hierarchy Patch hierarchy
     */
    void registerStructuresFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        const std::string& db_name_prefix,
        int num_structures,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set torque origins from center of mass positions.
     *
     * \param structure_COM Vector of COM positions for each structure [num_structures][3]
     */
    void setTorqueOriginsFromCOM(const std::vector<std::vector<double> >& structure_COM);

    /*!
     * \brief Register structures for visualization output.
     *
     * \param visit_data_writer VisIt data writer
     * \param patch_hierarchy Patch hierarchy
     */
    void registerVisualization(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_data_writer,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Update control volumes for all structures.
     *
     * Updates CV positions based on structure velocities. Handles adaptive
     * control volume logic to keep structures inside CVs.
     *
     * \param COM_velocities Center of mass velocities [num_structures][NDIM]
     * \param dt Timestep size
     * \param patch_hierarchy Patch hierarchy
     * \param coarse_grid_spacing Grid spacing on coarsest level (for adaptive CVs)
     */
    void updateAllControlVolumes(
        const std::vector<std::vector<double> >& COM_velocities,
        double dt,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
        const double* coarse_grid_spacing);

    /*!
     * \brief Compute hydrodynamic forces for all structures.
     *
     * \param u_idx Velocity patch data index
     * \param p_idx Pressure patch data index
     * \param patch_hierarchy Patch hierarchy
     * \param dt Timestep size
     * \param velocity_bc Velocity boundary conditions
     * \param pressure_bc Pressure boundary conditions
     */
    void computeAllForces(
        int u_idx,
        int p_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
        double dt,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* pressure_bc);

    /*!
     * \brief Get the underlying IBHydrodynamicForceEvaluator.
     *
     * Provides access to the IBAMR force evaluator for advanced usage.
     */
    SAMRAI::tbox::Pointer<IBAMR::IBHydrodynamicForceEvaluator> getForceEvaluator() const;

    /*!
     * \brief Get number of structures being tracked.
     */
    int getNumStructures() const;

private:
    /*!
     * \brief Object name.
     */
    std::string d_object_name;

    /*!
     * \brief Number of structures.
     */
    int d_num_structures = 0;

    /*!
     * \brief Underlying IBAMR force evaluator.
     */
    SAMRAI::tbox::Pointer<IBAMR::IBHydrodynamicForceEvaluator> d_hydro_force;

    /*!
     * \brief Track displacement for adaptive control volumes.
     */
    std::vector<double> d_box_disp;

    /*!
     * \brief Maximum box velocity (for adaptive CV logic).
     */
    IBTK::Vector3d d_max_box_vel;

    /*!
     * \brief Use adaptive control volumes.
     */
    bool d_use_adaptive_cvs;
};

} // namespace VINOD

#endif // #ifndef included_MultiStructureForceTracker
