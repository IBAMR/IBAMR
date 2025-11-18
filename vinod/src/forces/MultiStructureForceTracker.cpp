// ---------------------------------------------------------------------
//
// Copyright (c) 2025 by Vinod Thale
// Part of IBAMR Vinod Extensions
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "vinod/forces/MultiStructureForceTracker.h"

#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>
#include <tbox/Utilities.h>
#include <tbox/SAMRAI_MPI.h>

#include <cmath>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace VINOD
{

/////////////////////////////// PUBLIC ///////////////////////////////////////

MultiStructureForceTracker::MultiStructureForceTracker(const std::string& object_name,
                                                      double rho_fluid,
                                                      double mu_fluid,
                                                      double start_time,
                                                      bool use_adaptive_control_volumes)
    : d_object_name(object_name),
      d_use_adaptive_cvs(use_adaptive_control_volumes)
{
    // Create the underlying IBAMR force evaluator
    // Note: 5th parameter is register_for_restart (not adaptive CVs)
    // Adaptive CV logic is handled in updateAllControlVolumes()
    d_hydro_force = new IBAMR::IBHydrodynamicForceEvaluator(
        object_name, rho_fluid, mu_fluid, start_time, /*register_for_restart*/ true);

    // Initialize maximum box velocity
    d_max_box_vel.setZero();

    return;
} // Constructor

MultiStructureForceTracker::~MultiStructureForceTracker()
{
    // SAMRAI smart pointer handles cleanup
    return;
} // Destructor

void
MultiStructureForceTracker::registerStructuresFromDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    const std::string& db_name_prefix,
    int num_structures,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy)
{
    d_num_structures = num_structures;
    d_box_disp.resize(num_structures, 0.0);

    SAMRAI::tbox::plog << "MultiStructureForceTracker: Registering " << num_structures
                       << " structures from database prefix '" << db_name_prefix << "'" << std::endl;

    for (int struct_id = 0; struct_id < num_structures; ++struct_id)
    {
        // Construct database name
        const std::string box_db_name = db_name_prefix + "_" + std::to_string(struct_id);

        SAMRAI::tbox::plog << "  Structure " << struct_id << ": " << box_db_name << std::endl;

        // Input validation: check that database exists
        if (!input_db->isDatabase(box_db_name))
        {
            TBOX_ERROR("MultiStructureForceTracker::registerStructuresFromDatabase()\n"
                      << "  Missing control volume configuration: " << box_db_name << "\n"
                      << "  Each structure requires " << db_name_prefix << "_N database in input file.");
        }

        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> box_db = input_db->getDatabase(box_db_name);

        // Input validation: check that required keys exist
        if (!box_db->keyExists("lower_left_corner"))
        {
            TBOX_ERROR("MultiStructureForceTracker::registerStructuresFromDatabase()\n"
                      << "  Missing 'lower_left_corner' in " << box_db_name);
        }
        if (!box_db->keyExists("upper_right_corner"))
        {
            TBOX_ERROR("MultiStructureForceTracker::registerStructuresFromDatabase()\n"
                      << "  Missing 'upper_right_corner' in " << box_db_name);
        }
        if (!box_db->keyExists("init_velocity"))
        {
            TBOX_ERROR("MultiStructureForceTracker::registerStructuresFromDatabase()\n"
                      << "  Missing 'init_velocity' in " << box_db_name);
        }

        // Read control volume bounds and initial velocity
        IBTK::Vector3d box_X_lower, box_X_upper, box_init_vel;
        box_db->getDoubleArray("lower_left_corner", &box_X_lower[0], 3);
        box_db->getDoubleArray("upper_right_corner", &box_X_upper[0], 3);
        box_db->getDoubleArray("init_velocity", &box_init_vel[0], 3);

        // Input validation: check that lower < upper
        for (int d = 0; d < 3; ++d)
        {
            if (box_X_lower[d] >= box_X_upper[d])
            {
                TBOX_ERROR("MultiStructureForceTracker::registerStructuresFromDatabase()\n"
                          << "  Invalid control volume for " << box_db_name << "\n"
                          << "  lower_left_corner[" << d << "] = " << box_X_lower[d]
                          << " >= upper_right_corner[" << d << "] = " << box_X_upper[d] << "\n"
                          << "  Control volume bounds must satisfy lower < upper.");
            }
        }

        // Register structure with force evaluator
        d_hydro_force->registerStructure(box_X_lower, box_X_upper, patch_hierarchy, box_init_vel, struct_id);

        SAMRAI::tbox::plog << "    CV bounds: [" << box_X_lower[0] << ", " << box_X_lower[1] << ", " << box_X_lower[2]
                          << "] to [" << box_X_upper[0] << ", " << box_X_upper[1] << ", " << box_X_upper[2] << "]"
                          << std::endl;
    }

    SAMRAI::tbox::plog << "MultiStructureForceTracker: Successfully registered " << num_structures
                       << " structures" << std::endl;

    return;
} // registerStructuresFromDatabase

void
MultiStructureForceTracker::setTorqueOriginsFromCOM(const std::vector<std::vector<double> >& structure_COM)
{
    if (static_cast<int>(structure_COM.size()) != d_num_structures)
    {
        TBOX_ERROR("MultiStructureForceTracker::setTorqueOriginsFromCOM()\n"
                  << "  structure_COM size (" << structure_COM.size()
                  << ") does not match number of structures (" << d_num_structures << ")");
    }

    SAMRAI::tbox::plog << "MultiStructureForceTracker: Setting torque origins from COM" << std::endl;

    for (int struct_id = 0; struct_id < d_num_structures; ++struct_id)
    {
        IBTK::Vector3d torque_origin;
        for (int d = 0; d < 3; ++d)
        {
            torque_origin[d] = structure_COM[struct_id][d];
        }

        d_hydro_force->setTorqueOrigin(torque_origin, struct_id);

        SAMRAI::tbox::plog << "  Structure " << struct_id << " torque origin: ("
                          << torque_origin[0] << ", " << torque_origin[1] << ", " << torque_origin[2] << ")"
                          << std::endl;
    }

    return;
} // setTorqueOriginsFromCOM

void
MultiStructureForceTracker::registerVisualization(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_data_writer,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy)
{
    for (int struct_id = 0; struct_id < d_num_structures; ++struct_id)
    {
        d_hydro_force->registerStructurePlotData(visit_data_writer, patch_hierarchy, struct_id);
    }

    SAMRAI::tbox::plog << "MultiStructureForceTracker: Registered " << d_num_structures
                       << " structures for visualization" << std::endl;

    return;
} // registerVisualization

void
MultiStructureForceTracker::updateAllControlVolumes(
    const std::vector<std::vector<double> >& COM_velocities,
    double dt,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
    const double* coarse_grid_spacing)
{
    if (static_cast<int>(COM_velocities.size()) != d_num_structures)
    {
        TBOX_ERROR("MultiStructureForceTracker::updateAllControlVolumes()\n"
                  << "  COM_velocities size (" << COM_velocities.size()
                  << ") does not match number of structures (" << d_num_structures << ")");
    }

    // Update control volume for each structure
    for (int struct_id = 0; struct_id < d_num_structures; ++struct_id)
    {
        IBTK::Vector3d box_vel;
        box_vel.setZero();

        // Get velocity due to structure motion
        for (int d = 0; d < NDIM; ++d)
        {
            box_vel(d) = COM_velocities[struct_id][d];
        }

        if (d_use_adaptive_cvs && coarse_grid_spacing != nullptr)
        {
            // Adaptive control volume logic to keep structure inside CV
            // If COM has moved 0.9 coarse mesh widths, translate CV by 1 coarse mesh width

            d_box_disp[struct_id] += box_vel[0] * dt;

            if (std::abs(d_box_disp[struct_id]) >= std::abs(0.9 * coarse_grid_spacing[0]))
            {
                // Reset CV position
                box_vel.setZero();
                box_vel[0] = -coarse_grid_spacing[0] / dt;

                d_box_disp[struct_id] = 0.0;
            }
            else
            {
                // Keep CV stationary
                box_vel.setZero();
            }
        }

        // Update control volume position for this structure
        d_hydro_force->updateStructureDomain(box_vel, dt, patch_hierarchy, struct_id);
    }

    return;
} // updateAllControlVolumes

void
MultiStructureForceTracker::computeAllForces(
    int u_idx,
    int p_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
    double dt,
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& velocity_bc,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* pressure_bc)
{
    // Compute lagged momentum integral (MPI collective operation)
    d_hydro_force->computeLaggedMomentumIntegral(u_idx, patch_hierarchy, velocity_bc);

    // Compute hydrodynamic forces and torques for all structures
    d_hydro_force->computeHydrodynamicForce(
        u_idx, p_idx, /*f_idx*/ -1, patch_hierarchy, dt, velocity_bc, pressure_bc);

    return;
} // computeAllForces

SAMRAI::tbox::Pointer<IBAMR::IBHydrodynamicForceEvaluator>
MultiStructureForceTracker::getForceEvaluator() const
{
    return d_hydro_force;
} // getForceEvaluator

int
MultiStructureForceTracker::getNumStructures() const
{
    return d_num_structures;
} // getNumStructures

} // namespace VINOD

//////////////////////////////////////////////////////////////////////////////
