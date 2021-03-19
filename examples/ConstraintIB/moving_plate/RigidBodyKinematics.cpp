// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

//////////////////////////// INCLUDES /////////////////////////////////////////
#include "ibamr/namespaces.h"

#include "ibtk/IBTK_MPI.h"

#include "CartesianPatchGeometry.h"
#include "PatchLevel.h"
#include "RigidBodyKinematics.h"
#include "tbox/MathUtilities.h"

#include "muParser.h"

#include <cmath>
#include <fstream>
#include <iostream>

namespace IBAMR
{
//////////////////////////////////////////////////////////////////////////////

RigidBodyKinematics::RigidBodyKinematics(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         LDataManager* l_data_manager,
                                         Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                         bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_kinematics_data(),
      d_current_time(0.0),
      d_new_time(0.0),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_center_of_mass(3),
      d_incremented_angle_from_reference_axis(3),
      d_tagged_pt_position(3),
      d_mesh_width(NDIM)
{
    setImmersedBodyLayout(patch_hierarchy);
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();

        // At restart current and new states point to the same state.
        d_new_time = d_current_time;
        setRigidBodyVelocity(
            d_current_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);
        setShape(d_current_time, d_incremented_angle_from_reference_axis);
    }
    return;

} // RigidBodyKinematics

RigidBodyKinematics::~RigidBodyKinematics()
{
    return;
} // ~RigidBodyKinematics

void
RigidBodyKinematics::registerRigidBodyKinematics(RigidVelFcn trans_vel_fcn, RigidVelFcn rot_vel_fcn)
{
    d_kinematics_data.d_trans_vel_fcn = trans_vel_fcn;
    d_kinematics_data.d_rot_vel_fcn = rot_vel_fcn;

    return;
}

void
RigidBodyKinematics::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase

void RigidBodyKinematics::setImmersedBodyLayout(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/)
{
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    const int total_lag_pts = idx_range[0].second - idx_range[0].first;

    for (int d = 0; d < NDIM; ++d)
    {
        d_kinematics_vel[d].resize(total_lag_pts);
    }

    return;
}

void
RigidBodyKinematics::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }

    d_current_time = db->getDouble("d_current_time");
    db->getDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->getDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->getDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;
} // getFromRestart

void
RigidBodyKinematics::setRigidBodyVelocity(const double time,
                                          const std::vector<double>& /*incremented_angle_from_reference_axis*/,
                                          const std::vector<double>& /*center_of_mass*/,
                                          const std::vector<double>& /*tagged_pt_position*/)
{
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    const int total_lag_pts = idx_range[0].second - idx_range[0].first;

    Eigen::Vector3d trans_vel;
    trans_vel.setZero();
    if (d_kinematics_data.d_trans_vel_fcn)
    {
        d_kinematics_data.d_trans_vel_fcn(time, trans_vel);
    }

    for (int d = 0; d < NDIM; ++d)
    {
        for (int k = 0; k < total_lag_pts; ++k)
        {
            d_kinematics_vel[d][k] = trans_vel[d];
        }
    }

    return;
} // setRigidBodyVelocity

void
RigidBodyKinematics::setKinematicsVelocity(const double time,
                                           const std::vector<double>& incremented_angle_from_reference_axis,
                                           const std::vector<double>& center_of_mass,
                                           const std::vector<double>& tagged_pt_position)
{
    d_new_time = time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;

    setRigidBodyVelocity(d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

    return;

} // setNewKinematicsVelocity

const std::vector<std::vector<double> >&
RigidBodyKinematics::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getKinematicsVelocity

void
RigidBodyKinematics::setShape(const double time, const std::vector<double>& /*incremented_angle_from_reference_axis*/)
{
    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    if (position_update_method == "CONSTRAINT_VELOCITY")
    {
        d_current_time = d_new_time;
        return;
    }
    TBOX_ASSERT(d_new_time == time);
    d_current_time = d_new_time;

    return;
} // setShape

const std::vector<std::vector<double> >&
RigidBodyKinematics::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
