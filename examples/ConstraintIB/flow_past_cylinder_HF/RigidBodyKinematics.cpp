// Filename : RigidBodyKinematics.cpp
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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
// POSSIBILITY OF SUCH DAMAGE

//////////////////////////// INCLUDES /////////////////////////////////////////
#include <cmath>
#include <fstream>
#include <iostream>

#include "RigidBodyKinematics.h"
#include "PatchLevel.h"
#include "CartesianPatchGeometry.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include "ibamr/namespaces.h"
#include "muParser.h"

namespace IBAMR
{
//////////////////////////////////////////////////////////////////////////////

RigidBodyKinematics::RigidBodyKinematics(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         LDataManager* l_data_manager,
                                         Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                         bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_current_time(0.0),
      d_kinematics_data(),
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
                                 << " not found in restart file."
                                 << std::endl);
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
