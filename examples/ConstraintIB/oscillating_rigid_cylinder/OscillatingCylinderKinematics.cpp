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

#include "OscillatingCylinderKinematics.h"

/////////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include "ibtk/IBTK_MPI.h"

#include "CartesianPatchGeometry.h"
#include "PatchLevel.h"
#include "tbox/MathUtilities.h"
#include "tbox/PIO.h"
#include "tbox/Utilities.h"

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// IBTK INCLUDES

// IBTK THIRD PARTY INCLUDES
#include "muParser.h"

// C++ INCLUDES
#include <limits>
#include <string>

namespace IBAMR
{
OscillatingCylinderKinematics::OscillatingCylinderKinematics(const std::string& object_name,
                                                             Pointer<Database> input_db,
                                                             LDataManager* l_data_manager,
                                                             Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                                                             bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_prescribed_trans_vel(0.0),
      d_current_time(0.0)

{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name.

    // Get kinematic parameters from the input file

    if (input_db->keyExists("frequency"))
    {
        d_freq = input_db->getDouble("frequency");
        pout << "FREQUENCY ================= " << d_freq << "\n\n\n";
    }
    else
    {
        TBOX_ERROR(
            " OscillatingCylinderKinematics::OscillatingCylinderKinematics() :: frequency of cylinder does not exist "
            "in the InputDatabase \n\n"
            << std::endl);
    }

    if (input_db->keyExists("U_infinity"))
    {
        d_Uinf = input_db->getDouble("U_infinity");
        pout << "U_INFINITY ================= " << d_Uinf << "\n\n\n";
    }
    else
    {
        TBOX_ERROR(
            " OscillatingCylinderKinematics::OscillatingCylinderKinematics() :: maximum speed of cylinder does not "
            "exist in the InputDatabase \n\n"
            << std::endl);
    }

    // Set the size of vectors.
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    const int total_levels = finest_ln - coarsest_ln + 1;
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();

    if (total_levels > 1)
    {
        TBOX_ERROR(
            " OscillatingCylinderKinematics::OscillatingCylinderKinematics() :: Total levels > 1, OscillatingCylinder "
            "should be places on only one level \n\n"
            << std::endl);
    }

    // d_vec_coord.resize(number_lag_points);
    // d_vec_amp.resize(number_lag_points);
    d_new_kinematics_vel.resize(total_levels);
    d_current_kinematics_vel.resize(total_levels);

    for (int ln = 0; ln < total_levels; ++ln)
    {
        const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
        d_new_kinematics_vel[ln].resize(NDIM);
        d_current_kinematics_vel[ln].resize(NDIM);

        for (int d = 0; d < NDIM; ++d)
        {
            d_new_kinematics_vel[ln][d].resize(nodes_this_ln);
            d_current_kinematics_vel[ln][d].resize(nodes_this_ln);
        }
    }

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
    }

    // set current velocity using current time for initial and restarted runs.
    setOscillatingCylinderSpecificVelocity(d_current_time);

    return;

} // OscillatingCylinderKinematics

OscillatingCylinderKinematics::~OscillatingCylinderKinematics()
{
    // intentionally left blank

    return;

} //~OscillatingCylinderKinematics

void
OscillatingCylinderKinematics::getFromRestart()
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

    return;

} // getFromRestart

void
OscillatingCylinderKinematics::putToDatabase(Pointer<Database> db)
{
    IBAMR::ConstraintIBKinematics::putToDatabase(db);
    db->putDouble("d_current_time", d_current_time);
    return;

} // putToDatabase

void
OscillatingCylinderKinematics::setOscillatingCylinderSpecificVelocity(const double time)
{
    // Set the size of vectors.
    const StructureParameters& struct_param = getStructureParameters();
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    const int number_lag_points = idx_range[0].second - idx_range[0].first;

    for (int k = 0; k < number_lag_points; ++k)
    {
        d_current_kinematics_vel[0][0][k] = d_prescribed_trans_vel + d_Uinf * std::cos(2 * M_PI * d_freq * time);
    }

    return;

} // setOscillatingCylinderSpecificVelocity

void
OscillatingCylinderKinematics::setKinematicsVelocity(
    const double new_time,
    const std::vector<double>& /*incremented_angle_from_reference_axis*/,
    const std::vector<double>& /*center_of_mass*/,
    const std::vector<double>& /*tagged_pt_position*/)
{
    d_new_time = new_time;
    // fill current velocity at new time
    if (!MathUtilities<double>::equalEps(0.0, d_new_time)) setOscillatingCylinderSpecificVelocity(d_new_time);

    // swap current and new velocity
    if (MathUtilities<double>::equalEps(0.0, d_new_time))
        d_new_kinematics_vel = d_current_kinematics_vel;
    else
        d_new_kinematics_vel.swap(d_current_kinematics_vel);

    d_current_time = d_new_time;

    return;

} // setNewKinematicsVelocity

const std::vector<std::vector<double> >&
OscillatingCylinderKinematics::getKinematicsVelocity(const int level) const
{
    static const StructureParameters& struct_param = getStructureParameters();
    static const int coarsest_ln = struct_param.getCoarsestLevelNumber();

    const int offset = level - coarsest_ln;
    return d_new_kinematics_vel[offset];

} // getNewKinematicsVelocity

const std::vector<std::vector<double> >&
OscillatingCylinderKinematics::getCurrentKinematicsVelocity(const int level) const
{
    static const StructureParameters& struct_param = getStructureParameters();
    static const int coarsest_ln = struct_param.getCoarsestLevelNumber();

    const int offset = level - coarsest_ln;
    return d_current_kinematics_vel[offset];

} // getCurrentKinematicsVelocity

void
OscillatingCylinderKinematics::setShape(const double /*new_time*/, const std::vector<double>&)
{
    // intentionally left blank
    return;

} // setNewShape

const std::vector<std::vector<double> >&
OscillatingCylinderKinematics::getShape(const int /*level*/) const
{
    return d_new_shape;

} // getNewShape

} // namespace IBAMR
