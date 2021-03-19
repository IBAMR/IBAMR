// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include "RigidBodyKinematics.h"

/////////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include "ibtk/IBTK_MPI.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

// IBAMR INCLUDES
#include "ibamr/namespaces.h"

// IBTK INCLUDES

// IBTK THIRD-PARTY INCLUDES
#include "muParser.h"

// C++ INCLUDES
#include <string>

namespace IBAMR
{
namespace
{
static const double PII = 3.14159265358979323846264338327950288419716939937510;

} // namespace

RigidBodyKinematics::RigidBodyKinematics(const std::string& object_name,
                                         Pointer<Database> input_db,
                                         LDataManager* l_data_manager,
                                         Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                                         bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_parser_time(0.0),
      d_center_of_mass(3, 0.0),
      d_incremented_angle_from_reference_axis(3, 0.0),
      d_tagged_pt_position(3, 0.0)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name.

    // Read-in kinematics velocity functions
    for (int d = 0; d < NDIM; ++d)
    {
        const std::string postfix = "_function_" + std::to_string(d);
        std::string key_name = "kinematics_velocity" + postfix;

        if (input_db->isString(key_name))
        {
            d_kinematicsvel_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            d_kinematicsvel_function_strings.push_back("0.0");
            TBOX_WARNING("RigidBodyKinematics::RigidBodyKinematics() :\n"
                         << "  no function corresponding to key " << key_name << "found for dimension = " << d
                         << "; using kinematics_vel = 0.0. " << std::endl);
        }

        d_kinematicsvel_parsers.push_back(new mu::Parser());
        d_kinematicsvel_parsers.back()->SetExpr(d_kinematicsvel_function_strings.back());
        d_all_parsers.push_back(d_kinematicsvel_parsers.back());
    }

    // Define constants and variables for the parsers.
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {
        // Various names for pi.
        (*cit)->DefineConst("pi", PII);
        (*cit)->DefineConst("Pi", PII);
        (*cit)->DefineConst("PI", PII);

        // Variables
        (*cit)->DefineVar("T", &d_parser_time);
        (*cit)->DefineVar("t", &d_parser_time);
        for (int d = 0; d < NDIM; ++d)
        {
            const std::string postfix = std::to_string(d);
            (*cit)->DefineVar("X" + postfix, d_parser_posn.data() + d);
            (*cit)->DefineVar("x" + postfix, d_parser_posn.data() + d);
            (*cit)->DefineVar("X_" + postfix, d_parser_posn.data() + d);
            (*cit)->DefineVar("x_" + postfix, d_parser_posn.data() + d);
        }
    }

    // Set the size of vectors.
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    const int total_levels = finest_ln - coarsest_ln + 1;
    d_kinematics_vel.resize(total_levels);

    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    for (int ln = 0; ln < total_levels; ++ln)
    {
        const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
        d_kinematics_vel[ln].resize(NDIM);
        for (int d = 0; d < NDIM; ++d)
        {
            d_kinematics_vel[ln][d].resize(nodes_this_ln);
        }
    }

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();
        setRigidBodySpecificVelocity(
            d_current_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);
        setShape(d_current_time, d_incremented_angle_from_reference_axis);
    }

    return;

} // RigidBodyKinematics

RigidBodyKinematics::~RigidBodyKinematics()
{
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {
        delete (*cit);
    }

    return;

} //~RigidBodyKinematics

void
RigidBodyKinematics::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase

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
    d_new_time = d_current_time;
    db->getDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->getDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->getDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // getFromRestart

void
RigidBodyKinematics::setRigidBodySpecificVelocity(const double time,
                                                  const std::vector<double>& /*incremented_angle_from_reference_axis*/,
                                                  const std::vector<double>& center_of_mass,
                                                  const std::vector<double>& /*tagged_pt_position*/)
{
    std::vector<double> vel_parser(NDIM);
    d_parser_time = time;
    for (int d = 0; d < NDIM; ++d) d_parser_posn[d] = center_of_mass[d];
    for (int d = 0; d < NDIM; ++d) vel_parser[d] = d_kinematicsvel_parsers[d]->Eval();

    static const StructureParameters& struct_param = getStructureParameters();
    static const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    static const int finest_ln = struct_param.getFinestLevelNumber();
    static const int total_levels = finest_ln - coarsest_ln + 1;
    static const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();

    for (int ln = 0; ln < total_levels; ++ln)
    {
        const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
        for (int d = 0; d < NDIM; ++d)
        {
            for (int idx = 0; idx < nodes_this_ln; ++idx) d_kinematics_vel[ln][d][idx] = vel_parser[d];
        }
    }

    return;

} // setRigidBodySpecificVelocity

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

    setRigidBodySpecificVelocity(
        d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

    d_current_time = d_new_time;
    return;

} // setNewKinematicsVelocity

const std::vector<std::vector<double> >&
RigidBodyKinematics::getKinematicsVelocity(const int level) const
{
    static const StructureParameters& struct_param = getStructureParameters();
    static const int coarsest_ln = struct_param.getCoarsestLevelNumber();

#ifdef DEBUG_CHECK_ASSERTIONS
    static const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln <= level && level <= finest_ln);
#endif

    return d_kinematics_vel[level - coarsest_ln];

} // getKinematicsVelocity

void
RigidBodyKinematics::setShape(const double /*time*/,
                              const std::vector<double>& /*incremented_angle_from_reference_axis*/)
{
    // intentionally left blank
    return;

} // setShape

const std::vector<std::vector<double> >&
RigidBodyKinematics::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
