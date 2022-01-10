// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2022 by the IBAMR developers
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

// SAMRAI INCLUDES
#include "ibtk/IBTK_MPI.h"

#include <tbox/PIO.h>
#include <tbox/Utilities.h>

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

// C++ INCLUDES
#include <string>

// Application
#include "KnifeFishKinematics.h"

namespace IBAMR
{
namespace
{
static const bool DISCARD_COMMENTS = false;
static const double PII = 3.1415926535897932384626433832795;

inline std::string
discard_comments(const std::string& input_string)
{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    std::string output_string = input_string;
    std::istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    std::getline(string_stream, output_string, '%');
    string_stream.clear();
    return output_string;
} // discard_comments

} // namespace

KnifeFishKinematics::KnifeFishKinematics(const std::string& object_name,
                                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                         IBTK::LDataManager* l_data_manager,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                                         bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_fin_start_idx(0),
      d_current_time(0.0)
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name.

    // Get the angular frequency, wave number, level on which fin is placed & angle of excursion from the input data
    // base
    if (input_db->keyExists("length_fin"))
    {
        d_fin_length = input_db->getDouble("length_fin");
    }
    else
    {
        TBOX_ERROR(
            " KnifefishKinematics::KnifefishKinematics() :: Length of fin does not exist in the InputDatabase \n\n"
            << std::endl);
    }

    if (input_db->keyExists("numberofwaves"))
    {
        d_kappa = 2 * PII / (d_fin_length / (input_db->getDouble("numberofwaves")));
    }
    else
    {
        TBOX_ERROR(
            " KnifefishKinematics::KnifefishKinematics() :: number of waves of fin does not exist in the InputDatabase "
            "\n\n"
            << std::endl);
    }

    if (input_db->keyExists("theta_max"))
    {
        d_theta_max = input_db->getDouble("theta_max") * PII / 180.0;
    }
    else
    {
        TBOX_ERROR(
            " KnifefishKinematics::KnifefishKinematics() :: Theta_max for fin does not exist in the InputDatabase \n\n"
            << std::endl);
    }

    if (input_db->keyExists("frequency"))
    {
        d_omega = 2 * PII * input_db->getDouble("frequency");
    }
    else
    {
        TBOX_ERROR(
            " KnifefishKinematics::KnifefishKinematics() :: frequency of fin does not exist in the InputDatabase \n\n"
            << std::endl);
    }

    // Set the size of vectors.
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();

    if (input_db->keyExists("fin_starting_index")) d_fin_start_idx = input_db->getInteger("fin_starting_index");

    const int nodes_fin = idx_range[0].second - idx_range[0].first;
    const int number_fin_points = nodes_fin - d_fin_start_idx;

    // resize some vectors.
    d_vec_coord.resize(number_fin_points);
    d_vec_radius.resize(number_fin_points);
    d_vec_theta.resize(number_fin_points);
    for (int d = 0; d < NDIM; ++d)
    {
        d_kinematics_vel[d].resize(nodes_fin);
        d_shape[d].resize(nodes_fin);
    }

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();

    // get the name of file which has coordinates of the sphere.
    std::fstream coord_file_stream;
    std::string dir_name = input_db->getStringWithDefault("vertex_file_dir", "./");
    std::string filename = dir_name + "/" + input_db->getString("radius_amp");

    if (filename.empty())
    {
        TBOX_ERROR("KnifeFishKinemtics::KnifeFishKinematics() "
                   << "ERROR :: base file name containing coordinates of the rigid structure not found in input file"
                   << std::endl);
    }

    coord_file_stream.open(filename.c_str(), std::fstream::in);
    if (!coord_file_stream.is_open())
    {
        TBOX_ERROR("KnifeFishKinemtics::KnifeFishKinematics() "
                   << "could not open file" << filename << std::endl);
    }

    for (int k = 0; k < number_fin_points; ++k)
    {
        std::string line_string;

        if (!std::getline(coord_file_stream, line_string))
        {
            TBOX_ERROR("KnifeFishKinemtics::KnifeFishKinematics() "
                       << ":\n  Premature end to input file encountered before line " << k + 1 << " of file "
                       << filename << std::endl);
        }
        else
        {
            if (DISCARD_COMMENTS) line_string = discard_comments(line_string);
            std::istringstream line_stream(line_string);

            if (!(line_stream >> d_vec_coord[k]))
            {
                TBOX_ERROR("KnifeFishKinemtics::KnifeFishKinematics() "
                           << ":\n Invalid entry in input file encountered on line " << k + 2 << " of file " << filename
                           << " . Lagrangian X co-ordinate of material point at " << NDIM + 1 << "expected."
                           << std::endl);
            }

            if (!(line_stream >> d_vec_radius[k]))
            {
                TBOX_ERROR("KnifeFishKinemtics::KnifeFishKinematics() "
                           << ":\n Invalid entry in input file encountered on line " << k + 2 << " of file " << filename
                           << " .Radius of material point at " << NDIM + 1 << "expected." << std::endl);
            }

            if (!(line_stream >> d_vec_theta[k]))
            {
                TBOX_ERROR("KnifeFishKinemtics::KnifeFishKinematics() "
                           << ":\n Invalid entry in input file encountered on line " << k + 2 << " of file " << filename
                           << " .Angle of excursion at " << NDIM + 1 << "expected." << std::endl);
            }
        }
    } // end of k for loop
    coord_file_stream.close();

    return;

} // KnifeFishKinematics

KnifeFishKinematics::~KnifeFishKinematics()
{
    // intentionally left blank
    return;

} //~KnifeFishKinematics

void
KnifeFishKinematics::getFromRestart()
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
KnifeFishKinematics::putToDatabase(Pointer<Database> db)
{
    IBAMR::ConstraintIBKinematics::putToDatabase(db);
    db->putDouble("d_current_time", d_current_time);
    return;

} // putToDatabase

void
KnifeFishKinematics::setKnifefishSpecificVelocity(const double time)
{
    static const StructureParameters& struct_param = getStructureParameters();
    static const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    static const int nodes_fin = idx_range[0].second - idx_range[0].first;

    // Vector iterators
    std::vector<double>::iterator it_coord = d_vec_coord.begin(), it_theta = d_vec_theta.begin(),
                                  it_radius = d_vec_radius.begin();
    std::vector<double>::iterator /*it_vel_x = d_kinematics_vel[0].begin() + d_fin_start_idx,*/
        it_vel_y = d_kinematics_vel[1].begin() + d_fin_start_idx,
        it_vel_z = d_kinematics_vel[2].begin() + d_fin_start_idx;

    for (int k = d_fin_start_idx; k < nodes_fin;
         ++k, ++it_coord, ++it_theta, ++it_radius, /*++it_vel_x,*/ ++it_vel_y, ++it_vel_z)
    {
        const double X = (*it_coord);
        const double Theta = (*it_theta) * std::sin(d_kappa * X - d_omega * time);

        //*it_vel_x = 0.0;
        *it_vel_y = (*it_radius) * std::cos(Theta) * (*it_theta) * std::cos(d_kappa * X - d_omega * time) * (-d_omega);
        *it_vel_z = (*it_radius) * std::sin(Theta) * (*it_theta) * std::cos(d_kappa * X - d_omega * time) * (-d_omega);
    }

    return;

} // setKnifefishSpecificVelocity

void
KnifeFishKinematics::setKinematicsVelocity(const double new_time,
                                           const std::vector<double>& /*incremented_angle_from_reference_axis*/,
                                           const std::vector<double>& /*center_of_mass*/,
                                           const std::vector<double>& /*tagged_pt_position*/)
{
    d_new_time = new_time;
    // fill current velocity at new time
    if (!IBTK::abs_equal_eps(0.0, d_new_time)) setKnifefishSpecificVelocity(d_new_time);

    d_current_time = d_new_time;

    return;

} // setKinematicsVelocity

const std::vector<std::vector<double> >&
KnifeFishKinematics::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getKinematicsVelocity

void
KnifeFishKinematics::setShape(const double /*time*/,
                              const std::vector<double>& /*incremented_angle_from_reference_axis*/)
{
    // intentionally left blank
    return;

} // setShape

const std::vector<std::vector<double> >&
KnifeFishKinematics::getShape(const int /*level*/) const
{
    return d_shape;

} // getNewShape

} // namespace IBAMR
