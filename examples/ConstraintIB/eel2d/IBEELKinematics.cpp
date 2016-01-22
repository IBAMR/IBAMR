// Filename : IBEELKinematics.cpp
// Created by Amneet Bhalla on 6/16/2011.

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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

//////////////////////////// INCLUDES /////////////////////////////////////////
#include <cmath>
#include <fstream>
#include <iostream>

#include "IBEELKinematics.h"
#include "PatchLevel.h"
#include "CartesianPatchGeometry.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include "ibamr/namespaces.h"
#include "muParser.h"

namespace IBAMR
{
namespace
{
inline int
sign(const double X)
{
    return ((X > 0) ? 1 : ((X < 0) ? -1 : 0));
}

static const double PII = 3.1415926535897932384626433832795;
static const double __INFINITY = 1e9;

// Set Fish Related Parameters.
static const double LENGTH_FISH = 1.0;
static const double WIDTH_HEAD = 0.04 * LENGTH_FISH;
static const double LENGTH_HEAD = 0.04;

// prey capturing parameters.
static const double CUT_OFF_ANGLE = PII / 4;
static const double CUT_OFF_RADIUS = 0.7;
static const double LOWER_CUT_OFF_ANGLE = 7 * PII / 180;

} // namespace unknown

///////////////////////////////////////////////////////////////////////

IBEELKinematics::IBEELKinematics(const std::string& object_name,
                                 Pointer<Database> input_db,
                                 LDataManager* l_data_manager,
                                 Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                 bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_current_time(0.0),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_center_of_mass(3),
      d_incremented_angle_from_reference_axis(3),
      d_tagged_pt_position(3),
      d_mesh_width(NDIM),
      d_parser_time(new double),
      d_parser_posn(new double[NDIM]),
      d_parser_normal(new double[NDIM])
{
    // Read from inputdb
    d_initAngle_bodyAxis_x = input_db->getDoubleWithDefault("initial_angle_body_axis_0", 0.0);
    d_bodyIsManeuvering = input_db->getBoolWithDefault("body_is_maneuvering", false);
    d_maneuverAxisIsChangingShape = input_db->getBoolWithDefault("maneuvering_axis_is_changing_shape", false);

    // Read-in deformation velocity functions
    std::vector<std::string> deformationvel_function_strings;
    for (int d = 0; d < NDIM; ++d)
    {
        std::ostringstream stream;
        stream << "_function_" << d;
        const std::string postfix = stream.str();
        std::string key_name = "deformation_velocity" + postfix;

        if (input_db->isString(key_name))
        {
            deformationvel_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            deformationvel_function_strings.push_back("0.0");
            TBOX_WARNING("IBEELKinematics::IBEELKinematics() :\n"
                         << "  no function corresponding to key ``"
                         << key_name
                         << " '' found for dimension = "
                         << d
                         << "; using def_vel = 0.0. "
                         << std::endl);
        }

        d_deformationvel_parsers.push_back(new mu::Parser());
        d_deformationvel_parsers.back()->SetExpr(deformationvel_function_strings.back());
        d_all_parsers.push_back(d_deformationvel_parsers.back());
    }

    // Read-in the body shape parser
    {
        const std::string body_shape_equation = input_db->getString("body_shape_equation");
        d_body_shape_parser = new mu::Parser();
        d_body_shape_parser->SetExpr(body_shape_equation);
        d_all_parsers.push_back(d_body_shape_parser);
    }

    // Read in the maneuvering axis parser.
    if (d_bodyIsManeuvering)
    {
        const std::string maneuvering_axis_equation = input_db->getString("maneuvering_axis_equation");
        d_maneuvering_axis_parser = new mu::Parser();
        d_maneuvering_axis_parser->SetExpr(maneuvering_axis_equation);
        d_all_parsers.push_back(d_maneuvering_axis_parser);
    }

    // Define the default and the user-provided constants.
    const double pi = 3.1415926535897932384626433832795;
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {
        // Various names for pi.
        (*cit)->DefineConst("pi", pi);
        (*cit)->DefineConst("Pi", pi);
        (*cit)->DefineConst("PI", pi);

        // Variables
        (*cit)->DefineVar("T", d_parser_time);
        (*cit)->DefineVar("t", d_parser_time);
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();
            (*cit)->DefineVar("X" + postfix, &(d_parser_posn[d]));
            (*cit)->DefineVar("x" + postfix, &(d_parser_posn[d]));
            (*cit)->DefineVar("X_" + postfix, &(d_parser_posn[d]));
            (*cit)->DefineVar("x_" + postfix, &(d_parser_posn[d]));

            (*cit)->DefineVar("N" + postfix, &(d_parser_normal[d]));
            (*cit)->DefineVar("n" + postfix, &(d_parser_normal[d]));
            (*cit)->DefineVar("N_" + postfix, &(d_parser_normal[d]));
            (*cit)->DefineVar("n_" + postfix, &(d_parser_normal[d]));
        }
    }

    // set the location of the food particle from the input file.
    d_food_location.resizeArray(NDIM);
    for (int dim = 0; dim < NDIM; ++dim)
    {
        std::ostringstream stream;
        stream << "food_location_in_domain_" << dim;
        d_food_location[dim] = input_db->getDouble(stream.str());
    }

    // set how the immersed body is layout in reference frame.
    setImmersedBodyLayout(patch_hierarchy);

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();

        // At restart current and new states point to the same state.
        d_new_time = d_current_time;
        setEelSpecificVelocity(
            d_current_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);
        setShape(d_current_time, d_incremented_angle_from_reference_axis);
    }
    return;

} // IBEELKinematics

IBEELKinematics::~IBEELKinematics()
{
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {
        delete (*cit);
    }
    delete d_parser_time;
    delete[] d_parser_posn;
    delete[] d_parser_normal;

    return;

} // ~IBEELKinematics

void
IBEELKinematics::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase

void
IBEELKinematics::getFromRestart()
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
IBEELKinematics::setImmersedBodyLayout(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    // Set some vector sizes.
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    const int total_lag_pts = idx_range[0].second - idx_range[0].first;

    for (int d = 0; d < NDIM; ++d)
    {
        d_kinematics_vel[d].resize(total_lag_pts);
        d_shape[d].resize(total_lag_pts);
    }

    // Get Background mesh related data.
    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(finest_ln);
    PatchLevel<NDIM>::Iterator p(level);
    Pointer<Patch<NDIM> > patch = level->getPatch(p());
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (int dim = 0; dim < NDIM; ++dim)
    {
        d_mesh_width[dim] = dx[dim];
    }

    // No. of points on the backbone and till head.
    const int BodyNx = int(ceil(LENGTH_FISH / d_mesh_width[0]));
    const int HeadNx = int(ceil(LENGTH_HEAD / d_mesh_width[0]));

    d_ImmersedBodyData.clear();
    for (int i = 1; i <= HeadNx; ++i)
    {
        const double s = (i - 1) * d_mesh_width[0];
        const double section = sqrt(2 * WIDTH_HEAD * s - s * s);
        const int NumPtsInSection = 2 * int(ceil(section / d_mesh_width[1]));
        d_ImmersedBodyData.insert(std::make_pair(s, NumPtsInSection));
    }

    for (int i = HeadNx + 1; i <= BodyNx; ++i)
    {
        const double s = (i - 1) * d_mesh_width[0];
        const double section = WIDTH_HEAD * (LENGTH_FISH - s) / (LENGTH_FISH - LENGTH_HEAD);
        const int NumPtsInHeight = 2 * int(ceil(section / d_mesh_width[1]));
        d_ImmersedBodyData.insert(std::make_pair(s, NumPtsInHeight));
    }

    // Find the coordinates of the axis of maneuvering in the reference frame from the input file.
    if (d_bodyIsManeuvering)
    {
        d_maneuverAxisReferenceCoordinates_vec.clear();
        d_map_reference_tangent.clear();
        d_map_reference_sign.clear();

        std::vector<double> vec_axis_coord(2);
        for (std::map<double, int>::const_iterator mitr = d_ImmersedBodyData.begin(); mitr != d_ImmersedBodyData.end();
             ++mitr)
        {
            d_parser_posn[0] = mitr->first;
            vec_axis_coord[0] = mitr->first;
            vec_axis_coord[1] = d_maneuvering_axis_parser->Eval();
            d_maneuverAxisReferenceCoordinates_vec.push_back(vec_axis_coord);
        }

        // store the tangents to the reference maneuver axis in a reference_map.
        for (unsigned int i = 0; i <= (d_maneuverAxisReferenceCoordinates_vec.size() - 2); ++i)
        {
            std::vector<int> sign_vec(2);
            const double s = d_maneuverAxisReferenceCoordinates_vec[i][0];
            const double dX =
                (d_maneuverAxisReferenceCoordinates_vec[i + 1][0] - d_maneuverAxisReferenceCoordinates_vec[i][0]);
            const double dY =
                (d_maneuverAxisReferenceCoordinates_vec[i + 1][1] - d_maneuverAxisReferenceCoordinates_vec[i][1]);
            sign_vec[0] = sign(dX);
            sign_vec[1] = sign(dY);
            const double theta = std::atan(std::abs(dY / dX));
            d_map_reference_tangent.insert(std::make_pair(s, theta));
            d_map_reference_sign.insert(std::make_pair(s, sign_vec));
        }

        // Fill in the last point in the map.
        d_map_reference_tangent.insert(std::make_pair((d_maneuverAxisReferenceCoordinates_vec.back())[0],
                                                      (d_map_reference_tangent.rbegin())->second));
        d_map_reference_sign.insert(std::make_pair((d_maneuverAxisReferenceCoordinates_vec.back())[0],
                                                   (d_map_reference_sign.rbegin())->second));

        // Find the COM of the maneuver axis.
        double maneuverAxis_x_cm = 0.0;
        double maneuverAxis_y_cm = 0.0;
        for (unsigned int i = 0; i < d_maneuverAxisReferenceCoordinates_vec.size(); ++i)
        {
            maneuverAxis_x_cm += d_maneuverAxisReferenceCoordinates_vec[i][0];
            maneuverAxis_y_cm += d_maneuverAxisReferenceCoordinates_vec[i][1];
        }
        maneuverAxis_x_cm /= d_maneuverAxisReferenceCoordinates_vec.size();
        maneuverAxis_y_cm /= d_maneuverAxisReferenceCoordinates_vec.size();

        // Shift the reference so that maneuver Axis coordinate COM coincides with the origin.
        for (unsigned int i = 0; i < d_maneuverAxisReferenceCoordinates_vec.size(); ++i)
        {
            d_maneuverAxisReferenceCoordinates_vec[i][0] -= maneuverAxis_x_cm;
            d_maneuverAxisReferenceCoordinates_vec[i][1] -= maneuverAxis_y_cm;
        }
    } // body is maneuvering

    return;

} // setImmersedBodyLayout

void
IBEELKinematics::transformManeuverAxisAndCalculateTangents(const double angleFromHorizontal)
{
    d_maneuverAxisTransformedCoordinates_vec.clear();
    d_map_transformed_tangent.clear();
    d_map_transformed_sign.clear();

    const int BodyNx = int(ceil(LENGTH_FISH / d_mesh_width[0]));
    std::vector<double> transformed_coord(2);
    for (int i = 0; i <= (BodyNx - 1); ++i)
    {
        transformed_coord[0] = d_maneuverAxisReferenceCoordinates_vec[i][0] * cos(angleFromHorizontal) -
                               d_maneuverAxisReferenceCoordinates_vec[i][1] * sin(angleFromHorizontal);
        transformed_coord[1] = d_maneuverAxisReferenceCoordinates_vec[i][0] * sin(angleFromHorizontal) +
                               d_maneuverAxisReferenceCoordinates_vec[i][1] * cos(angleFromHorizontal);
        d_maneuverAxisTransformedCoordinates_vec.push_back(transformed_coord);
    }

    for (int i = 0; i <= (BodyNx - 2); ++i)
    {
        std::vector<int> sign_vec(2);
        const double s = i * d_mesh_width[0];
        const double dX =
            (d_maneuverAxisTransformedCoordinates_vec[i + 1][0] - d_maneuverAxisTransformedCoordinates_vec[i][0]);
        const double dY =
            (d_maneuverAxisTransformedCoordinates_vec[i + 1][1] - d_maneuverAxisTransformedCoordinates_vec[i][1]);
        sign_vec[0] = sign(dX);
        sign_vec[1] = sign(dY);
        const double theta = std::atan(std::abs(dY / dX));
        d_map_transformed_tangent.insert(std::make_pair(s, theta));
        d_map_transformed_sign.insert(std::make_pair(s, sign_vec));
    }

    // Fill in the last point in the map.
    d_map_transformed_tangent.insert(
        std::make_pair((BodyNx - 1) * d_mesh_width[0], (d_map_transformed_tangent.rbegin())->second));
    d_map_transformed_sign.insert(
        std::make_pair((BodyNx - 1) * d_mesh_width[0], (d_map_transformed_sign.rbegin())->second));

    return;

} // transformManeuverAxisAndCalculateTangents

void
IBEELKinematics::setEelSpecificVelocity(const double time,
                                        const std::vector<double>& incremented_angle_from_reference_axis,
                                        const std::vector<double>& center_of_mass,
                                        const std::vector<double>& tagged_pt_position)
{
    *d_parser_time = time;
    const double angleFromHorizontal = d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];

    if (d_bodyIsManeuvering)
    {
        if (d_maneuverAxisIsChangingShape)
        {
            // calculate the radius of the circular path on which the fish will have its backbone.
            double radius_circular_path;
            std::vector<double> bodyline_vector(NDIM), foodline_vector(NDIM);
            double mag_bodyline_vector = 0.0, mag_foodline_vector = 0.0;

            for (int dim = 0; dim < NDIM; ++dim)
            {
                bodyline_vector[dim] = tagged_pt_position[dim] - center_of_mass[dim];
                foodline_vector[dim] = d_food_location[dim] - tagged_pt_position[dim];
                mag_bodyline_vector += std::pow(bodyline_vector[dim], 2);
                mag_foodline_vector += std::pow(foodline_vector[dim], 2);
            }

            // Normalize the vectors.
            for (int dim = 0; dim < NDIM; ++dim)
            {
                bodyline_vector[dim] /= sqrt(mag_bodyline_vector);
                foodline_vector[dim] /= sqrt(mag_foodline_vector);
            }

            // Find the angle between bodyline_axis and foodline_axis
            // angle = sign(aXb)* acos(a.b/|a||b|)
            const double angle_bw_target_vision =
                sign(bodyline_vector[0] * foodline_vector[1] - bodyline_vector[1] * foodline_vector[0]) *
                std::acos(bodyline_vector[0] * foodline_vector[0] + bodyline_vector[1] * foodline_vector[1]);

            if (angle_bw_target_vision >= CUT_OFF_ANGLE)
            {
                radius_circular_path = CUT_OFF_RADIUS;
            }
            else if (angle_bw_target_vision <= -CUT_OFF_ANGLE)
            {
                radius_circular_path = CUT_OFF_RADIUS;
            }
            else if (MathUtilities<double>::equalEps(MathUtilities<double>::Abs(angle_bw_target_vision), 0.0))
            {
                radius_circular_path = __INFINITY;
            }
            else if (angle_bw_target_vision >= -LOWER_CUT_OFF_ANGLE && angle_bw_target_vision <= LOWER_CUT_OFF_ANGLE)
            {
                radius_circular_path = std::abs(CUT_OFF_RADIUS * std::pow((CUT_OFF_ANGLE / LOWER_CUT_OFF_ANGLE), 1));
            }
            else
            {
                radius_circular_path = std::abs(CUT_OFF_RADIUS * std::pow((CUT_OFF_ANGLE / angle_bw_target_vision), 1));
            }
            // set the reference maneuver axis coordinates.
            const int BodyNx = int(ceil(LENGTH_FISH / d_mesh_width[0]));
            if (radius_circular_path != __INFINITY)
            {
                const double angle_sector = LENGTH_FISH / radius_circular_path;
                const double dtheta = angle_sector / (BodyNx - 1);

                d_maneuverAxisReferenceCoordinates_vec.clear();
                std::vector<double> vec_axis_coord(2);
                for (int i = 1; i <= BodyNx; ++i)
                {
                    const double angleFromVertical = -angle_sector / 2 + (i - 1) * dtheta;
                    vec_axis_coord[0] = radius_circular_path * sin(angleFromVertical);
                    vec_axis_coord[1] = radius_circular_path * cos(angleFromVertical);
                    d_maneuverAxisReferenceCoordinates_vec.push_back(vec_axis_coord);
                }
            }
            else
            {
                d_maneuverAxisReferenceCoordinates_vec.clear();
                std::vector<double> vec_axis_coord(2);
                for (int i = 1; i <= BodyNx; ++i)
                {
                    vec_axis_coord[0] = (i - 1) * d_mesh_width[0];
                    vec_axis_coord[1] = 0.0;
                    d_maneuverAxisReferenceCoordinates_vec.push_back(vec_axis_coord);
                }
            }

            // Find the COM of the maneuver axis.
            double maneuverAxis_x_cm = 0.0;
            double maneuverAxis_y_cm = 0.0;
            for (unsigned int i = 0; i < d_maneuverAxisReferenceCoordinates_vec.size(); ++i)
            {
                maneuverAxis_x_cm += d_maneuverAxisReferenceCoordinates_vec[i][0];
                maneuverAxis_y_cm += d_maneuverAxisReferenceCoordinates_vec[i][1];
            }
            maneuverAxis_x_cm /= d_maneuverAxisReferenceCoordinates_vec.size();
            maneuverAxis_y_cm /= d_maneuverAxisReferenceCoordinates_vec.size();

            // Shift the reference so that maneuver Axis coordinate COM coincides with the origin.
            for (unsigned int i = 0; i < d_maneuverAxisReferenceCoordinates_vec.size(); ++i)
            {
                d_maneuverAxisReferenceCoordinates_vec[i][0] -= maneuverAxis_x_cm;
                d_maneuverAxisReferenceCoordinates_vec[i][1] -= maneuverAxis_y_cm;
            }

            // Find the tangents on this reference axis for shape update.
            d_map_reference_tangent.clear();
            d_map_reference_sign.clear();
            for (int i = 0; i <= (BodyNx - 2); ++i)
            {
                std::vector<int> sign_vec(2);
                const double s = i * d_mesh_width[0];
                const double dX =
                    (d_maneuverAxisReferenceCoordinates_vec[i + 1][0] - d_maneuverAxisReferenceCoordinates_vec[i][0]);
                const double dY =
                    (d_maneuverAxisReferenceCoordinates_vec[i + 1][1] - d_maneuverAxisReferenceCoordinates_vec[i][1]);
                sign_vec[0] = sign(dX);
                sign_vec[1] = sign(dY);
                const double theta = std::atan(std::abs(dY / dX));
                d_map_reference_tangent.insert(std::make_pair(s, theta));
                d_map_reference_sign.insert(std::make_pair(s, sign_vec));
            }
            // Fill in the last point in the map.
            d_map_reference_tangent.insert(
                std::make_pair((BodyNx - 1) * d_mesh_width[0], (d_map_reference_tangent.rbegin())->second));
            d_map_reference_sign.insert(
                std::make_pair((BodyNx - 1) * d_mesh_width[0], (d_map_reference_sign.rbegin())->second));
        } // maneuverAxisIsChangingShape

        // Rotate the reference axis and calculate tangents in the rotated frame.
        transformManeuverAxisAndCalculateTangents(angleFromHorizontal);
    } // bodyIsManeuvering

    // Set the deformation velocity in the body frame.
    std::vector<double> vec_vel(NDIM);
    int lag_idx = 0;
    for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    {
        d_parser_posn[0] = itr->first;
        const int NumPtsInSection = itr->second;

        if (d_bodyIsManeuvering)
        {
            d_parser_normal[0] =
                -sin(d_map_transformed_tangent[d_parser_posn[0]]) * d_map_transformed_sign[d_parser_posn[0]][1];
            d_parser_normal[1] =
                cos(d_map_transformed_tangent[d_parser_posn[0]]) * d_map_transformed_sign[d_parser_posn[0]][0];
        }
        else
        {
            d_parser_normal[0] = -sin(angleFromHorizontal);
            d_parser_normal[1] = cos(angleFromHorizontal);
        }

        vec_vel[0] = d_deformationvel_parsers[0]->Eval();
        vec_vel[1] = d_deformationvel_parsers[1]->Eval();

        const int lowerlimit = lag_idx;
        const int upperlimit = lag_idx + NumPtsInSection;
        for (int d = 0; d < NDIM; ++d)
        {
            for (int i = lowerlimit; i < upperlimit; ++i) d_kinematics_vel[d][i] = vec_vel[d];
        }

        lag_idx = upperlimit;
    }

    return;
} // setEelSpecificVelocity

void
IBEELKinematics::setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position)
{
    d_new_time = time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;

    setEelSpecificVelocity(d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

    return;

} // setNewKinematicsVelocity

const std::vector<std::vector<double> >&
IBEELKinematics::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getKinematicsVelocity

void
IBEELKinematics::setShape(const double time, const std::vector<double>& /*incremented_angle_from_reference_axis*/)
{
    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    if (position_update_method == "CONSTRAINT_VELOCITY") return;

    // Find the deformed shape. Rotate the shape about center of mass.
    TBOX_ASSERT(d_new_time == time);
    *d_parser_time = time;
    std::vector<double> shape_new(NDIM);

    int lag_idx = -1;
    int reference_axis_idx = -1;
    for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    {
        const int NumPtsInSection = itr->second;
        d_parser_posn[0] = itr->first;
        const double y_shape_base = d_body_shape_parser->Eval();

        if (d_bodyIsManeuvering)
        {
            ++reference_axis_idx;
            const double x_maneuver_base = d_maneuverAxisReferenceCoordinates_vec[reference_axis_idx][0];
            const double y_maneuver_base = d_maneuverAxisReferenceCoordinates_vec[reference_axis_idx][1];

            for (int j = 1; j <= NumPtsInSection / 2; ++j)
            {
                const double nx = (-1 * sin(d_map_reference_tangent[itr->first]) * d_map_reference_sign[itr->first][1]);
                const double ny = (cos(d_map_reference_tangent[itr->first]) * d_map_reference_sign[itr->first][0]);

                shape_new[0] = x_maneuver_base + (y_shape_base + (j - 1) * d_mesh_width[1]) * nx;
                shape_new[1] = y_maneuver_base + (y_shape_base + (j - 1) * d_mesh_width[1]) * ny;

                d_shape[0][++lag_idx] = shape_new[0];
                d_shape[1][lag_idx] = shape_new[1];
            }

            for (int j = 1; j <= NumPtsInSection / 2; ++j)
            {
                const double nx = (-1 * sin(d_map_reference_tangent[itr->first]) * d_map_reference_sign[itr->first][1]);
                const double ny = (cos(d_map_reference_tangent[itr->first]) * d_map_reference_sign[itr->first][0]);

                shape_new[0] = x_maneuver_base + (y_shape_base - (j)*d_mesh_width[1]) * nx;
                shape_new[1] = y_maneuver_base + (y_shape_base - (j)*d_mesh_width[1]) * ny;

                d_shape[0][++lag_idx] = shape_new[0];
                d_shape[1][lag_idx] = shape_new[1];
            }
        } // bodyIsManeuvering.
        else
        {
            for (int j = 1; j <= NumPtsInSection / 2; ++j)
            {
                d_shape[0][++lag_idx] = itr->first;
                d_shape[1][lag_idx] = y_shape_base + (j - 1) * d_mesh_width[1];
            }

            for (int j = 1; j <= NumPtsInSection / 2; ++j)
            {
                d_shape[0][++lag_idx] = itr->first;
                d_shape[1][lag_idx] = y_shape_base - j * d_mesh_width[1];
            }
        }
    }

    // Find the c.m of this new shape.
    std::vector<double> center_of_mass(NDIM, 0.0);
    const int total_lag_pts = d_shape[0].size();
    for (int d = 0; d < NDIM; ++d)
    {
        for (std::vector<double>::const_iterator citr = d_shape[d].begin(); citr != d_shape[d].end(); ++citr)
        {
            center_of_mass[d] += *citr;
        }
    }

    for (int d = 0; d < NDIM; ++d) center_of_mass[d] /= total_lag_pts;

    // Shift the c.m to the origin to apply the rotation
    for (int d = 0; d < NDIM; ++d)
    {
        for (std::vector<double>::iterator itr = d_shape[d].begin(); itr != d_shape[d].end(); ++itr)
        {
            *itr -= center_of_mass[d];
        }
    }

    // Now rotate the shape about origin or center of mass.
    const double angleFromHorizontal = d_initAngle_bodyAxis_x + d_incremented_angle_from_reference_axis[2];
    for (int i = 0; i < total_lag_pts; ++i)
    {
        const double x_rotated = d_shape[0][i] * cos(angleFromHorizontal) - d_shape[1][i] * sin(angleFromHorizontal);
        const double y_rotated = d_shape[0][i] * sin(angleFromHorizontal) + d_shape[1][i] * cos(angleFromHorizontal);
        d_shape[0][i] = x_rotated;
        d_shape[1][i] = y_rotated;
    }

    d_current_time = d_new_time;

    return;
} // setShape

const std::vector<std::vector<double> >&
IBEELKinematics::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
