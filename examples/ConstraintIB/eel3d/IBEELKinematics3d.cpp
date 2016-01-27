// Filename : IBEELKinematics3d.cpp
// Created by Amneet Bhalla on 05/19/2012.

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
// POSSIBILITY OF SUCH DAMAGE.

//////////////////////////////////// INCLUDES ////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include <sstream>

#include "IBEELKinematics3d.h"
#include "PatchLevel.h"
#include "CartesianPatchGeometry.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/MathUtilities.h"
#include "tbox/Utilities.h"
#include "ibamr/namespaces.h"
#include "muParser.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_integration.h"

namespace IBAMR
{
namespace
{
static const double PII = 3.1415926535897932384626433832795;

// Set Fish Related Parameters.
static const double LENGTH_FISH = 1.0;
static const double WIDTH_HEAD = 0.04 * LENGTH_FISH;
static const double LENGTH_TILLHEAD = 0.04 * LENGTH_FISH;
static const double MAJOR_AXIS = 0.51 * LENGTH_FISH;
static const double MINOR_AXIS = 0.08 * LENGTH_FISH;

// If kappa is the curvature of backbone,i.e, kappa = d/ds (alpha); where s
// is the arclength of the backbone:
// then x(s,t) = integral_0^s {cos(alpha(s,t)) ds }
// and  y(s,t) = integral_0^s {sin(alpha(s,t)) ds }
//
// Deformational velocity in frame of fish:
// d/dt[x(s,t)] = integral_0^s { d/dt[cos(alpha(s,t))] ds }
// d/dt[y(s,t)] = integral_0^s { d/dt[sin(alpha(s,t))] ds }

double
yVelocity(double s, void* params)
{
    double* input = static_cast<double*>(params);
    double a0 = input[0];
    double a1 = input[1];
    double a2 = input[2];
    double a3 = input[3];
    double tau = input[4];
    double t = input[5];
    double T = input[6];

    double dydt =
        std::cos((1 / (8 * std::pow(PII * tau, 4))) * (std::tanh(PII * t / T)) *
                 (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
                  2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                      cos((2 * PII * (t - s * T * tau)) / T) +
                  (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
                  (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                      std::sin((2 * PII * (t - s * T * tau)) / T))) *
        ((1 / (8 * std::pow(PII * tau, 4))) * ((2 * PII) / T) * (std::tanh(PII * t / T)) *
             (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::sin((2 * PII * t) / T) * (-1) +
              2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                  std::sin((2 * PII * (t - s * T * tau)) / T) * (-1) +
              (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
              (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                  std::cos((2 * PII * (t - s * T * tau)) / T)) +
         (1 / (8 * std::pow(PII * tau, 4))) * (PII / T) * (std::pow(std::cosh(PII * t / T), -2)) *
             (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
              2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                  cos((2 * PII * (t - s * T * tau)) / T) +
              (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
              (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                  std::sin((2 * PII * (t - s * T * tau)) / T))

             );

    return dydt;

} // yVelocity

double
xVelocity(double s, void* params)
{
    double* input = static_cast<double*>(params);
    double a0 = input[0];
    double a1 = input[1];
    double a2 = input[2];
    double a3 = input[3];
    double tau = input[4];
    double t = input[5];
    double T = input[6];

    double dxdt =
        -std::sin((1 / (8 * std::pow(PII * tau, 4))) * (std::tanh(PII * t / T)) *
                  (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
                   2 * PII * tau *
                       (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                       cos((2 * PII * (t - s * T * tau)) / T) +
                   (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
                   (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                       std::sin((2 * PII * (t - s * T * tau)) / T))) *
        ((1 / (8 * std::pow(PII * tau, 4))) * ((2 * PII) / T) * (std::tanh(PII * t / T)) *
             (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::sin((2 * PII * t) / T) * (-1) +
              2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                  std::sin((2 * PII * (t - s * T * tau)) / T) * (-1) +
              (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
              (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                  std::cos((2 * PII * (t - s * T * tau)) / T)) +
         (1 / (8 * std::pow(PII * tau, 4))) * (PII / T) * (std::pow(std::cosh(PII * t / T), -2)) *
             (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
              2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                  cos((2 * PII * (t - s * T * tau)) / T) +
              (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
              (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                  std::sin((2 * PII * (t - s * T * tau)) / T))

             );

    return dxdt;

} // xVelocity

double
xPosition(double s, void* params)
{
    double* input = static_cast<double*>(params);
    double a0 = input[0];
    double a1 = input[1];
    double a2 = input[2];
    double a3 = input[3];
    double tau = input[4];
    double t = input[5];
    double T = input[6];

    double x =
        std::cos((1 / (8 * std::pow(PII * tau, 4))) * (std::tanh(PII * t / T)) *
                 (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
                  2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                      cos((2 * PII * (t - s * T * tau)) / T) +
                  (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
                  (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                      std::sin((2 * PII * (t - s * T * tau)) / T)));

    return x;

} // xposition

double
yPosition(double s, void* params)
{
    double* input = static_cast<double*>(params);
    double a0 = input[0];
    double a1 = input[1];
    double a2 = input[2];
    double a3 = input[3];
    double tau = input[4];
    double t = input[5];
    double T = input[6];

    double y =
        std::sin((1 / (8 * std::pow(PII * tau, 4))) * (std::tanh(PII * t / T)) *
                 (2 * PII * tau * (a2 - 2 * a0 * std::pow(PII * tau, 2)) * std::cos((2 * PII * t) / T) +
                  2 * PII * tau * (-a2 - 3 * a3 * s + 2 * PII * PII * (a0 + s * (a1 + s * (a2 + a3 * s))) * tau * tau) *
                      cos((2 * PII * (t - s * T * tau)) / T) +
                  (3 * a3 - 2 * a1 * std::pow(PII * tau, 2)) * sin((2 * PII * t) / T) +
                  (-3 * a3 + 2 * PII * PII * (a1 + 2 * a2 * s + 3 * a3 * s * s) * tau * tau) *
                      std::sin((2 * PII * (t - s * T * tau)) / T)));

    return y;

} // yPosition

} // namespace anonymous

IBEELKinematics3d::IBEELKinematics3d(const std::string& object_name,
                                     Pointer<Database> input_db,
                                     LDataManager* l_data_manager,
                                     Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                     bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_mesh_width(NDIM),
      d_current_time(0.0),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_center_of_mass(3),
      d_incremented_angle_from_reference_axis(3),
      d_tagged_pt_position(3),
      d_interp_coefs()
{
    // NOTE: Parent class constructor registers class with the restart manager, sets object name.

    // Get the curvature of the body from the input file.
    d_interp_coefs = input_db->getDoubleArray("interp_coefs");
    d_tau_tail = input_db->getDouble("tau_tail");
    d_time_period = input_db->getDouble("time_period");
    d_initAngle_bodyAxis_x = input_db->getDoubleWithDefault("initial_angle_horizontal", 0.0);

    // set how the immersed body is layout in reference frame.
    setImmersedBodyLayout(patch_hierarchy);

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart)
    {
        getFromRestart();

        // Set the new state equal to current state.
        d_new_time = d_current_time;
        setEelSpecificVelocity(
            d_current_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);
        setShape(d_current_time, d_incremented_angle_from_reference_axis);
    }

    return;

} // IBEELKinematics3d

IBEELKinematics3d::~IBEELKinematics3d()
{
    return;

} // ~IBEELKinematics3d

void
IBEELKinematics3d::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase

void
IBEELKinematics3d::getFromRestart()
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
IBEELKinematics3d::setImmersedBodyLayout(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
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
        d_shape[d].resize(total_lag_pts);
    }

    // Get Background mesh related data.
    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(coarsest_ln);
    PatchLevel<NDIM>::Iterator p(level);
    Pointer<Patch<NDIM> > patch = level->getPatch(p());
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (int dim = 0; dim < NDIM; ++dim)
    {
        d_mesh_width[dim] = dx[dim];
    }

    // No. of points on the backbone and till head.
    d_HeadNs = int(ceil(LENGTH_TILLHEAD / d_mesh_width[0]));
    d_TailNs = int(ceil((LENGTH_FISH - LENGTH_TILLHEAD) / d_mesh_width[0]));
    d_BodyNs = d_HeadNs + d_TailNs + 1;

    d_IBPts.resize(d_BodyNs);
    d_IBWidthHeight.resize(d_BodyNs);

    for (int i = 1; i <= d_HeadNs + 1; ++i)
    {
        const double s = (i - 1) * d_mesh_width[0];
        const double section = sqrt(2 * WIDTH_HEAD * s - s * s);
        const double height = MINOR_AXIS * std::sqrt(1 - pow((s - MAJOR_AXIS) / MAJOR_AXIS, 2));
        const int NumPtsInSection = int(ceil(section / d_mesh_width[1]));
        const int NumPtsInHeight = int(ceil(height / d_mesh_width[2]));
        d_IBPts[i - 1] = std::make_pair(NumPtsInSection, NumPtsInHeight);
        d_IBWidthHeight[i - 1] = std::make_pair(section, height);
    }

    for (int i = d_HeadNs + 2; i <= d_BodyNs; ++i)
    {
        const double s = (i - 1) * d_mesh_width[0];
        const double section = WIDTH_HEAD * (LENGTH_FISH - s) / (LENGTH_FISH - LENGTH_TILLHEAD);
        const double height = MINOR_AXIS * std::sqrt(1 - pow((s - MAJOR_AXIS) / MAJOR_AXIS, 2));
        const int NumPtsInSection = int(ceil(section / d_mesh_width[1]));
        const int NumPtsInHeight = int(ceil(height / d_mesh_width[2]));
        d_IBPts[i - 1] = std::make_pair(NumPtsInSection, NumPtsInHeight);
        d_IBWidthHeight[i - 1] = std::make_pair(section, height);
    }

    return;

} // setImmersedBodyLayout

void
IBEELKinematics3d::setEelSpecificVelocity(const double time,
                                          const std::vector<double>& incremented_angle_from_reference_axis,
                                          const std::vector<double>& /*center_of_mass*/,
                                          const std::vector<double>& /*tagged_pt_position*/)
{
    const double angleFromHorizontal = d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];
    std::vector<double> vec_vel(NDIM, 0.0);
    int lag_idx = 0;

    double input[7];
    double dydt, dxdt, errory, errorx;
    size_t nevalsy, nevalsx;
    gsl_function Fx, Fy;
    Fx.function = xVelocity;
    Fx.params = input;
    Fy.function = yVelocity;
    Fy.params = input;

    input[0] = d_interp_coefs[0];
    input[1] = d_interp_coefs[1];
    input[2] = d_interp_coefs[2];
    input[3] = d_interp_coefs[3];
    input[4] = d_tau_tail;
    input[5] = time;
    input[6] = d_time_period;

    for (int i = 1; i <= d_BodyNs; ++i)
    {
        const int NumPtsInSection = d_IBPts[i - 1].first;
        const int NumPtsInHeight = d_IBPts[i - 1].second;

        if (NumPtsInSection && NumPtsInHeight)
        {
            const double width = d_IBWidthHeight[i - 1].first;
            const double depth = d_IBWidthHeight[i - 1].second;
            const double S = (i - 1) * d_mesh_width[0];

            gsl_integration_qng(&Fx, 0, S, 1e-8, 0.0, &dxdt, &errorx, &nevalsx);
            gsl_integration_qng(&Fy, 0, S, 1e-8, 0.0, &dydt, &errory, &nevalsy);

            vec_vel[0] = dxdt * (std::cos(angleFromHorizontal)) + dydt * (-std::sin(angleFromHorizontal));
            vec_vel[1] = dydt * (std::cos(angleFromHorizontal)) + dxdt * (std::sin(angleFromHorizontal));

            int pts_this_xsection = 2 * NumPtsInHeight + 1;
            for (int j = 1; j <= NumPtsInSection; ++j)
            {
                const double y = j * d_mesh_width[1];
                for (int k = -NumPtsInHeight; k <= NumPtsInHeight; ++k)
                {
                    const double z = k * d_mesh_width[2];
                    if ((std::pow(y / width, 2) + std::pow(z / depth, 2)) <= 1) // use elliptical cross sections
                        pts_this_xsection += 2;
                }
            } // cross section filled

            const int lowerlimit = lag_idx;
            const int upperlimit = lag_idx + pts_this_xsection;
            for (int d = 0; d < NDIM; ++d)
                for (int j = lowerlimit; j < upperlimit; ++j) d_kinematics_vel[d][j] = vec_vel[d];

            lag_idx = upperlimit;
        }
    }

    return;

} // setEelSpecificvelocity

void
IBEELKinematics3d::setKinematicsVelocity(const double new_time,
                                         const std::vector<double>& incremented_angle_from_reference_axis,
                                         const std::vector<double>& center_of_mass,
                                         const std::vector<double>& tagged_pt_position)
{
    d_new_time = new_time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;

    setEelSpecificVelocity(d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

    return;

} // setKinematicsVelocity

const std::vector<std::vector<double> >&
IBEELKinematics3d::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getNewKinematicsVelocity

void
IBEELKinematics3d::setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis)
{
    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    if (position_update_method == "CONSTRAINT_VELOCITY")
    {
        return;
    }
    else
    {
        TBOX_ASSERT(d_new_time == time);
        double input[7];
        input[0] = d_interp_coefs[0];
        input[1] = d_interp_coefs[1];
        input[2] = d_interp_coefs[2];
        input[3] = d_interp_coefs[3];
        input[4] = d_tau_tail;
        input[5] = d_new_time;
        input[6] = d_time_period;

        gsl_function Fx, Fy;
        Fx.function = xPosition;
        Fx.params = input;
        Fy.function = yPosition;
        Fy.params = input;

        double ybase, xbase, errory, errorx;
        size_t nevalsy, nevalsx;

        // Find the deformed shape. Rotate the shape about center of mass.
        std::vector<double> shape_new(NDIM);
        int lag_idx = -1;
        for (int i = 1; i <= d_BodyNs; ++i)
        {
            const int NumPtsInSection = d_IBPts[i - 1].first;
            const int NumPtsInHeight = d_IBPts[i - 1].second;

            if (NumPtsInSection && NumPtsInHeight)
            {
                const double width = d_IBWidthHeight[i - 1].first;
                const double depth = d_IBWidthHeight[i - 1].second;
                const double S = (i - 1) * d_mesh_width[0];

                gsl_integration_qng(&Fx, 0, S, 1e-8, 0.0, &xbase, &errorx, &nevalsx);
                gsl_integration_qng(&Fy, 0, S, 1e-8, 0.0, &ybase, &errory, &nevalsy);

                // Fill the middle line first.
                for (int k = -NumPtsInHeight; k <= NumPtsInHeight; ++k)
                {
                    d_shape[0][++lag_idx] = xbase;
                    d_shape[1][lag_idx] = ybase;
                    d_shape[2][lag_idx] = k * d_mesh_width[2];
                }

                // Fill the rest of the cross section next.
                for (int j = 1; j <= NumPtsInSection; ++j)
                {
                    const double y = j * d_mesh_width[1];
                    for (int k = -NumPtsInHeight; k <= NumPtsInHeight; ++k)
                    {
                        const double z = k * d_mesh_width[2];
                        if ((std::pow(y / width, 2) + std::pow(z / depth, 2)) <= 1) // use elliptical cross sections
                        {
                            d_shape[0][++lag_idx] = xbase; // right side.
                            d_shape[1][lag_idx] = ybase + y;
                            d_shape[2][lag_idx] = z;

                            d_shape[0][++lag_idx] = xbase; // left side.
                            d_shape[1][lag_idx] = ybase - y;
                            d_shape[2][lag_idx] = z;
                        }
                    }
                } // cross section filled
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

        // Now rotate the shape about the origin (or the shifted center of mass).
        const double angleFromHorizontal = d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];
        for (int i = 0; i < total_lag_pts; ++i)
        {
            const double x_rotated =
                d_shape[0][i] * cos(angleFromHorizontal) - d_shape[1][i] * sin(angleFromHorizontal);
            const double y_rotated =
                d_shape[0][i] * sin(angleFromHorizontal) + d_shape[1][i] * cos(angleFromHorizontal);
            d_shape[0][i] = x_rotated;
            d_shape[1][i] = y_rotated;
        }

        return;
    }

    d_current_time = d_new_time;

} // setShape

const std::vector<std::vector<double> >&
IBEELKinematics3d::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
