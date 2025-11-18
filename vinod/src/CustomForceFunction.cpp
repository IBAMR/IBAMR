// Filename: CustomForceFunction.cpp
// Created by: Vinod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "CustomForceFunction.h"

#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>

#include <CartesianPatchGeometry.h>
#include <PatchLevel.h>
#include <tbox/Database.h>
#include <tbox/Utilities.h>

#include <cmath>

/////////////////////////////// PUBLIC ///////////////////////////////////////

CustomForceFunction::CustomForceFunction(const std::string& object_name,
                                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : d_object_name(object_name)
{
    // Read parameters from input database if provided
    if (input_db)
    {
        if (input_db->keyExists("force_magnitude"))
        {
            d_force_magnitude = input_db->getDouble("force_magnitude");
        }

        if (input_db->keyExists("force_direction"))
        {
            input_db->getDoubleArray("force_direction", d_force_direction, NDIM);
        }

        if (input_db->keyExists("start_time"))
        {
            d_start_time = input_db->getDouble("start_time");
        }

        if (input_db->keyExists("ramp_time"))
        {
            d_ramp_time = input_db->getDouble("ramp_time");
        }
    }

    return;
} // CustomForceFunction

CustomForceFunction::~CustomForceFunction()
{
    // Intentionally blank
    return;
} // ~CustomForceFunction

void
CustomForceFunction::setForceMagnitude(double force_mag)
{
    d_force_magnitude = force_mag;
    return;
} // setForceMagnitude

void
CustomForceFunction::setForceDirection(double fx, double fy)
{
    d_force_direction[0] = fx;
    d_force_direction[1] = fy;
#if (NDIM == 3)
    d_force_direction[2] = 0.0;
#endif
    return;
} // setForceDirection

void
CustomForceFunction::setLagrangianForce(SAMRAI::tbox::Pointer<IBTK::LData> F_data,
                                       SAMRAI::tbox::Pointer<IBTK::LData> X_data,
                                       SAMRAI::tbox::Pointer<IBTK::LData> U_data,
                                       const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > /*hierarchy*/,
                                       const int /*level_number*/,
                                       const double data_time,
                                       IBTK::LDataManager* const /*lag_manager*/)
{
    // Compute time-dependent ramping factor
    double time_factor = 1.0;
    if (data_time < d_start_time)
    {
        time_factor = 0.0;
    }
    else if (d_ramp_time > 0.0 && data_time < d_start_time + d_ramp_time)
    {
        // Smooth ramp using cosine function
        const double t = (data_time - d_start_time) / d_ramp_time;
        time_factor = 0.5 * (1.0 - std::cos(M_PI * t));
    }

    // Get force vector
    double* const F_vals = F_data->getLocalFormVecArray()->data();
    const int F_depth = F_data->getDepth();

    // Get position vector (if needed for spatially-varying forces)
    const double* const X_vals = X_data->getLocalFormVecArray()->data();
    const int X_depth = X_data->getDepth();

    // Get number of local nodes
    const int num_local_nodes = F_data->getLocalNodeCount();

    // Apply force to all Lagrangian points
    for (int k = 0; k < num_local_nodes; ++k)
    {
        // Current position
        const double* const X = &X_vals[k * X_depth];

        // Apply force (can be made position-dependent here)
        double* const F = &F_vals[k * F_depth];
        for (int d = 0; d < NDIM; ++d)
        {
            F[d] += time_factor * d_force_magnitude * d_force_direction[d];
        }

        // Example: Position-dependent force (uncomment to use)
        // if (X[0] > 0.5) // Apply force only to points with x > 0.5
        // {
        //     F[0] += time_factor * d_force_magnitude * d_force_direction[0];
        //     F[1] += time_factor * d_force_magnitude * d_force_direction[1];
        // }
    }

    // Restore arrays
    F_data->restoreArrays();
    X_data->restoreArrays();
    if (U_data) U_data->restoreArrays();

    return;
} // setLagrangianForce

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
