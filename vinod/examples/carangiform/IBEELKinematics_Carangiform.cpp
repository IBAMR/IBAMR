// IBEELKinematics_Carangiform.cpp
// Implementation of CARANGIFORM (tuna-like) swimming kinematics
// Amplitude envelope: A(X) = a₀ + a₁X + a₂X²

#include "IBEELKinematics.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBEELKinematics::IBEELKinematics(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, register_for_restart)
{
    // Read swimming parameters from input database
    d_frequency = input_db->getDoubleWithDefault("frequency", 3.0);  // Hz
    d_wavelength = input_db->getDoubleWithDefault("wavelength", 1.5);  // fraction of body length

    // Carangiform quadratic amplitude envelope coefficients
    d_a0 = input_db->getDoubleWithDefault("a0", 0.01);   // Constant term
    d_a1 = input_db->getDoubleWithDefault("a1", -0.05);  // Linear term (negative!)
    d_a2 = input_db->getDoubleWithDefault("a2", 0.14);   // Quadratic term

    // Fish geometry
    d_length = input_db->getDoubleWithDefault("length", 1.0);  // Body length

    // Calculate derived parameters
    d_omega = 2.0 * M_PI * d_frequency;  // Angular frequency
    d_wavenumber = 2.0 * M_PI / (d_wavelength * d_length);  // Wave number k = 2π/λ

    // Initial center of mass
    d_initCenterOfMass.resize(NDIM);
    input_db->getDoubleArray("initCenterOfMass", d_initCenterOfMass.data(), NDIM);

    return;
} // IBEELKinematics

IBEELKinematics::~IBEELKinematics()
{
    // intentionally blank
    return;
} // ~IBEELKinematics

void
IBEELKinematics::setKinematicsVelocity(
    const double time,
    const std::vector<double>& incremented_angle_from_reference_axis,
    const std::vector<double>& center_of_mass,
    const std::vector<double>& tagged_pt_position)
{
    d_current_time = time;
    d_center_of_mass = center_of_mass;

    return;
} // setKinematicsVelocity

void
IBEELKinematics::getVelocity(
    const Eigen::Vector3d& point,
    Eigen::Vector3d& velocity)
{
    // Calculate position along body (x-coordinate relative to center)
    double x = point[0] - d_initCenterOfMass[0];

    // Normalize position: X = x/L, where X ∈ [0, 1]
    double X = (x / d_length) + 0.5;  // Assuming fish centered at origin, X=0 is head, X=1 is tail

    // Clamp to [0, 1]
    if (X < 0.0) X = 0.0;
    if (X > 1.0) X = 1.0;

    // Calculate CARANGIFORM amplitude envelope: A(X) = a₀ + a₁X + a₂X²
    // This creates a quadratic increase with most amplitude at the tail
    // The negative a₁ term creates minimal amplitude in anterior body
    double amplitude = d_a0 + d_a1 * X + d_a2 * X * X;

    // Ensure amplitude is non-negative (important near head where a₁ is negative)
    if (amplitude < 0.0) amplitude = 0.0;

    // Calculate phase: φ = kx - ωt
    double phase = d_wavenumber * x - d_omega * d_current_time;

    // Calculate lateral velocity: v_y = -A(x) * ω * cos(kx - ωt)
    // The negative sign comes from ∂/∂t[sin(kx - ωt)] = -ω*cos(kx - ωt)
    double velocity_y = -amplitude * d_omega * cos(phase);

    // Set velocity components
    velocity[0] = 0.0;  // No forward swimming velocity (let IBAMR compute from reaction)
    velocity[1] = velocity_y;  // Lateral (vertical) velocity

    if (NDIM == 3)
    {
        velocity[2] = 0.0;  // No z-velocity in 2D simulation
    }

    return;
} // getVelocity

void
IBEELKinematics::setShape(
    const double time,
    const std::vector<double>& /*incremented_angle_from_reference_axis*/)
{
    // Update time for velocity calculations
    d_current_time = time;

    return;
} // setShape

/////////////////////////////// PRIVATE //////////////////////////////////////

} // namespace IBAMR
