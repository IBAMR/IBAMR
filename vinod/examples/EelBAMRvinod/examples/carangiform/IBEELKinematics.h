// IBEELKinematics.h
// Header file for CARANGIFORM swimming kinematics

#ifndef included_IBEELKinematics
#define included_IBEELKinematics

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/ConstraintIBKinematics.h>
#include <ibtk/LDataManager.h>

namespace IBTK
{
class LDataManager;
} // namespace IBTK

/*!
 * \brief Class IBEELKinematics provides kinematic velocity for carangiform
 * (tuna-like) swimming.
 *
 * The swimming kinematics follows the equation:
 *   y(x,t) = A(x) * sin(kx - ωt)
 *
 * where the CARANGIFORM amplitude envelope is:
 *   A(x) = a₀ + a₁(x/L) + a₂(x/L)²
 *
 * This creates a quadratic amplitude distribution that concentrates
 * motion at the posterior body and tail, characteristic of tuna-like
 * swimming with high efficiency.
 */
class IBEELKinematics : public IBAMR::ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    IBEELKinematics(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBEELKinematics();

    /*!
     * \brief Set kinematics velocity at specified time.
     */
    virtual void setKinematicsVelocity(
        const double time,
        const std::vector<double>& incremented_angle_from_reference_axis,
        const std::vector<double>& center_of_mass,
        const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Get velocity at a point.
     */
    virtual void getVelocity(
        const Eigen::Vector3d& point,
        Eigen::Vector3d& velocity);

    /*!
     * \brief Set shape at specified time.
     */
    virtual void setShape(
        const double time,
        const std::vector<double>& incremented_angle_from_reference_axis);

private:
    /*!
     * \brief Swimming parameters
     */
    double d_frequency;      // Swimming frequency (Hz)
    double d_wavelength;     // Wavelength (fraction of body length)
    double d_a0;             // Quadratic amplitude coefficient (constant)
    double d_a1;             // Quadratic amplitude coefficient (linear)
    double d_a2;             // Quadratic amplitude coefficient (quadratic)
    double d_length;         // Fish body length

    /*!
     * \brief Derived parameters
     */
    double d_omega;          // Angular frequency ω = 2πf
    double d_wavenumber;     // Wave number k = 2π/λ

    /*!
     * \brief Current state
     */
    double d_current_time;
    std::vector<double> d_center_of_mass;
    std::vector<double> d_initCenterOfMass;
};

#endif // included_IBEELKinematics
