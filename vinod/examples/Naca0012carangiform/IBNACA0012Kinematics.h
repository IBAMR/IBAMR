// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBNACA0012Kinematics
#define included_IBNACA0012Kinematics

/////////////////////////////////////// INCLUDES ////////////////////////////////
#include "ibamr/ConstraintIBKinematics.h"

#include "PatchHierarchy.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <iostream>
#include <map>
#include <vector>

namespace mu
{
class Parser;
} // namespace mu

///////////////////////////////////////////////////////////////// CLASS DEFORMATIONAL KINEMATICS //////////////////

namespace IBAMR
{
/*!
 * \brief IBNACA0012Kinematics is a concrete class which calculates the deformation velocity and updated shape
 * for a NACA0012 airfoil representing a swimmer's body.
 *
 * ==================================================================================
 * GEOMETRIC MODELING:
 * ==================================================================================
 * We employ a NACA0012 airfoil profile to model the swimmers' bodies, where the chord
 * represents the spine of a swimmer at the time of their static equilibrium. The symmetric
 * thickness distribution of NACA0012 provides a realistic representation of the body
 * cross-section, while the chord line serves as the centerline/spine that undergoes
 * undulatory deformation during swimming.
 *
 * NACA0012 Thickness Distribution:
 *    y_t(x) = (t/0.2) * [0.2969√(x/c) - 0.1260(x/c) - 0.3516(x/c)² + 0.2843(x/c)³ - 0.1015(x/c)⁴]
 *    where t = 0.12 (12% thickness), c = chord length = body length L
 *
 * ==================================================================================
 * SWIMMING MODES:
 * ==================================================================================
 * Two types of wavy kinematic modes are considered to model different swimming strategies:
 *
 * 1. ANGUILLIFORM (Eel-like Swimming):
 *    --------------------------------------------------------------------------
 *    Biological Examples: Eels, lampreys, sea snakes
 *
 *    Characteristics:
 *    - Wave amplitude increases gradually and continuously from head to tail
 *    - Large amplitude undulations propagate along the entire body length
 *    - Wavelength typically equals or exceeds body length (λ ≥ L)
 *    - Entire body participates in thrust generation
 *    - More than one complete wavelength may be present on the body
 *
 *    Amplitude Envelope:
 *    A(x/L) = c₀ + c₁*(x/L) + c₂*(x/L)²
 *    where coefficients increase posteriorly (e.g., c₀=0.0367, c₁=0.0323, c₂=0.0310)
 *
 *    Physical Interpretation:
 *    - c₀: baseline amplitude at head (small but non-zero)
 *    - c₁: linear growth rate
 *    - c₂: quadratic acceleration toward tail
 *
 *    Performance: High maneuverability, moderate efficiency, effective in confined spaces
 *
 * 2. CARANGIFORM (Fish-like Swimming):
 *    --------------------------------------------------------------------------
 *    Biological Examples: Tunas, mackerels, jacks, most fast-swimming teleost fish
 *
 *    Characteristics:
 *    - Wave amplitude concentrated primarily in posterior half of body
 *    - Anterior body (head region) remains relatively rigid
 *    - Posterior body and caudal fin provide majority of thrust
 *    - Wavelength approximately equals body length (λ ≈ L)
 *    - Typically less than one complete wavelength on the body
 *
 *    Amplitude Envelope (Khalid et al. 2016):
 *    A(x/L) = 0.02 - 0.0825*(x/L) + 0.1625*(x/L)²
 *
 *    Physical Interpretation:
 *    - Minimum amplitude at x/L ≈ 0.25 (anterior rigid body region)
 *    - Rapid amplitude increase in posterior half (x/L > 0.5)
 *    - Maximum amplitude at tail (x/L = 1.0)
 *    - The quadratic coefficient (0.1625) dominates, creating the characteristic
 *      posterior-concentrated motion
 *
 *    Performance: High swimming speed, high efficiency, streamlined for cruising
 *
 * ==================================================================================
 * MATHEMATICAL FORMULATION:
 * ==================================================================================
 * The undulatory swimming kinematics are modeled using a traveling wave formulation:
 *
 *    Centerline Displacement:  y(x,t) = A(x/L) * cos[2π(x/λ - ft)]
 *
 *    Deformation Velocity:     ∂y/∂t = A(x/L) * 2πf * sin[2π(x/λ - ft)]
 *
 * where:
 *    x   = chordwise position along spine [0, L]
 *    L   = body length (chord length)
 *    λ   = wavelength of undulation
 *    f   = swimming frequency (Hz)
 *    t   = time (s)
 *    A(x/L) = amplitude envelope function (mode-dependent)
 *
 * The wave travels posteriorly with phase velocity c_p = λf, and the
 * deformation velocity is normal to the local body surface.
 *
 * ==================================================================================
 * IMPLEMENTATION DETAILS:
 * ==================================================================================
 * - Swimming mode is selected by choosing appropriate amplitude envelope A(x/L) in input file
 * - Deformation velocities computed using mu::Parser for flexible function definitions
 * - Supports optional maneuvering (curved swimming paths for food tracking/turning)
 * - Compatible with IBAMR's Constraint IB method for fluid-structure interaction
 *
 * ==================================================================================
 * REFERENCES:
 * ==================================================================================
 * - Khalid, M.S.U., et al. "A bio-inspired study on tuna-mimetic soft robot with a
 *   compliant caudal fin." Journal of Fluids and Structures, 66:19-35 (2016).
 * - Lighthill, M.J. "Note on the swimming of slender fish." Journal of Fluid
 *   Mechanics, 9(2):305-317 (1960).
 * - Sfakiotakis, M., et al. "Review of fish swimming modes for aquatic locomotion."
 *   IEEE Journal of Oceanic Engineering, 24(2):237-252 (1999).
 */

class IBNACA0012Kinematics : public ConstraintIBKinematics

{
public:
    /*!
     * \brief ctor. This is the only ctor for this object.
     */
    IBNACA0012Kinematics(const std::string& object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         IBTK::LDataManager* l_data_manager,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                         bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBNACA0012Kinematics();

    /*!
     * \brief Set kinematics velocity for NACA0012 carangiform fish.
     * \see IBAMR::ConstraintIBKinematics::setKinematicsVelocity
     */
    virtual void setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Get the kinematics velocity on the specified level.
     * \see IBAMR::ConstraintIBKinematics::getKinematicsVelocity
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const;

    /*!
     * \brief Set the shape of NACA0012 fish at the required time.
     * \see IBAMR::ConstraintIBKinematics::setShape
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * \brief Get the shape of NACA0012 fish at the required level.
     * \see IBAMR::ConstraintIBKinematics::getShape
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const;

    /*!
     * \brief Override the ConstraintIBkinematics base class method.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief The default constructor is not implemented and should not be used.
     */
    IBNACA0012Kinematics();

    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    IBNACA0012Kinematics(const IBNACA0012Kinematics& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    IBNACA0012Kinematics& operator=(const IBNACA0012Kinematics& that);

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief set NACA0012 fish body shape related data.
     */
    void setImmersedBodyLayout(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set deformation kinematics velocity of the NACA0012 fish.
     */
    void setNACA0012SpecificVelocity(const double time,
                                     const std::vector<double>& incremented_angle_from_reference_axis,
                                     const std::vector<double>& center_of_mass,
                                     const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Rotate the maneuver axis and calculate tangents in this orientation.
     */
    void transformManeuverAxisAndCalculateTangents(const double angleFromHorizontal);

    /*!
     * Current time (t) and new time (t+dt).
     */
    double d_current_time, d_new_time;

    /*!
     * Deformational velocity and shape vectors.
     */
    std::vector<std::vector<double> > d_kinematics_vel;
    std::vector<std::vector<double> > d_shape;

    /*!
     * Save COM, tagged point position and incremented angle from reference axis for restarted runs.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;

    /*!
     * Eulerian Mesh width parameters.
     */
    std::vector<double> d_mesh_width;

    /*!
     * The following map is used to store NACA0012 fish body shape specific data.
     * The arc length 's' varies from 0 - 1. In the std::map the arc length 's' is used as a key.
     * d_ImmersedBodyData is used to store the no of material points which represents a cross section. The
     * width of cross section of NACA0012 varies with arc length.
     */
    std::map<double, int> d_ImmersedBodyData;

    /*!
     * Initial orientation of the body axis
     */
    double d_initAngle_bodyAxis_x;

    /*!
     * Boolean value indicating if NACA0012 fish is maneuvering or not.
     * If fish is maneuvering then the traveling wave will be along a curved axis, otherwise it will be on a straight
     * line
     */
    bool d_bodyIsManeuvering;

    /*!
     * Boolean to indicate if shape of maneuver axis is changing. The maneuver axis will change shape in food tracking
     * cases.
     */
    bool d_maneuverAxisIsChangingShape;

    /*!
     * Vector of coordinates defining the axis of maneuvering. The reference axis will rotate with body omega.
     * The coordinates of the rotated maneuver axis is stored separately.
     */
    std::vector<std::vector<double> > d_maneuverAxisReferenceCoordinates_vec, d_maneuverAxisTransformedCoordinates_vec;

    /*!
     * map of tangents along the body/maneuver axis in rotated frame. The key used is arc length 's' and it stores only
     * the abs(theta).
     * Sign of tangent is stored separately. This is done to avoid a lot of if conditions needed to determine the
     * quadrant of the
     * angle.
     */
    std::map<double, double> d_map_transformed_tangent;

    /*!
     * map of tangents along the body/maneuver axis in reference/unrotated frame.The key used is arc length 's' and it
     * stores only the abs(theta).
     * Sign of tangent is stored separately. This is done to avoid a lot of if conditions needed to determine the
     * quadrant of the
     * angle.
     */
    std::map<double, double> d_map_reference_tangent;

    /*!
     * Sign of tangent vector in rotated frame. The key used is arc length 's'. 'mapped_value' is a vector which has
     * sign of t_x and t_y
     * respectively.
     */
    std::map<double, std::vector<int> > d_map_transformed_sign;

    /*!
     * Sign of tangent vector in reference/unrotated frame. The key used is arc length 's'. 'mapped_value' is a vector
     * which has sign of t_x and t_y
     * respectively.
     */
    std::map<double, std::vector<int> > d_map_reference_sign;

    /*!
     * mu::Parser object which evaluates the maneuvering axis equation.
     */
    mu::Parser* d_maneuvering_axis_parser;

    /*!
     * mu::Parser object which evaluates the shape of the body.
     */
    mu::Parser* d_body_shape_parser;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser*> d_deformationvel_parsers;
    std::vector<mu::Parser*> d_all_parsers;

    /*!
     * Time and position variables.
     */
    mutable double d_parser_time;
    mutable IBTK::Point d_parser_posn;
    mutable IBTK::Point d_parser_normal;

    /*!
     * Array containing initial coordinates of the food location.
     */
    SAMRAI::tbox::Array<double> d_food_location;

}; // IBNACA0012Kinematics

} // namespace IBAMR
#endif // #ifndef included_IBNACA0012Kinematics