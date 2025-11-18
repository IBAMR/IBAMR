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

/////////////////////////////// OVERVIEW //////////////////////////////////////
//
// ==================================================================================
// NACA0012 Foil Kinematics for Swimmer Modeling
// ==================================================================================
//
// GEOMETRIC MODELING APPROACH:
// ----------------------------
// We employ a NACA0012 airfoil profile to model the swimmers' bodies, where the
// chord represents the spine of a swimmer at the time of their static equilibrium.
// The symmetric thickness distribution of NACA0012 (12% thickness-to-chord ratio)
// provides a realistic representation of streamlined aquatic organisms' body
// cross-sections. During swimming, the chord (spine) undergoes undulatory
// deformations while maintaining the NACA0012 thickness profile perpendicular
// to the local spine position.
//
// ==================================================================================
// TWO TYPES OF WAVY KINEMATIC MODES:
// ==================================================================================
//
// 1. ANGUILLIFORM MODE (Eel-like Swimming):
//    -----------------------------------------------------------------------
//    Biological Examples: Eels (Anguilla), lampreys, sea snakes
//
//    Physical Characteristics:
//    - Wave amplitude increases gradually and continuously from head to tail
//    - Large amplitude undulations propagate along the entire body length
//    - Entire body participates actively in thrust generation
//    - Wavelength typically equals or exceeds body length (λ ≥ L)
//    - More than one complete wavelength may be present along the body
//
//    Mathematical Model - Amplitude Envelope:
//    A(x/L) = c₀ + c₁*(x/L) + c₂*(x/L)²
//    Example coefficients: c₀ = 0.0367, c₁ = 0.0323, c₂ = 0.0310
//
//    Interpretation:
//    - c₀: Baseline amplitude at head (small but non-zero for gradual onset)
//    - c₁: Linear growth rate along body
//    - c₂: Quadratic acceleration toward tail region
//    - All coefficients positive → monotonically increasing amplitude
//
//    Swimming Performance:
//    - High maneuverability in confined spaces
//    - Moderate swimming efficiency
//    - Effective for navigating complex environments
//    - Lower maximum speed compared to carangiform
//
// 2. CARANGIFORM MODE (Fish-like Swimming):
//    -----------------------------------------------------------------------
//    Biological Examples: Tunas (Thunnus), mackerels, jacks, sailfish
//
//    Physical Characteristics:
//    - Wave amplitude concentrated primarily in posterior half of body
//    - Anterior body (head and thorax) remains relatively rigid
//    - Posterior body and caudal peduncle provide majority of thrust
//    - Wavelength approximately equals body length (λ ≈ L)
//    - Typically less than one complete wavelength on the body
//
//    Mathematical Model - Amplitude Envelope (Khalid et al. 2016):
//    A(x/L) = 0.02 - 0.0825*(x/L) + 0.1625*(x/L)²
//
//    Interpretation:
//    - At x/L = 0 (head): A = 0.02 (very small amplitude)
//    - Minimum amplitude occurs at x/L ≈ 0.25 (rigid anterior region)
//    - Amplitude increases rapidly for x/L > 0.5 (flexible posterior)
//    - At x/L = 1.0 (tail): A ≈ 0.10 (maximum amplitude, 5× head amplitude)
//    - Negative linear term + dominant positive quadratic → characteristic
//      posterior-concentrated motion profile
//
//    Swimming Performance:
//    - High cruising speed and efficiency
//    - Streamlined for sustained fast swimming
//    - Reduced drag from rigid anterior body
//    - Optimized for open-water locomotion
//
// ==================================================================================
// MATHEMATICAL FORMULATION OF UNDULATORY SWIMMING:
// ==================================================================================
//
//    Centerline Displacement:  y(x,t) = A(x/L) * cos[2π(x/λ - ft)]
//
//    Deformation Velocity:     ∂y/∂t = A(x/L) * 2πf * sin[2π(x/λ - ft)]
//
// Parameters:
//    x      = chordwise position along spine [0, L]
//    L      = body length (chord length)
//    λ      = wavelength of undulation
//    f      = swimming frequency (Hz)
//    t      = time (s)
//    A(x/L) = amplitude envelope function (mode-dependent)
//
// Wave Propagation:
//    Phase velocity: c_p = λf (wave travels posteriorly)
//    Wave number:    k = 2π/λ
//    Angular freq:   ω = 2πf
//
// ==================================================================================
// IMPLEMENTATION AND USAGE:
// ==================================================================================
// The swimming mode is selected by choosing the appropriate amplitude envelope
// function A(x/L) in the input2d configuration file. Users specify:
//   - body_shape_equation: Centerline shape y(x,t)
//   - deformation_velocity_function_0/1: Components of ∂y/∂t
//
// The implementation supports:
//   - Constraint-based IB method (ConstraintIBKinematics)
//   - Optional maneuvering along curved paths (food tracking, turning)
//   - Adaptive mesh refinement (AMR) with IBAMR
//   - Multiple swimming frequencies and wavelengths
//
// ==================================================================================
// REFERENCES:
// ==================================================================================
// - Khalid, M.S.U., et al. "A bio-inspired study on tuna-mimetic soft robot with a
//   compliant caudal fin." J. Fluids Structures, 66:19-35 (2016).
// - Lighthill, M.J. "Note on the swimming of slender fish." J. Fluid Mech., 9(2):305 (1960).
// - Sfakiotakis, M., et al. "Review of fish swimming modes for aquatic locomotion."
//   IEEE J. Oceanic Eng., 24(2):237-252 (1999).
//
//////////////////////////// INCLUDES /////////////////////////////////////////
#include "ibtk/IBTK_MPI.h"

#include "CartesianPatchGeometry.h"
#include "IBNACA0012Kinematics.h"
#include "PatchLevel.h"
#include "tbox/MathUtilities.h"

#include "muParser.h"

#include <cmath>
#include <fstream>
#include <iostream>

#include "ibamr/namespaces.h"

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

// NACA0012 Foil Parameters for Swimmer Modeling
// The chord represents the spine of the swimmer at static equilibrium
static const double CHORD_LENGTH = 1.0;              // Chord length (body length)
static const double MAX_THICKNESS = 0.12;            // NACA0012: 12% thickness-to-chord ratio

// NACA 4-digit airfoil thickness distribution coefficients
// y_t/c = (t/c) * [a0*sqrt(x/c) + a1*(x/c) + a2*(x/c)^2 + a3*(x/c)^3 + a4*(x/c)^4]
static const double NACA_A0 =  0.2969;
static const double NACA_A1 = -0.1260;
static const double NACA_A2 = -0.3516;
static const double NACA_A3 =  0.2843;
static const double NACA_A4 = -0.1015;  // Closed trailing edge (original: -0.1036 for open TE)

// Maneuvering parameters (for food tracking/turning behavior)
static const double CUT_OFF_ANGLE = PII / 4;
static const double CUT_OFF_RADIUS = 0.7;
static const double LOWER_CUT_OFF_ANGLE = 7 * PII / 180;

/*!
 * \brief Calculate NACA0012 half-thickness at given chordwise position
 *
 * For NACA 4-digit series: y_t = (t/0.2) * c * [a0*sqrt(x/c) + a1*(x/c) + ... + a4*(x/c)^4]
 * where t = maximum thickness ratio (0.12 for NACA0012)
 *
 * \param x_c Chordwise position normalized by chord (x/c), range [0, 1]
 * \return Half-thickness at position x_c
 */
inline double
naca0012_thickness(const double x_c)
{
    if (x_c < 0.0 || x_c > 1.0) return 0.0;

    const double sqrt_x_c = std::sqrt(x_c);
    const double x_c2 = x_c * x_c;
    const double x_c3 = x_c2 * x_c;
    const double x_c4 = x_c3 * x_c;

    // NACA 4-digit thickness distribution
    const double y_t = (MAX_THICKNESS / 0.2) * (
        NACA_A0 * sqrt_x_c +
        NACA_A1 * x_c +
        NACA_A2 * x_c2 +
        NACA_A3 * x_c3 +
        NACA_A4 * x_c4
    );

    return y_t * CHORD_LENGTH;  // Return dimensional thickness
}

} // namespace

///////////////////////////////////////////////////////////////////////

IBNACA0012Kinematics::IBNACA0012Kinematics(const std::string& object_name,
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
      d_parser_time(0.0)
{
    // Read from inputdb
    d_initAngle_bodyAxis_x = input_db->getDoubleWithDefault("initial_angle_body_axis_0", 0.0);
    d_bodyIsManeuvering = input_db->getBoolWithDefault("body_is_maneuvering", false);
    d_maneuverAxisIsChangingShape = input_db->getBoolWithDefault("maneuvering_axis_is_changing_shape", false);

    // Read-in deformation velocity functions
    std::vector<std::string> deformationvel_function_strings;
    for (int d = 0; d < NDIM; ++d)
    {
        const std::string postfix = "_function_" + std::to_string(d);
        std::string key_name = "deformation_velocity" + postfix;

        if (input_db->isString(key_name))
        {
            deformationvel_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            deformationvel_function_strings.push_back("0.0");
            TBOX_WARNING("IBNACA0012Kinematics::IBNACA0012Kinematics() :\n"
                         << "  no function corresponding to key ``" << key_name << " '' found for dimension = " << d
                         << "; using def_vel = 0.0. " << std::endl);
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
        (*cit)->DefineVar("T", &d_parser_time);
        (*cit)->DefineVar("t", &d_parser_time);
        for (int d = 0; d < NDIM; ++d)
        {
            const std::string postfix = std::to_string(d);
            (*cit)->DefineVar("X" + postfix, d_parser_posn.data() + d);
            (*cit)->DefineVar("x" + postfix, d_parser_posn.data() + d);
            (*cit)->DefineVar("X_" + postfix, d_parser_posn.data() + d);
            (*cit)->DefineVar("x_" + postfix, d_parser_posn.data() + d);

            (*cit)->DefineVar("N" + postfix, d_parser_normal.data() + d);
            (*cit)->DefineVar("n" + postfix, d_parser_normal.data() + d);
            (*cit)->DefineVar("N_" + postfix, d_parser_normal.data() + d);
            (*cit)->DefineVar("n_" + postfix, d_parser_normal.data() + d);
        }
    }

    // Set the location of the food/target particle (only needed if maneuvering with changing axis)
    if (d_bodyIsManeuvering && d_maneuverAxisIsChangingShape)
    {
        d_food_location.resizeArray(NDIM);
        for (int dim = 0; dim < NDIM; ++dim)
        {
            if (input_db->keyExists("food_location_in_domain_" + std::to_string(dim)))
            {
                d_food_location[dim] = input_db->getDouble("food_location_in_domain_" + std::to_string(dim));
            }
            else
            {
                TBOX_WARNING("IBNACA0012Kinematics::IBNACA0012Kinematics() :\n"
                             << "  body_is_maneuvering and maneuvering_axis_is_changing_shape are TRUE,\n"
                             << "  but food_location_in_domain_" << dim << " not found. Using default 0.0." << std::endl);
                d_food_location[dim] = 0.0;
            }
        }
    }
    else
    {
        // Not needed for non-maneuvering cases
        d_food_location.resizeArray(NDIM);
        for (int dim = 0; dim < NDIM; ++dim)
        {
            d_food_location[dim] = 0.0;
        }
    }

    // set how the immersed body is layout in reference frame.
    setImmersedBodyLayout(patch_hierarchy);

    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();

    return;

} // IBNACA0012Kinematics

IBNACA0012Kinematics::~IBNACA0012Kinematics()
{
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {
        delete (*cit);
    }
    return;

} // ~IBNACA0012Kinematics

void
IBNACA0012Kinematics::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase

void
IBNACA0012Kinematics::getFromRestart()
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
    db->getDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->getDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->getDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;
} // getFromRestart

void
IBNACA0012Kinematics::setImmersedBodyLayout(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
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

    // Number of points along the chord (backbone)
    // For NACA0012: discretize the chord from leading edge (x=0) to trailing edge (x=CHORD_LENGTH)
    const int ChordNx = static_cast<int>(ceil(CHORD_LENGTH / d_mesh_width[0]));

    d_ImmersedBodyData.clear();

    // Build the NACA0012 cross-sectional data using the thickness distribution
    // At each chordwise station, we need to discretize the thickness
    for (int i = 0; i < ChordNx; ++i)
    {
        const double x = i * d_mesh_width[0];              // Dimensional chordwise position
        const double x_c = x / CHORD_LENGTH;                // Normalized chordwise position (x/c)

        // Get NACA0012 half-thickness at this station
        const double half_thickness = naca0012_thickness(x_c);

        // Number of points needed to discretize this cross-section
        // We need points on both upper and lower surfaces
        int NumPtsInSection;
        if (half_thickness > 1.0e-10)  // Avoid division by zero at leading/trailing edges
        {
            NumPtsInSection = 2 * static_cast<int>(ceil(half_thickness / d_mesh_width[1]));
            // Ensure at least 2 points (one on each surface)
            if (NumPtsInSection < 2) NumPtsInSection = 2;
        }
        else
        {
            NumPtsInSection = 2;  // Minimal points at sharp edges
        }

        d_ImmersedBodyData.insert(std::make_pair(x, NumPtsInSection));
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
IBNACA0012Kinematics::transformManeuverAxisAndCalculateTangents(const double angleFromHorizontal)
{
    d_maneuverAxisTransformedCoordinates_vec.clear();
    d_map_transformed_tangent.clear();
    d_map_transformed_sign.clear();

    const int ChordNx = static_cast<int>(ceil(CHORD_LENGTH / d_mesh_width[0]));
    std::vector<double> transformed_coord(2);
    for (int i = 0; i < ChordNx; ++i)
    {
        transformed_coord[0] = d_maneuverAxisReferenceCoordinates_vec[i][0] * cos(angleFromHorizontal) -
                               d_maneuverAxisReferenceCoordinates_vec[i][1] * sin(angleFromHorizontal);
        transformed_coord[1] = d_maneuverAxisReferenceCoordinates_vec[i][0] * sin(angleFromHorizontal) +
                               d_maneuverAxisReferenceCoordinates_vec[i][1] * cos(angleFromHorizontal);
        d_maneuverAxisTransformedCoordinates_vec.push_back(transformed_coord);
    }

    for (int i = 0; i < (ChordNx - 1); ++i)
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
    if (ChordNx > 0)
    {
        d_map_transformed_tangent.insert(
            std::make_pair((ChordNx - 1) * d_mesh_width[0], (d_map_transformed_tangent.rbegin())->second));
        d_map_transformed_sign.insert(
            std::make_pair((ChordNx - 1) * d_mesh_width[0], (d_map_transformed_sign.rbegin())->second));
    }

    return;

} // transformManeuverAxisAndCalculateTangents

void
IBNACA0012Kinematics::setNACA0012SpecificVelocity(const double time,
                                                   const std::vector<double>& incremented_angle_from_reference_axis,
                                                   const std::vector<double>& center_of_mass,
                                                   const std::vector<double>& tagged_pt_position)
{
    d_parser_time = time;
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
            else if (IBTK::abs_equal_eps(MathUtilities<double>::Abs(angle_bw_target_vision), 0.0))
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
            const int ChordNx = static_cast<int>(ceil(CHORD_LENGTH / d_mesh_width[0]));
            if (radius_circular_path != __INFINITY)
            {
                const double angle_sector = CHORD_LENGTH / radius_circular_path;
                const double dtheta = angle_sector / (ChordNx - 1);

                d_maneuverAxisReferenceCoordinates_vec.clear();
                std::vector<double> vec_axis_coord(2);
                for (int i = 0; i < ChordNx; ++i)
                {
                    const double angleFromVertical = -angle_sector / 2 + i * dtheta;
                    vec_axis_coord[0] = radius_circular_path * sin(angleFromVertical);
                    vec_axis_coord[1] = radius_circular_path * cos(angleFromVertical);
                    d_maneuverAxisReferenceCoordinates_vec.push_back(vec_axis_coord);
                }
            }
            else
            {
                d_maneuverAxisReferenceCoordinates_vec.clear();
                std::vector<double> vec_axis_coord(2);
                for (int i = 0; i < ChordNx; ++i)
                {
                    vec_axis_coord[0] = i * d_mesh_width[0];
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
            for (int i = 0; i < (ChordNx - 1); ++i)
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
            if (ChordNx > 0)
            {
                d_map_reference_tangent.insert(
                    std::make_pair((ChordNx - 1) * d_mesh_width[0], (d_map_reference_tangent.rbegin())->second));
                d_map_reference_sign.insert(
                    std::make_pair((ChordNx - 1) * d_mesh_width[0], (d_map_reference_sign.rbegin())->second));
            }
        } // maneuverAxisIsChangingShape

        // Rotate the reference axis and calculate tangents in the rotated frame.
        transformManeuverAxisAndCalculateTangents(angleFromHorizontal);
    } // bodyIsManeuvering

    // Set the deformation velocity in the body frame.
    // The velocity is calculated using mu::Parser with equations from input2d file:
    // deformation_velocity_function_0 and deformation_velocity_function_1
    // For NACA0012 carangiform: velocity = A(x/L) * 2πf * sin[2π(x/L - ft)] * Normal
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

        // Evaluate the deformation velocity using mu::Parser
        // The equations from input2d will be evaluated here
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
} // setNACA0012SpecificVelocity

void
IBNACA0012Kinematics::setKinematicsVelocity(const double time,
                                            const std::vector<double>& incremented_angle_from_reference_axis,
                                            const std::vector<double>& center_of_mass,
                                            const std::vector<double>& tagged_pt_position)
{
    d_new_time = time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;

    setNACA0012SpecificVelocity(d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

    return;

} // setKinematicsVelocity

const std::vector<std::vector<double> >&
IBNACA0012Kinematics::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getKinematicsVelocity

void
IBNACA0012Kinematics::setShape(const double time, const std::vector<double>& /*incremented_angle_from_reference_axis*/)
{
    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    if (position_update_method == "CONSTRAINT_VELOCITY") return;

    // Find the deformed shape. Rotate the shape about center of mass.
    TBOX_ASSERT(d_new_time == time);
    d_parser_time = time;
    std::vector<double> shape_new(NDIM);

    int lag_idx = -1;
    int reference_axis_idx = -1;
    for (std::map<double, int>::const_iterator itr = d_ImmersedBodyData.begin(); itr != d_ImmersedBodyData.end(); itr++)
    {
        const int NumPtsInSection = itr->second;
        d_parser_posn[0] = itr->first;
        
        // Evaluate the body shape at this streamwise position using mu::Parser
        // For NACA0012: shape = A(x/L) * cos[2π(x/L - ft)]
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
IBNACA0012Kinematics::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR