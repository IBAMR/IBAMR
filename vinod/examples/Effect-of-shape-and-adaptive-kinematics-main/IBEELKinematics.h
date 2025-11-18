// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBEELKinematics
#define included_IBEELKinematics

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/ConstraintIBKinematics.h>

#include <ibtk/LDataManager.h>
#include <ibtk/ibtk_utilities.h>

#include <CartesianGridGeometry.h>
#include <PatchHierarchy.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace mu
{
class Parser;
}

namespace IBAMR
{
/*!
 * \brief Class IBEELKinematics provides support for kinematics of an eel-like undulatory foil
 * with Reynolds number dependent adaptive kinematics and variable thickness.
 *
 * This class implements:
 * - Reynolds number dependent swimming kinematics
 * - Variable foil thickness (shape parameter)
 * - Adaptive amplitude and frequency modulation
 * - Anguilliform and carangiform swimming modes
 */
class IBEELKinematics : public ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    IBEELKinematics(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    IBTK::LDataManager* l_data_manager,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                    bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBEELKinematics();

    /*!
     * \brief Set kinematics velocity.
     */
    virtual void setKinematicsVelocity(const double time,
                                      const std::vector<double>& incremented_angle_from_reference_axis,
                                      const std::vector<double>& center_of_mass,
                                      const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Get kinematics velocity.
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const;

    /*!
     * \brief Set shape of the body.
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * \brief Get shape of the body.
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const;

    /*!
     * \brief Write state to restart database.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Copy constructor (not implemented).
     */
    IBEELKinematics(const IBEELKinematics& from);

    /*!
     * \brief Assignment operator (not implemented).
     */
    IBEELKinematics& operator=(const IBEELKinematics& that);

    /*!
     * \brief Read object state from restart file.
     */
    void getFromRestart();

    /*!
     * \brief Set immersed body layout.
     */
    void setImmersedBodyLayout(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set eel specific velocity.
     */
    void setEelSpecificVelocity(const double time,
                               const std::vector<double>& incremented_angle_from_reference_axis,
                               const std::vector<double>& center_of_mass,
                               const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Transform maneuver axis and calculate tangents.
     */
    void transformManeuverAxisAndCalculateTangents(const double angleFromHorizontal);

    /*!
     * \brief Calculate adaptive kinematics parameters based on Reynolds number.
     */
    void calculateAdaptiveKinematics(const double time);

    /*!
     * \brief Write performance metrics to file.
     */
    void writePerformanceMetrics(const double time);

    /*!
     * Current time (t) and new time (t+dt).
     */
    double d_current_time, d_new_time;

    /*!
     * Deformation velocity and shape of the body.
     */
    std::vector<std::vector<double> > d_kinematics_vel;
    std::vector<std::vector<double> > d_shape;

    /*!
     * Center of mass, tagged point position, and incremented rotation angle.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;

    /*!
     * Background mesh data.
     */
    SAMRAI::tbox::Array<double> d_mesh_width;

    /*!
     * Parser variables.
     */
    std::array<double, NDIM> d_parser_posn;
    std::array<double, NDIM> d_parser_normal;
    double d_parser_time;

    /*!
     * Parsers for deformation velocity, body shape, and maneuvering axis.
     */
    std::vector<mu::Parser*> d_deformationvel_parsers;
    std::vector<mu::Parser*> d_all_parsers;
    mu::Parser* d_body_shape_parser;
    mu::Parser* d_maneuvering_axis_parser;

    /*!
     * Body kinematics flags.
     */
    bool d_bodyIsManeuvering;
    bool d_maneuverAxisIsChangingShape;
    double d_initAngle_bodyAxis_x;

    /*!
     * Immersed body data: map from s-coordinate to number of points.
     */
    std::map<double, int> d_ImmersedBodyData;

    /*!
     * Maneuvering axis coordinates and tangent data.
     */
    std::vector<std::vector<double> > d_maneuverAxisReferenceCoordinates_vec;
    std::vector<std::vector<double> > d_maneuverAxisTransformedCoordinates_vec;
    std::map<double, double> d_map_reference_tangent;
    std::map<double, double> d_map_transformed_tangent;
    std::map<double, std::vector<int> > d_map_reference_sign;
    std::map<double, std::vector<int> > d_map_transformed_sign;

    /*!
     * Food location (for adaptive maneuvering).
     */
    SAMRAI::tbox::Array<double> d_food_location;

    /*!
     * Reynolds number and thickness parameters for adaptive kinematics.
     */
    double d_reynolds_number;
    double d_thickness_ratio;           // h/L ratio
    double d_base_amplitude;            // Base amplitude
    double d_base_frequency;            // Base frequency
    double d_swimming_mode;             // 0 = anguilliform, 1 = carangiform, intermediate values = mixed

    /*!
     * Adaptive kinematics parameters (computed based on Re and thickness).
     */
    double d_adapted_amplitude;
    double d_adapted_frequency;
    double d_adapted_wavelength;
    double d_envelope_power;           // Power for amplitude envelope

    /*!
     * Performance metrics tracking.
     */
    bool d_track_performance;
    std::string d_performance_log_file;
    double d_instantaneous_thrust;
    double d_instantaneous_power;
    double d_swimming_speed;

    /*!
     * Shape adaptation parameters.
     */
    bool d_enable_shape_adaptation;
    double d_head_width_ratio;         // Ratio of head width to body length
    double d_tail_width_ratio;         // Ratio of tail width (for carangiform)

}; // IBEELKinematics

} // namespace IBAMR

#endif // #ifndef included_IBEELKinematics
