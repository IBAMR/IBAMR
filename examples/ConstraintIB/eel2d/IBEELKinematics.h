// Filename : IBEELKinematics.h
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

#ifndef included_IBEELKinematics
#define included_IBEELKinematics

/////////////////////////////////////// INCLUDES ////////////////////////////////
#include <iostream>
#include <vector>
#include <map>

#include "ibamr/ConstraintIBKinematics.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "PatchHierarchy.h"

namespace mu
{
class Parser;
} // namespace mu

///////////////////////////////////////////////////////////////// CLASS DEFORMATIONAL KINEMATICS //////////////////

namespace IBAMR
{
/*!
 * \brief IBEELKinematics is a concrete class which calculates the deformation velocity and updated shape
 * for 2D eel. It also provides routines for maneuvering and food tracking cases. Example taken from:
 *
 *  Bhalla et al. A unified mathematical framework and an adaptive numerical method for
 *  fluid-structure interaction with rigid, deforming, and elastic bodies. J Comput Phys, 250:446-476 (2013).
 */

class IBEELKinematics : public ConstraintIBKinematics

{
public:
    /*!
     * \brief ctor. This is the only ctor for this object.
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
     * \brief Set kinematics velocity for eel.
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
     * \brief Set the shape of eel at the required time.
     * \see IBAMR::ConstraintIBKinematics::setShape
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * \brief Get the shape of eel at the required level.
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
    IBEELKinematics();

    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    IBEELKinematics(const IBEELKinematics& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    IBEELKinematics& operator=(const IBEELKinematics& that);

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief set eel body shape related data.
     */
    void setImmersedBodyLayout(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set deformation kinematics velocity of the eel.
     */
    void setEelSpecificVelocity(const double time,
                                const std::vector<double>& incremented_angle_from_reference_axis,
                                const std::vector<double>& center_of_mass,
                                const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Rotate the maneuver axis and caluclate tangents in this orientation.
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
     * The following map is used to store eel body shape specific data.
     * The arc length 's' varies from 0 - 1. In the std::map  the arc length 's' is used as a key.
     * d_ImmersedBodyData is used to store the no of material points which represents a cross section. The
     * width of cross section of eel varies with arc length.
     */
    std::map<double, int> d_ImmersedBodyData;

    /*!
     * Initial orientation of the body axis
     */
    double d_initAngle_bodyAxis_x;

    /*!
     * Boolean value indicating if eel is maneuvering or not.
     * If eel is maneuvering then the traveling wave will be along a curved axis, otherwise it will be on a straight
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
    double* d_parser_time;
    double* d_parser_posn;
    double* d_parser_normal;

    /*!
     * Array containing initial coordinates of the food location.
     */
    SAMRAI::tbox::Array<double> d_food_location;

}; // IBEELKinematics

} // IBAMR
#endif //#ifndef included_IBEELKinematics
