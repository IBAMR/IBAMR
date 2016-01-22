// Filename : IBEELKinematics3d.h
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

#ifndef included_IBEELKinematics3d
#define included_IBEELKinematics3d

/////////////////////////////////////// INCLUDES ///////////////////////////////////////////
#include <iostream>
#include <vector>

#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Array.h"
#include "PatchHierarchy.h"
#include "ibamr/ConstraintIBKinematics.h"

namespace mu
{
class Parser;
}

namespace IBAMR
{
/*!
 * \brief IBEELKinematics3d Class.
 *
 * IBEELKinematics3d is a concrete class which calculates the deformation velocity and updated shape
 * for a 3d eel.
 */

class IBEELKinematics3d : public IBAMR::ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    IBEELKinematics3d(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      IBTK::LDataManager* l_data_manager,
                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                      bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBEELKinematics3d();

    /*!
     * \brief Set kinematics velocity for the ed eel at specified time.
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
     * \brief Set the shape of eel at specified time. The shape should have
     * its center of mass at the origin, with appropriate rigid body rotation applied
     * to it.
     * \see IBAMR::ConstraintIBKinematics::setShape
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * \brief Get the shape of eel on the specified level. The shape should have
     * its center of mass at the origin, with appropriate rigid body rotation applied
     * to it.
     * \see IBAMR::ConstraintIBKinematics::getShape
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const;

    /*!
     * \brief Write state necessary for restarted runs.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief The default constructor is not implemented and should not be used.
     */
    IBEELKinematics3d();

    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    IBEELKinematics3d(const IBEELKinematics3d& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    IBEELKinematics3d& operator=(const IBEELKinematics3d& that);

    /*!
     * \brief Calculate the coefficients of the cubic polynomial
     * from the curvature coefficients.
     */
    void getInterpCoefs();

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief Set eel body shape related data.
     */
    void setImmersedBodyLayout(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set eel kinematics velocity.
     */
    void setEelSpecificVelocity(const double time,
                                const std::vector<double>& incremented_angle_from_reference_axis,
                                const std::vector<double>& center_of_mass,
                                const std::vector<double>& tagged_pt_position);

    /*!
     * The following vector is used to store eel body shape specific data.
     * d_IBPts is used to store the no of material points which represents a cross section.
     * d_IBWidthHeight stores width and height of a cross section.
     * width of cross section of eel varies with arc length.
     */
    std::vector<std::pair<int, int> > d_IBPts;
    std::vector<std::pair<double, double> > d_IBWidthHeight;

    /*!
     * Eulerian Mesh width parameters.
     */
    std::vector<double> d_mesh_width;

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
     * Curvatures from input file.
     */
    SAMRAI::tbox::Array<double> d_interp_coefs;

    /*!
     * Tau tail from the input file.
     */
    double d_tau_tail;

    /*!
     * Time period of traveling wave.
     */
    double d_time_period;

    /*!
     * Initial angle of the body axis from the horizontal.
     */
    double d_initAngle_bodyAxis_x;

    /*!
     * No. of points on eel.
     */
    int d_HeadNs, d_HeadTailNs, d_TailNs, d_BodyNs;

}; // IBEELKinematics3d

} // IBAMR

#endif //#ifndef included_IBEELKinematics3d
