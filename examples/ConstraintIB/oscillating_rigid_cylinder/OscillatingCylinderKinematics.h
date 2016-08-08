// Filename: OscillatingCylinderKinematics.h
//
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
// POSSIBILITY OF SUCH DAMAGE

#ifndef included_OscillatingCylinderKinematics
#define included_OscillatingCylinderKinematics

///////////////////////////////////////// INCLUDES //////////////////////////////////////////

// IBAMR INCLUDES
#include "ibamr/ConstraintIBKinematics.h"

// C++ INCLUDES
#include <iostream>
#include <map>
#include <vector>

// SAMRAI INCLUDES
#include "PatchHierarchy.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

/////////////////////////////////////// FORWARD DECLARATION ////////////////////////////////

namespace mu
{
class Parser;
} // namespace mu

namespace IBAMR
{
/*!
 * \brief Class OscillatingCylinderKinematics provides definition for base ConstraintIBKinematics class.
 */
class OscillatingCylinderKinematics : public ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    OscillatingCylinderKinematics(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                  IBTK::LDataManager* l_data_manager,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                  bool register_for_restart = true);

    /*!
     * Destructor.
     */
    virtual ~OscillatingCylinderKinematics();

    /*!
     * Set kinematics velocity at new time for OscillatingCylinder.
     */
    virtual void setKinematicsVelocity(const double new_time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position);

    /*!
     * Get the kinematics velocity at new time for OscillatingCylinder on the specified level.
     *
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const;

    /*!
     * Get the kinematics velocity at current time for OscillatingCylinder on the specified level.
     *
     */
    virtual const std::vector<std::vector<double> >& getCurrentKinematicsVelocity(const int level) const;

    /*!
     * Set the shape of OscillatingCylinder at new time for the structure on all levels.
     */
    virtual void setShape(const double new_time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * Get the shape of structure at new time  on the specified level.
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const;

    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * Deleted default ctor.
     */
    OscillatingCylinderKinematics();

    /*!
     * Deleted default copy ctor.
     */
    OscillatingCylinderKinematics(const OscillatingCylinderKinematics& from);

    /*!
     * Deleted default assignment.
     */
    OscillatingCylinderKinematics& operator=(const OscillatingCylinderKinematics& that);

    /*!
     *Intialize object for Restarted runs
     */
    void getFromRestart();

    /*!
     * Calculate the kinematic velocity of the OscillatingCylinder
     */
    void setOscillatingCylinderSpecificVelocity(const double time);

    /*!
    * New and current kinematics velocity. New shape of the body.
    */
    std::vector<std::vector<std::vector<double> > > d_new_kinematics_vel, d_current_kinematics_vel;
    std::vector<std::vector<double> > d_new_shape;

    /*!
     * Max speed and frequency of oscillation
     */
    double d_Uinf, d_freq;

    /*!
     * Prescribed translational velocity
     */
    double d_prescribed_trans_vel;

    /*!
     * current(t) and new time(t+dt)
     */
    double d_current_time, d_new_time;

}; // OscillatingCylinderKinematics

} // IBAMR

#endif // included_OscillatingCylinderKinematics
