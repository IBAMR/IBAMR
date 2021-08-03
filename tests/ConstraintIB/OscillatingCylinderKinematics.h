// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

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

} // namespace IBAMR

#endif // included_OscillatingCylinderKinematics
