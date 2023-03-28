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

#ifndef included_RigidBodyKinematics
#define included_RigidBodyKinematics

///////////////////////////// INCLUDES ///////////////////////////////////////
#include "ibamr/ConstraintIBKinematics.h"

#include "PatchHierarchy.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "Eigen/Core"
IBTK_ENABLE_EXTRA_WARNINGS

#include <iostream>
#include <map>
#include <vector>

//////////////////////////// CLASS DEFORMATIONAL KINEMATICS //////////////////

namespace IBAMR
{
/*!
 * \brief RigidBodyKinematics is a concrete class which provides rigid body velocity
 * at the marker points
 */

class RigidBodyKinematics : public ConstraintIBKinematics

{
public:
    /*!
     * \brief ctor. This is the only ctor for this object.
     */
    RigidBodyKinematics(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        IBTK::LDataManager* l_data_manager,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                        bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~RigidBodyKinematics();

    /*!
     * \brief Functions providing time-dependent translational and rotational velocity.
     */
    typedef void (*RigidVelFcn)(const double time, Eigen::Vector3d& rigid_vel);
    struct RigidKinematicsData
    {
        RigidKinematicsData() : d_trans_vel_fcn(NULL), d_rot_vel_fcn(NULL)
        {
            return;
        } // RigidKinematicsData

        RigidKinematicsData(RigidVelFcn trans_vel_fcn, RigidVelFcn rot_vel_fcn = NULL)
            : d_trans_vel_fcn(trans_vel_fcn), d_rot_vel_fcn(rot_vel_fcn)
        {
            return;
        } // RigidKinematicsData

        RigidVelFcn d_trans_vel_fcn, d_rot_vel_fcn;
    };

    /*!
     * \brief Set the COM translational and rotational velocity.
     */
    void registerRigidBodyKinematics(RigidVelFcn tran_vel_fcn = NULL, RigidVelFcn rot_vel_fcn = NULL);

    /*!
     * \brief Set kinematics velocity for the rigid body.
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
     * \brief The copy constructor is not implemented and should not be used.
     */
    RigidBodyKinematics(const RigidBodyKinematics& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    RigidBodyKinematics& operator=(const RigidBodyKinematics& that);

    /*!
     * \brief Set the rigid body related data.
     */
    void setImmersedBodyLayout(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief Set kinematics velocity at the marker points.
     */
    void setRigidBodyVelocity(const double time,
                              const std::vector<double>& incremented_angle_from_reference_axis,
                              const std::vector<double>& center_of_mass,
                              const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Kinematics data.
     */
    RigidKinematicsData d_kinematics_data;

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

}; // RigidBodyKinematics

} // namespace IBAMR
#endif //#ifndef included_RigidBodyKinematics
