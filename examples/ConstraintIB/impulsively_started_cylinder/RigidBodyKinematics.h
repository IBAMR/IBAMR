// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2019 by the IBAMR developers
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

///////////////////////////////////////// INCLUDES //////////////////////////////////////////

// IBAMR INCLUDES
#include "ibamr/ConstraintIBKinematics.h"

// IBTK INCLUDES
#include "ibtk/LDataManager.h"

// SAMRAI INCLUDES
#include "PatchHierarchy.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

// C++ INCLUDES
#include "string"
#include "vector"

/////////////////////////////////////// FORWARD DECLARATION ////////////////////////////////

namespace mu
{
class Parser;
} // namespace mu

///////////////////////////////////////////////////////////////// CLASS DEFORMATIONAL KINEMATICS //////////////////

namespace IBAMR
{
/*!
 * \brief Class RigidBodyKinematics provides definition for base ConstraintIBKinematics class.
 */
class RigidBodyKinematics : public ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
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
     * \brief Set kinematics velocity at new time for rigid body.
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
     * \brief Set the shape of rigid body at new time on all levels.
     * \see IBAMR::ConstraintIBKinematics::setShape
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * \brief Get the shape of rigid body at new time on the specified level.
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
    RigidBodyKinematics();

    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    RigidBodyKinematics(const RigidBodyKinematics& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    RigidBodyKinematics& operator=(const RigidBodyKinematics& that);

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief Set rigid body velocity.
     */
    void setRigidBodySpecificVelocity(const double time,
                                      const std::vector<double>& incremented_angle_from_reference_axis,
                                      const std::vector<double>& center_of_mass,
                                      const std::vector<double>& tagged_pt_position);

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser*> d_kinematicsvel_parsers;
    std::vector<mu::Parser*> d_all_parsers;

    /*!
     * Input kiematics velocity functions.
     */
    std::vector<std::string> d_kinematicsvel_function_strings;

    /*!
     * Parser variables.
     */
    mutable double d_parser_time;
    mutable IBTK::Point d_parser_posn;

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

}; // RigidBodyKinematics

} // namespace IBAMR

#endif // included_RgidBodyKinematics
