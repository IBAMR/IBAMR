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

#ifndef included_KnifeFishKinematics
#define included_KnifeFishKinematics

///////////////////////////////////////// INCLUDES //////////////////////////////////////////

#include "ibamr/ConstraintIBKinematics.h"

#include "PatchHierarchy.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <map>
#include <vector>

namespace IBAMR
{
/*!
 * \brief KnifeFishKinematics is a concrete class which calculates the deformation velocity and updated shape
 * for 3D knifefish. Example taken from:
 *
 *  Bhalla et al. A unified mathematical framework and an adaptive numerical method for
 *  fluid-structure interaction with rigid, deforming, and elastic bodies. J Comput Phys, 250:446-476 (2013).
 */
class KnifeFishKinematics : public ConstraintIBKinematics
{
public:
    /*!
     * \brief ctor. This is the only ctor for this object.
     */
    KnifeFishKinematics(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        IBTK::LDataManager* l_data_manager,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                        bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~KnifeFishKinematics();

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
    KnifeFishKinematics();

    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    KnifeFishKinematics(const KnifeFishKinematics& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    KnifeFishKinematics& operator=(const KnifeFishKinematics& that);

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief Set the knifefish velocity.
     */
    void setKnifefishSpecificVelocity(const double time);

    /*!
     * Deformational velocity and shape vectors.
     */
    std::vector<std::vector<double> > d_kinematics_vel;
    std::vector<std::vector<double> > d_shape;

    /*!
     * Radius and angle of excursion of the fin
     */
    std::vector<double> d_vec_radius, d_vec_theta, d_vec_coord;

    /*!
     * Name of the object.
     */
    std::string d_body_name;

    /*!
     * fin length, wave number, angular frequency, Theta_max for the fin.
     */
    double d_fin_length, d_kappa, d_omega, d_theta_max;

    /*!
     * Global index of starting point of fin.
     */
    int d_fin_start_idx;

    /*!
     * current(t) and new time(t+dt)
     */
    double d_current_time, d_new_time;

}; // KnifeFishKinematics

} // namespace IBAMR

#endif // included_KnifeFishKinematics
