// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_ConstraintIBKinematics
#define included_IBAMR_ConstraintIBKinematics

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/LData.h"
#include "ibtk/LDataManager.h"

#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

#include <string>
#include <utility>
#include <vector>

namespace IBTK
{
class LDataManager;
} // namespace IBTK
namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Class ConstraintIBKinematics encapsulates structure information and provides abstraction to get
 * kinematics (deformational or imposed) of immersed structure to ConstraintIBMethod class.
 */
class ConstraintIBKinematics : public virtual SAMRAI::tbox::DescribedClass, public SAMRAI::tbox::Serializable
{
public:
    class StructureParameters
    {
    public:
        /*!
         * \brief Constructor.
         */
        StructureParameters(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, IBTK::LDataManager* l_data_manager);

        /*!
         * \brief Lagrangian point to tag on this structure.
         */
        inline int getTaggedPtIdx() const
        {
            return d_tagged_pt_idx;
        } // getTaggedPtIdx

        /*!
         * \brief Get the unlocked components of translational momentum.
         */
        inline SAMRAI::tbox::Array<int> getCalculateTranslationalMomentum() const
        {
            return d_calculate_trans_mom;
        } // getCalculateTranslationalMomentum

        /*!
         * \brief Get the unlocked components of rotational momentum.
         */
        inline SAMRAI::tbox::Array<int> getCalculateRotationalMomentum() const
        {
            return d_calculate_rot_mom;
        } // getCalculateRotationalMomentum

        /*!
         * \brief Check if the structure has translational degree unlocked.
         */
        inline bool getStructureIsSelfTranslating() const
        {
            return d_struct_is_self_translating;
        } // getStructureIsSelfTranslating

        /*!
         * \brief Check if the structure has rotational degree unlocked.
         */
        inline bool getStructureIsSelfRotating() const
        {
            return d_struct_is_self_rotating;
        } // getStructureIsSelfRotating

        /*!
         * \brief The coarsest level on which the structure resides.
         */
        inline int getCoarsestLevelNumber() const
        {
            return d_coarsest_ln;
        } // getCoarsestLevelNumber

        /*!
         * \brief The finest level on which the structure resides.
         */
        inline int getFinestLevelNumber() const
        {
            return d_finest_ln;
        } // getFinestLevelNumber

        /*!
         * \brief Global Lagrangian indices managed for this structure.
         */
        inline const std::vector<std::pair<int, int> >& getLagIdxRange() const
        {
            return d_idx_range;
        } // getLagIdxRange

        /*!
         * \brief Total number of Lagrangian nodes managed for this structure.
         */
        inline int getTotalNodes() const
        {
            return d_total_nodes;
        } // getTotalNodes

        /*!
         * \brief Lagrangian nodes update method for this structure.
         */
        inline std::string getPositionUpdateMethod() const
        {
            return d_lag_position_update_method;
        } // getPositionUpdateMethod

    private:
        std::string d_lag_position_update_method;
        int d_coarsest_ln, d_finest_ln;
        std::vector<std::pair<int, int> > d_idx_range;
        int d_total_nodes;
        int d_tagged_pt_idx;
        SAMRAI::tbox::Array<int> d_calculate_trans_mom, d_calculate_rot_mom;
        bool d_struct_is_self_translating, d_struct_is_self_rotating;
    };

    /*!
     * \brief Constructor.
     */
    ConstraintIBKinematics(std::string object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           IBTK::LDataManager* l_data_manager,
                           bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~ConstraintIBKinematics();

    /*!
     * \brief Get the object enclosing this structure's parameters.
     */
    inline const StructureParameters& getStructureParameters() const
    {
        return d_struct_param;
    } // getStructureParameters

    /*!
     * \brief Set the kinematics velocity (deformational or imposed) at Lagrangian points managed by this object.
     *
     * \param time Time at which kinematics velocity is to be set.
     *
     * \param incremented_angle_from_reference_axis Angle made with x,y & z axis due to rigid rotational velocity.
     * \f$ \theta_n = \theta_{n-1} + \omega_{n-1} \triangle t \f$
     *
     * \param center_of_mass COM of the structure at the given time.
     *
     * \param tagged_pt_position Coordinates of the tagged point of the structure.
     */
    virtual void setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position) = 0;

    /*!
     * \brief Get the kinematics velocity for the structure on the specified level.
     *
     * \param level Kinematics velocity of the structure on this level.
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const = 0;

    /*!
     * \brief Set the shape of structure at the required time. The shape should have its center of
     * mass at origin, with appropriate rigid body ratation applied to it.
     *
     * \param time Time at which kinematics velocity is to be set.
     *
     * \param incremented_angle_from_reference_axis Angle made with x,y & z axis due to rigid rotational velocity.
     * \f$ \theta_n = \theta_{n-1} + \omega_{n-1} \triangle t \f$.
     *
     * \note The call to setShape() is made after setKinematicsVelocity( ) with the same value of the common arguments.
     *
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis) = 0;

    /*!
     * \brief Get the shape of structure on the specified level.
     *
     * \param level
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const = 0;

    /*!
     * \brief Write out object state to the given database.
     *
     * \note An empty default implementation is provided.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

protected:
    /*!
     * Name of the object.
     */
    std::string d_object_name;

    /*!
     * If the object is registred for restart.
     */
    bool d_registered_for_restart = false;

private:
    /*!
     * \brief Deleted default ctor.
     */
    ConstraintIBKinematics() = delete;

    /*!
     * \brief Deleted default copy ctor.
     */
    ConstraintIBKinematics(const ConstraintIBKinematics& from) = delete;

    /*!
     * \brief Deleted default assignment.
     */
    ConstraintIBKinematics& operator=(const ConstraintIBKinematics& that) = delete;

    /*!
     * \brief Object enclosing all the parameters of the structure.
     */
    StructureParameters d_struct_param;
};
} // namespace IBAMR

/////////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_ConstraintIBKinematics
