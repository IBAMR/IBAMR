// Filename: RigidBodyKinematics.h
//
// Copyright (c) 2002-2017, Amneet Bhalla and Nishant Nangia
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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic SAMRAI objects
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/ConstraintIBKinematics.h>
#include <ibtk/CartGridFunctionSet.h>

#include <muParser.h>

#ifndef included_IBAMR_multiphase_flow_RigidBodyKinematics
#define included_IBAMR_multiphase_flow_RigidBodyKinematics

class RigidBodyKinematics : public IBAMR::ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    RigidBodyKinematics(const std::string& object_name,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        IBTK::LDataManager* l_data_manager,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                        bool register_for_restart = true)
        : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
          d_parser_time(0.0),
          d_center_of_mass(3, 0.0),
          d_incremented_angle_from_reference_axis(3, 0.0),
          d_tagged_pt_position(3, 0.0)
    {
        // NOTE: Parent class constructor registers class with the restart manager, sets object name.

        // Read-in kinematics velocity functions
        for (int d = 0; d < NDIM; ++d)
        {
            const std::string postfix = "_function_" + std::to_string(d);
            std::string key_name = "kinematics_velocity" + postfix;

            if (input_db->isString(key_name))
            {
                d_kinematicsvel_function_strings.push_back(input_db->getString(key_name));
            }
            else
            {
                d_kinematicsvel_function_strings.push_back("0.0");
                TBOX_WARNING("RigidBodyKinematics::RigidBodyKinematics() :\n"
                             << "  no function corresponding to key " << key_name << "found for dimension = " << d
                             << "; using kinematics_vel = 0.0. " << std::endl);
            }

            d_kinematicsvel_parsers.push_back(new mu::Parser());
            d_kinematicsvel_parsers.back()->SetExpr(d_kinematicsvel_function_strings.back());
            d_all_parsers.push_back(d_kinematicsvel_parsers.back());
        }

        // Define constants and variables for the parsers.
        for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
        {
            // Various names for pi.
            (*cit)->DefineConst("pi", M_PI);
            (*cit)->DefineConst("Pi", M_PI);
            (*cit)->DefineConst("PI", M_PI);

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
            }
        }

        // Set the size of vectors.
        const StructureParameters& struct_param = getStructureParameters();
        const int coarsest_ln = struct_param.getCoarsestLevelNumber();
        const int finest_ln = struct_param.getFinestLevelNumber();
        const int total_levels = finest_ln - coarsest_ln + 1;
        d_kinematics_vel.resize(total_levels);

        const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
        for (int ln = 0; ln < total_levels; ++ln)
        {
            const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
            d_kinematics_vel[ln].resize(NDIM);
            for (int d = 0; d < NDIM; ++d)
            {
                d_kinematics_vel[ln][d].resize(nodes_this_ln);
            }
        }

        // No need for restart code in this example.
    } // RigidBodyKinematics

    /*!
     * \brief Destructor.
     */
    virtual ~RigidBodyKinematics()
    {
        for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
        {
            delete (*cit);
        }

        return;

    } //~RigidBodyKinematics

    /*!
     * \brief Set kinematics velocity at new time for rigid body.
     */
    virtual void setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position)
    {
        d_new_time = time;
        d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
        d_center_of_mass = center_of_mass;
        d_tagged_pt_position = tagged_pt_position;

        setRigidBodySpecificVelocity(
            d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);

        d_current_time = d_new_time;
        return;

    } // setNewKinematicsVelocity

    /*!
     * \brief Get the kinematics velocity at new time for rigid body on the specified level.
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const
    {
        static const StructureParameters& struct_param = getStructureParameters();
        static const int coarsest_ln = struct_param.getCoarsestLevelNumber();

#ifdef DEBUG_CHECK_ASSERTIONS
        static const int finest_ln = struct_param.getFinestLevelNumber();
        TBOX_ASSERT(coarsest_ln <= level && level <= finest_ln);
#endif

        return d_kinematics_vel[level - coarsest_ln];

    } // getKinematicsVelocity

    /*!
     * \brief Set the shape of rigid body at new time on all levels.
     */
    virtual void setShape(const double /*time*/, const std::vector<double>& /*incremented_angle_from_reference_axis*/)
    {
        // intentionally left blank
        return;

    } // setShape

    /*!
     * \brief Get the shape of rigid body at new time on the specified level.
     */
    virtual const std::vector<std::vector<double> >& getShape(const int /*level*/) const
    {
        return d_shape;
    } // getShape

    // No need for database code in this example
private:
    /*!
     * \brief Deleted default ctor.
     */
    RigidBodyKinematics();

    /*!
     * \brief Deleted default copy ctor.
     */
    RigidBodyKinematics(const RigidBodyKinematics& from);

    /*!
     * \brief Deleted default assignment.
     */
    RigidBodyKinematics& operator=(const RigidBodyKinematics& that);

    // No need for restart code in this example

    /*!
     * \brief Set rigid body velocity.
     */
    void setRigidBodySpecificVelocity(const double time,
                                      const std::vector<double>& /*incremented_angle_from_reference_axis*/,
                                      const std::vector<double>& center_of_mass,
                                      const std::vector<double>& /*tagged_pt_position*/)
    {
        std::vector<double> vel_parser(NDIM);
        d_parser_time = time;
        for (int d = 0; d < NDIM; ++d) d_parser_posn[d] = center_of_mass[d];
        for (int d = 0; d < NDIM; ++d) vel_parser[d] = d_kinematicsvel_parsers[d]->Eval();

        static const StructureParameters& struct_param = getStructureParameters();
        static const int coarsest_ln = struct_param.getCoarsestLevelNumber();
        static const int finest_ln = struct_param.getFinestLevelNumber();
        static const int total_levels = finest_ln - coarsest_ln + 1;
        static const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();

        for (int ln = 0; ln < total_levels; ++ln)
        {
            const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
            for (int d = 0; d < NDIM; ++d)
            {
                for (int idx = 0; idx < nodes_this_ln; ++idx) d_kinematics_vel[ln][d][idx] = vel_parser[d];
            }
        }

        return;

    } // setRigidBodySpecificVelocity

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
     * New kinematics velocity. New shape of the body.
     *
     * \NOTE Current velocity is always equal to new velocity. Position is
     * updated via CONSTRAINT_VELOCITY method, so new shape is not filled in.
     */
    std::vector<std::vector<std::vector<double> > > d_kinematics_vel;
    std::vector<std::vector<double> > d_shape;

    /*!
     * Save COM, tagged point position and incremented angle from reference axis for restarted runs.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;
}; // RigidBodyKinematics

#endif
