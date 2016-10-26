// Filename: IBHydrodynamicForceEvaluator.h
// Created on 22 Oct 2016 by Amneet Bhalla
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

#ifndef included_IBHydrodynamicForceEvaluator
#define included_IBHydrodynamicForceEvaluator

#include <map>
#include <vector>

#include "Box.h"
#include "Eigen/Core"
#include "tbox/Serializable.h"

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace tbox
{
template <class TYPE>
class Pointer;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

namespace IBAMR
{
class IBMethod;
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHydrodynamicForceEvaluator computes hydrodynamic force and
 * torque on immersed bodies.
 */
class IBHydrodynamicForceEvaluator : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Default constructor.
     *
     * \param rho Fluid/structure density.
     * \param mu Fluid Viscosity.
     */
    IBHydrodynamicForceEvaluator(const std::string& object_name,
                                 double rho,
                                 double mu,
                                 bool register_for_restart = true);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBHydrodynamicForceEvaluator();

    struct IBHydrodynamicForceObject
    {
        // Structure details.
        int strct_id, strct_ln;

        // Force, torque, and momentum of the body.
        Eigen::Vector3d F_current, T_current, P_current, L_current;
        Eigen::Vector3d F_new, T_new, P_new, L_new;

        // Integration domain.
        SAMRAI::hier::Box<NDIM> box_current, box_new;
        Eigen::Vector3d P_box_current, L_box_current;
        Eigen::Vector3d P_box_new, L_box_new;
        Eigen::Vector3d box_u_current, box_u_new;
        Eigen::Vector3d box_X_lower_current, box_X_upper_current;
        Eigen::Vector3d box_X_lower_new, box_X_upper_new;

    }; // IBHydrodynamicForceObject

    /*!
     * \brief Register structure ID, level number, and integration domain
     * with the class.
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     *
     * \param strct_ln Integer representing the level number on which the structure resides.
     *
     * \param box_vel Initial (typically at time = 0) velocity of the integration
     *  domain in three Cartesian directions.
     *
     * \param box_X_lower Initial (typically at time = 0) position of lower left
     * corner of the integration domain.
     *
     * \param box_X_upper Initial (typically at time = 0) position of upper right
     * corner of the integration domain.
     */
    void registerStructure(int strct_id,
                           int strct_ln,
                           const Eigen::Vector3d& box_vel,
                           const Eigen::Vector3d& box_X_lower,
                           const Eigen::Vector3d& box_X_upper);

    /*!
     * \brief Update the domain of integration as a result of structure motion.
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     *
     * \param strct_ln Integer representing the level number on which the structure resides.
     *
     * \param current_time Current time of integration.
     *
     * \param new_time Time after integration.
     *
     * \param box_vel_new Velocity of the integration domain in three Cartesian directions
     * during the time step.
     *
     * \param P_strct_new Linear momentum of the structure after the integration.
     *
     * \param L_strct_new Angular momentum of the structure after the integration.
     */
    void updateStructureDomain(int strct_id,
                               int strct_ln,
                               double current_time,
                               double new_time,
                               const Eigen::Vector3d& box_vel_new,
                               const Eigen::Vector3d& P_strct_new,
                               const Eigen::Vector3d& L_strct_new);

    /*!
     * \brief Preprocess data for the current timestep.
     */
    virtual void preprocessIntegrateData(double current_time, double new_time);

    /*!
     * \brief Get access to hydrodynamic data of the given structure id.
     */
    const IBHydrodynamicForceObject& getHydrodynamicForceObject(int strct_id);

    /*!
     * \brief Compute hydrodynamic force.
     */
    virtual void computeHydrodynamicForce(int u_idx,
                                          int p_idx,
                                          int f_idx,
                                          int wgt_sc_idx,
                                          int wgt_cc_idx,
                                          const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& F_data,
                                          const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data,
                                          const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& U_data,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                          int coarsest_ln,
                                          int finest_ln,
                                          double current_time,
                                          double new_time,
                                          IBMethod* ib_method);

    /*!
     * \brief Postprocess data for the next timestep.
     */
    virtual void postprocessIntegrateData(double current_time, double new_time);

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHydrodynamicForceEvaluator(const IBHydrodynamicForceEvaluator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHydrodynamicForceEvaluator& operator=(const IBHydrodynamicForceEvaluator& that);

    /*!
     * \brief Compute weight of the cell face.
     */
    void computeFaceWeight(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Object name.
     */
    std::string d_object_name;

    /*!
     * \brief Fluid density and viscosity.
     */
    double d_rho, d_mu;

    /*!
     * \brief Patch data index for face weights.
     */
    int d_face_wgt_sc_idx;

    /*!
     * \brief Data structure encapsulating hydrodynamic force on an object.
     */
    std::map<int, IBHydrodynamicForceObject> d_hydro_objs;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHydrodynamicForceEvaluator
