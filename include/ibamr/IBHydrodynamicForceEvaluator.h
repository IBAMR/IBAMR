// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

#ifndef included_IBHydrodynamicForceEvaluator
#define included_IBHydrodynamicForceEvaluator

#include <ibamr/config.h>

#include <ibtk/ibtk_utilities.h>

#include "Box.h"
#include "CellVariable.h"
#include "RobinBcCoefStrategy.h"
#include "SideVariable.h"
#include "tbox/Pointer.h"
#include "tbox/Serializable.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include <Eigen/Core>
#include <Eigen/Geometry>
IBTK_ENABLE_EXTRA_WARNINGS

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace appu
{
template <int DIM>
class VisItDataWriter;
} // namespace appu
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class PatchLevel;
} // namespace hier
namespace pdat
{
template <int DIM>
class SideIndex;
} // namespace pdat
namespace tbox
{
class Database;
template <class TYPE>
class Pointer;
} // namespace tbox
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHydrodynamicForceEvaluator computes hydrodynamic force and
 * torque on immersed bodies. The class uses Reynolds transport theorem to integrate
 * momentum over a Cartesian box region that moves with an arbitrary rigid body
 * translation velocity.
 *
 * References
 * Flavio Noca, <A HREF="http://thesis.library.caltech.edu/3081/1/Noca_f_1997.pdf">On the evaluation
 * of time-dependent fluid-dynamic forces on bluff bodies.</A>
 *
 * Nangia et al., <A HREF="https://www.sciencedirect.com/science/article/pii/S0021999117305016">A moving control volume
 * approach to computing
 * hydrodynamic forces and
 * torques on immersed bodies.
 *
 * \note  The Cartesian box should enclose the body entirely.
 * \note  Various IB methods need to provide linear and angular momentum of the
 *  enclosed body to the class.
 */
class IBHydrodynamicForceEvaluator : public SAMRAI::tbox::Serializable
{
public:
    /*!
     * \brief Default constructor.
     *
     * \param rho Fluid/structure density.
     * \param mu Fluid Viscosity.
     * \param current_time Current integration time.
     */
    IBHydrodynamicForceEvaluator(std::string object_name,
                                 double rho,
                                 double mu,
                                 double current_time,
                                 bool register_for_restart = true);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~IBHydrodynamicForceEvaluator();

    struct IBHydrodynamicForceObject
    {
        // Since this structure has Eigen members we must explicitly override
        // its new and delete operators:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        // Force, torque, and momentum of the body.
        IBTK::Vector3d F_current, T_current, P_current, L_current;
        IBTK::Vector3d F_new, T_new, P_new, L_new;

        // Origin of the r vector from which torque is calculated
        IBTK::Vector3d r0;

        // Momentum of the box
        IBTK::Vector3d P_box_current, L_box_current;
        IBTK::Vector3d P_box_new, L_box_new;

        // Velocity of the box
        IBTK::Vector3d box_u_current, box_u_new;

        // Integration domain.
        IBTK::Vector3d box_X_lower_current, box_X_upper_current;
        IBTK::Vector3d box_X_lower_new, box_X_upper_new;

        // Box volume (area in 2D)
        double box_vol_current, box_vol_new;

        // Indicator variable index of the control volume for plotting
        int inside_strct_idx;

        // Structure details.
        int strct_id;

        // File streams associated for the output.
        std::ofstream *drag_CV_stream, *torque_CV_stream;

    }; // IBHydrodynamicForceObject

    /*!
     * \brief Register structure ID, level number, and integration domain
     * with the class.
     *
     * \param box_X_lower Initial (typically at time = 0) position of lower left
     * corner of the integration domain.
     *
     * \param box_X_upper Initial (typically at time = 0) position of upper right
     * corner of the integration domain.
     *
     * \param box_vel Initial (typically at time = 0) velocity of the integration
     *  domain in three Cartesian directions.
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     */
    void registerStructure(IBTK::Vector3d& box_X_lower,
                           IBTK::Vector3d& box_X_upper,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                           const IBTK::Vector3d& box_vel = IBTK::Vector3d::Zero(),
                           int strct_id = 0);

    /*!
     * \brief Update the domain of integration as a result of structure motion.
     * This should be called before computeLaggedMomentumIntegral to insure that the correct box is being used.
     *
     * \param box_vel_new Velocity of the integration domain in three Cartesian directions
     * during the time step.
     *
     * \param dt Time step from the integrator. Pass IBHierarchyIntegrator::getMaximumTimeStepSize() instead of
     * new_time - old_time to avoid floating point errors from subtraction
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     *
     */
    void updateStructureDomain(const IBTK::Vector3d& box_vel_new,
                               double dt,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                               int strct_id = 0);

    /*!
     * \brief Set the origin of the position vector used to compute torques
     *
     * \param X0 A 3D vector corresponding to the origin of the position vector used to compute torques
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     */
    void setTorqueOrigin(const IBTK::Vector3d& X0, int strct_id = 0);

    /*!
     * \brief Preprocess data for the current timestep.
     */
    virtual void preprocessIntegrateData(double current_time, double new_time);

    /*!
     * \brief Get access to hydrodynamic data of the given structure id.
     */
    const IBHydrodynamicForceObject& getHydrodynamicForceObject(int strct_id);

    /*!
     * \brief Compute the momentum integral for the velocity variable at the previous time step over the new control
     * volume.
     * This should be called before advancing the hierarchy but after updateStructureDomain
     *
     * \param u_old_idx Patch index of velocity variable (before advancing the hierarchy) with appropriate ghost cell
     * width.
     *
     * \param wgt_sc_idx Patch index of volume weights associated with faces.
     *
     * \param u_src_bc_coef Velocity boundary condition object maintained by the integrator.
     *
     * \param p_src_bc_coef Pressure boundary condition object maintained by the integrator.
     *
     */

    void computeLaggedMomentumIntegral(int u_old_idx,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                       const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_src_bc_coef =
                                           std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>());

    /*!
     * \brief Update the new momenta of the bodies within the structure.
     * This should be called after advancing the hierarchy.
     *
     * \param P_strct_new Linear momentum of the structure after the integration.
     *
     * \param L_strct_new Angular momentum of the structure after the integration.
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     */

    void
    updateStructureMomentum(const IBTK::Vector3d& P_strct_new, const IBTK::Vector3d& L_strct_new, int strct_id = 0);

    /*!
     * \brief Compute hydrodynamic force.
     * This should be called after advancing the hierarchy.
     *
     * \param u_idx Patch index of velocity variable (after advancing the hierarchy) with appropriate ghost cell width.
     *
     * \param p_idx Patch index of pressure variable with appropriate ghost cell width.
     *
     * \param f_idx Patch index of body force variable with appropriate ghost cell width.
     *
     * \param dt Time step from the integrator. Pass IBHierarchyIntegrator::getMaximumTimeStepSize() instead of
     * new_time - old_time to avoid floating point errors from subtraction
     *
     * \param u_src_bc_coef Velocity boundary condition object maintained by the integrator.
     *
     * \param p_src_bc_coef Pressure boundary condition object maintained by the integrator.
     *
     */
    virtual void computeHydrodynamicForce(int u_idx,
                                          int p_idx,
                                          int f_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                          double dt,
                                          const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_src_bc_coef =
                                              std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(),
                                          SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_src_bc_coef = NULL);

    /*!
     * \brief Postprocess data for the next timestep.
     */
    virtual void postprocessIntegrateData(double current_time, double new_time);

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \brief Create the control volume plot data and register it with the VisIt data writer
     *
     * \param strct_id A unique integer id to associate with an integration domain.
     */
    void registerStructurePlotData(SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_data_writer,
                                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                   int strct_id = 0);

    /*!
     * \brief Update the plot variable for the new location of the control volume box
     *
     * \param strct_id A unique integer id to associate with an integration domain
     */
    void updateStructurePlotData(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                 int strct_id = 0);

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHydrodynamicForceEvaluator(const IBHydrodynamicForceEvaluator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHydrodynamicForceEvaluator& operator=(const IBHydrodynamicForceEvaluator& that) = delete;

    /*!
     * \brief Reset weight of the cell face to face area.
     */
    void resetFaceAreaWeight(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Reset weight of the cell face to cell volume.
     */
    void resetFaceVolWeight(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Allocate and fill velocity and pressure patch data.
     */
    void fillPatchData(const int u_src_idx,
                       const int p_src_idx,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                       const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_src_bc_coef,
                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_src_bc_coef,
                       const double fill_time);

    /*!
     * \brief Compute the physical coordinate of a given side index
     */
    void getPhysicalCoordinateFromSideIndex(IBTK::Vector3d& side_coord,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                            const SAMRAI::pdat::SideIndex<NDIM> side_idx,
                                            const int axis);

    /*!
     * \brief Object name.
     */
    std::string d_object_name;

    /*!
     * \brief Fluid density and viscosity.
     */
    double d_rho, d_mu;

    /*!
     * \brief Current integrator time.
     */
    double d_current_time;

    /*!
     * \brief Fluid velocity and pressure with appropriate ghost width.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_p_var;
    int d_u_idx, d_p_idx;

    /*!
     * \brief Patch data index for face weights.
     */
    int d_face_wgt_sc_idx, d_vol_wgt_sc_idx;

    /*!
     * \brief Data structure encapsulating hydrodynamic force on an object.
     */
    std::map<int, IBHydrodynamicForceObject> d_hydro_objs;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHydrodynamicForceEvaluator
