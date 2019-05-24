// Filename: ConstraintIBMethod.h
// Created on 1 Dec 2011 by Amneet Bhalla
//
// Copyright (c) 2002-2017, Amneet Bhalla and Boyce Griffith
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

#ifndef included_IBAMR_ConstraintIBMethod
#define included_IBAMR_ConstraintIBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <fstream>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "LocationIndexRobinBcCoefs.h"
#include "PoissonSpecifications.h"
#include "VariableContext.h"
#include "ibamr/ConstraintIBKinematics.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBMethod.h"
#include "ibtk/CCLaplaceOperator.h"
#include "ibtk/CCPoissonPointRelaxationFACOperator.h"
#include "ibtk/FACPreconditioner.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/PETScKrylovPoissonSolver.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
/*!
 * \brief Class ConstraintIBMethod implements the rigidity constraint for rigid and deforming bodies
 * using the constraint based IB method.
 *
 * References
 *  Bhalla et al. A unified mathematical framework and an adaptive numerical method for
 *  fluid-structure interaction with rigid, deforming, and elastic bodies. J Comput Phys, 250:446-476 (2013).
 */
class ConstraintIBMethod : public IBAMR::IBMethod
{
public:
    /*!
     * \brief Constructor
     */
    ConstraintIBMethod(std::string object_name,
                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                       const int no_structures,
                       bool register_for_restart = true);

    /*!
     * \brief Destructor
     */
    ~ConstraintIBMethod();

    /*!
     *  \brief Initialize Hierarchy operators and data at initial time.
     */
    void initializeHierarchyOperatorsandData();

    /*!
     * \brief Register Eulerian variables with base IBStrategy class.
     */
    virtual void registerEulerianVariables() override;

    /*!
     *  \brief Create Lagrangian workspace.
     */
    virtual void preprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     *  \brief Destroy Lagrangian workspace.
     */
    virtual void postprocessIntegrateData(double current_time, double new_time, int num_cycles) override;

    /*!
     * \brief Register kinematics of the immersed structure(s) with this class.
     */
    void registerConstraintIBKinematics(
        const std::vector<SAMRAI::tbox::Pointer<IBAMR::ConstraintIBKinematics> >& ib_kinematics_op);

    /*!
     * \brief Register any preprocess fluid solve callback functions.
     */
    inline void registerPreProcessSolveFluidEquationsCallBackFunction(
        void (*ptr_preprocess_callbackfnc)(const double, const double, const int, void*),
        void* ctx)
    {
        d_prefluidsolve_callback_fns.push_back(ptr_preprocess_callbackfnc);
        d_prefluidsolve_callback_fns_ctx.push_back(ctx);
        return;
    }

    /*!
     * \brief Calculate any body forces for INS solver over here.
     */
    virtual void preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num) override;

    /*!
     * \brief Register any postprocess fluid solve callback functions.
     */
    inline void registerPostProcessSolveFluidEquationsCallBackFunction(
        void (*ptr_postprocess_callbackfnc)(const double, const double, const int, void*),
        void* ctx)
    {
        d_postfluidsolve_callback_fns.push_back(ptr_postprocess_callbackfnc);
        d_postfluidsolve_callback_fns_ctx.push_back(ctx);
        return;
    }

    /*!
     * \brief Apply the FuRMoRP algorithm in the postprocessSolveFluidEquations method.
     */
    virtual void postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num) override;

    /*!
     * \brief Override the forwardEulerStep method of the base IBMethod class.
     */
    virtual void forwardEulerStep(double current_time, double new_time) override;

    /*!
     * \brief Override the midpointStep method of the base IBMethod class.
     */
    virtual void midpointStep(double current_time, double new_time) override;

    /*!
     * \brief Override the putToDatabase method of the base Serializable class.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * \brief Get the volume element associated with material points
     * of all structures.
     */
    inline const std::vector<double>& getVolumeElement()
    {
        return d_vol_element;
    }

    /*!
     * \brief Get the current COM velocity associated with Lagrangian structures.
     */
    inline const std::vector<std::vector<double> >& getCurrentCOMVelocity()
    {
        return d_rigid_trans_vel_current;
    }

    /*!
     * \brief Get LData associated with Lagrange multiplier force field.
     */
    inline const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& getLagrangeMultiplier()
    {
        return d_l_data_U_correction;
    }

    /*!
     * \brief Get the current center of mass for all Lagrangian structures
     */
    inline const std::vector<std::vector<double> >& getCurrentStructureCOM()
    {
        return d_center_of_mass_current;
    }

    /*
     * Set velocity physical boundary options
     */
    inline void setVelocityPhysBdryOp(IBTK::RobinPhysBdryPatchStrategy* u_phys_bdry_op)
    {
        d_u_phys_bdry_op = u_phys_bdry_op;
        return;
    }

    /*
     * Set the volume element for each Lagrangian node for each individual structure
     */
    inline void setVolumeElement(double vol_element, int struct_no)
    {
        d_vol_element[struct_no] = vol_element;
        d_vol_element_is_set[struct_no] = true;
        return;
    }

    /*
     * Set the volume element for each Lagrangian node for all the structures
     */
    inline void setVolumeElement(std::vector<double> vol_element)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(vol_element.size() == ((size_t)d_no_structures));
#endif
        d_vol_element = vol_element;
        d_vol_element_is_set = std::vector<bool>(d_no_structures, true);
    }

    /*
     * \brief Get the total volume for all the Lagrangian structures
     */
    inline const std::vector<double>& getStructureVolume()
    {
        return d_structure_vol;
    }

    /*!
     * \brief Get the total linear momentum for all the Lagrangian structures
     */
    inline const std::vector<std::vector<double> >& getStructureMomentum()
    {
        if (!d_calculate_structure_linear_mom)
        {
            TBOX_ERROR("ConstraintIBMethod::getStructureMomentum() called with calculate_structure_linear_mom = FALSE");
        }
        return d_structure_mom;
    }

    /*!
     * \brief Get the total rotational momentum for all the Lagrangian structures with respect to their COM
     */
    inline const std::vector<std::vector<double> >& getStructureRotationalMomentum()
    {
        if (!d_calculate_structure_rotational_mom)
        {
            TBOX_ERROR(
                "ConstraintIBMethod::getStructureRotationalMomentum() called with calculate_structure_rotational_mom = "
                "FALSE");
        }
        return d_structure_rotational_mom;
    }

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    ConstraintIBMethod() = delete;

    /*!
     * \brief Default copy constructor.
     *
     * \note This copy constructor is not implemented and should not be used.
     */
    ConstraintIBMethod(const ConstraintIBMethod& from) = delete;

    /*!
     * \brief Default assignment operator.
     *
     * \note This assignment operator is not implemented and should not be used.
     */
    ConstraintIBMethod& operator=(const ConstraintIBMethod& that) = delete;

    /*!
     * \brief Get values from input file.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db, const bool from_restart);

    /*!
     * \brief Get values from restart file.
     */
    void getFromRestart();

    /*!
     * \brief Set initial Lagrangian velocity on material points.
     */
    void setInitialLagrangianVelocity();

    /*!
     * \brief Calculate center of mass and moment of inertia of immersed
     * structures.
     */
    void calculateCOMandMOIOfStructures();

    /*!
     * \brief Calculate the kinematics velocity for all structures handled by this class.
     */
    void calculateKinematicsVelocity();

    /*!
     * \brief Calculate momentum of kinematics velocity. This is extraneous momentum
     * that needs to be subtracted from the kinematics velocity.
     */
    void calculateMomentumOfKinematicsVelocity(const int position_handle);

    /*!
     * \brief Calculate volume element associated with material points.
     */
    void calculateVolumeElement();

    /*!
     * \brief Set the counter for this method.
     */
    inline void setCounter()
    {
        ++d_timestep_counter;
        return;
    }

    /*!
     * \brief Set the time at which FuRMoRP is applied.
     */
    void setFuRMoRPTime(const double current_time, const double new_time)
    {
        d_FuRMoRP_current_time = current_time;
        d_FuRMoRP_new_time = new_time;
        return;
    }

    /*!
     * \brief Copy vector.
     */
    void copyFluidVariable(int copy_from_idx, int copy_to);

    /*!
     * \brief Copy density patch data.
     */
    void copyDensityVariable(int copy_from_idx, int copy_to);

    /*!
     * \brief Interpolate fluid solve velocity from Eulerian grid onto the Lagrangian mesh.
     */
    void interpolateFluidSolveVelocity();

    /*!
     * \brief Calculate the rigid translational velocity.
     */
    void calculateRigidTranslationalMomentum();

    /*!
     * \brief Calculate the rigid rotational velocity.
     */
    void calculateRigidRotationalMomentum();

    /*!
     * \brief Calculate current velocity on the material points.
     */
    void calculateCurrentLagrangianVelocity();

    /*!
     * \brief Correct velocity on Lagrangian mesh. Set the velocity on Lagrangian
     * mesh as U_lag_corr = U_trans + Omega X r + U_def - U_interpolated.
     */
    void correctVelocityOnLagrangianMesh();

    /*!
     * \brief Spread the corrected velocity at the Lagrangian mesh to the Eulerian Grid.
     */
    void spreadCorrectedLagrangianVelocity();

    /*!
     * \brief The correction on Eulerian grid can lead to non-divergence free velocity field.
     * We project the corrected velocity field onto a divergence free field.
     */
    void applyProjection();

    /*!
     * \brief Predict the position of structures according to forward Euler step method.
     */
    void updateStructurePositionEulerStep();

    /*!
     * \brief Update the position of structures according to mid point step method.
     */
    void updateStructurePositionMidPointStep();

    /*!
     * \brief Compute U_half = 0.5(U_current + U_new);
     */
    void calculateMidPointVelocity();

    /*!
     * \brief Calculate hydrodynamic drag on the immersed structures.
     */
    void calculateDrag();

    /*!
     * \brief Calculate hydrodynamic torque on the immersed structures.
     */
    void calculateTorque();

    /*!
     * \brief Calculate power spent during swimming.
     */
    void calculatePower();

    /*!
     * \brief Calculate Eulerian Momentum.
     */
    void calculateEulerianMomentum();

    /*!
     * \brief Calculate total translational momentum of all Lagrangian structures
     */
    void calculateStructureMomentum();

    /*!
     * \brief Calculate the total rotational momentum of all Lagrangian structures with respect to their COM
     */
    void calculateStructureRotationalMomentum();

    /*!
     * No of immersed structures.
     */
    const int d_no_structures;

    /*!
     * Pointer to the kinematics of the immersed structures.
     */
    std::vector<SAMRAI::tbox::Pointer<IBAMR::ConstraintIBKinematics> > d_ib_kinematics;

    /*!
     * FuRMoRP apply time.
     */
    double d_FuRMoRP_current_time = 0.0, d_FuRMoRP_new_time = 0.0;

    /*!
     * Volume element associated with material points.
     */
    std::vector<double> d_vol_element;

    /*
     * Whether or not the volume has been set for the material structure
     */
    std::vector<bool> d_vol_element_is_set;

    /*!
     * Volume associated with each immersed structure
     */
    std::vector<double> d_structure_vol;

    /*!
     * Linear momentum associated with each immersed structure
     */
    std::vector<std::vector<double> > d_structure_mom;

    /*!
     * Rotational momentum associated with each immersed structure with respect to their COM
     */
    std::vector<std::vector<double> > d_structure_rotational_mom;

    /*!
     * If divergence free projection is needed after FuRMoRP algorithm?
     */
    bool d_needs_div_free_projection = false;

    /*!
     * Rigid translational velocity of the structures.
     */
    std::vector<std::vector<double> > d_rigid_trans_vel_current, d_rigid_trans_vel_new;

    /*!
     * Rigid rotational velocity of the structures.
     */
    std::vector<std::vector<double> > d_rigid_rot_vel_current, d_rigid_rot_vel_new;

    /*!
     * Incremented angle from x, y and z axis when the body is rotating.
     */
    std::vector<std::vector<double> > d_incremented_angle_from_reference_axis;

    /*!
     * Translational velocity of the structures due to deformational kinematics.
     */
    std::vector<std::vector<double> > d_vel_com_def_current, d_vel_com_def_new;

    /*!
     * Rotational velocity of the structures due to deformational kinematics.
     */
    std::vector<std::vector<double> > d_omega_com_def_current, d_omega_com_def_new;

    /*!
     * Center of mass of the immersed structures.
     */
    std::vector<std::vector<double> > d_center_of_mass_current, d_center_of_mass_new;

    /*!
     * Moment of inertia of the structures.
     */
    IBTK::EigenAlignedVector<Eigen::Matrix3d> d_moment_of_inertia_current, d_moment_of_inertia_new;

    /*!
     * Tag a Lagrangian point (generally eye of the fish) of the immersed structures.
     */
    std::vector<int> d_tagged_pt_lag_idx;

    /*!
     * Coordinates of the tagged points of different structures.
     */
    std::vector<std::vector<double> > d_tagged_pt_position;

    /*!
     * Density of the structures.
     */
    std::vector<double> d_rho_solid;

    /*!
     * Density of the fluid in constant coefficient case.
     */
    double d_rho_fluid = std::numeric_limits<double>::quiet_NaN();

    /*!
     * Whether or not the density from the integrator is constant
     */
    bool d_rho_is_const = true;

    /*!
     * Bools for computing linear and rotational momentums of the body
     */
    bool d_calculate_structure_linear_mom = false, d_calculate_structure_rotational_mom = false;

    /*!
     * Iteration_counter for printing stuff.
     */
    int d_timestep_counter = 0, d_output_interval = 1;

    /*!
     * Bools for outputing stuff which is calculated on the fly.
     */
    bool d_print_output = false, d_output_drag = false, d_output_torque = false, d_output_power = false,
         d_output_trans_vel = false, d_output_rot_vel = false, d_output_COM_coordinates = false, d_output_MOI = false,
         d_output_eul_mom = false;

    /*!
     * output file name string.
     */
    std::string d_dir_name = "./ConstraintIBMethodDump", d_base_output_filename = "ImmersedStructure";

    /*!
     * Store LData for only those levels which contain immersed structures.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_l_data_U_interp, d_l_data_U_correction, d_l_data_U_new,
        d_l_data_U_current, d_l_data_U_half, d_l_data_X_half_Euler, d_l_data_X_new_MidPoint;

    /*!
     * Hierarchy operations object. Needed for projection step.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_no_fill_op;
    int d_wgt_cc_idx, d_wgt_sc_idx;
    double d_volume;

    /*!
     *  Variables and variable contexts associated with calculating divergence free projection.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_u_fluidSolve_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_u_var;

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_scratch_context;
    int d_u_fluidSolve_idx, d_u_fluidSolve_cib_idx, d_phi_idx, d_Div_u_scratch_idx;

    /*!
     * Variables associated with the spatially varying density field, which is maintained by an integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_rho_var;
    int d_rho_ins_idx = IBTK::invalid_index, d_rho_scratch_idx = IBTK::invalid_index;

    /*!
     * The following variables are needed to solve cell centered poison equation for \f$ \phi \f$ ,which is
     * used to project the corrected background fluid velocity on divergence free field to remove kinkiness
     * introduced via FuRMoRP algorithm.
     */
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM> d_velcorrection_projection_bc_coef;
    std::unique_ptr<SAMRAI::solv::PoissonSpecifications> d_velcorrection_projection_spec;
    SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> d_velcorrection_projection_op;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovPoissonSolver> d_velcorrection_projection_solver;
    SAMRAI::tbox::Pointer<IBTK::CCPoissonPointRelaxationFACOperator> d_velcorrection_projection_fac_op;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_velcorrection_projection_fac_pc_db;
    SAMRAI::tbox::Pointer<IBTK::FACPreconditioner> d_velcorrection_projection_fac_pc;

    /*!
     * File streams associated for the output.
     */
    std::vector<std::unique_ptr<std::ofstream> > d_trans_vel_stream, d_rot_vel_stream, d_drag_force_stream,
        d_moment_of_inertia_stream, d_torque_stream, d_position_COM_stream, d_power_spent_stream;

    /*!
     * Stream for calculating Eulerian momentum.
     */
    std::fstream d_eulerian_mom_stream;

    /*!
     * Pre and post fluid solve call back functions and contexts.
     *
     * TODO: Update these to use std::function.
     */
    std::vector<void (*)(const double, const double, const int, void*)> d_prefluidsolve_callback_fns,
        d_postfluidsolve_callback_fns;
    std::vector<void*> d_prefluidsolve_callback_fns_ctx, d_postfluidsolve_callback_fns_ctx;

    // Velocity boundary operator.
    IBTK::RobinPhysBdryPatchStrategy* d_u_phys_bdry_op = nullptr;
};
} // namespace IBAMR

/////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_ConstraintIBMethod
