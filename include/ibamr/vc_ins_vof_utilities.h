// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_vc_ins_vof_utilities
#define included_IBAMR_vc_ins_vof_utilities

#include <ibamr/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/CartGridFunction.h>

#include <tbox/Pointer.h>

namespace IBTK
{
class HierarchyMathOps;
}

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////
namespace IBAMR
{

namespace VCINSVOFUtilities
{
/*!
 * \brief VOF initial condition derived from an existing level set initial condition.
 *
 * The goal is to keep examples general: they specify only a level set initial condition
 * (geometry/physics), and the VOF field is generated implicitly from that same object.
 *
 * Implementation strategy:
 *   1) allocate a scratch cell-centered patch data index,
 *   2) call the provided LS CartGridFunction to write \phi into that scratch storage, and
 *   3) compute the cell volume fraction using the same vof_alpha() reconstruction used
 *      by the runtime LS->VOF update callback.
 */
class VOFInitialConditionFromLevelSet : public IBTK::CartGridFunction
{
public:
    VOFInitialConditionFromLevelSet(const std::string& object_name,
                                    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> ls_ic);

    ~VOFInitialConditionFromLevelSet() override = default;

    bool isTimeDependent() const override;

    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        const double data_time,
                        const bool initial_time,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level) override;

private:
    std::string d_object_name;

    // Keep the LS initial condition alive for the duration of the run.
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_ls_ic;

    // Scratch storage used to ask the LS IC to write \phi values on each patch.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_phi_scratch_var;
    int d_phi_scratch_idx = IBTK::invalid_index;
};

/*!\brief Helper that recomputes VOF from the transported level set field each time step.
 *
 * This is the "VOF initialization class" you were using as an advection-diffusion
 * integrate callback in the example. It lives here so the example stays clean.
 */
class VOFFromLevelSetInitializer
{
public:
    VOFFromLevelSetInitializer(const std::string& object_name,
                               SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_integrator,
                               SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> ls_var,
                               SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vof_var);

    ~VOFFromLevelSetInitializer() = default;

    /*!\brief Register the integrateHierarchy callback on the provided integrator. */
    void registerIntegrateHierarchyCallback();

    /*!\brief Manually compute VOF from LS at a given time.
     *
     * This is useful right after hierarchy initialization (t0) if you want VOF
     * populated before the first time step.
     */
    void computeVOFFromLevelSet(const double time, const bool use_new_context);

    /*!\brief Static trampoline with the signature expected by AdvDiffHierarchyIntegrator. */
    static void integrateHierarchyCallback(double current_time, double new_time, int cycle_num, void* ctx);

private:
    void computeVOFInternal(const double time, const int ls_idx, const int vof_idx);

    std::string d_object_name;
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_integrator;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_ls_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_vof_var;
};

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetVOFBasedDensity(int rho_idx,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> rho_var,
                            SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                            const int cycle_num,
                            const double time,
                            const double current_time,
                            const double new_time,
                            void* ctx);

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetVOFBasedViscosity(int mu_idx,
                              SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> mu_var,
                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                              const int cycle_num,
                              const double time,
                              const double current_time,
                              const double new_time,
                              void* ctx);

/*!
 * \brief Class SetVOFBasedFluidProperties is a utility class which sets the fluid
 * Eulerian density and viscosity based on the current level set and VOF information.
 */
class SetVOFBasedFluidProperties
{
public:
    /*!
     * Constructor for this class.
     */
    SetVOFBasedFluidProperties(const std::string& object_name,
                               SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                               SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vof_var,
                               const double rho_liquid,
                               const double rho_gas,
                               const double mu_liquid,
                               const double mu_gas);

    // define multiple for liquid, solid, gas and other things

    /*!
     * Destructor for this class.
     */
    ~SetVOFBasedFluidProperties() = default;

    /*!
     * Set the density based on the current level set information.
     */
    inline void setVOFBasedDensityPatchData(int rho_idx,
                                            SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> rho_var,
                                            SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                            const int cycle_num,
                                            const double time,
                                            const double current_time,
                                            const double new_time)
    {
        setVOFBasedDensityPatchData2PhaseFlows(
            rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

        return;
    } // setVOFBasedDensityPatchData

    /*!
     * Set the viscosity based on the current level set information.
     */
    inline void setVOFBasedViscosityPatchData(int mu_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> mu_var,
                                              SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                              const int cycle_num,
                                              const double time,
                                              const double current_time,
                                              const double new_time)
    {
        setVOFBasedViscosityPatchData2PhaseFlows(
            mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

        return;
    } // setVOFBasedViscosityPatchData

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetVOFBasedFluidProperties();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetVOFBasedFluidProperties& operator=(const SetVOFBasedFluidProperties& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetVOFBasedFluidProperties(const SetVOFBasedFluidProperties& from);

    /*!
     * Set the density based on the current level set information for two-phase flows.
     */
    void setVOFBasedDensityPatchData2PhaseFlows(int rho_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> rho_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                const int cycle_num,
                                                const double time,
                                                const double current_time,
                                                const double new_time);

    /*!
     * Set the viscosity based on the current level set information for two-phase flows.
     */
    void setVOFBasedViscosityPatchData2PhaseFlows(int mu_idx,
                                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> mu_var,
                                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                  const int cycle_num,
                                                  const double time,
                                                  const double current_time,
                                                  const double new_time);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set and VOF variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_vof_var;

    /*!
     * Density.
     */
    double d_rho_liquid, d_rho_gas;

    /*!
     * Viscosity.
     */
    double d_mu_liquid, d_mu_gas;

    /*!
     * Number of phases. Valid options are 2 and 3.
     */
    int d_num_phases = 2;

}; // SetFluidProperties

class VOFBasedGravityForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor for this class. Applies the gravitational force throughout the
     * computational domain with the density retrieved from the provided <code>ins_hierarchy_integrator<\code>.
     *
     * grav_const stores the acceleration due to gravity.
     *
     */
    VOFBasedGravityForcing(const std::string& object_name,
                           SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                           std::vector<double> grav_const);

    /*!
     * \brief Constructor for this class. Applies the gravitational force using the
     * flow density field, which is computed from the fluid VOF function.
     *
     * @param input_db provides parameters such as rho_neg, and rho_pos.
     * grav_const stores the acceleration due to gravity.
     *
     */
    VOFBasedGravityForcing(const std::string& object_name,
                           SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> vof_var,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           std::vector<double> grav_const);

    /*!
     * \brief Empty destructor.
     */
    ~VOFBasedGravityForcing() = default;

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete GravityForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>(nullptr)) override;

    //\}

private:
    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to INSVC solver.
     */
    SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_hierarchy_integrator;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hierarchy_integrator;

    /*!
     * VOF variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double>> d_vof_var;

    /*!
     * Vector to store the acceleration due to gravity.
     */
    const std::vector<double> d_grav_const;

    /*!
     * String to specify the type of gravity force. Valid options are: "FULL" and "FLOW".
     * "FULL" - compute volumetric gravitational source term \f$ \rho g \f$.
     * "FLOW" - sets the \f$ \rho^\text{flow} g \f$ where \f$ \rho^\text{flow} \f$ is the density of the flow phase
     * (excludes solid).
     */
    std::string d_grav_type;

    /*!
     * rho_neg - density where the fluid level set (vof) takes a negative (zero) value.
     * rho_pos - density where the fluid level set (vof) takes a positive (one) value.
     */
    double d_rho_neg, d_rho_pos;

    /*!
     * VOF scratch data.
     */
    int d_vof_scratch_idx = IBTK::invalid_index;
};

} // namespace VCINSVOFUtilities

} // namespace IBAMR

#endif // #ifndef included_IBAMR_vc_ins_vof_utilities
