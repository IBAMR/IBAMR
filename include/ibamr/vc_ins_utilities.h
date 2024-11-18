// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_vc_ins_utilities
#define included_IBAMR_vc_ins_utilities

#include <ibtk/config.h>

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/ibamr_enums.h>

#include "ibtk/CartGridFunction.h"

#include "tbox/Pointer.h"

namespace IBTK
{
class HierarchyMathOps;
}

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////
namespace IBAMR
{
/*!
 * \brief The vc_ins_utilities class can be utilized to set fluid properties such as density and viscosity
 * for both two-phase and three-phase flows throughout the entire domain. Additionally, this class provides
 * implementations for gravity force calculations in two-phase and three-phase flows.
 *
 * \note Various options are available for computing the side-centered density within this class.
 */
namespace VCINSUtilities
{

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetDensityCallbackFunction(int rho_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
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

void callSetViscosityCallbackFunction(int mu_idx,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                      SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                      const int cycle_num,
                                      const double time,
                                      const double current_time,
                                      const double new_time,
                                      void* ctx);

/*!
 * \brief Class SetFluidProperties is a utility class which sets the fluid and
 * solid Eulerian density based on the current level set information.
 */
class SetFluidProperties
{
public:
    /*!
     * Constructor for this class.
     *
     * num_interface_cells - number of cells over which the Heaviside function is smoothed on
     * either side of the interface.
     *
     */
    SetFluidProperties(const std::string& object_name,
                       SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                       const double rho_liquid,
                       const double rho_gas,
                       const double mu_liquid,
                       const double mu_gas,
                       const double num_interface_cells);

    /*!
     * Constructor for this class.
     *
     * num_gas_interface_cells - number of cells over which the Heaviside function is smoothed on
     * either side of the gas interface.
     *
     * num_solid_interface_cells - number of cells over which the Heaviside function is smoothed on
     * either side of the solid interface.
     *
     * set_mu_solid - If it is true, then visocity is set in the entire domain. If it is false, the viscosity is set
     * in the fluid and gas phases but not the solid.
     *
     */
    SetFluidProperties(const std::string& object_name,
                       SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                       const double rho_liquid,
                       const double rho_gas,
                       const double rho_solid,
                       const double mu_fluid,
                       const double mu_gas,
                       const double mu_solid,
                       const double num_gas_interface_cells,
                       const double num_solid_interface_cells,
                       const bool set_mu_solid);

    /*!
     * Destructor for this class.
     */
    ~SetFluidProperties() = default;

    /*!
     * Set the density based on the current level set information.
     */
    inline void setDensityPatchData(int rho_idx,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int cycle_num,
                                    const double time,
                                    const double current_time,
                                    const double new_time)
    {
        if (d_num_phases == 2)
        {
            setDensityPatchData2PhaseFlows(rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);
        }
        else if (d_num_phases == 3)
        {
            setDensityPatchData3PhaseFlows(rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);
        }

        return;
    } // setDensityPatchData

    /*!
     * Set the viscosity based on the current level set information.
     */
    inline void setViscosityPatchData(int mu_idx,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                      SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                      const int cycle_num,
                                      const double time,
                                      const double current_time,
                                      const double new_time)
    {
        if (d_num_phases == 2)
        {
            setViscosityPatchData2PhaseFlows(mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);
        }
        else if (d_num_phases == 3)
        {
            setViscosityPatchData3PhaseFlows(mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);
        }

        return;
    } // setViscosityPatchData

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetFluidProperties();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetFluidProperties& operator=(const SetFluidProperties& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetFluidProperties(const SetFluidProperties& from);

    /*!
     * Set the density based on the current level set information for two-phase flows.
     */
    void setDensityPatchData2PhaseFlows(int rho_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        const int cycle_num,
                                        const double time,
                                        const double current_time,
                                        const double new_time);

    /*!
     * Set the density based on the current level set information for three-phase flows.
     */
    void setDensityPatchData3PhaseFlows(int rho_idx,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                        const int cycle_num,
                                        const double time,
                                        const double current_time,
                                        const double new_time);

    /*!
     * Set the viscosity based on the current level set information for two-phase flows.
     */
    void setViscosityPatchData2PhaseFlows(int mu_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                          const int cycle_num,
                                          const double time,
                                          const double current_time,
                                          const double new_time);

    /*!
     * Set the viscosity based on the current level set information for three-phase flows.
     */
    void setViscosityPatchData3PhaseFlows(int mu_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
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
     * Level set variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var, d_ls_solid_var;

    /*!
     * Density.
     */
    double d_rho_liquid, d_rho_gas, d_rho_solid;

    /*!
     * Viscosity.
     */
    double d_mu_liquid, d_mu_gas, d_mu_solid;

    /*!
     * Number of interface cells over which to smooth the material properties.
     */
    double d_num_gas_interface_cells, d_num_solid_interface_cells;

    /*!
     * Whether or not to set viscosity in the solid region.
     */
    bool d_set_mu_solid;

    /*!
     * Number of phases. Valid options are 2 and 3.
     */
    int d_num_phases = 2;

}; // SetFluidProperties

/*!
 * \brief The GravityForcing class provides an implementation of gravity force.
 * This class can be utilized to apply the gravitational force \f$ \rho g \f$ using the density field,
 * which includes all three phases: liquid, gas, and solid; or using the flow density field, which includes
 * only liquid and gas phases and excludes the solid phase.
 */
class GravityForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor for this class. Applies the gravitational force throughout the
     * computational domain with the density retrieved from the provided <code>ins_hierarchy_integrator<\code>.
     *
     * grav_const stores the acceleration due to gravity.
     *
     */
    GravityForcing(const std::string& object_name,
                   SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                   std::vector<double> grav_const);

    /*!
     * \brief Constructor for this class. Applies the gravitational force using the
     * flow density field, which is computed from the fluid level set function.
     *
     * @param input_db provides parameters such as rho_neg, rho_pos, and num_gas_interface_cells.
     * grav_const stores the acceleration due to gravity.
     *
     */
    GravityForcing(const std::string& object_name,
                   SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator,
                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                   std::vector<double> grav_const);

    /*!
     * \brief Empty destructor.
     */
    ~GravityForcing() = default;

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
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1) override;

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(nullptr)) override;

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
     * Level set variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var;

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
     * rho_neg - density where the fluid level set takes a negative value.
     * rho_pos - density where the fluid level set takes a positive value.
     */
    double d_rho_neg, d_rho_pos;

    /*!
     * Number of interface cells over which to smooth the material properties.
     */
    int d_num_gas_interface_cells;

    /*!
     * Level set scratch data.
     */
    int d_ls_gas_scratch_idx = IBTK::invalid_index;
};

} // namespace VCINSUtilities

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_vc_ins_utilities
