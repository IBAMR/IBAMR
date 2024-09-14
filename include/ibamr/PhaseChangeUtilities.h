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

#ifndef included_IBAMR_PhaseChangeUtilities
#define included_IBAMR_PhaseChangeUtilities

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include "tbox/Pointer.h"

namespace IBTK
{
class HierarchyMathOps;
}

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////
namespace IBAMR
{
/*!
 * \brief The PhaseChangeUtilities class can be utilized to set fluid properties such as density and viscosity
 * for both two-phase and three-phase flows throughout the entire domain.
 *
 * \note Various options are available for computing the side-centered density within this class.
 */
namespace PhaseChangeUtilities
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
                                    int cycle_num,
                                    double time,
                                    double current_time,
                                    double new_time,
                                    void* ctx);

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param kappa_idx a patch data index for the current thermal conductivity variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetThermalConductivityCallbackFunction(int kappa_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > kappa_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                int cycle_num,
                                                double time,
                                                double current_time,
                                                double new_time,
                                                void* ctx);

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param specific_heat_idx a patch data index for the current specific heat variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetSpecificHeatCallbackFunction(int specific_heat_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > specific_heat_var,
                                         SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         int cycle_num,
                                         double time,
                                         double current_time,
                                         double new_time,
                                         void* ctx);

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param mu_idx a patch data index for the current viscosity variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetViscosityCallbackFunction(int mu_idx,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                      SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                      int cycle_num,
                                      double time,
                                      double current_time,
                                      double new_time,
                                      void* ctx);

class SetFluidProperties
{
public:
    /*!
     * Constructor for this class.
     */
    SetFluidProperties(const std::string& object_name,
                       SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var,
                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* H_bc_coef,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lf_bc_coef,
                       double rho_liquid,
                       double rho_solid,
                       double rho_gas,
                       double kappa_liquid,
                       double kappa_solid,
                       double kappa_gas,
                       double Cp_liquid,
                       double Cp_solid,
                       double Cp_gas,
                       double mu_liquid,
                       double mu_solid,
                       double mu_gas);

    /*!
     * Constructor for this class.
     *
     * This constructor can be used when only the energy equation is solved.
     */
    SetFluidProperties(const std::string& object_name,
                       SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > H_var,
                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* H_bc_coef,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lf_bc_coef,
                       double rho_liquid,
                       double rho_solid,
                       double rho_gas,
                       double kappa_liquid,
                       double kappa_solid,
                       double kappa_gas,
                       double Cp_liquid,
                       double Cp_solid,
                       double Cp_gas);

    /*!
     * Destructor for this class.
     */
    ~SetFluidProperties() = default;

    /*!
     * Set the density based on the current Heaviside and liquid fraction information.
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             int cycle_num,
                             double time,
                             double current_time,
                             double new_time);

    /*!
     * Set the thermal conductivity based on the current Heaviside and liquid fraction information.
     */
    void setThermalConductivityPatchData(int kappa_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > kappa_var,
                                         SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         int cycle_num,
                                         double time,
                                         double current_time,
                                         double new_time);

    /*!
     * Set the specific heat based on the current Heaviside and liquid fraction information.
     */
    void setSpecificHeatPatchData(int specific_heat_idx,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > specific_heat_var,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                  int cycle_num,
                                  double time,
                                  double current_time,
                                  double new_time);

    /*!
     * Set the viscosity based on the current Heaviside and liquid fraction information.
     */
    void setViscosityPatchData(int mu_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               int cycle_num,
                               double time,
                               double current_time,
                               double new_time);

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
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Heaviside variable and its boundary condition pointer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_H_var;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_H_bc_coef = nullptr;

    /*!
     * Liquid fraction variable and its boundary condition pointer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_lf_bc_coef = nullptr;

    /*!
     * Density.
     */
    double d_rho_liquid, d_rho_solid, d_rho_gas;

    /*!
     * Thermal conductivity.
     */
    double d_kappa_liquid, d_kappa_solid, d_kappa_gas;

    /*!
     * Specific heat.
     */
    double d_specific_heat_liquid, d_specific_heat_solid, d_specific_heat_gas;

    /*!
     * Viscosity.
     */
    double d_mu_liquid, d_mu_solid, d_mu_gas;

}; // SetFluidProperties

/*!
 * \brief A lightweight class to tag grid cells based on the liquid fraction gradient value for grid refinement.
 */
class TagLiquidFractionRefinementCells
{
public:
    /*!
     * \brief Constructor of the class.
     */
    TagLiquidFractionRefinementCells(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_grad_var,
                                     double tag_min_value = 0.0,
                                     double tag_max_value = 0.0)
        : d_adv_diff_solver(adv_diff_integrator),
          d_lf_var(lf_var),
          d_lf_grad_var(lf_grad_var),
          d_tag_min_value(tag_min_value),
          d_tag_max_value(tag_max_value)
    {
        // intentionally left blank
    } // TagLevelSetRefinementCells

    /*!
     * Tag the liquid-solid interface cells based on the liquid fraction and/or liquid fraction gradient value.
     */
    void tagLiquidFractionCells(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                int level_number,
                                double error_data_time,
                                int tag_index,
                                bool initial_time,
                                bool uses_richardson_extrapolation_too,
                                void* ctx);

private:
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_grad_var;
    double d_tag_min_value = 0.0;
    double d_tag_max_value = 0.0;
}; // TagLiquidFractionRefinementCells

/*!
 * \brief Preprocessing call back function to be hooked into IBAMR::HierarchyIntegrator class
 * to tag the cells for grid refinement based on the given tagging criteria.
 *
 * This static member should be registered with an appropriate hierarchy integrator
 * via registerApplyGradientDetectorCallback().
 *
 * \param ctx is the pointer to the TagLiquidFractionRefinementCells class object.
 */
void
callTagLiquidFractionCellsCallbackFunction(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int level_number,
                                           double error_data_time,
                                           int tag_index,
                                           bool initial_time,
                                           bool uses_richardson_extrapolation_too,
                                           void* ctx);

} // namespace PhaseChangeUtilities

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_PhaseChangeUtilities
