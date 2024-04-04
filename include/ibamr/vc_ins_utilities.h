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
namespace VcINSUtilities
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

class SetFluidProperties
{
    /*!
     * \brief Class SetFluidProperties is a utility class which sets the fluid and
     * solid Eulerian density based on the current level set information
     */
public:
    /*!
     * Constructor for this class.
     */
    SetFluidProperties(const std::string& object_name,
                       SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_var,
                       const double rho_liquid,
                       const double rho_gas,
                       const double mu_liquid,
                       const double mu_gas,
                       const int ls_reinit_interval,
                       const double num_interface_cells,
                       const std::string& num_phases);

    /*!
     * Constructor for this class.
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
                       const int ls_reinit_interval,
                       const double num_gas_interface_cells,
                       const double num_solid_interface_cells,
                       const bool set_mu_solid,
                       const std::string& num_phases);

    /*!
     * Destructor for this class.
     */
    ~SetFluidProperties();

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
     * Return the number of phases.
     */
    const std::string& getNumberOfPhases() const
    {
        return d_num_phases;
    } // getNumberOfPhases

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
     * Level set reinitialization interval.
     */
    int d_ls_reinit_interval;

    /*!
     * Number of interface cells over which to smooth the material properties.
     */
    double d_num_gas_interface_cells, d_num_solid_interface_cells;

    /*!
     * Whether or not to set viscosity in the solid region.
     */
    bool d_set_mu_solid;

    /*!
     * Number of phases. Valid options are "TWO_PHASE" and "THREE_PHASE".
     */
    const std::string d_num_phases = "TWO_PHASE";

}; // SetFluidProperties

class GravityForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Class constructor.
     */
    GravityForcing(const std::string& object_name,
                   SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> ins_hierarchy_integrator,
                   std::vector<double> grav_const,
                   std::string grav_type = "FULL",
                   SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hierarchy_integrator = nullptr,
                   SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_gas_var = nullptr,
                   SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr);

    /*!
     * \brief Empty destructor.
     */
    ~GravityForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete GravityForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

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
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        const double data_time,
                        const bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

private:
    GravityForcing();

    GravityForcing(const GravityForcing& from);

    GravityForcing& operator=(const GravityForcing& that);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to INSVC solver.
     */
    SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_hierarchy_integrator = nullptr;

    /*!
     * Vector to store the acceleration due to gravity.
     */
    std::vector<double> d_grav_const;

    /*!
     * String to specify the type of gravity force. Valid options are: "FULL" and "FLOW".
     */
    std::string d_grav_type = "FULL";

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hierarchy_integrator = nullptr;

    /*!
     * Level set variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_gas_var = nullptr;

    /*!
     * Density.
     */
    double d_rho_neg, d_rho_pos;

    /*!
     * Number of interface cells over which to smooth the material properties.
     */
    int d_num_gas_interface_cells;
};

} // namespace VcINSUtilities

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_vc_ins_utilities
