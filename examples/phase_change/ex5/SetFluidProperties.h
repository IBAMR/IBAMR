// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetFluidProperties
#define included_SetFluidProperties

///////////////////////////// INCLUDES ///////////////////////////////////
#include <ibamr/AdvDiffHierarchyIntegrator.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back function to be hooked into
 * IBAMR::VCINSStaggeredHierarchyIntegratorclass.
 *
 * \param rho_idx a patch data index for the current density variable maintained
 * by the integrator. \param ctx is the pointer to SetFluidProperties class
 * object.
 */

void callSetLiquidSolidDensityCallbackFunction(int rho_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                               const int cycle_num,
                                               const double time,
                                               const double current_time,
                                               const double new_time,
                                               void* ctx);

/*!
 * Pre processing call back function to be hooked into
 * IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param Cp_idx a patch data index for the current Cp variable maintained by
 * the integrator. \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetLiquidSolidSpecificHeatCallbackFunction(int Cp_idx,
                                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > Cp_var,
                                                    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                    const int cycle_num,
                                                    const double time,
                                                    const double current_time,
                                                    const double new_time,
                                                    void* ctx);

/*!
 * Pre processing call back function to be hooked into
 * IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param k_idx a patch data index for the current k variable maintained by
 * the integrator. \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetLiquidSolidConductivityCallbackFunction(int D_idx,
                                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > D_var,
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
     * solid density.
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidProperties(const std::string& object_name,
                       const SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                       const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > lf_var,
                       SAMRAI::solv::RobinBcCoefStrategy<NDIM>* lf_bc_coef,
                       const double rho_liquid,
                       const double rho_solid,
                       const double kappa_liquid,
                       const double kappa_solid,
                       const double Cp_liquid,
                       const double Cp_solid);

    /*!
     * Destructor for this class.
     */
    ~SetFluidProperties();

    /*!
     * Set the diffusion coefficient based on the current level set information
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int cycle_num,
                             const double time,
                             const double current_time,
                             const double new_time);

    /*!
     * Set the diffusion coefficient based on the current level set information
     */
    void setDiffusionCoefficientPatchData(int D_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > D_var,
                                          SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                          const int cycle_num,
                                          const double time,
                                          const double current_time,
                                          const double new_time);

    /*!
     * Set the specific heat based on the current level set information
     */
    void setSpecificHeatPatchData(int Cp_idx,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > Cp_var,
                                  SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                  const int cycle_num,
                                  const double time,
                                  const double current_time,
                                  const double new_time);

private:
    /*!
     * Deleted default constructor.
     */
    SetFluidProperties() = delete;

    /*!
     * Deleted copy constructor.
     */
    SetFluidProperties(const SetFluidProperties& from) = delete;

    /*!
     * Deleted assignment operator.
     */
    SetFluidProperties& operator=(const SetFluidProperties& that) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Liquid fraction variable and its boundary condition pointer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_lf_var;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_lf_bc_coef = nullptr;

    /*!
     * Density.
     */
    double d_rho_liquid, d_rho_solid;

    /*!
     * Diffusion coefficient.
     */
    double d_kappa_liquid, d_kappa_solid;

    /*!
     * Specific heat.
     */
    double d_Cp_liquid, d_Cp_solid;

    /*!
     * Level set reinitialization interval
     */
    int d_ls_reinit_interval;

    /*!
     * Number of interface cells over which to smooth the material properties
     */
    double d_num_interface_cells;
}; // SetFluidProperties

#endif // #ifndef included_SetFluidProperties