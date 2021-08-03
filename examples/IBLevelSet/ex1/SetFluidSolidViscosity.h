// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
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

#ifndef included_SetFluidSolidViscosity
#define included_SetFluidSolidViscosity

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <Variable.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back function to be hooked into IBAMR::VCINSStaggeredHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidSolidViscosity class object.
 *
 * \TODO: Let's move this out of the global namespace and use "snake case" for static function names.
 */
void callSetFluidSolidViscosityCallbackFunction(int mu_idx,
                                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                                SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                                const int cycle_num,
                                                const double time,
                                                const double current_time,
                                                const double new_time,
                                                void* ctx);

/*!
 * \brief Class SetFluidSolidViscosity is a utility class which sets the fluid and
 * solid Eulerian density based on the current level set information
 */
class SetFluidSolidViscosity
{
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidSolidViscosity(const std::string& object_name,
                           SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> adv_diff_solver,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > ls_solid_var,
                           const double mu_fluid,
                           const double mu_solid,
                           const int ls_reinit_interval,
                           const double num_solid_interface_cells,
                           const bool set_mu_solid);

    /*!
     * Destructor for this class.
     */
    ~SetFluidSolidViscosity() = default;

    /*!
     * Set the density based on the current level set information
     */
    void setViscosityPatchData(int mu_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int cycle_num,
                               const double time,
                               const double current_time,
                               const double new_time);

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetFluidSolidViscosity() = delete;

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetFluidSolidViscosity& operator=(const SetFluidSolidViscosity& that) = delete;

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetFluidSolidViscosity(const SetFluidSolidViscosity& from) = delete;

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver.
     */
    SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Level set variables
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_ls_solid_var, d_ls_gas_var;

    /*!
     * Density of the fluid and solid.
     */
    double d_mu_fluid, d_mu_solid;

    /*!
     * Number of cells over which to transition between values
     */
    double d_num_solid_interface_cells;

    /*!
     * Whether or not to set viscosity in the solid region
     */
    bool d_set_mu_solid;

}; // SetFluidSolidViscosity

#endif // #ifndef included_SetFluidSolidViscosity
