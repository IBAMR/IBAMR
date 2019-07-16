// Filename: SetFluidProperties.h
// Created on March 6, 2018 by Nishant Nangia
//
// Copyright (c) 2002-2018, Amneet Bhalla and Boyce Griffith
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


/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetFluidProperties
#define included_SetFluidProperties

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>

#include <Variable.h>

namespace IBTK
{
class HierarchyMathOps;
class CartGridFunction;
} // namespace IBTK

/*!
 * Pre processing call back function to be hooked into IBAMR::INSVCStaggeredHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetFluidDensityCallbackFunction(int rho_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                         SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         const int cycle_num,
                                         const double time,
                                         const double current_time,
                                         const double new_time,
                                         void* ctx);

/*!
 * Pre processing call back function to be hooked into IBAMR::INSVCStaggeredHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetFluidViscosityCallbackFunction(int mu_idx,
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
     * solid Eulerian density based on the input file
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidProperties(const std::string& object_name,
                       SAMRAI::tbox::Pointer<IBTK::CartGridFunction> rho_fcn,
                       SAMRAI::tbox::Pointer<IBTK::CartGridFunction> mu_fcn);

    /*!
     * Destructor for this class.
     */
    ~SetFluidProperties();

    /*!
     * Set the density
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int cycle_num,
                             const double time,
                             const double current_time,
                             const double new_time);

    /*!
     * Set the viscosity
     */
    void setViscosityPatchData(int mu_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int cycle_num,
                               const double time,
                               const double current_time,
                               const double new_time);

    //////////////// PRIVATE /////////////////////////////

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
     * Pointer to density and viscosity functions
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_rho_fcn, d_mu_fcn;

}; // SetFluidProperties

#endif // #ifndef included_SetFluidProperties
