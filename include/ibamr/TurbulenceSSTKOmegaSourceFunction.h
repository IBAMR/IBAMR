// Filename: TurbulenceSSTKOmegaSourceFunction.h
// Created on 23 Sep 2019 by Ramakrishnan Thirumalaisamy
//
// Copyright (c) 2002-2017, Amneet Bhalla
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

#ifndef included_IBAMR_TurbulenceSSTKOmegaSourceFunction
#define included_IBAMR_TurbulenceSSTKOmegaSourceFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

#include <string>

// IBAMR INCLUDES
#include <ibamr/app_namespaces.h>

namespace IBAMR
{
class INSVCStaggeredHierarchyIntegrator;
class TwoEquationTurbulenceHierarchyIntegrator;
} // namespace IBAMR

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class SideData;
template <int DIM, class TYPE>
class CellData;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class TurbulenceSSTKOmegaSourceFunction describes the implementation of source terms in k and omega equations
 *
 *
 */
class TurbulenceSSTKOmegaSourceFunction : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    TurbulenceSSTKOmegaSourceFunction(const std::string& object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                      TwoEquationTurbulenceHierarchyIntegrator* turb_kw_integrator,
                                      INSVCStaggeredHierarchyIntegrator* ins_hier_integrator);

    /*!
     * \brief Destructor.
     */
    ~TurbulenceSSTKOmegaSourceFunction() = default;

    /*!
     * \note This concrete IBTK::CartGridFunction is time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     *
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    void setDataOnPatchHierarchy(int data_idx,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                 double data_time,
                                 bool initial_time = false,
                                 int coarsest_ln = -1,
                                 int finest_ln = -1) override;
    /*!
     * Set the data on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

    //\}
private:
    /*!
     * Get the F1 variable index registered with the TwoEquationTurbulenceHierarchyIntegrator hierarchy
     * integrator.
     */
    int getBlendingFunctionVariableIndex(
        SAMRAI::tbox::Pointer<TwoEquationTurbulenceHierarchyIntegrator> turb_hier_integrator);

    /*!
     * Get the production variable index registered with the TwoEquationTurbulenceHierarchyIntegrator
     * hierarchy integrator.
     */

    int
    getProductionVariableIndex(SAMRAI::tbox::Pointer<TwoEquationTurbulenceHierarchyIntegrator> turb_hier_integrator);

    /* set the data on patch cells for k*/
    void setDataOnPatchCellForK(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > k_f_data,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                const double data_time,
                                const bool initial_time,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level);

    /* set the data on patch cells for omega*/
    void setDataOnPatchCellForOmega(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > w_f_data,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                    const double data_time,
                                    const bool initial_time,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level);

    std::string d_object_name;
    SAMRAI::tbox::Pointer<TwoEquationTurbulenceHierarchyIntegrator> d_turb_hier_integrator;
    SAMRAI::tbox::Pointer<INSVCStaggeredHierarchyIntegrator> d_ins_hier_integrator;

    IBTK::Vector2d d_gravity;

    int d_mu_t_new_idx, d_rho_new_idx;
    int d_rho_scratch_idx;
    int d_k_new_idx, d_k_scratch_idx, d_w_new_idx, d_w_scratch_idx, d_f1_scratch_idx, d_p_scratch_idx;
};
} // namespace IBAMR
#endif //#ifndef included_IBAMR_SurfaceTensionForceFunction
