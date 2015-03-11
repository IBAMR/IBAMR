// Filename: AdvDiffStochasticForcing.h
// Created on 29 Apr 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#ifndef included_AdvDiffStochasticForcing
#define included_AdvDiffStochasticForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "ibtk/CartGridFunction.h"
#include "muParser.h"
#include "SAMRAI/tbox/Array.h"

namespace IBAMR
{
class AdvDiffSemiImplicitHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
class Patch;
class PatchHierarchy;
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffStochasticForcing provides an interface for specifying a
 * stochastic forcing term for cell-centered advection-diffusion solver solver.
 */
class AdvDiffStochasticForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic fluxes at the faces of the Cartesian grid.
     */
    AdvDiffStochasticForcing(const std::string& object_name,
                             boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > C_var,
                             const AdvDiffSemiImplicitHierarchyIntegrator* adv_diff_solver);

    /*!
     * \brief Empty destructor.
     */
    ~AdvDiffStochasticForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete AdvDiffStochasticForcing object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy.
     */
    void setDataOnPatchHierarchy(const int data_idx,
                                 boost::shared_ptr<SAMRAI::hier::Variable> var,
                                 boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
                                 const double data_time,
                                 const bool initial_time = false,
                                 const int coarsest_ln = -1,
                                 const int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interior.
     */
    void setDataOnPatch(const int data_idx,
                        boost::shared_ptr<SAMRAI::hier::Variable> var,
                        boost::shared_ptr<SAMRAI::hier::Patch> patch,
                        const double data_time,
                        const bool initial_time = false,
                        boost::shared_ptr<SAMRAI::hier::PatchLevel> patch_level = NULL);

    //\}

protected:
    /*!
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffStochasticForcing();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffStochasticForcing(const AdvDiffStochasticForcing& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffStochasticForcing& operator=(const AdvDiffStochasticForcing& that);

    /*!
     * boost::shared_ptr to the concentration variable associated with this source term
     * generator.
     */
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > d_C_var;

    /*!
     * Concentration-dependent flux scaling function.
     */
    mu::Parser d_f_parser;

    /*!
     * boost::shared_ptr to the advection-diffusion solver object that is using this
     * stochastic source term generator.
     */
    const AdvDiffSemiImplicitHierarchyIntegrator* const d_adv_diff_solver;

    /*!
     * Weighting data.
     */
    double d_std;
    int d_num_rand_vals;
    std::vector<std::vector<double> > d_weights;

    /*!
     * Boundary condition scalings.
     */
    double d_dirichlet_bc_scaling, d_neumann_bc_scaling;

    /*!
     * VariableContext and Variable objects for storing the components of the
     * stochastic fluxes.
     */
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > d_C_cc_var;
    int d_C_current_cc_idx, d_C_half_cc_idx, d_C_new_cc_idx;
    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > d_F_sc_var;
    int d_F_sc_idx;
    std::vector<int> d_F_sc_idxs;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffStochasticForcing
