// Filename: INSStaggeredStochasticForcing.h
// Created on 02 Feb 2011 by Boyce Griffith
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

#ifndef included_INSStaggeredStochasticForcing
#define included_INSStaggeredStochasticForcing

#include <stddef.h>
#include <string>
#include <vector>

#include "CellVariable.h"
#include "EdgeVariable.h" // IWYU pragma: keep
#include "IntVector.h"
#include "NodeVariable.h" // IWYU pragma: keep
#include "PatchLevel.h"
#include "VariableContext.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "tbox/Array.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI
/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredStochasticForcing provides an interface for
 * specifying a stochastic forcing term for a staggered-grid incompressible
 * Navier-Stokes solver.
 */
class INSStaggeredStochasticForcing : public IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic stresses at the centers and nodes of the Cartesian
     * grid.
     */
    INSStaggeredStochasticForcing(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                  const INSStaggeredHierarchyIntegrator* fluid_solver);

    /*!
     * \brief Empty destructor.
     */
    ~INSStaggeredStochasticForcing();

    /*!
     * \name Methods to set patch data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete INSStaggeredStochasticForcing object is
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
    INSStaggeredStochasticForcing();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredStochasticForcing(const INSStaggeredStochasticForcing& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredStochasticForcing& operator=(const INSStaggeredStochasticForcing& that);

    /*!
     * Pointer to the fluid solver object that is using this stochastic force
     * generator.
     */
    const INSStaggeredHierarchyIntegrator* const d_fluid_solver;

    /*!
     * Type of stress tensor (correlated or uncorrelated).
     */
    StochasticStressTensorType d_stress_tensor_type;

    /*!
     * Weighting data.
     */
    double d_std;
    int d_num_rand_vals;
    std::vector<SAMRAI::tbox::Array<double> > d_weights;

    /*!
     * Boundary condition scalings.
     */
    double d_velocity_bc_scaling, d_traction_bc_scaling;

    /*!
     * VariableContext and Variable objects for storing the components of the
     * stochastic stresses.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_W_cc_var;
    int d_W_cc_idx;
    std::vector<int> d_W_cc_idxs;
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_W_nc_var;
    int d_W_nc_idx;
    std::vector<int> d_W_nc_idxs;
#endif
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_W_ec_var;
    int d_W_ec_idx;
    std::vector<int> d_W_ec_idxs;
#endif
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredStochasticForcing
