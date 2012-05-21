// Filename: AdvDiffStochasticForcing.h
// Created on 29 Apr 2011 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

// IBAMR INCLUDES
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

// IBTK THIRD-PARTY INCLUDES
#include <ibtk/muParser.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <SideVariable.h>
#include <VariableContext.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffStochasticForcing provides an interface for specifying a
 * stochastic forcing term for a concentration that is coupled to the
 * staggered-grid incompressible Navier-Stokes solver.
 */
class AdvDiffStochasticForcing
    : public virtual IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic fluxes at the faces of the Cartesian grid.
     *
     * The optional parameter \p f_expression specifies the functional form of a
     * concentration-dependent scaling factor.  The optional parameter \p g may
     * be used to setup a pseudo concentration gradient for periodic problems.
     */
    AdvDiffStochasticForcing(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > C_var,
        const INSStaggeredHierarchyIntegrator* const fluid_solver,
        const AdvDiffSourceTermType& eval_type,
        const std::string& f_expression="1.0",
        const SAMRAI::tbox::Array<double>& g=SAMRAI::tbox::Array<double>());

    /*!
     * \brief Empty virtual destructor.
     */
    virtual
    ~AdvDiffStochasticForcing();

    /*!
     * Set an object to provide boundary conditions for a scalar-valued
     * concentration.
     */
    void
    setPhysicalBcCoefs(
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* C_bc_coef);

    /*!
     * Set objects to provide boundary conditions for a vector-valued
     * concentration.
     */
    void
    setPhysicalBcCoefs(
        std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> C_bc_coef);

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete AdvDiffStochasticForcing object is
     * time-dependent.
     */
    virtual bool
    isTimeDependent() const;

    /*!
     * \brief Set the diffusion coefficient.
     */
    void
    setDiffusionCoefficient(
        const double kappa);

    /*!
     * \brief Set the timestep size.
     */
    void
    setDt(
        const double dt);

    /*!
     * \brief Set the scale factor std.
     */
    void
    setStd(
        const double std);

    /*!
     * \brief Set a vector of integer values that specify whether to evaluate
     * the stochastic fluxes for a particular Runge-Kutta cycle.
     *
     * For each cycle k, if regen_rand_cycle[k] == 1, then a new set of random
     * values is generated; if regen_rand_cycle[k] == 0, then the random values
     * from the previous cycle are used; and if regen_rand_cycle[k] == -1, then
     * the stochastic forcing term is set to zero.
     */
    void
    setRegenRandCycle(
        const SAMRAI::tbox::Array<int>& regen_rand_cycle);

    /*!
     * \brief Set the value of f(C), concentration-dependent flux scaling
     * function.
     */
    void
    setFExpression(
        const std::string& f_expression);

    /*!
     * \brief Set the value of g, the pseudo concentration gradient vector
     * coefficient.
     */
    void
    setG(
        const SAMRAI::tbox::Array<double>& g);

    /*!
     * \brief Set the manner in which the source term should be evaluated.
     */
    void
    setAdvDiffSourceTermType(
        const AdvDiffSourceTermType& eval_type);

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatchLevel().
     *
     * \see setDataOnPatch
     */
    virtual void
    setDataOnPatchHierarchy(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const double data_time,
        const bool initial_time=false,
        const int coarsest_ln=-1,
        const int finest_ln=-1);

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * level of the patch hierarchy using the virtual function setDataOnPatch().
     *
     * \note This function also allocates data for storing the stochastic stress
     * components, evaluates those components, and synchronizes the values of
     * those components before calling setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    virtual void
    setDataOnPatchLevel(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        const double data_time,
        const bool initial_time=false);

    /*!
     * \brief Pure virtual function to evaluate the function on the patch
     * interior.
     */
    virtual void
    setDataOnPatch(
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        const double data_time,
        const bool initial_time=false,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level=SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

protected:
    /*
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
    AdvDiffStochasticForcing(
        const AdvDiffStochasticForcing& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffStochasticForcing&
    operator=(
        const AdvDiffStochasticForcing& that);

    /*!
     * Pointer to the concentration variable registered with the class.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_C_var;

    /*!
     * Boundary condition specification objects.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_C_bc_coef;

    /*!
     * Concentration-dependent scaling factor function.
     */
    mu::Parser d_f_parser;

    /*!
     * Pesudo-concentration gradient coefficient.
     */
    SAMRAI::tbox::Array<double> d_g;

    /*!
     * Pointer to the fluid solver object that is using this stochastic force
     * generator.
     */
    const INSStaggeredHierarchyIntegrator* const d_fluid_solver;

    /*!
     * Forcing type (midpoint rule or trapezoidal rule).
     */
    AdvDiffSourceTermType d_eval_type;

    /*!
     * VariableContext and Variable objects for storing the components of the
     * stochastic stresses.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_W_sc_var;
    int d_W_sc_idx, d_C_cc_idx;
    double d_kappa, d_dt, d_std;
    SAMRAI::tbox::Array<int> d_regen_rand_cycle;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/AdvDiffStochasticForcing.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffStochasticForcing
