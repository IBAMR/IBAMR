// Filename: INSStaggeredStochasticForcing.h
// Created on 02 Feb 2011 by Boyce Griffith
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

#ifndef included_INSStaggeredStochasticForcing
#define included_INSStaggeredStochasticForcing

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

// IBTK INCLUDES
#include <ibtk/CartGridFunction.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <NodeVariable.h>
#include <SideVariable.h>
#include <VariableContext.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredStochasticForcing provides an interface for
 * specifying a stochastic forcing term for the staggered-grid incompressible
 * Navier-Stokes solver.
 */
class INSStaggeredStochasticForcing
    : public virtual IBTK::CartGridFunction
{
public:
    /*!
     * \brief This constructor creates Variable and VariableContext objects for
     * storing the stochastic stresses at the centers and nodes of the Cartesian
     * grid.
     */
    INSStaggeredStochasticForcing(
        const INSStaggeredHierarchyIntegrator* const fluid_solver);

    /*!
     * \brief Empty virtual destructor.
     */
    virtual
    ~INSStaggeredStochasticForcing();

    /*!
     * Set objects to provide boundary conditions for the fluid velocity field.
     */
    void
    setPhysicalBcCoefs(
        blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> u_bc_coef);

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete INSStaggeredStochasticForcing object is
     * time-dependent.
     */
    virtual bool
    isTimeDependent() const;

    /*!
     * \brief Set the form of the stress tensor.
     */
    void
    setStochasticStressTensorType(
        const StochasticStressTensorType stress_tensor_type);

    /*!
     * \brief Set the density of the fluid.
     */
    void
    setFluidDensity(
        const double rho);

    /*!
     * \brief Set the dynamic viscosity of the fluid.
     */
    void
    setFluidViscosity(
        const double mu);

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
     * \brief Set the stochastic flux scaling factor to use at Dirichlet
     * boundaries.
     */
    void
    setDirichletBcScaling(
        const double dirichlet_scaling);

    /*!
     * \brief Set the stochastic flux scaling factor to use at Neumann
     * boundaries.
     */
    void
    setNeumannBcScaling(
        const double neumann_scaling);

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
    INSStaggeredStochasticForcing();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredStochasticForcing(
        const INSStaggeredStochasticForcing& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredStochasticForcing&
    operator=(
        const INSStaggeredStochasticForcing& that);

    /*!
     * Boundary condition specification objects.
     */
    blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> d_u_bc_coef;
    double d_dirichlet_scaling, d_neumann_scaling;

    /*!
     * Pointer to the fluid solver object that is using this stochastic force
     * generator.
     */
    const INSStaggeredHierarchyIntegrator* const d_fluid_solver;

    /*!
     * VariableContext and Variable objects for storing the components of the
     * stochastic stresses.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_W_cc_var;
    int d_W_cc_idx;
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> > d_W_nc_var;
    int d_W_nc_idx;
#endif
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM,double> > d_W_ec_var;
    int d_W_ec_idx;
#endif
    StochasticStressTensorType d_stress_tensor_type;
    double d_rho, d_mu, d_dt, d_std;
    SAMRAI::tbox::Array<int> d_regen_rand_cycle;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/INSStaggeredStochasticForcing.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredStochasticForcing
