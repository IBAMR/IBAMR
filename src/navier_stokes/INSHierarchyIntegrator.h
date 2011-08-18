// Filename: INSHierarchyIntegrator.h
// Created on 10 Aug 2011 by Boyce Griffith
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

#ifndef included_INSHierarchyIntegrator
#define included_INSHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/HierarchyIntegrator.h>

// SAMRAI INCLUDES
#include <LocationIndexRobinBcCoefs.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSHierarchyIntegrator provides an abstract interface for a time
 * integrator for the incompressible Navier-Stokes equations on an AMR grid
 * hierarchy, along with basic data management for variables defined on that
 * hierarchy.
 */
class INSHierarchyIntegrator
    : public HierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart=true);

    /*!
     * The destructor for class INSHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSHierarchyIntegrator();

    /*!
     * Supply initial conditions for the velocity field.
     */
    void
    registerVelocityInitialConditions(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> U_init);

    /*!
     * Supply a physical boundary conditions specificaion for the velocity
     * field.
     */
    void
    registerPhysicalBoundaryConditions(
        const blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs);

    /*!
     * Supply initial conditions for the pressure.
     *
     * \note These initial conditions are used for output purposes only.  They
     * are not actually used in the computation.
     */
    void
    registerPressureInitialConditions(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> P_init);

    /*!
     * Supply a body force (optional).
     */
    void
    registerBodyForceFunction(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Supply a fluid source/sink distribution (optional).
     */
    void
    registerFluidSourceFunction(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q_fcn);

    /*!
     * Return a pointer to the fluid velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >
    getVelocityVariable() const;

    /*!
     * Return a pointer to the fluid pressure state variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >
    getPressureVariable() const;

    /*!
     * Return a pointer to the body force variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >
    getForceVariable() const;

    /*!
     * Return a pointer to the source strength variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >
    getSourceVariable() const;

protected:
    /*
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_U_var, d_P_var, d_F_var, d_Q_var;

    /*
     * Objects to set initial conditions, boundary conditions, body forces, and
     * fluid source/sink distributions.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_U_init, d_P_init;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM> d_default_bc_coefs;
    blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> d_bc_coefs;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_F_fcn, d_Q_fcn;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSHierarchyIntegrator(
        const INSHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSHierarchyIntegrator&
    operator=(
        const INSHierarchyIntegrator& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSHierarchyIntegrator
