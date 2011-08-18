// Filename: INSStaggeredHierarchyIntegrator.C
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

#include "INSHierarchyIntegrator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSHierarchyIntegrator::INSHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart),
      d_U_var(NULL),
      d_P_var(NULL),
      d_F_var(NULL),
      d_Q_var(NULL),
      d_U_init(NULL),
      d_P_init(NULL),
      d_default_bc_coefs(d_object_name+"::default_bc_coefs", Pointer<Database>(NULL)),
      d_bc_coefs(static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_F_fcn(NULL),
      d_Q_fcn(NULL)
{
    // Setup default boundary condition objects that specify homogeneous
    // Dirichlet (solid-wall) boundary conditions for the velocity.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coefs.setBoundaryValue(2*d  ,0.0);
        d_default_bc_coefs.setBoundaryValue(2*d+1,0.0);
    }
    registerPhysicalBoundaryConditions(blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>(&d_default_bc_coefs));
    return;
}// INSHierarchyIntegrator

INSHierarchyIntegrator::~INSHierarchyIntegrator()
{
    // intentionally blank
    return;
}// ~INSHierarchyIntegrator

void
INSHierarchyIntegrator::registerVelocityInitialConditions(
    Pointer<CartGridFunction> U_init)
{
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
INSHierarchyIntegrator::registerPhysicalBoundaryConditions(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    d_bc_coefs = bc_coefs;
    return;
}// registerPhysicalBoundaryConditions

void
INSHierarchyIntegrator::registerPressureInitialConditions(
    Pointer<CartGridFunction> P_init)
{
    d_P_init = P_init;
    return;
}// registerPressureInitialConditions

void
INSHierarchyIntegrator::registerBodyForceFunction(
    Pointer<CartGridFunction> F_fcn)
{
    d_F_fcn = F_fcn;
    return;
}// registerBodyForceFunction

void
INSHierarchyIntegrator::registerFluidSourceFunction(
    Pointer<CartGridFunction> Q_fcn)
{
    d_Q_fcn = Q_fcn;
    return;
}// registerFluidSourceFunction

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getVelocityVariable() const
{
    return d_U_var;
}// getVelocityVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getPressureVariable() const
{
    return d_P_var;
}// getPressureVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getForceVariable() const
{
    return d_F_var;
}// getForceVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getSourceVariable() const
{
    return d_Q_var;
}// getSourceVariable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::INSHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
