// Filename: PoissonSolver.cpp
// Created on 07 Apr 2012 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "ibtk/PoissonSolver.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Database.h"
#include "tbox/Pointer.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

PoissonSolver::PoissonSolver() : d_poisson_spec(""), d_default_bc_coef(NULL), d_bc_coefs()
{
    // intentionally blank
    return;
} // PoissonSolver()

PoissonSolver::~PoissonSolver()
{
    delete d_default_bc_coef;
    d_default_bc_coef = NULL;
    return;
} // ~PoissonSolver()

void
PoissonSolver::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

void
PoissonSolver::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
PoissonSolver::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l])
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
} // setPhysicalBcCoefs

void
PoissonSolver::initSpecialized(const std::string& object_name, const bool /*homogeneous_bc*/)
{
    // Initialize the Poisson specifications.
    PoissonSpecifications poisson_spec(object_name + "::poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);
    setPoissonSpecifications(poisson_spec);

    // Initialize the boundary conditions.
    d_default_bc_coef = new LocationIndexRobinBcCoefs<NDIM>(object_name + "::default_bc_coef", Pointer<Database>(NULL));
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        LocationIndexRobinBcCoefs<NDIM>* p_default_bc_coef =
            dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_bc_coef);
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }
    setPhysicalBcCoef(d_default_bc_coef);
} // initSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
