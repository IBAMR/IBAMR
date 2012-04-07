// Filename: LaplaceOperator.C
// Created on 07 Apr 2012 by Boyce Griffith
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

#include "LaplaceOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LaplaceOperator::LaplaceOperator(
    PoissonSpecifications poisson_spec,
    RobinBcCoefStrategy<NDIM>* default_bc_coef,
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
    bool is_symmetric,
    bool homogeneous_bc)
    : LinearOperator(is_symmetric, homogeneous_bc),
      d_apply_time(0.0),
      d_poisson_spec("LaplaceOperator::poisson_spec"),
      d_default_bc_coef(default_bc_coef),
      d_bc_coefs()
{
    setPoissonSpecifications(poisson_spec);
    setPhysicalBcCoefs(bc_coefs);
    return;
}// LaplaceOperator()

LaplaceOperator::~LaplaceOperator()
{
    deallocateOperatorState();
    if (d_default_bc_coef != NULL) delete d_default_bc_coef;
    d_default_bc_coef = NULL;
    return;
}// ~LaplaceOperator()

void
LaplaceOperator::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
}// setPoissonSpecifications

void
LaplaceOperator::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1,bc_coef));
    return;
}// setPhysicalBcCoef

void
LaplaceOperator::setPhysicalBcCoefs(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    d_bc_coefs.resize(bc_coefs.size());
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        if (bc_coefs[l] != NULL)
        {
            d_bc_coefs[l] = bc_coefs[l];
        }
        else
        {
            d_bc_coefs[l] = d_default_bc_coef;
        }
    }
    return;
}// setPhysicalBcCoefs

void
LaplaceOperator::setPhysicalBcCoefs(
    const blitz::TinyVector<RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(&bc_coefs[0],&bc_coefs[0]+NDIM));
    return;
}// setPhysicalBcCoefs

void
LaplaceOperator::setTime(
    const double time)
{
    d_apply_time = time;
    return;
}// setTime

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
