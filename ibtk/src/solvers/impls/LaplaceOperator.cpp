// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/LaplaceOperator.h>
#include <ibtk/LinearOperator.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBox.h>
#include <SAMRAIDatabase.h>
#include <SAMRAILocationIndexRobinBcCoefs.h>
#include <SAMRAIPointer.h>
#include <SAMRAIPoissonSpecifications.h>
#include <SAMRAIRobinBcCoefStrategy.h>

#include <string>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

LaplaceOperator::LaplaceOperator(std::string object_name, bool homogeneous_bc)
    : LinearOperator(std::move(object_name), homogeneous_bc),
      d_poisson_spec(d_object_name + "::poisson_spec"),
      d_default_bc_coef(new SAMRAILocationIndexRobinBcCoefs(d_object_name + "::default_bc_coef",
                                                            SAMRAIPointer<SAMRAIDatabase>(nullptr))),
      d_bc_coefs(1, d_default_bc_coef.get())
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    auto p_default_bc_coef = dynamic_cast<SAMRAILocationIndexRobinBcCoefs*>(d_default_bc_coef.get());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }
    return;
} // LaplaceOperator()

void
LaplaceOperator::setPoissonSpecifications(const SAMRAIPoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

const SAMRAIPoissonSpecifications&
LaplaceOperator::getPoissonSpecifications() const
{
    return d_poisson_spec;
} // getPoissonSpecifications

void
LaplaceOperator::setPhysicalBcCoef(SAMRAIRobinBcCoefStrategy* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<SAMRAIRobinBcCoefStrategy*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
LaplaceOperator::setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs)
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
            d_bc_coefs[l] = d_default_bc_coef.get();
        }
    }
    return;
} // setPhysicalBcCoefs

const std::vector<SAMRAIRobinBcCoefStrategy*>&
LaplaceOperator::getPhysicalBcCoefs() const
{
    return d_bc_coefs;
} // getPhysicalBcCoefs

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
