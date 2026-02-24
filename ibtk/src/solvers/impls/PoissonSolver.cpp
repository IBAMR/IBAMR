// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2024 by the IBAMR developers
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

#include <ibtk/PoissonSolver.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIBox.h>
#include <SAMRAIDatabase.h>
#include <SAMRAILocationIndexRobinBcCoefs.h>
#include <SAMRAIPointer.h>
#include <SAMRAIPoissonSpecifications.h>
#include <SAMRAIRobinBcCoefStrategy.h>

#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PoissonSolver::setPoissonSpecifications(const SAMRAIPoissonSpecifications& poisson_spec)
{
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

void
PoissonSolver::setPhysicalBcCoef(SAMRAIRobinBcCoefStrategy* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<SAMRAIRobinBcCoefStrategy*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
PoissonSolver::setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs)
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

void
PoissonSolver::initSpecialized(const std::string& object_name, const bool /*homogeneous_bc*/)
{
    // Initialize the Poisson specifications.
    SAMRAIPoissonSpecifications poisson_spec(object_name + "::poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);
    setPoissonSpecifications(poisson_spec);

    // Initialize the boundary conditions.
    d_default_bc_coef = std::make_unique<SAMRAILocationIndexRobinBcCoefs>(object_name + "::default_bc_coef",
                                                                          SAMRAIPointer<SAMRAIDatabase>(nullptr));
    auto p_default_bc_coef = dynamic_cast<SAMRAILocationIndexRobinBcCoefs*>(d_default_bc_coef.get());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }
    setPhysicalBcCoef(d_default_bc_coef.get());
} // initSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
