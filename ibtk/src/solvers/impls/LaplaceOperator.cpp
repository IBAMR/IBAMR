// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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
#include <ibtk/PoissonUtilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>

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
      d_default_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_bc_coef", Pointer<Database>(nullptr))),
      d_bc_coefs(1, d_default_bc_coef.get())
{
    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions.
    auto p_default_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_bc_coef.get());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        p_default_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
    }
    return;
} // LaplaceOperator()

void
LaplaceOperator::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    PoissonUtilities::checkCoefficientDepths(poisson_spec, d_object_name);
    d_poisson_spec = poisson_spec;
    return;
} // setPoissonSpecifications

const PoissonSpecifications&
LaplaceOperator::getPoissonSpecifications() const
{
    return d_poisson_spec;
} // getPoissonSpecifications

void
LaplaceOperator::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
LaplaceOperator::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
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

const std::vector<RobinBcCoefStrategy<NDIM>*>&
LaplaceOperator::getPhysicalBcCoefs() const
{
    return d_bc_coefs;
} // getPhysicalBcCoefs

bool
LaplaceOperator::supportsBdryInterpType(const std::string& bdry_interp_type) const
{
    return bdry_interp_type == "LINEAR";
} // supportsBdryInterpType

void
LaplaceOperator::setBdryInterpType(const std::string& bdry_interp_type)
{
    if (!supportsBdryInterpType(bdry_interp_type))
    {
        TBOX_ERROR(d_object_name << "::setBdryInterpType():\n"
                                 << "  this operator does not support bdry_interp_type = " << bdry_interp_type << "\n"
                                 << "  use bdry_interp_type = LINEAR" << std::endl);
    }
    d_bdry_interp_type = bdry_interp_type;
    return;
} // setBdryInterpType

void
LaplaceOperator::setBdryExtrapType(const std::string& /*bdry_extrap_type*/)
{
    return;
} // setBdryExtrapType

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
