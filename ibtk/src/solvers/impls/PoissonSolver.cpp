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

#include <ibtk/PoissonSolver.h>
#include <ibtk/PoissonUtilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>

#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
PoissonSolver::setPoissonSpecifications(const PoissonSpecifications& poisson_spec)
{
    PoissonUtilities::checkCoefficientDepths(poisson_spec, d_object_name);
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
            d_bc_coefs[l] = d_default_bc_coef.get();
        }
    }
    return;
} // setPhysicalBcCoefs

bool
PoissonSolver::supportsBdryInterpType(const BdryInterpType bdry_interp_type) const
{
    return bdry_interp_type == BdryInterpType::LINEAR;
} // supportsBdryInterpType

void
PoissonSolver::setBdryInterpType(const BdryInterpType bdry_interp_type)
{
    TBOX_ASSERT(bdry_interp_type != BdryInterpType::UNKNOWN_BDRY_INTERP_TYPE);
    if (!supportsBdryInterpType(bdry_interp_type))
    {
        TBOX_ERROR(d_object_name << "::setBdryInterpType():\n"
                                 << "  this solver does not support bdry_interp_type = "
                                 << enum_to_string(bdry_interp_type) << "\n"
                                 << "  use bdry_interp_type = LINEAR" << std::endl);
    }
    d_bdry_interp_type = bdry_interp_type;
    return;
} // setBdryInterpType

void
PoissonSolver::setBdryInterpTypeFromInput(const Pointer<Database>& input_db)
{
    if (input_db && input_db->isString("bdry_interp_type"))
    {
        const std::string type_name = input_db->getString("bdry_interp_type");
        const BdryInterpType bdry_interp_type = string_to_enum<BdryInterpType>(type_name);
        if (bdry_interp_type == BdryInterpType::UNKNOWN_BDRY_INTERP_TYPE)
        {
            TBOX_ERROR(d_object_name << "::setBdryInterpTypeFromInput():\n"
                                     << "  unknown bdry_interp_type \"" << type_name << "\"\n"
                                     << "  valid values are \"LINEAR\" and \"QUADRATIC\"\n");
        }
        setBdryInterpType(bdry_interp_type);
    }
    return;
} // setBdryInterpTypeFromInput

void
PoissonSolver::initSpecialized(const std::string& object_name, const bool /*homogeneous_bc*/)
{
    // Initialize the Poisson specifications.
    PoissonSpecifications poisson_spec(object_name + "::poisson_spec");
    poisson_spec.setCZero();
    poisson_spec.setDConstant(-1.0);
    setPoissonSpecifications(poisson_spec);

    // Initialize the boundary conditions.
    d_default_bc_coef = std::make_unique<LocationIndexRobinBcCoefs<NDIM>>(object_name + "::default_bc_coef",
                                                                          Pointer<Database>(nullptr));
    auto p_default_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_bc_coef.get());
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
