// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_ViscosityModel
#define included_IBAMR_ViscosityModel

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class ViscosityModel is an interface for a viscosity that depends on the velocity, such as that of a
 * generalized Newtonian fluid.
 *
 * INSVCStaggeredHierarchyIntegrator::setViscosityModel() gives a model to the integrator, which then creates the
 * cell-centered viscosity variable and calls setDataOnPatchHierarchy() from the velocity of the initial conditions and
 * wherever it would call the viscosity reset functions: for the viscosity at the current time at the start of each time
 * step, and for the viscosity at the new time at the start of each cycle and after the last one.
 */
class ViscosityModel
{
public:
    virtual ~ViscosityModel() = default;

    /*!
     * Return the number of ghost cells of the velocity that setDataOnPatchHierarchy() reads.
     */
    virtual int getVelocityGhostWidth() const = 0;

    /*!
     * Set the interior values of the cell-centered viscosity with patch data index \p mu_idx on every level of the
     * patch hierarchy of \p hier_math_ops from the side-centered velocity with patch data index \p U_idx at time
     * \p time. The velocity has at least getVelocityGhostWidth() ghost cells, which are filled.
     */
    virtual void setDataOnPatchHierarchy(int mu_idx, int U_idx, IBTK::HierarchyMathOps& hier_math_ops, double time) = 0;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ViscosityModel
