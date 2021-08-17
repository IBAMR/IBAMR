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

#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/INSStaggeredPressureBcCoef.h"
#include "ibamr/StokesBcCoefStrategy.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/ExtendedRobinBcCoefStrategy.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "RobinBcCoefStrategy.h"
#include "SideData.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredPressureBcCoef::INSStaggeredPressureBcCoef(const INSStaggeredHierarchyIntegrator* fluid_solver,
                                                       const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                       const TractionBcType traction_bc_type,
                                                       const bool homogeneous_bc)
    : d_fluid_solver(fluid_solver), d_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(nullptr))
{
    setStokesSpecifications(d_fluid_solver->getStokesSpecifications());
    setPhysicalBcCoefs(bc_coefs);
    setTractionBcType(traction_bc_type);
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSStaggeredPressureBcCoef

void
INSStaggeredPressureBcCoef::setStokesSpecifications(const StokesSpecifications* problem_coefs)
{
    StokesBcCoefStrategy::setStokesSpecifications(problem_coefs);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setStokesSpecifications(problem_coefs);
    }
    return;
} // setStokesSpecifications

void
INSStaggeredPressureBcCoef::setTargetVelocityPatchDataIndex(int u_target_data_idx)
{
    StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(u_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
    }
    return;
} // setTargetVelocityPatchDataIndex

void
INSStaggeredPressureBcCoef::clearTargetVelocityPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetVelocityPatchDataIndex();
    }
    return;
} // clearTargetVelocityPatchDataIndex

void
INSStaggeredPressureBcCoef::setTargetPressurePatchDataIndex(int p_target_data_idx)
{
    StokesBcCoefStrategy::setTargetPressurePatchDataIndex(p_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
    }
    return;
} // setTargetPressurePatchDataIndex

void
INSStaggeredPressureBcCoef::clearTargetPressurePatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPressurePatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPressurePatchDataIndex();
    }
    return;
} // clearTargetPressurePatchDataIndex

void
INSStaggeredPressureBcCoef::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    TBOX_ASSERT(bc_coefs.size() == NDIM);
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSStaggeredPressureBcCoef::setSolutionTime(const double /*solution_time*/)
{
    // intentionally blank
    return;
} // setSolutionTime

void
INSStaggeredPressureBcCoef::setTimeInterval(const double /*current_time*/, const double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSStaggeredPressureBcCoef::setTargetPatchDataIndex(int target_idx)
{
    StokesBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
} // setTargetPatchDataIndex

void
INSStaggeredPressureBcCoef::clearTargetPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // clearTargetPatchDataIndex

void
INSStaggeredPressureBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
} // setHomogeneousBc

void
INSStaggeredPressureBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                       Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                       Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                       const Pointer<Variable<NDIM> >& variable,
                                       const Patch<NDIM>& patch,
                                       const BoundaryBox<NDIM>& bdry_box,
                                       double /*fill_time*/) const
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);

    Box<NDIM> bc_coef_box = acoef_data->getBox();
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(!gcoef_data || (bc_coef_box == gcoef_data->getBox()));

    // Set the unmodified velocity bc coefs.
    const int location_index = bdry_box.getLocationIndex();
    const int bdry_normal_axis = location_index / 2;
    const bool at_lower_bdry = location_index % 2 == 0;
    const double half_time = d_fluid_solver->getIntegratorTime() + 0.5 * d_fluid_solver->getCurrentTimeStepSize();
    d_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, half_time);

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc && gcoef_data) gcoef_data->fillAll(0.0);

    // Get the target velocity data.
    Pointer<SideData<NDIM, double> > u_current_data, u_target_data;
    if (gcoef_data)
    {
        u_current_data = patch.getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
        TBOX_ASSERT(u_current_data);

        if (d_u_target_data_idx >= 0)
            u_target_data = patch.getPatchData(d_u_target_data_idx);
        else if (d_target_data_idx >= 0)
            u_target_data = patch.getPatchData(d_target_data_idx);
        TBOX_ASSERT(u_target_data);

        const Box<NDIM> ghost_box = u_current_data->getGhostBox() * u_target_data->getGhostBox();
        for (int d = 0; d < NDIM; ++d)
        {
            if (d != bdry_normal_axis)
            {
                bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), ghost_box.lower(d));
                bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), ghost_box.upper(d));
            }
        }
    }

    // Update the boundary condition coefficients.  Normal velocity boundary
    // conditions are converted into Neumann conditions for the pressure, and
    // normal traction boundary conditions are converted into Dirichlet
    // conditions for the pressure.
    const double mu = d_fluid_solver->getStokesSpecifications()->getMu();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        double dummy_val = std::numeric_limits<double>::quiet_NaN();
        double& a = acoef_data ? (*acoef_data)(i, 0) : dummy_val;
        double& b = bcoef_data ? (*bcoef_data)(i, 0) : dummy_val;
        double& g = gcoef_data ? (*gcoef_data)(i, 0) : dummy_val;
        const bool velocity_bc = abs(b) < sqrt(std::numeric_limits<double>::epsilon());
        if (velocity_bc)
        {
            a = 0.0;
            b = 1.0;
            g = 0.0;
        }
        else
        {
            // Place i_i in the interior cell abutting the boundary.
            pdat::CellIndex<NDIM> i_i(i);
            if (!at_lower_bdry)
            {
                i_i(bdry_normal_axis) -= 1;
            }

            // Extract the first NVALS interior values of the normal velocity, starting at the boundary and moving
            // towards the interior. We use these values to approximate du_norm/dx_norm.
            static const int NVALS = 3;
            std::array<double, NVALS> u_current, u_new;
            SideIndex<NDIM> i_s(i_i, bdry_normal_axis, at_lower_bdry ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
            for (int k = 0; k < NVALS; ++k, i_s(bdry_normal_axis) += (at_lower_bdry ? 1 : -1))
            {
                u_current[k] = u_current_data ? 0.0 : (*u_current_data)(i_s);
                u_new[k] = u_target_data ? 0.0 : (*u_target_data)(i_s);
            }

            const double u_norm_current = (at_lower_bdry ? -1.0 : +1.0) * u_current[0];
            const double u_norm_new = (at_lower_bdry ? -1.0 : +1.0) * u_new[0];

            const double h = dx[bdry_normal_axis];
            const double du_norm_dx_norm_current =
                (at_lower_bdry ? +1.0 : -1.0) * (2.0 * u_current[1] - 1.5 * u_current[0] - 0.5 * u_current[2]) / h;
            const double du_norm_dx_norm_new =
                (at_lower_bdry ? +1.0 : -1.0) * (2.0 * u_new[1] - 1.5 * u_new[0] - 0.5 * u_new[2]) / h;

            // Here we reinterpret the velocity boundary condition coefficients to provide boundary conditions for the
            // pressure.
            double a_new = std::numeric_limits<double>::quiet_NaN();
            double b_new = std::numeric_limits<double>::quiet_NaN();
            double g_new = std::numeric_limits<double>::quiet_NaN();

            switch (d_traction_bc_type)
            {
            case TRACTION: // a u_norm + b (-p + 2*mu*du_norm/dx_norm) = g
                           // ===> -p + 2 mu du_norm/dx_norm = (g - a u_norm) / b
            {
                // The boundary condition is -p + 2*mu*du_norm/dx_norm = (g - a u_norm) / b.
                //
                // Because p is centered about t^{n+1/2}, we compute this
                // as:
                //
                // p^{n+1/2} = mu*du_norm/dx_norm^{n} + mu*du_norm/dx_norm^{n+1}
                //             - (g^{n+1/2} - 0.5 a (u_norm^{n} + u_norm^{n+1})) / b.
                a_new = 1.0;
                b_new = 0.0;
                g_new = (d_homogeneous_bc ? 0.0 : mu * du_norm_dx_norm_current + 0.5 * a * u_norm_current / b) +
                        mu * du_norm_dx_norm_new - (g - 0.5 * a * u_norm_new) / b;
                break;
            }
            case PSEUDO_TRACTION: // a u_norm + b (-p) = g
                                  // ===> -p = (g - a u_norm) / b
            {
                a_new = 1.0;
                b_new = 0.0;
                g_new = (d_homogeneous_bc ? 0.0 : 0.5 * a * u_norm_current / b) - (g - 0.5 * a * u_norm_new) / b;
                break;
            }
            default:
            {
                TBOX_ERROR(
                    "INSStaggeredPressureBcCoef::setBcCoefs(): unrecognized or "
                    "unsupported traction boundary condition type: "
                    << enum_to_string<TractionBcType>(d_traction_bc_type) << "\n");
            }
            }

            a = a_new;
            b = b_new;
            g = g_new;
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
INSStaggeredPressureBcCoef::numberOfExtensionsFillable() const
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
