// Filename: INSVCStaggeredPressureBcCoef.cpp
// Created on 25 Sep 2017 by Nishant Nangia
//
// Copyright (c) 2002-2014, Nishant Nangia and Amneet Bhalla
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

#include <algorithm>
#include <limits>
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

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
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredPressureBcCoef.h"
#include "ibamr/StokesBcCoefStrategy.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/ExtendedRobinBcCoefStrategy.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredPressureBcCoef::INSVCStaggeredPressureBcCoef(const INSVCStaggeredHierarchyIntegrator* fluid_solver,
                                                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                           const TractionBcType traction_bc_type,
                                                           const bool homogeneous_bc)
    : d_fluid_solver(fluid_solver), d_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL))
{
    setStokesSpecifications(d_fluid_solver->getStokesSpecifications());
    setPhysicalBcCoefs(bc_coefs);
    setTractionBcType(traction_bc_type);
    setHomogeneousBc(homogeneous_bc);

    // Set some default values.
    d_mu_interp_type = VC_HARMONIC_INTERP;
    return;
} // INSVCStaggeredPressureBcCoef

INSVCStaggeredPressureBcCoef::~INSVCStaggeredPressureBcCoef()
{
    // intentionally blank
    return;
} // ~INSVCStaggeredPressureBcCoef

void
INSVCStaggeredPressureBcCoef::setStokesSpecifications(const StokesSpecifications* problem_coefs)
{
    StokesBcCoefStrategy::setStokesSpecifications(problem_coefs);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        StokesBcCoefStrategy* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setStokesSpecifications(problem_coefs);
    }
    return;
} // setStokesSpecifications

void
INSVCStaggeredPressureBcCoef::setTargetVelocityPatchDataIndex(int u_target_data_idx)
{
    StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(u_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        StokesBcCoefStrategy* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
    }
    return;
} // setTargetVelocityPatchDataIndex

void
INSVCStaggeredPressureBcCoef::clearTargetVelocityPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        StokesBcCoefStrategy* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetVelocityPatchDataIndex();
    }
    return;
} // clearTargetVelocityPatchDataIndex

void
INSVCStaggeredPressureBcCoef::setTargetPressurePatchDataIndex(int p_target_data_idx)
{
    StokesBcCoefStrategy::setTargetPressurePatchDataIndex(p_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        StokesBcCoefStrategy* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
    }
    return;
} // setTargetPressurePatchDataIndex

void
INSVCStaggeredPressureBcCoef::clearTargetPressurePatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPressurePatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        StokesBcCoefStrategy* p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPressurePatchDataIndex();
    }
    return;
} // clearTargetPressurePatchDataIndex

void
INSVCStaggeredPressureBcCoef::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSVCStaggeredPressureBcCoef::setSolutionTime(const double /*solution_time*/)
{
    // intentionally blank
    return;
} // setSolutionTime

void
INSVCStaggeredPressureBcCoef::setTimeInterval(const double /*current_time*/, const double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSVCStaggeredPressureBcCoef::setTargetPatchDataIndex(int target_idx)
{
    StokesBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
} // setTargetPatchDataIndex

void
INSVCStaggeredPressureBcCoef::clearTargetPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // clearTargetPatchDataIndex

void
INSVCStaggeredPressureBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
} // setHomogeneousBc

void
INSVCStaggeredPressureBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                         Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                         Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                         const Pointer<Variable<NDIM> >& variable,
                                         const Patch<NDIM>& patch,
                                         const BoundaryBox<NDIM>& bdry_box,
                                         double /*fill_time*/) const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);
#endif
    Box<NDIM> bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(!gcoef_data || (bc_coef_box == gcoef_data->getBox()));
#endif
    // Set the unmodified velocity bc coefs.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    const bool is_lower = location_index % 2 == 0;
    const double half_time = d_fluid_solver->getIntegratorTime() + 0.5 * d_fluid_solver->getCurrentTimeStepSize();
    d_bc_coefs[bdry_normal_axis]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, half_time);

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc && gcoef_data) gcoef_data->fillAll(0.0);

    // Get the target velocity data.
    Pointer<SideData<NDIM, double> > u_target_data;
    if (d_u_target_data_idx >= 0)
        u_target_data = patch.getPatchData(d_u_target_data_idx);
    else if (d_target_data_idx >= 0)
        u_target_data = patch.getPatchData(d_target_data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(u_target_data);
#endif
    Pointer<SideData<NDIM, double> > u_current_data =
        patch.getPatchData(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
#if !defined(NDEBUG)
    TBOX_ASSERT(u_current_data);
#endif
    const Box<NDIM> ghost_box = u_target_data->getGhostBox() * u_current_data->getGhostBox();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d != bdry_normal_axis)
        {
            bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), ghost_box.lower(d));
            bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), ghost_box.upper(d));
        }
    }

    // Update the boundary condition coefficients.  Normal velocity boundary
    // conditions are converted into Neumann conditions for the pressure, and
    // normal traction boundary conditions are converted into Dirichlet
    // conditions for the pressure.
    Array<double> A_scale = d_fluid_solver->getScalingFactor();
    const double p_scale = A_scale[patch.getPatchLevelNumber()];
    double mu = d_fluid_solver->muIsConstant() ? d_problem_coefs->getMu() : -1;
    int mu_idx = -1;
    Pointer<CellData<NDIM, double> > mu_data;
    if (!d_fluid_solver->muIsConstant())
    {
        mu_idx = d_fluid_solver->getLinearOperatorMuPatchDataIndex();
#if !defined(NDEBUG)
        TBOX_ASSERT(mu_idx >= 0);
#endif
        mu_data = patch.getPatchData(mu_idx);
    }

    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const Index<NDIM>& i = it();
        double dummy_val;
        double& alpha = acoef_data ? (*acoef_data)(i, 0) : dummy_val;
        double& beta = bcoef_data ? (*bcoef_data)(i, 0) : dummy_val;
        double& gamma = gcoef_data ? (*gcoef_data)(i, 0) : dummy_val;
        const bool velocity_bc = MathUtilities<double>::equalEps(alpha, 1.0);
        const bool traction_bc = MathUtilities<double>::equalEps(beta, 1.0);
#if !defined(NDEBUG)
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (velocity_bc)
        {
            alpha = 0.0;
            beta = 1.0;
            gamma = 0.0;
        }
        else if (traction_bc)
        {
#if !defined(NDEBUG)
            if (!d_fluid_solver->muIsConstant()) TBOX_ASSERT(mu_data);
#endif
            switch (d_traction_bc_type)
            {
            case TRACTION: // -p + 2*mu*du_n/dx_n = g.
            {
#if !defined(NDEBUG)
                if (!d_fluid_solver->muIsConstant()) TBOX_ASSERT(mu_data);
#endif
                // Place i_i in the interior cell abutting the boundary, and
                // place i_g in the ghost cell abutting the boundary.
                Index<NDIM> i_i(i), i_g(i);
                if (is_lower)
                {
                    i_g(bdry_normal_axis) -= 1;
                }
                else
                {
                    i_i(bdry_normal_axis) -= 1;
                }
                if (!d_fluid_solver->muIsConstant())
                {
                    // Average mu onto side centers
                    // In certain use cases with traction boundary conditions, this class will attempt to fill
                    // Robin BC coefficient values along an extended physical boundary outside of the physical domain
                    // that will never be used. However, viscosity values will not be available at those locations, so
                    // we
                    // need to ensure we don't access those unallocated data. The unphysical value should result in bad
                    // stuff if
                    // it gets used for whatever reason.
                    const bool mu_contains_cells =
                        (mu_data->getGhostBox().contains(i_g) && mu_data->getGhostBox().contains(i_i));
                    if (d_mu_interp_type == VC_AVERAGE_INTERP)
                    {
                        mu = mu_contains_cells ? 0.5 * ((*mu_data)(i_g) + (*mu_data)(i_i)) : -1.0e305;
                    }
                    else if (d_mu_interp_type == VC_HARMONIC_INTERP)
                    {
                        mu = mu_contains_cells ?
                                 2.0 * ((*mu_data)(i_g) * (*mu_data)(i_i)) / ((*mu_data)(i_g) + (*mu_data)(i_i)) :
                                 -1.0e305;
                    }
                    else
                    {
                        TBOX_ERROR("this statement should not be reached");
                    }
                }

                // The boundary condition is -p + 2*mu*du_n/dx_n = g.
                //
                // Because p is centered about t^{n+1/2}, we compute this
                // as:
                //
                // p^{n+1/2} = mu*du_n/dx_n^{n} + mu*du_n/dx_n^{n+1} - g^{n+1/2}.
                static const int NVALS = 3;
                double u_current[NVALS], u_new[NVALS];
                SideIndex<NDIM> i_s(i_i, bdry_normal_axis, is_lower ? SideIndex<NDIM>::Lower : SideIndex<NDIM>::Upper);
                for (int k = 0; k < NVALS; ++k, i_s(bdry_normal_axis) += (is_lower ? 1 : -1))
                {
                    u_current[k] = (*u_current_data)(i_s);
                    u_new[k] = (*u_target_data)(i_s);
                }
                const double h = dx[bdry_normal_axis];
                const double du_norm_current_dx_norm =
                    (is_lower ? +1.0 : -1.0) * (2.0 * u_current[1] - 1.5 * u_current[0] - 0.5 * u_current[2]) / h;
                const double du_norm_new_dx_norm =
                    (is_lower ? +1.0 : -1.0) * (2.0 * u_new[1] - 1.5 * u_new[0] - 0.5 * u_new[2]) / h;
                alpha = 1.0;
                beta = 0.0;
                gamma = ((d_homogeneous_bc ? 0.0 : mu * du_norm_current_dx_norm) + mu * du_norm_new_dx_norm - gamma) *
                        p_scale;
                break;
            }
            case PSEUDO_TRACTION: // -p = g.
            {
                alpha = 1.0;
                beta = 0.0;
                gamma = -gamma * p_scale;
                break;
            }
            default:
            {
                TBOX_ERROR(
                    "INSVCStaggeredPressureBcCoef::setBcCoefs(): unrecognized or unsupported "
                    "traction boundary condition type: "
                    << enum_to_string<TractionBcType>(d_traction_bc_type) << "\n");
            }
            }
        }
        else
        {
            TBOX_ERROR("this statement should not be reached!\n");
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
INSVCStaggeredPressureBcCoef::numberOfExtensionsFillable() const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    IntVector<NDIM> ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = IntVector<NDIM>::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
} // numberOfExtensionsFillable

void
INSVCStaggeredPressureBcCoef::setViscosityInterpolationType(const IBTK::VCInterpType mu_interp_type)
{
    d_mu_interp_type = mu_interp_type;
    return;
} // setViscosityInterpolationType

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
