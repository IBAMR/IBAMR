// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2026 by the IBAMR developers
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

#include <ibamr/BoundaryTreatment.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredVelocityBcCoef.h>
#include <ibamr/StokesBcCoefStrategy.h>
#include <ibamr/StokesSpecifications.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>

#include <tbox/MathUtilities.h>
#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianPatchGeometry.h>
#include <EdgeData.h>
#include <Index.h>
#include <IntVector.h>
#include <NodeData.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>

#include <limits>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "./ins_staggered_traction_stencil.h"
#include "./ins_velocity_bc.h"

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredVelocityBcCoef::INSVCStaggeredVelocityBcCoef(const unsigned int comp_idx,
                                                           const INSVCStaggeredHierarchyIntegrator* fluid_solver,
                                                           const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                           const TractionBcType traction_bc_type,
                                                           const bool homogeneous_bc)
    : d_comp_idx(comp_idx), d_fluid_solver(fluid_solver), d_bc_coefs(NDIM, nullptr)
{
    setStokesSpecifications(d_fluid_solver->getStokesSpecifications());
    setPhysicalBcCoefs(bc_coefs);
    setTractionBcType(traction_bc_type);
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSVCStaggeredVelocityBcCoef

void
INSVCStaggeredVelocityBcCoef::setStokesSpecifications(const StokesSpecifications* problem_coefs)
{
    StokesBcCoefStrategy::setStokesSpecifications(problem_coefs);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setStokesSpecifications(problem_coefs);
    }
    return;
} // setStokesSpecifications

void
INSVCStaggeredVelocityBcCoef::setTargetVelocityPatchDataIndex(int u_target_data_idx)
{
    StokesBcCoefStrategy::setTargetVelocityPatchDataIndex(u_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetVelocityPatchDataIndex(u_target_data_idx);
    }
    return;
} // setTargetVelocityPatchDataIndex

void
INSVCStaggeredVelocityBcCoef::clearTargetVelocityPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetVelocityPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetVelocityPatchDataIndex();
    }
    return;
} // clearTargetVelocityPatchDataIndex

void
INSVCStaggeredVelocityBcCoef::setTargetPressurePatchDataIndex(int p_target_data_idx)
{
    StokesBcCoefStrategy::setTargetPressurePatchDataIndex(p_target_data_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPressurePatchDataIndex(p_target_data_idx);
    }
    return;
} // setTargetPressurePatchDataIndex

void
INSVCStaggeredVelocityBcCoef::clearTargetPressurePatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPressurePatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<StokesBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPressurePatchDataIndex();
    }
    return;
} // clearTargetPressurePatchDataIndex

void
INSVCStaggeredVelocityBcCoef::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSVCStaggeredVelocityBcCoef::setSolutionTime(const double /*solution_time*/)
{
    // intentionally blank
    return;
} // setSolutionTime

void
INSVCStaggeredVelocityBcCoef::setTimeInterval(const double /*current_time*/, const double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSVCStaggeredVelocityBcCoef::setBoundaryTreatment(const BoundaryTreatment& bdry_treatment)
{
    d_bdry_interp_type = bdry_treatment.bdry_interp_type;
    return;
} // setBoundaryTreatment

void
INSVCStaggeredVelocityBcCoef::setTargetPatchDataIndex(int target_idx)
{
    StokesBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
} // setTargetPatchDataIndex

void
INSVCStaggeredVelocityBcCoef::clearTargetPatchDataIndex()
{
    StokesBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // clearTargetPatchDataIndex

void
INSVCStaggeredVelocityBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
} // setHomogeneousBc

bool
INSVCStaggeredVelocityBcCoef::dependsOnTargetPatchData(const Patch<NDIM>& patch,
                                                       const BoundaryBox<NDIM>& bdry_box,
                                                       const double fill_time) const
{
    return d_traction_bc_type == TRACTION && velocity_bc::traction_data_depend_on_velocity(
                                                 d_bc_coefs[d_comp_idx], d_comp_idx, patch, bdry_box, fill_time);
} // dependsOnTargetPatchData

void
INSVCStaggeredVelocityBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double>>& acoef_data,
                                         Pointer<ArrayData<NDIM, double>>& bcoef_data,
                                         Pointer<ArrayData<NDIM, double>>& gcoef_data,
                                         const Pointer<Variable<NDIM>>& variable,
                                         const Patch<NDIM>& patch,
                                         const BoundaryBox<NDIM>& bdry_box,
                                         double fill_time) const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    // Set the unmodified velocity bc coefs.
    d_bc_coefs[d_comp_idx]->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // We do not make any further modifications to the values of acoef_data and
    // bcoef_data beyond this point.
    if (!gcoef_data) return;
#if !defined(NDEBUG)
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);
#endif

    // Ensure homogeneous boundary conditions are enforced.
    if (d_homogeneous_bc) gcoef_data->fillAll(0.0);

    // Get the target velocity data.
    Pointer<SideData<NDIM, double>> u_target_data;
    if (d_u_target_data_idx >= 0)
        u_target_data = patch.getPatchData(d_u_target_data_idx);
    else if (d_target_data_idx >= 0)
        u_target_data = patch.getPatchData(d_target_data_idx);
#if !defined(NDEBUG)
    TBOX_ASSERT(u_target_data);
#endif

    // Where appropriate, update boundary condition coefficients.
    //
    // Dirichlet boundary conditions are not modified.
    //
    // Neumann boundary conditions on the normal component of the velocity are
    // interpreted as "open" boundary conditions, and we set du/dn = 0.
    //
    // Neumann boundary conditions on the tangential component of the velocity
    // are interpreted as traction (stress) boundary conditions, and we update
    // the boundary condition coefficients accordingly.
    const unsigned int location_index = bdry_box.getLocationIndex();
    const unsigned int bdry_normal_axis = location_index / 2;
    const bool is_lower = location_index % 2 == 0;
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch.getPatchGeometry();
    const double* const dx = pgeom->getDx();
    std::optional<BoxArray<NDIM>> domain;
    std::optional<traction_stencil::NormalVelocityPatch> normal_velocity_patch;

    double mu = d_fluid_solver->muIsConstant() ? d_problem_coefs->getMu() : -1;
#if (NDIM == 2)
    Pointer<NodeData<NDIM, double>> mu_data;
#elif (NDIM == 3)
    Pointer<EdgeData<NDIM, double>> mu_data;
#endif
    if (!d_fluid_solver->muIsConstant())
    {
        const int mu_idx = d_fluid_solver->getInterpolatedLinearOperatorMuPatchDataIndex();

#if !defined(NDEBUG)
        TBOX_ASSERT(mu_idx >= 0);
#endif
        mu_data = patch.getPatchData(mu_idx);
    }
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const hier::Index<NDIM>& i = it();
        double& alpha = (*acoef_data)(i, 0);
        double& beta = (*bcoef_data)(i, 0);
        double& gamma = (*gcoef_data)(i, 0);
        if (velocity_bc::get_type(alpha, beta, "INSVCStaggeredVelocityBcCoef::setBcCoefs()") ==
            velocity_bc::Type::TRACTION)
        {
            const bool mu_data_exists = mu_data;
            if (d_comp_idx == bdry_normal_axis)
            {
                // Set du/dn = 0.
                //
                // NOTE: We would prefer to determine the ghost cell value of
                // the normal velocity so that div u = 0 in the ghost cell.
                // This could be done here, but it is more convenient to do so
                // as a post-processing step after the tangential velocity ghost
                // cell values have all been set.
                alpha = 0.0;
                beta = 1.0;
                gamma = 0.0;
            }
            else
            {
                // Use mu = 1 on a patch without viscosity data (a temporary coarse patch, or a new patch of a regrid),
                // which is exact for g = 0.
                if (!d_fluid_solver->muIsConstant() && !mu_data_exists)
                {
                    if (gamma != 0.0)
                    {
                        TBOX_ERROR("INSVCStaggeredVelocityBcCoef::setBcCoefs():\n"
                                   << "  a nonzero tangential traction at the face " << i << " of axis " << d_comp_idx
                                   << " is not supported on a patch without viscosity data (a temporary coarse patch, "
                                      "or a new patch filled during a regrid)\n"
                                   << "  use a zero tangential traction or a constant viscosity\n");
                    }
                    mu = 1.0;
                }
                switch (d_traction_bc_type)
                {
                case TRACTION: // mu*(du_tan/dx_norm + du_norm/dx_tan) = g.
                {
                    // Compute the tangential derivative of the normal
                    // component of the velocity at the boundary.
                    if (!domain)
                    {
                        domain = traction_stencil::get_periodic_physical_domain(d_fluid_solver->getPatchHierarchy(),
                                                                                pgeom->getRatio());
                        normal_velocity_patch.emplace(patch, bdry_normal_axis, d_comp_idx);
                    }
                    const double du_norm_dx_tan =
                        traction_stencil::get_normal_velocity_difference(*u_target_data,
                                                                         i,
                                                                         bdry_normal_axis,
                                                                         is_lower,
                                                                         d_comp_idx,
                                                                         d_bc_coefs[bdry_normal_axis],
                                                                         d_homogeneous_bc,
                                                                         *normal_velocity_patch,
                                                                         *domain,
                                                                         d_bdry_interp_type == "QUADRATIC",
                                                                         fill_time) /
                        dx[d_comp_idx];

                    if (!d_fluid_solver->muIsConstant())
                    {
                        // The fill can set faces beyond the viscosity data. There the viscosity is NaN, so that
                        // a ghost value computed from it is NaN.
                        if (mu_data_exists)
                        {
#if (NDIM == 2)
                            const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData();

#elif (NDIM == 3)
                            const int perp =
                                2 * (bdry_normal_axis + d_comp_idx) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                            const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData(perp);
#endif
                            mu = (mu_array_data.getBox().contains(i)) ? mu_array_data(i, 0) :
                                                                        std::numeric_limits<double>::quiet_NaN();
                        }
                    }

                    // Correct the boundary condition value.
                    alpha = 0.0;
                    beta = 1.0;
                    gamma = (is_lower ? -1.0 : +1.0) * (gamma / mu - du_norm_dx_tan);
                    break;
                }
                case PSEUDO_TRACTION: // mu*du_tan/dx_norm = g.
                {
                    if (!d_fluid_solver->muIsConstant())
                    {
                        // The fill can set faces beyond the viscosity data. There the viscosity is NaN, so that
                        // a ghost value computed from it is NaN.
                        if (mu_data_exists)
                        {
#if (NDIM == 2)
                            const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData();

#elif (NDIM == 3)
                            const int perp =
                                2 * (bdry_normal_axis + d_comp_idx) % 3; // 2 if {0,1}, 1 if {0,2} and 0 if {1,2}
                            const ArrayData<NDIM, double>& mu_array_data = mu_data->getArrayData(perp);
#endif
                            mu = (mu_array_data.getBox().contains(i)) ? mu_array_data(i, 0) :
                                                                        std::numeric_limits<double>::quiet_NaN();
                        }
                    }
                    alpha = 0.0;
                    beta = 1.0;
                    gamma = (is_lower ? -1.0 : +1.0) * (gamma / mu);
                    break;
                }
                default:
                {
                    TBOX_ERROR(
                        "INSVCStaggeredVelocityBcCoef::setBcCoefs(): unrecognized or "
                        "unsupported "
                        "traction boundary condition type: "
                        << enum_to_string<TractionBcType>(d_traction_bc_type) << "\n");
                }
                }
            }
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
INSVCStaggeredVelocityBcCoef::numberOfExtensionsFillable() const
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

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
