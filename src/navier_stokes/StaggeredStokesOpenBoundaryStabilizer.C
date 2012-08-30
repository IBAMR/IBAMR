// Filename: StaggeredStokesOpenBoundaryStabilizer.C
// Created on 29 Aug 2012 by Boyce Griffith
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

#include "StaggeredStokesOpenBoundaryStabilizer.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesOpenBoundaryStabilizer::StaggeredStokesOpenBoundaryStabilizer(
    unsigned int comp_idx,
    RobinBcCoefStrategy<NDIM>* comp_bc_coef,
    Pointer<Database> input_db,
    const INSStaggeredHierarchyIntegrator* fluid_solver)
    : d_alpha(std::numeric_limits<double>::quiet_NaN()),
      d_beta(std::numeric_limits<double>::quiet_NaN()),
      d_comp_idx(comp_idx),
      d_comp_bc_coef(comp_bc_coef),
      d_open_bdry(false),
      d_inflow_bdry(false),
      d_outflow_bdry(false),
      d_fluid_solver(fluid_solver)
{
    if (input_db)
    {
        if (input_db->keyExists("alpha")) d_alpha = input_db->getDouble("alpha");
        if (input_db->keyExists("beta")) d_beta = input_db->getDouble("beta");
        for (int location_index = 0; location_index < 2*NDIM; ++location_index)
        {
            std::ostringstream stabilization_type_stream;
            stabilization_type_stream << "stabilization_type_" << location_index;
            const std::string stabilization_type_key = stabilization_type_stream.str();
            if (input_db->keyExists(stabilization_type_key))
            {
                const std::string stabilization_type = input_db->getString(stabilization_type_key);
                if (stabilization_type == "INFLOW")
                {
                    d_open_bdry   [location_index] = true;
                    d_inflow_bdry [location_index] = true;
                    d_outflow_bdry[location_index] = false;
                }
                else if (stabilization_type == "OUTFLOW")
                {
                    d_open_bdry   [location_index] = true;
                    d_inflow_bdry [location_index] = false;
                    d_outflow_bdry[location_index] = true;
                }
                else if (stabilization_type != "NONE")
                {
                    TBOX_ERROR("StaggeredStokesOpenBoundaryStabilizer::StaggeredStokesOpenBoundaryStabilizer():\n"
                               << "  unsupported stabilization type: ``" << stabilization_type << "''\n"
                               << "  supported values are: ``INFLOW'', ``OUTFLOW'', or ``NONE''\n");
                }
            }
        }
    }
    return;
}// StaggeredStokesOpenBoundaryStabilizer

StaggeredStokesOpenBoundaryStabilizer::~StaggeredStokesOpenBoundaryStabilizer()
{
    // intentionally blank
    return;
}// ~StaggeredStokesOpenBoundaryStabilizer

void
StaggeredStokesOpenBoundaryStabilizer::setTargetPatchDataIndex(
    int target_idx)
{
    ExtendedRobinBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_comp_bc_coef);
    if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    return;
}// setTargetPatchDataIndex

void
StaggeredStokesOpenBoundaryStabilizer::setHomogeneousBc(
    bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    ExtendedRobinBcCoefStrategy* p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_comp_bc_coef);
    if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    return;
}// setHomogeneousBc

void
StaggeredStokesOpenBoundaryStabilizer::setBcCoefs(
    Pointer<ArrayData<NDIM,double> >& acoef_data,
    Pointer<ArrayData<NDIM,double> >& bcoef_data,
    Pointer<ArrayData<NDIM,double> >& gcoef_data,
    const Pointer<Variable<NDIM> >& variable,
    const Patch<NDIM>& patch,
    const BoundaryBox<NDIM>& bdry_box,
    double fill_time) const
{
    // Set the unmodified velocity bc coefs.
    d_comp_bc_coef->setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);

    // If we are not setting inhomogeneous coefficients, at an open boundary, or
    // operating on the correct velocity component, then there is nothing else
    // to do.
    if (!gcoef_data || d_homogeneous_bc) return;
    const unsigned int location_index = bdry_box.getLocationIndex();
    if (!d_open_bdry[location_index]) return;
    const unsigned int bdry_normal_axis = location_index/2;
    if (bdry_normal_axis != d_comp_idx) return;

    // Where appropriate, update normal traction boundary conditions to penalize
    // flow reversal.
    const bool is_lower = location_index%2 == 0;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(acoef_data);
    TBOX_ASSERT(bcoef_data);
    TBOX_ASSERT(gcoef_data);
#endif
    Box<NDIM> bc_coef_box = acoef_data->getBox();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(bc_coef_box == acoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(bc_coef_box == gcoef_data->getBox());
#endif
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getCurrentContext());
    const int u_new_idx     = var_db->mapVariableAndContextToIndex(d_fluid_solver->getVelocityVariable(), d_fluid_solver->getNewContext());
    Pointer<SideData<NDIM,double> > u_current_data = patch.getPatchData(u_current_idx);
    Pointer<SideData<NDIM,double> > u_new_data     = patch.getPatchData(u_new_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(u_current_data);
    TBOX_ASSERT(u_new_data);
#endif
    const Box<NDIM>& ghost_box = u_current_data->getGhostBox() * u_new_data->getGhostBox();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (d != bdry_normal_axis)
        {
            bc_coef_box.lower(d) = std::max(bc_coef_box.lower(d), ghost_box.lower(d));
            bc_coef_box.upper(d) = std::min(bc_coef_box.upper(d), ghost_box.upper(d));
        }
    }
    const int cycle_num = d_fluid_solver->getCurrentCycleNumber();
    for (Box<NDIM>::Iterator it(bc_coef_box); it; it++)
    {
        const Index<NDIM>& i = it();
        const double& alpha = (*acoef_data)(i,0);
        const double& beta  = (*bcoef_data)(i,0);
        double& gamma = (*gcoef_data)(i,0);
        const bool traction_bc = MathUtilities<double>::equalEps(beta ,1.0);
#ifdef DEBUG_CHECK_ASSERTIONS
        const bool velocity_bc = MathUtilities<double>::equalEps(alpha,1.0);
        TBOX_ASSERT((velocity_bc || traction_bc) && !(velocity_bc && traction_bc));
#endif
        if (traction_bc)
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_comp_idx == bdry_normal_axis);
#endif
            const SideIndex<NDIM> i_s(i, bdry_normal_axis, SideIndex<NDIM>::Lower);
            const double u_n = (is_lower ? -1.0 : 1.0) * (cycle_num == 0 ? (*u_current_data)(i_s) : 0.5*((*u_current_data)(i_s) + (*u_new_data)(i_s)));
            if (d_inflow_bdry[location_index] && u_n > 0.0)
            {
                gamma -= d_alpha*pow(std::abs(u_n),d_beta);
            }
            else if (d_outflow_bdry[location_index] && u_n < 0.0)
            {
                gamma += d_alpha*pow(std::abs(u_n),d_beta);
            }
        }
    }
    return;
}// setBcCoefs

SAMRAI::hier::IntVector<NDIM>
StaggeredStokesOpenBoundaryStabilizer::numberOfExtensionsFillable() const
{
    return d_comp_bc_coef->numberOfExtensionsFillable();
}// numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
