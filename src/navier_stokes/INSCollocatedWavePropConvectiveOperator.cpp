// Filename: INSCollocatedWavePropConvectiveOperator.cpp
// Created on 23 Feb 2018 by Aaron Barrett
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include <utility>

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "FaceData.h"
#include "IBAMR_config.h"
#include "Index.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSCollocatedWavePropConvectiveOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

extern "C"
{
#if (NDIM == 2)
    void adv_diff_wp_convective_op2d_(const double*,
                                      const int&,
                                      const double*,
                                      const double*,
                                      const int&,
                                      const double*,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const int&,
                                      const double*,
                                      const int&);
#endif
#if (NDIM == 3)
    void adv_diff_wp_convective_op3d_(const double* q_data,
                                      const int& q_gcw,
                                      const double* u_data_0,
                                      const double* u_data_1,
                                      const double* u_data_2,
                                      const int& u_gcw,
                                      const double* r_data,
                                      const int& r_gcw,
                                      const int& depth,
                                      const int& ilower0,
                                      const int& ilower1,
                                      const int& ilower2,
                                      const int& iupper0,
                                      const int& iupper1,
                                      const int& iupper2,
                                      const double* dx,
                                      const int& k);
#endif
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

INSCollocatedWavePropConvectiveOperator::INSCollocatedWavePropConvectiveOperator(
    std::string object_name,
    Pointer<Database> input_db,
    const ConvectiveDifferencingType difference_form,
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs)
    : ConvectiveOperator(std::move(object_name), difference_form), d_bc_coefs(std::move(bc_coefs))
{
    if (d_difference_form != ADVECTIVE /* && d_difference_form != CONSERVATIVE && d_difference_form != SKEW_SYMMETRIC*/)
    {
        TBOX_ERROR("INSCollocatedWavePropConvectiveOperator::INSCollocatedWavePropConvectiveOperator():\n"
                   << "  unsupported differencing form: "
                   << enum_to_string<ConvectiveDifferencingType>(d_difference_form) << " \n"
                   << "  valid choices are: ADVECTIVE\n");
    }

    if (input_db)
    {
        if (input_db->keyExists("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
    }

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = var_db->getContext("INSCollocatedWavePropConvectiveOperator::CONTEXT");

    const std::string U_var_name = "INSCollocatedWavePropConvectiveOperator::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var)
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
    else
    {
        d_U_var = new CellVariable<NDIM, double>(U_var_name, NDIM);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, IntVector<NDIM>(d_k + 1));
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    return;
} // INSCollocatedWavePropConvectiveOperator

INSCollocatedWavePropConvectiveOperator::~INSCollocatedWavePropConvectiveOperator()
{
    deallocateOperatorState();
    return;
} // ~INSCollocatedWavePropConvectiveOperator

void
INSCollocatedWavePropConvectiveOperator::applyConvectiveOperator(const int U_idx, const int N_idx)
{
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSCollocatedWavePropConvectiveOperator::applyConvectiveOperator():\n"
                   << "  operator must be initialized prior to call to applyConvectiveOperator\n");
    }
#endif

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_U_scratch_idx);
    }

    // Fill ghost cell values for all components.
    HierarchyMathOps hier_math_ops("HierarchyMathOps", d_hierarchy);
    static const bool homogeneous_bc = false;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(1);
    transaction_comps[0] = InterpolationTransactionComponent(d_U_scratch_idx,
                                                             U_idx,
                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                             true,
                                                             "CONSERVATIVE_COARSEN",
                                                             d_bdry_extrap_type,
                                                             false,
                                                             d_bc_coefs);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the convective derivative.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());

            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();

            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();
            const IntVector<NDIM>& patch_upper = patch_box.upper();

            Pointer<CellData<NDIM, double> > N_data = patch->getPatchData(N_idx);
            const IntVector<NDIM> N_gcw = N_data->getGhostCellWidth();
            Pointer<CellData<NDIM, double> > U_data = patch->getPatchData(d_U_scratch_idx);
            const IntVector<NDIM> U_gcw = U_data->getGhostCellWidth();
            Pointer<FaceData<NDIM, double> > U_sp_data = patch->getPatchData(d_u_idx);
            const IntVector<NDIM> U_sp_gcw = U_sp_data->getGhostCellWidth();

            // Do differencing
            for (int d = 0; d < NDIM; ++d)
            {
#if (NDIM == 2)
                adv_diff_wp_convective_op2d_(U_data->getPointer(d),
                                             U_gcw.max(),
                                             U_sp_data->getPointer(0),
                                             U_sp_data->getPointer(1),
                                             U_sp_gcw.max(),
                                             N_data->getPointer(d),
                                             N_gcw.max(),
                                             1,
                                             patch_lower(0),
                                             patch_lower(1),
                                             patch_upper(0),
                                             patch_upper(1),
                                             dx,
                                             d_k);
#endif
#if (NDIM == 3)
                adv_diff_wp_convective_op3d_(U_data->getPointer(d),
                                             U_gcw.max(),
                                             U_sp_data->getPointer(0),
                                             U_sp_data->getPointer(1),
                                             U_sp_data->getPointer(2),
                                             U_sp_gcw.max(),
                                             N_data->getPointer(d),
                                             N_gcw.max(),
                                             1,
                                             patch_lower(0),
                                             patch_lower(1),
                                             patch_lower(2),
                                             patch_upper(0),
                                             patch_upper(1),
                                             patch_upper(2),
                                             dx,
                                             d_k);
#endif
            }
        }
    }
    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_U_scratch_idx);
    }

    return;
} // applyConvectiveOperator

void
INSCollocatedWavePropConvectiveOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                                 const SAMRAIVectorReal<NDIM, double>& out)
{
    if (d_is_initialized) deallocateOperatorState();

    // Get the hierarchy configuration.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();
#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#else
    NULL_USE(out);
#endif

    // Setup the interpolation transaction information.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.resize(1);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_U_scratch_idx,
                                                               in.getComponentDescriptorIndex(0),
                                                               "CONSERVATIVE_LINEAR_REFINE",
                                                               false,
                                                               "CONSERVATIVE_COARSEN",
                                                               d_bdry_extrap_type,
                                                               false,
                                                               d_bc_coefs);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy);

    d_is_initialized = true;

    return;
} // initializeOperatorState

void
INSCollocatedWavePropConvectiveOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    // Deallocate the communications operators and BC helpers.
    d_hier_bdry_fill.setNull();

    d_is_initialized = false;

    return;
} // deallocateOperatorState
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
