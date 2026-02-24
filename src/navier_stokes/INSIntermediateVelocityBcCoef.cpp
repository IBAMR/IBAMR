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

#include <ibamr/INSIntermediateVelocityBcCoef.h>

#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/samrai_compatibility_names.h>

#include <SAMRAIArrayData.h>
#include <SAMRAIBoundaryBox.h>
#include <SAMRAIIntVector.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPointer.h>
#include <SAMRAIRobinBcCoefStrategy.h>
#include <SAMRAIVariable.h>

#include <algorithm>
#include <limits>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Variable;
template <int DIM>
class Patch;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSIntermediateVelocityBcCoef::INSIntermediateVelocityBcCoef(const int comp_idx,
                                                             const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs,
                                                             const bool homogeneous_bc)
    : d_comp_idx(comp_idx), d_bc_coefs(NDIM, nullptr)
{
    setPhysicalBcCoefs(bc_coefs);
    setHomogeneousBc(homogeneous_bc);
    return;
} // INSIntermediateVelocityBcCoef

void
INSIntermediateVelocityBcCoef::setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
INSIntermediateVelocityBcCoef::setSolutionTime(double /*solution_time*/)
{
    // intentionally blank
    return;
} // setSolutionTime

void
INSIntermediateVelocityBcCoef::setTimeInterval(double /*current_time*/, double /*new_time*/)
{
    // intentionally blank
    return;
} // setTimeInterval

void
INSIntermediateVelocityBcCoef::setTargetPatchDataIndex(int target_idx)
{
    ExtendedRobinBcCoefStrategy::setTargetPatchDataIndex(target_idx);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setTargetPatchDataIndex(target_idx);
    }
    return;
} // setTargetPatchDataIndex

void
INSIntermediateVelocityBcCoef::clearTargetPatchDataIndex()
{
    ExtendedRobinBcCoefStrategy::clearTargetPatchDataIndex();
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->clearTargetPatchDataIndex();
    }
    return;
} // clearTargetPatchDataIndex

void
INSIntermediateVelocityBcCoef::setHomogeneousBc(bool homogeneous_bc)
{
    ExtendedRobinBcCoefStrategy::setHomogeneousBc(homogeneous_bc);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_comp_bc_coef = dynamic_cast<ExtendedRobinBcCoefStrategy*>(d_bc_coefs[d]);
        if (p_comp_bc_coef) p_comp_bc_coef->setHomogeneousBc(homogeneous_bc);
    }
    return;
} // setHomogeneousBc

void
INSIntermediateVelocityBcCoef::setBcCoefs(SAMRAIPointer<SAMRAIArrayData<double>>& acoef_data,
                                          SAMRAIPointer<SAMRAIArrayData<double>>& bcoef_data,
                                          SAMRAIPointer<SAMRAIArrayData<double>>& gcoef_data,
                                          const SAMRAIPointer<SAMRAIVariable>& variable,
                                          const SAMRAIPatch& patch,
                                          const SAMRAIBoundaryBox& bdry_box,
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
    if (d_homogeneous_bc && gcoef_data) gcoef_data->fillAll(0.0);
    return;
} // setBcCoefs

SAMRAIIntVector
INSIntermediateVelocityBcCoef::numberOfExtensionsFillable() const
{
#if !defined(NDEBUG)
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        TBOX_ASSERT(d_bc_coefs[d]);
    }
#endif
    SAMRAIIntVector ret_val(std::numeric_limits<int>::max());
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ret_val = SAMRAIIntVector::min(ret_val, d_bc_coefs[d]->numberOfExtensionsFillable());
    }
    return ret_val;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
