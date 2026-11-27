// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2026 by the IBAMR developers
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

#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleRT0Coarsen.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <Box.h>
#include <BoxGeometry.h>
#include <BoxOverlap.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <CoarsenOperator.h>
#include <Geometry.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <PatchLevel.h>
#include <Variable.h>

#include <algorithm>
#include <cmath>
#include <string>

#include <ibtk/app_namespaces.h>

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// The variable and operator name that identify a geometry that has received the transfer operators. A geometry has an
// operator for CONSERVATIVE_COARSEN of cell data, so looking it up for this variable never fails, and it returns the
// RegistrationSentinel if the geometry has one.
const std::string REGISTRATION_VARIABLE_NAME = "IBTK_TRANSFER_OPERATORS_REGISTRATION";
const std::string REGISTRATION_OPERATOR_NAME = "CONSERVATIVE_COARSEN";

// Coarsen operator that matches only that variable and name. The geometry owns it, so a geometry that is destroyed
// takes its record of the registration along.
class RegistrationSentinel : public CoarsenOperator<NDIM>
{
public:
    bool findCoarsenOperator(const Pointer<Variable<NDIM>>& var, const std::string& op_name) const override
    {
        return !var.isNull() && var->getName() == REGISTRATION_VARIABLE_NAME && op_name == REGISTRATION_OPERATOR_NAME;
    }

    const std::string& getOperatorName() const override
    {
        return d_name;
    }

    int getOperatorPriority() const override
    {
        return 0;
    }

    IntVector<NDIM> getStencilWidth() const override
    {
        return IntVector<NDIM>(0);
    }

    void coarsen(Patch<NDIM>& /*coarse*/,
                 const Patch<NDIM>& /*fine*/,
                 int /*dst_component*/,
                 int /*src_component*/,
                 const Box<NDIM>& /*coarse_box*/,
                 const IntVector<NDIM>& /*ratio*/) const override
    {
    }

private:
    const std::string d_name = "REGISTRATION_SENTINEL";
};

// Whether the geometry has received the transfer operators.
bool
has_transfer_operators(const xfer::Geometry<NDIM>& geometry)
{
    Pointer<Variable<NDIM>> variable = new CellVariable<NDIM, double>(REGISTRATION_VARIABLE_NAME);
    Pointer<CoarsenOperator<NDIM>> op = geometry.lookupCoarsenOperator(variable, REGISTRATION_OPERATOR_NAME);
    return dynamic_cast<const RegistrationSentinel*>(op.getPointer()) != nullptr;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

void
register_ibtk_transfer_operators(xfer::Geometry<NDIM>& geometry)
{
    if (has_transfer_operators(geometry)) return;

    geometry.addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
    geometry.addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen());
    geometry.addSpatialCoarsenOperator(new CartSideDoubleRT0Coarsen());
    geometry.addSpatialRefineOperator(new CartSideDoubleRT0Refine());
    // The geometry searches operators from the most recently added, so the sentinel is found before the operators that
    // the geometry has by default. Add it last.
    geometry.addSpatialCoarsenOperator(new RegistrationSentinel());
} // register_ibtk_transfer_operators

double
get_min_patch_dx(const PatchLevel<NDIM>& patch_level)
{
    double result = std::numeric_limits<double>::max();

    // Some processors might not have any patches so its easier to just quit
    // after one loop operation than to check
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level.getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double patch_dx_min = *std::min_element(patch_dx, patch_dx + NDIM);
        result = std::min(result, patch_dx_min);
        break; // all patches on the same level have the same dx values
    }

    result = IBTK_MPI::minReduction(result);

    return result;
} // get_min_patch_dx

void
copy_ghost_region(const PatchLevel<NDIM>& patch_level, const int dst_idx, const int src_idx)
{
    for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = patch_level.getPatch(p());
        Pointer<PatchData<NDIM>> dst_data = patch->getPatchData(dst_idx);
        Pointer<PatchData<NDIM>> src_data = patch->getPatchData(src_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(dst_data);
        TBOX_ASSERT(src_data);
#endif
        Pointer<PatchDescriptor<NDIM>> descriptor = patch->getPatchDescriptor();
        Pointer<BoxGeometry<NDIM>> dst_geometry =
            descriptor->getPatchDataFactory(dst_idx)->getBoxGeometry(patch->getBox());
        Pointer<BoxGeometry<NDIM>> src_geometry =
            descriptor->getPatchDataFactory(src_idx)->getBoxGeometry(patch->getBox());
        // Without overwriting the interior, the overlap is the ghost region of dst.
        Pointer<BoxOverlap<NDIM>> overlap = dst_geometry->calculateOverlap(
            *src_geometry, src_data->getGhostBox(), /*overwrite_interior*/ false, IntVector<NDIM>(0));
        if (!overlap)
        {
            TBOX_ERROR("copy_ghost_region():\n"
                       << "  patch data indices " << dst_idx << " and " << src_idx
                       << " have geometries that do not support an overlap" << std::endl);
        }
        dst_data->copy(*src_data, *overlap);
    }
    return;
} // copy_ghost_region
} // namespace IBTK
