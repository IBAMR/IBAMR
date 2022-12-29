// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2022 by the IBAMR developers
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

#include "ibamr/RelaxationLSBcCoefs.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CellData.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"

#include <utility>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RelaxationLSBcCoefs::RelaxationLSBcCoefs(std::string object_name) : d_object_name(std::move(object_name))
{
    // intentionally blank
    return;
} // RelaxationLSBcCoefs

void
RelaxationLSBcCoefs::setLSPatchDataIndex(int phi_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(phi_idx >= 0);
#endif
    d_phi_idx = phi_idx;
    return;
} // setLSPatchDataIndex

void
RelaxationLSBcCoefs::resetLSPatchDataIndex()
{
    d_phi_idx = invalid_index;
    return;
} // resetLSPatchDataIndex

void
RelaxationLSBcCoefs::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                const Pointer<Variable<NDIM> >& /*variable*/,
                                const Patch<NDIM>& patch,
                                const BoundaryBox<NDIM>& bdry_box,
                                double /*fill_time*/) const
{
    Pointer<CellData<NDIM, double> > phi_data = patch.getPatchData(d_phi_idx);

    const int location_index = bdry_box.getLocationIndex();
    const int axis = location_index / 2;
    const bool is_upper = (location_index % 2 == 1);
#if !defined(NDEBUG)
    TBOX_ASSERT(!acoef_data.isNull());
#endif
    const Box<NDIM>& bc_coef_box = acoef_data->getBox();
#if !defined(NDEBUG)
    TBOX_ASSERT(bcoef_data.isNull() || bc_coef_box == bcoef_data->getBox());
    TBOX_ASSERT(gcoef_data.isNull() || bc_coef_box == gcoef_data->getBox());
#endif
    for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
    {
        const hier::Index<NDIM>& i = bc();
        hier::Index<NDIM> i_bdry = bc();
        hier::Index<NDIM> i_intr = bc();

        if (is_upper)
        {
            i_bdry(axis) -= 1;
            i_intr(axis) -= 2;
        }
        else
        {
            i_intr(axis) += 1;
        }
        double dummy;
        double& a = (!acoef_data.isNull() ? (*acoef_data)(i, 0) : dummy);
        double& b = (!bcoef_data.isNull() ? (*bcoef_data)(i, 0) : dummy);
        double& g = (!gcoef_data.isNull() ? (*gcoef_data)(i, 0) : dummy);

        a = 1.0;
        b = 0.0;
        g = 1.5 * (*phi_data)(i_bdry, 0) - 0.5 * (*phi_data)(i_intr, 0);
    }
    return;
} // setBcCoefs

IntVector<NDIM>
RelaxationLSBcCoefs::numberOfExtensionsFillable() const
{
    return 128;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
