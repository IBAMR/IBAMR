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

#include "ibtk/CopyToRootSchedule.h"
#include "ibtk/CopyToRootTransaction.h"
#include "ibtk/IBTK_MPI.h"

#include "BoxArray.h"
#include "GridGeometry.h"
#include "IntVector.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "tbox/Pointer.h"
#include "tbox/Schedule.h"
#include "tbox/Transaction.h"

#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

CopyToRootSchedule::CopyToRootSchedule(const int root_proc,
                                       const Pointer<PatchLevel<NDIM> > patch_level,
                                       const int src_patch_data_idx)
    : d_root_proc(root_proc), d_patch_level(patch_level), d_src_patch_data_idxs(1, src_patch_data_idx)
{
    commonClassCtor();
    return;
} // CopyToRootSchedule

CopyToRootSchedule::CopyToRootSchedule(const int root_proc,
                                       const Pointer<PatchLevel<NDIM> > patch_level,
                                       std::vector<int> src_patch_data_idxs)
    : d_root_proc(root_proc), d_patch_level(patch_level), d_src_patch_data_idxs(std::move(src_patch_data_idxs))
{
    commonClassCtor();
    return;
} // CopyToRootSchedule

void
CopyToRootSchedule::communicate()
{
    d_schedule.communicate();
    return;
} // communicate

const std::vector<Pointer<PatchData<NDIM> > >&
CopyToRootSchedule::getRootPatchData() const
{
    return d_root_patch_data;
} // getRootPatchData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

void
CopyToRootSchedule::commonClassCtor()
{
    Pointer<GridGeometry<NDIM> > grid_geom = d_patch_level->getGridGeometry();
#if !defined(NDEBUG)
    TBOX_ASSERT(grid_geom->getDomainIsSingleBox());
#endif
    const Box<NDIM>& domain_box = grid_geom->getPhysicalDomain()[0];

    const size_t num_vars = d_src_patch_data_idxs.size();

    d_root_patch_data.resize(num_vars, Pointer<PatchData<NDIM> >(nullptr));
    if (IBTK_MPI::getRank() == d_root_proc)
    {
        for (unsigned int k = 0; k < num_vars; ++k)
        {
            Pointer<PatchDataFactory<NDIM> > pdat_factory =
                d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idxs[k]);
            d_root_patch_data[k] = pdat_factory->allocate(domain_box);
        }
    }

    const int mpi_nodes = IBTK_MPI::getNodes();
    for (int src_proc = 0; src_proc < mpi_nodes; ++src_proc)
    {
        for (unsigned int k = 0; k < num_vars; ++k)
        {
            d_schedule.appendTransaction(new CopyToRootTransaction(
                src_proc, d_root_proc, d_patch_level, d_src_patch_data_idxs[k], d_root_patch_data[k]));
        }
    }
    return;
} // commonClassCtor

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
