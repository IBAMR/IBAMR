// Filename: CopyToRootSchedule.cpp
// Created on 04 May 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#include <stddef.h>
#include <ostream>
#include <vector>

#include "BoxArray.h"
#include "GridGeometry.h"
#include "IntVector.h"
#include "PatchData.h"
#include "PatchDataFactory.h"
#include "PatchDescriptor.h"
#include "PatchLevel.h"
#include "ibtk/CopyToRootSchedule.h"
#include "ibtk/CopyToRootTransaction.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Pointer.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Schedule.h"
#include "tbox/Transaction.h"
#include "tbox/Utilities.h"

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
    : d_root_proc(root_proc),
      d_patch_level(patch_level),
      d_src_patch_data_idxs(1, src_patch_data_idx),
      d_root_patch_data(),
      d_schedule()
{
    commonClassCtor();
    return;
} // CopyToRootSchedule

CopyToRootSchedule::CopyToRootSchedule(const int root_proc,
                                       const Pointer<PatchLevel<NDIM> > patch_level,
                                       const std::vector<int>& src_patch_data_idxs)
    : d_root_proc(root_proc),
      d_patch_level(patch_level),
      d_src_patch_data_idxs(src_patch_data_idxs),
      d_root_patch_data(),
      d_schedule()
{
    commonClassCtor();
    return;
} // CopyToRootSchedule

CopyToRootSchedule::~CopyToRootSchedule()
{
    // intentionally blank
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

    d_root_patch_data.resize(num_vars, Pointer<PatchData<NDIM> >(NULL));
    if (SAMRAI_MPI::getRank() == d_root_proc)
    {
        for (unsigned int k = 0; k < num_vars; ++k)
        {
            Pointer<PatchDataFactory<NDIM> > pdat_factory =
                d_patch_level->getPatchDescriptor()->getPatchDataFactory(d_src_patch_data_idxs[k]);
            d_root_patch_data[k] = pdat_factory->allocate(domain_box);
        }
    }

    const int mpi_nodes = SAMRAI_MPI::getNodes();
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
