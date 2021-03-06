// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include "ibtk/IBTK_MPI.h"
#include "ibtk/MergingLoadBalancer.h"
#include "ibtk/box_utilities.h"

#include "Box.h"
#include "PatchHierarchy.h"
#include "ProcessorMapping.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////
void
MergingLoadBalancer::loadBalanceBoxes(hier::BoxArray<NDIM>& out_boxes,
                                      hier::ProcessorMapping& mapping,
                                      const hier::BoxList<NDIM>& in_boxes,
                                      const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
                                      int level_number,
                                      const hier::BoxArray<NDIM>& physical_domain,
                                      const hier::IntVector<NDIM>& ratio_to_hierarchy_level_zero,
                                      const hier::IntVector<NDIM>& min_size,
                                      const hier::IntVector<NDIM>& max_size,
                                      const hier::IntVector<NDIM>& cut_factor,
                                      const hier::IntVector<NDIM>& bad_interval) const
{
    mesh::LoadBalancer<NDIM>::loadBalanceBoxes(out_boxes,
                                               mapping,
                                               in_boxes,
                                               hierarchy,
                                               level_number,
                                               physical_domain,
                                               ratio_to_hierarchy_level_zero,
                                               min_size,
                                               max_size,
                                               cut_factor,
                                               bad_interval);

    // pairs of processors and boxes
    std::vector<std::pair<int, hier::Box<NDIM> > > new_boxes;

    const int n_nodes = IBTK_MPI::getNodes();
    for (int r = 0; r < n_nodes; ++r)
    {
        // get all boxes on processor r.
        std::vector<hier::Box<NDIM> > boxes;
        for (int i = 0; i < out_boxes.size(); ++i)
            if (mapping.getProcessorAssignment(i) == r) boxes.push_back(out_boxes[i]);

        // populate new_boxes with merged boxes.
        const auto merged_boxes = IBTK::merge_boxes_by_longest_edge(boxes);
        for (const hier::Box<NDIM>& box : merged_boxes) new_boxes.emplace_back(r, box);
    }

    // Overwrite what the parent class did with the merged boxes.
    mapping.setMappingSize(new_boxes.size());
    out_boxes.resizeBoxArray(new_boxes.size());
    for (unsigned int i = 0; i < new_boxes.size(); ++i)
    {
        mapping.setProcessorAssignment(i, new_boxes[i].first);
        out_boxes[i] = new_boxes[i].second;
    }
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
