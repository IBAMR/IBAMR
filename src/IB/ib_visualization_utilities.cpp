// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/ib_visualization_utilities.h>

#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>

#include <ibamr/app_namespaces.h>

namespace IBAMR
{
void
regenerate_structure_viz(SAMRAI::tbox::Pointer<LSiloDataWriter> silo_writer,
                         LDataManager* data_manager,
                         const std::string& struct_name,
                         int ln)
{
    std::vector<std::string> struct_names = { struct_name };
    regenerate_structure_viz(silo_writer, data_manager, struct_names, ln);
}

void
regenerate_structure_viz(SAMRAI::tbox::Pointer<LSiloDataWriter> silo_writer,
                         LDataManager* data_manager,
                         int struct_id,
                         int ln)
{
    regenerate_structure_viz(silo_writer, data_manager, { struct_id }, ln);
}

void
regenerate_structure_viz(SAMRAI::tbox::Pointer<LSiloDataWriter> silo_writer,
                         LDataManager* data_manager,
                         const std::set<int>& struct_ids,
                         int ln)
{
    // Map structure ids to names
    std::vector<std::string> struct_names;
    for (const auto& struct_id : struct_ids)
        struct_names.push_back(data_manager->getLagrangianStructureName(struct_id, ln));
    regenerate_structure_viz(silo_writer, data_manager, struct_names, ln);
}

void
regenerate_structure_viz(Pointer<LSiloDataWriter> silo_writer,
                         LDataManager* data_manager,
                         const std::vector<std::string>& struct_names,
                         int ln)
{
    // Get Lagrangian index ranges.
    const size_t num_structs = struct_names.size();
    std::vector<std::pair<int, int> > lag_idx_ranges;
    for (const auto& struct_name : struct_names)
    {
        const int struct_id = data_manager->getLagrangianStructureID(struct_name, ln);
        lag_idx_ranges.push_back(data_manager->getLagrangianStructureIndexRange(struct_id, ln));
    }

    // We need to build the edge map so the silo_writer knows the updated connections. The edge_map needs to be on the
    // root process and contain all links on other processes. We first build a vector of all the links, then send them
    // to the root node.
    Pointer<LMesh> lag_mesh = data_manager->getLMesh(ln);
    const std::vector<LNode*>& local_nodes = lag_mesh->getLocalNodes();
    std::vector<std::vector<std::pair<int, int> > > struct_edge_vec(num_structs);

    // Note we assume the "master index" for a spring is the same as the Lagrangian index. If this were not true, we
    // would need edge_vec to consist of a triple of int values.
    for (const auto& node : local_nodes)
    {
        const int lag_idx = node->getLagrangianIndex();
        auto force_spec = node->getNodeDataItem<IBSpringForceSpec>();
        if (force_spec)
        {
            // We need to accumulate all links to regenerate the edge map.
            for (size_t struct_num = 0; struct_num < num_structs; ++struct_num)
            {
                // Is the Lagrangian point within the range for this structure?
                if (lag_idx >= lag_idx_ranges[struct_num].first && lag_idx < lag_idx_ranges[struct_num].second)
                {
                    const std::vector<int>& slave_idxs = force_spec->getSlaveNodeIndices();
                    for (const auto& slave_idx : slave_idxs)
                        struct_edge_vec[struct_num].push_back({ lag_idx, slave_idx });
                }
            }
        }
    }

    // Now rebuild the edge map for each structure
    for (size_t struct_num = 0; struct_num < num_structs; ++struct_num)
    {
        std::multimap<int, std::pair<int, int> > edge_map;
        build_edge_map(struct_edge_vec[struct_num], edge_map);
        // Register the map with the silo writer. This should only be done on the root process.
        if (IBTK_MPI::getRank() == 0)
            silo_writer->registerUnstructuredMesh(struct_names[struct_num] + "_mesh", edge_map, ln);
    }
}

void
build_edge_map(const std::vector<std::pair<int, int> >& edge_vec, std::multimap<int, std::pair<int, int> >& edge_map)
{
    const int rank = IBTK_MPI::getRank();
    const int nodes = IBTK_MPI::getNodes();
    // Allocate data on the root process
    std::vector<int> num_links_per_proc(nodes);
    std::vector<int> strides_per_proc(nodes);
    num_links_per_proc[rank] = edge_vec.size();
    MPI_Gather(
        &num_links_per_proc[rank], 1, MPI_INT, num_links_per_proc.data(), 1, MPI_INT, 0, IBTK_MPI::getCommunicator());
    std::vector<std::pair<int, int> > root_edge_vec;
    if (rank == 0)
    {
        const int tot_idxs = std::accumulate(num_links_per_proc.begin(), num_links_per_proc.end(), 0);
        root_edge_vec.resize(tot_idxs);
        strides_per_proc[0] = 0;
        for (int i = 1; i < nodes; ++i) strides_per_proc[i] = strides_per_proc[i - 1] + num_links_per_proc[i - 1];
    }

    // Now gather to the root process
    MPI_Gatherv(edge_vec.data(),
                num_links_per_proc[rank],
                MPI_2INT,
                root_edge_vec.data(),
                num_links_per_proc.data(),
                strides_per_proc.data(),
                MPI_2INT,
                0,
                IBTK_MPI::getCommunicator());

    // Now generate the map on the root process.
    if (rank == 0)
    {
        for (const auto& pair : root_edge_vec)
        {
            edge_map.insert(std::make_pair(pair.first, pair));
        }
    }
}
} // namespace IBAMR
