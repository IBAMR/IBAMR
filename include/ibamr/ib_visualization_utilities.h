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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////
#ifndef included_IBAMR_visualization_utilities
#define included_IBAMR_visualization_utilities

#include <ibtk/LDataManager.h>
#include <ibtk/LSiloDataWriter.h>

namespace IBAMR
{
/*!
 * Inform the Silo data writer that the unstructured mesh corresponding to the provided list of structures has changed.
 * For each structure (either by name or id), we rebuild the edge_map by checking for IBSpringForceSpec objects at each
 * Lagrangian node. The structure names or ids need to be valid structures in the LDataManager object.
 *
 * @note This can potentially be expensive because only the root processor writes Silo data. All edges need to be
 * communicated to the root processor.
 *
 * @note This assumes the "master index" for a spring is the same as the Lagrangian index.
 *
 * @note This will not generate the correct meshes for structures with XSprings.
 */
//@{
void regenerate_structure_viz(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer,
                              IBTK::LDataManager* data_manager,
                              const std::string& struct_names,
                              int ln);
void regenerate_structure_viz(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer,
                              IBTK::LDataManager* data_manager,
                              int struct_ids,
                              int ln);
void regenerate_structure_viz(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer,
                              IBTK::LDataManager* data_manager,
                              const std::set<int>& struct_ids,
                              int ln);
void regenerate_structure_viz(SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer,
                              IBTK::LDataManager* data_manager,
                              const std::vector<std::string>& struct_names,
                              int ln);
//@}

/*!
 * Given a distributed vector of edges, build the edge map on the root process. Only the root process has a valid
 * edge_map.
 *
 * @param[in] edge_vec The list of edges on each process.
 * @param[out] edge_map The edge map between points constructed from @p edge_vec. Only valid on the root process.
 */
void build_edge_map(const std::vector<std::pair<int, int> >& edge_vec,
                    std::multimap<int, std::pair<int, int> >& edge_map);
} // namespace IBAMR
#endif
