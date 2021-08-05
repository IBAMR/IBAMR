// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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

#include "ibamr/FESurfaceDistanceEvaluator.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellIndex.h"
#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "ProcessorMapping.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefineSchedule.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/MathUtilities.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

#include "libmesh/boundary_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/point.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/type_vector.h"
#include "libmesh/variant_filter_iterator.h"
#include "libmesh/vector_value.h"

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <set>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "ibamr/app_namespaces.h"

namespace libMesh
{
class Node;
} // namespace libMesh
/////////////////////////////// NAMESPACE ////////////////////////////////////

// FORTRAN ROUTINES
#if (NDIM == 2)
#define SIGN_SWEEP_FC IBAMR_FC_FUNC(signsweep2d, SIGNSWEEP2D)
#endif

#if (NDIM == 3)
#define SIGN_SWEEP_FC IBAMR_FC_FUNC(signsweep3d, SIGNSWEEP3D)
#endif

extern "C"
{
    void SIGN_SWEEP_FC(double* U,
                       const int& U_gcw,
                       const int& ilower0,
                       const int& iupper0,
                       const int& ilower1,
                       const int& iupper1,
#if (NDIM == 3)
                       const int& ilower2,
                       const int& iupper2,
#endif
                       const double& large_dist,
                       int& n_updates);
}

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
// Timers.
static Pointer<Timer> t_collectNeighboringPatchElements;
static Pointer<Timer> t_buildIntersectionMap;

template <class T1, class T2>
inline T1
make_edge(T2 n0, T2 n1)
{
    T1 e0(n1(0) - n0(0), n1(1) - n0(1), n1(2) - n0(2));
    return (e0);
} // make_edge
} // namespace

const double FESurfaceDistanceEvaluator::s_large_distance = 1234567.0;

/////////////////////////////// PUBLIC //////////////////////////////////////
FESurfaceDistanceEvaluator::FESurfaceDistanceEvaluator(std::string object_name,
                                                       Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                       const libMesh::Mesh& mesh,
                                                       const BoundaryMesh& bdry_mesh,
                                                       const int gcw,
                                                       bool use_extracted_bdry_mesh)
    : d_object_name(std::move(object_name)),
      d_patch_hierarchy(patch_hierarchy),
      d_mesh(mesh),
      d_bdry_mesh(bdry_mesh),
      d_gcw(gcw),
      d_use_vol_extracted_bdry_mesh(use_extracted_bdry_mesh)
{
// The only supported element type for this class
#if (NDIM == 2)
    d_supported_elem_type = EDGE2;
#endif
#if (NDIM == 3)
    d_supported_elem_type = TRI3;
#endif

    // Note that this class is specialized to work on a boundary mesh with
    //
    //    dim = NDIM - 1
    //
    // derived from a volumetric mesh.
    if (d_use_vol_extracted_bdry_mesh && d_bdry_mesh.mesh_dimension() != NDIM - 1)
    {
        TBOX_ERROR(
            "FESurfaceDistanceEvaluator presently requires a boundary mesh "
            "with dim = NDIM - 1 to be registered");
    }

    // Set up timers
    IBTK_DO_ONCE(t_collectNeighboringPatchElements = TimerManager::getManager()->getTimer(
                     "FESurfaceDistanceEvaluator::collectNeighboringPatchElements()", true);
                 t_buildIntersectionMap =
                     TimerManager::getManager()->getTimer("FESurfaceDistanceEvaluator::buildIntersectionMap()", true););

    return;
} // FESurfaceDistanceEvaluator

void
FESurfaceDistanceEvaluator::mapIntersections()
{
    // Loop over each cell and a set of triangles and map cell-triangle
    // intersections.
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();

    IBTK_TIMER_START(t_collectNeighboringPatchElements);
    collectNeighboringPatchElements(finest_ln);
    IBTK_TIMER_STOP(t_collectNeighboringPatchElements);
    IBTK_TIMER_START(t_buildIntersectionMap);

    // Clear out the data structure.
    d_cell_elem_neighbor_map.clear();

    // Loop over patches on finest level, while keeping track of the local patch
    // indexing.
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(finest_ln);
    int local_patch_num = 0;

    // Desired ghost cell width.
    IntVector<NDIM> ghost_width = d_gcw;

    // Compute a bounding box for the entire structure, relying on LibMesh
    // parallel decomposition.
    MeshBase::const_element_iterator el_it = d_mesh.active_local_elements_begin();
    MeshBase::const_element_iterator el_end = d_mesh.active_local_elements_end();
    if (d_use_vol_extracted_bdry_mesh)
    {
        el_it = d_bdry_mesh.active_local_elements_begin();
        el_end = d_bdry_mesh.active_local_elements_end();
    }
    IBTK::Vector3d elem_bl, elem_tr;
    elem_bl << std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max();
    elem_tr << std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min();
    for (; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;

        // Get the coordinates of the nodes.
        const libMesh::Point& n0 = elem->point(0);
        const libMesh::Point& n1 = elem->point(1);

#if (NDIM == 3)
        const libMesh::Point& n2 = elem->point(2);
#endif

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            elem_bl[d] = std::min(elem_bl[d], n0(d));
            elem_bl[d] = std::min(elem_bl[d], n1(d));
            elem_tr[d] = std::max(elem_tr[d], n0(d));
            elem_tr[d] = std::max(elem_tr[d], n1(d));

#if (NDIM == 3)
            elem_bl[d] = std::min(elem_bl[d], n2(d));
            elem_tr[d] = std::max(elem_tr[d], n2(d));
#endif
        }
    }
    IBTK_MPI::minReduction(elem_bl.data(), 3);
    IBTK_MPI::maxReduction(elem_tr.data(), 3);

    // Structure bounding box, taking into account ghost cell width.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const double* const dx0 = grid_geom->getDx();
    const IntVector<NDIM>& level_ratio = level->getRatio();
    double level_dx[NDIM] = { 0.0 };
    for (int d = 0; d < NDIM; ++d) level_dx[d] = dx0[d] / level_ratio(d);

    IBTK::Vector3d struct_bl, struct_tr, large_struct_bl, large_struct_tr;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        struct_bl[d] = elem_bl[d] - (ghost_width[d] * level_dx[d]);
        struct_tr[d] = elem_tr[d] + (ghost_width[d] * level_dx[d]);
        large_struct_bl[d] = elem_bl[d] - ((ghost_width[d] + 1) * level_dx[d]);
        large_struct_tr[d] = elem_tr[d] + ((ghost_width[d] + 1) * level_dx[d]);
    }
    Box<NDIM> struct_box(IndexUtilities::getCellIndex(struct_bl.data(), grid_geom, level_ratio),
                         IndexUtilities::getCellIndex(struct_tr.data(), grid_geom, level_ratio));

    // Create a bounding box with some additional ghost cell width than the
    // given one to eliminate some corner cases that may show up in the sign
    // update sweeping algorithm.
    d_large_struct_box = Box<NDIM>(IndexUtilities::getCellIndex(large_struct_bl.data(), grid_geom, level_ratio),
                                   IndexUtilities::getCellIndex(large_struct_tr.data(), grid_geom, level_ratio));

    // Map the neighbor intersections.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        // The relevant collection of elements.
        const std::vector<Elem*>& patch_elems = d_active_neighbor_patch_bdry_elem_map[local_patch_num];
        const size_t num_active_patch_elems = patch_elems.size();
        if (!num_active_patch_elems) continue;

        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_X_lower = patch_geom->getXLower();
        const SAMRAI::hier::Index<NDIM>& patch_lower_index = patch_box.lower();
        const double* const patch_dx = patch_geom->getDx();

        // If the patch box doesn't intersect the structure box, no need to do
        // computations.
        if (!patch_box.intersects(struct_box)) continue;

        // Loop over cells
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            // Get the coordinates/dimensions a box grown out of the cell
            // center of ci by the ghost cell width. This will ensure that we
            // capture not only the elements intersecting the cell, but also
            // the elements intersecting the neighbor.
            CellIndex<NDIM> ci(it());
            IBTK::Vector3d r_bl, r_tr;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                r_bl[d] = patch_X_lower[d] +
                          patch_dx[d] * (static_cast<double>(ci(d) - ghost_width[d] - patch_lower_index(d)));
                r_tr[d] = patch_X_lower[d] +
                          patch_dx[d] * (static_cast<double>(ci(d) + 1 + ghost_width[d] - patch_lower_index(d)));
            }

            // Do a simple bounding box check first, before carrying out expensive
            // checkIntersection routines.
            Box<NDIM> ghost_box(IndexUtilities::getCellIndex(r_bl.data(), grid_geom, level_ratio),
                                IndexUtilities::getCellIndex(r_tr.data(), grid_geom, level_ratio));
            if (!ghost_box.intersects(struct_box)) continue;

                // Prepare the required vectors.
#if (NDIM == 2)
            IBTK::Vector3d r_br, r_tl;
            r_br(1) = r_bl(1);
            r_br(0) = r_tr(0);
            r_tl(1) = r_tr(1);
            r_tl(0) = r_bl(0);
#endif
#if (NDIM == 3)
            // Get the coordinates of the cell center.
            IBTK::Vector3d grown_box_center, grown_box_half_dx;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                grown_box_center[d] =
                    patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_index(d)) + 0.5);
                grown_box_half_dx[d] = (0.5 + ghost_width[d]) * patch_dx[d];
            }
#endif

            // Loop over elements in the patch.
            for (const auto& elem : patch_elems)
            {
                // Get the coordinates of the nodes.
                const libMesh::Point& n0 = elem->point(0);
                const libMesh::Point& n1 = elem->point(1);

#if (NDIM == 3)
                const libMesh::Point& n2 = elem->point(2);
#endif

                // Intersection detection routines.
#if (NDIM == 2)
                const bool found_intersection = checkIntersection2D(r_bl, r_tr, r_br, r_tl, n0, n1);
#endif
#if (NDIM == 3)
                // In 3D, it is more convenient to pass vertices as IBTK::Vector3d.
                IBTK::Vector3d vert0, vert1, vert2;
                vert0 << n0(0), n0(1), n0(2);
                vert1 << n1(0), n1(1), n1(2);
                vert2 << n2(0), n2(1), n2(2);
                const bool found_intersection =
                    checkIntersection3D(grown_box_center, grown_box_half_dx, vert0, vert1, vert2);
#endif
                if (found_intersection)
                {
                    d_cell_elem_neighbor_map[ci].insert(elem);
                }
            }
        }
    }

    IBTK_TIMER_STOP(t_buildIntersectionMap);

    return;
} // mapIntersections

const std::map<CellIndex<NDIM>, std::set<Elem*>, CellIndexFortranOrder>&
FESurfaceDistanceEvaluator::getNeighborIntersectionsMap()
{
    return d_cell_elem_neighbor_map;
} // getNeighborIntersectionsMap

void
FESurfaceDistanceEvaluator::computeSignedDistance(int n_idx, int d_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(n_idx >= 0);
    TBOX_ASSERT(d_idx >= 0);
#endif

    // Loop over patches on finest level.
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const SAMRAI::hier::Index<NDIM>& patch_lower_index = patch_box.lower();
        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* patch_X_lower = patch_geom->getXLower();
        const double* const patch_dx = patch_geom->getDx();
        Pointer<CellData<NDIM, double> > n_data = patch->getPatchData(n_idx);
        Pointer<CellData<NDIM, double> > d_data = patch->getPatchData(d_idx);

        // Note that we only work with cells that satisfy the intersecting criteria.
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            CellIndex<NDIM> ci(it());
            const bool contains_key = (d_cell_elem_neighbor_map.find(ci) != d_cell_elem_neighbor_map.end());
            if (contains_key)
            {
                std::set<Elem*>& elem_set = d_cell_elem_neighbor_map[ci];
                const int num_elements = static_cast<int>(elem_set.size());
                (*n_data)(ci) = num_elements;

                // Loop over the cutting elements and find the minimum distance.
                IBTK::VectorNd P;
                for (int d = 0; d < NDIM; ++d)
                {
                    P[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_index(d)) + 0.5);
                }
                double min_dist = std::numeric_limits<double>::max();

                // Create a pair to take care of projection and pseudo-normal for cells
                // equidistant to multiple elements.
                std::vector<std::pair<IBTK::VectorNd, IBTK::VectorNd> > vec_equidistant_pair;

                for (auto& elem : elem_set)
                {
                    IBTK::VectorNd v, w, proj;
                    double dist = std::numeric_limits<double>::max();
#if (NDIM == 2)
                    // Get the nodes.
                    const libMesh::Point& n0 = elem->point(0);
                    const libMesh::Point& n1 = elem->point(1);

                    v << n0(0), n0(1);
                    w << n1(0), n1(1);

                    const double L2 = (v - w).squaredNorm();
                    if (MathUtilities<double>::equalEps(L2, 0.0))
                    {
                        // Special case where line element collapses to a
                        // point. Shouldn't happen.
                        proj = v;
                        dist = (proj - P).norm();
                    }
                    else
                    {
                        // Parameterize and project Note that this will
                        // take care of the edge case where the projection
                        // does not fall on the line.
                        const double t = std::max(0.0, std::min(1.0, (P - v).dot(w - v) / L2));
                        proj = v + t * (w - v);
                        dist = (P - proj).norm();
                    }

                    // If the distance is the same as the minimal
                    // distance, then the cell is equidistant to multiple
                    // elements, so add it to the set.
                    if (MathUtilities<double>::equalEps(dist, min_dist))
                    {
                        vec_equidistant_pair.push_back(std::make_pair(proj, d_elem_face_normal[elem]));
                    }
                    else if (dist < min_dist)
                    {
                        // If a new minimal element is found, clear the
                        // previous vector of pairs and simply keep this
                        // one.
                        min_dist = dist;
                        vec_equidistant_pair.clear();
                        vec_equidistant_pair.push_back(std::make_pair(proj, d_elem_face_normal[elem]));
                    }
#endif

#if (NDIM == 3)
                    // Get the closest point on the element from the cell and the angle-weighted
                    // pseudo-normal of the element
                    const auto proj_and_normal = getClosestPointandAngleWeightedNormal3D(P, elem);
                    proj = proj_and_normal.first;
                    dist = (P - proj).norm();

                    // If the distance is the same as the minimal
                    // distance, then the cell is equidistant to multiple
                    // elements, so add it to the set.
                    if (MathUtilities<double>::equalEps(dist, min_dist))
                    {
                        vec_equidistant_pair.push_back(std::make_pair(proj, proj_and_normal.second));
                    }
                    else if (dist < min_dist)
                    {
                        // If a new minimal element is found, clear the
                        // previous vector of pairs and simply keep this
                        // one.
                        min_dist = dist;
                        vec_equidistant_pair.clear();
                        vec_equidistant_pair.push_back(std::make_pair(proj, proj_and_normal.second));
                    }
#endif
                }

                // Determine the sign based on the angle-weighted pseudo-normal vectors
                // If the cell is equidistant to multiple elements, take
                // the average pseudo-normal from those elements.
                double sgn = 0.0;
                const size_t vec_length = vec_equidistant_pair.size();
                TBOX_ASSERT(vec_length > 0);
                IBTK::VectorNd avg_unit_normal, avg_proj;
                avg_unit_normal.setZero();
                avg_proj.setZero();
                for (const auto& elem_vec_pair : vec_equidistant_pair)
                {
                    const IBTK::VectorNd proj = elem_vec_pair.first;
                    avg_proj += proj;

                    const IBTK::VectorNd N = elem_vec_pair.second;
                    avg_unit_normal += N;
                }

                // Take the average normal and normalize it.
                avg_unit_normal /= static_cast<double>(vec_length);
                avg_unit_normal /= avg_unit_normal.norm();

                // Average the proj point.
                avg_proj /= static_cast<double>(vec_length);

                // Compute the signed distance function.
                sgn = avg_unit_normal.dot(P - avg_proj) <= 0.0 ? -1.0 : 1.0;
                (*d_data)(ci) = sgn * min_dist;
            }
        }
    }
} // computeSignedDistance

void
FESurfaceDistanceEvaluator::calculateSurfaceNormals()
{
    // Clear out the data structure
    d_elem_face_normal.clear();
    d_node_to_elem.clear();
    d_edge_to_elem.clear();

    // Loop over the active elements
    MeshBase::const_element_iterator el_it = d_mesh.active_elements_begin();
    MeshBase::const_element_iterator el_end = d_mesh.active_elements_end();
    if (d_use_vol_extracted_bdry_mesh)
    {
        el_it = d_bdry_mesh.active_elements_begin();
        el_end = d_bdry_mesh.active_elements_end();
    }
    for (; el_it != el_end; ++el_it)
    {
        Elem* elem = *el_it;

        std::array<VectorValue<double>, 2> Edge;
        const libMesh::Point& n0 = elem->point(0);
        const libMesh::Point& n1 = elem->point(1);
        Edge[0] = make_edge<VectorValue<double>, libMesh::Point>(n0, n1);

        if (NDIM == 2)
        {
            Edge[1] = VectorValue<double>(0.0, 0.0, 1.0);
        }

#if (NDIM == 3)
        const libMesh::Point& n2 = elem->point(2);
        Edge[1] = make_edge<VectorValue<double>, libMesh::Point>(n0, n2);
#endif
        // cross product of edge 0 and edge 1 results the outward face normal.
        VectorValue<double> N = Edge[0].cross(Edge[1]);
        IBTK::VectorNd normal;

#if (NDIM == 2)
        normal << N(0), N(1);
        normal /= normal.norm();
#endif

#if (NDIM == 3)
        normal << N(0), N(1), N(2);
        normal /= normal.norm();
#endif

        d_elem_face_normal[elem] = normal;

#if (NDIM == 3)
        // for each node, store the current element in the node to element map.
        const std::size_t n_nodes = elem->n_nodes();
        for (std::size_t node_n = 0; node_n < n_nodes; ++node_n) d_node_to_elem[elem->node_ptr(node_n)].insert(elem);

        // for each edge, store the current element in the edge to element map.
        std::vector<Node*> edge_nodes(2);
        constexpr unsigned int edge_pairs[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
        for (std::size_t edge_n = 0; edge_n < 3; ++edge_n)
        {
            for (std::size_t i = 0; i < 2; ++i)
            {
                std::size_t node_n = edge_pairs[edge_n][i];
                edge_nodes[i] = elem->node_ptr(node_n);
            }
            std::sort(edge_nodes.begin(), edge_nodes.end());
            std::pair<Node*, Node*> edge(edge_nodes[0], edge_nodes[1]);
            d_edge_to_elem[edge].insert(elem);
        }
#endif
    }
    return;
} // calculateSurfaceNormal

/////////////////////////////// STATIC ///////////////////////////////////////
void
FESurfaceDistanceEvaluator::updateSignAwayFromInterface(int D_idx,
                                                        Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                        double large_distance)
{
    // Allocate scratch variable on the finest level.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SAMRAI::hier::Variable<NDIM> > data_var;
    var_db->mapIndexToVariable(D_idx, data_var);
    Pointer<CellVariable<NDIM, double> > D_var = data_var;
#if !defined(NDEBUG)
    TBOX_ASSERT(!D_var.isNull());
#endif
    static const IntVector<NDIM> cell_ghosts = 1;
    const int D_iter_idx = var_db->registerVariableAndContext(
        D_var, var_db->getContext("FESurfaceDistanceEvaluator::updateSignAwayFromInterface::ITER"), cell_ghosts);
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(finest_ln);
    level->allocatePatchData(D_iter_idx, /*time*/ 0.0);

    // Copy d_idx to D_iter_idx.
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy, finest_ln, finest_ln);
    hier_cc_data_ops.setToScalar(D_iter_idx,
                                 large_distance,
                                 /*interior_only*/ false);
    hier_cc_data_ops.copyData(D_iter_idx, D_idx);

    // Fill ghost cells.
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(D_iter_idx, D_iter_idx, D_iter_idx, nullptr);
    Pointer<RefineSchedule<NDIM> > ghost_fill_sched = ghost_fill_alg.createSchedule(level);

    int n_global_updates = 1;
    while (n_global_updates > 0)
    {
        ghost_fill_sched->fillData(/*time*/ 0.0);
        int n_local_updates = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            // If the patch box doesn't intersect the structure box, no need to do
            // computations.
            if (!patch_box.intersects(d_large_struct_box)) continue;

            const SAMRAI::hier::Index<NDIM>& patch_lower_index = patch_box.lower();
            const SAMRAI::hier::Index<NDIM>& patch_upper_index = patch_box.upper();
            Pointer<CellData<NDIM, double> > D_iter_data = patch->getPatchData(D_iter_idx);
            double* const D = D_iter_data->getPointer(0);
            const int D_ghosts = (D_iter_data->getGhostCellWidth()).max();

            SIGN_SWEEP_FC(D,
                          D_ghosts,
                          patch_lower_index(0),
                          patch_upper_index(0),
                          patch_lower_index(1),
                          patch_upper_index(1),
#if (NDIM == 3)
                          patch_lower_index(2),
                          patch_upper_index(2),
#endif
                          large_distance,
                          n_local_updates);
        }
        n_global_updates = IBTK_MPI::sumReduction(n_local_updates);
    }

    // Copy D_iter_idx to D_idx and deallocate D_iter_idx.
    hier_cc_data_ops.copyData(D_idx, D_iter_idx);
    level->deallocatePatchData(D_iter_idx);
    var_db->removePatchDataIndex(D_iter_idx);
    return;
} // updateSignAwayFromInterface

#if (NDIM == 2)
bool
FESurfaceDistanceEvaluator::checkIntersection2D(const IBTK::Vector3d& box_bl,
                                                const IBTK::Vector3d& box_tr,
                                                const IBTK::Vector3d& box_br,
                                                const IBTK::Vector3d& box_tl,
                                                const libMesh::Point& n0,
                                                const libMesh::Point& n1)
{
    // If the line element is entirely contained within the box, then count as
    // "intersected".
    if ((box_bl(0) <= n0(0) && n0(0) <= box_tr(0)) && (box_bl(1) <= n0(1) && n0(1) <= box_tr(1)) &&
        (box_bl(0) <= n1(0) && n1(0) <= box_tr(0)) && (box_bl(1) <= n1(1) && n1(1) <= box_tr(1)))
    {
        return true;
    }

    auto line_equation = [](const IBTK::Vector3d& coord, const libMesh::Point& n0, const libMesh::Point& n1)
    {
        const double x0 = n0(0);
        const double y0 = n0(1);
        const double x1 = n1(0);
        const double y1 = n1(1);

        return (y1 - y0) * coord(0) + (x0 - x1) * coord(1) + (x1 * y0 - x0 * y1);
    }; // line_equation

    const double F_br = line_equation(box_br, n0, n1);
    const double F_bl = line_equation(box_bl, n0, n1);
    const double F_tr = line_equation(box_tr, n0, n1);
    const double F_tl = line_equation(box_tl, n0, n1);

    if ((F_br <= 0 && F_bl <= 0 && F_tr <= 0 && F_tl <= 0) || (F_br >= 0 && F_bl >= 0 && F_tr >= 0 && F_tl >= 0))
    {
        return false;
    }
    const double x0 = n0(0);
    const double y0 = n0(1);
    const double x1 = n1(0);
    const double y1 = n1(1);
    const double x_tr = box_tr(0);
    const double x_bl = box_bl(0);
    const double y_tr = box_tr(1);
    const double y_bl = box_bl(1);

    if (x0 > x_tr && x1 > x_tr) return false;
    if (x0 < x_bl && x1 < x_bl) return false;
    if (y0 > y_tr && y1 > y_tr) return false;
    if (y0 < y_bl && y1 < y_bl) return false;

    return true;
} // checkIntersection2D
#endif

#if (NDIM == 3)
bool
FESurfaceDistanceEvaluator::checkIntersection3D(const IBTK::Vector3d& box_center,
                                                const IBTK::Vector3d& box_half_dx,
                                                const IBTK::Vector3d& vert0,
                                                const IBTK::Vector3d& vert1,
                                                const IBTK::Vector3d& vert2)
{
    // If the tri element is entirely contained within the box, then count as
    // "intersected".
    bool contains_vert0 = true;
    bool contains_vert1 = true;
    bool contains_vert2 = true;

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        const double c = box_center(d);
        const double dx = box_half_dx(d);
        contains_vert0 = contains_vert0 && (c - dx <= vert0(d) && vert0(d) <= c + dx);
        contains_vert1 = contains_vert1 && (c - dx <= vert1(d) && vert1(d) <= c + dx);
        contains_vert2 = contains_vert2 && (c - dx <= vert1(d) && vert2(d) <= c + dx);
    }
    if (contains_vert0 && contains_vert1 && contains_vert2)
    {
        return true;
    }

    // Use the separating axis theorem to test the overlap between the tri element
    // and the box.
    double pmin, pmax, p0, p1, p2, rad, fex, fey, fez, a, b, fa, fb;
    IBTK::Vector3d v0, v1, v2, normal, e0, e1, e2;
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    // Shift everything so that box_center is the origin.
    v0 = vert0 - box_center;
    v1 = vert1 - box_center;
    v2 = vert2 - box_center;

    // Compute the triangle edges.
    e0 = v1 - v0;
    e1 = v2 - v1;
    e2 = v0 - v2;

    // Outer test 3 (9 tests in total)
    fex = std::abs(e0(X));
    fey = std::abs(e0(Y));
    fez = std::abs(e0(Z));

    // test 1
    a = e0[Z];
    b = e0[Y];
    fa = fez;
    fb = fey;
    p0 = a * v0[Y] - b * v0[Z];
    p2 = a * v2[Y] - b * v2[Z];
    pmin = std::min(p0, p2);
    pmax = std::max(p0, p2);
    rad = fa * box_half_dx[Y] + fb * box_half_dx[Z];
    if (pmin > rad || pmax < -rad) return false;

    // test 2
    a = e0[Z];
    b = e0[X];
    fa = fez;
    fb = fex;
    p0 = -a * v0[X] + b * v0[Z];
    p2 = -a * v2[X] + b * v2[Z];
    pmin = std::min(p0, p2);
    pmax = std::max(p0, p2);
    rad = fa * box_half_dx[X] + fb * box_half_dx[Z];
    if (pmin > rad || pmax < -rad) return false;

    // test 3
    a = e0[Y];
    b = e0[X];
    fa = fey;
    fb = fex;
    p1 = a * v1[X] - b * v1[Y];
    p2 = a * v2[X] - b * v2[Y];
    pmin = std::min(p1, p2);
    pmax = std::max(p1, p2);
    rad = fa * box_half_dx[X] + fb * box_half_dx[Y];
    if (pmin > rad || pmax < -rad) return false;

    fex = std::abs(e1(X));
    fey = std::abs(e1(Y));
    fez = std::abs(e1(Z));

    // test 4
    a = e1[Z];
    b = e1[Y];
    fa = fez;
    fb = fey;
    p0 = a * v0[Y] - b * v0[Z];
    p2 = a * v2[Y] - b * v2[Z];
    pmin = std::min(p0, p2);
    pmax = std::max(p0, p2);
    rad = fa * box_half_dx[Y] + fb * box_half_dx[Z];
    if (pmin > rad || pmax < -rad) return false;

    // test 5
    a = e1[Z];
    b = e1[X];
    fa = fez;
    fb = fex;
    p0 = -a * v0[X] + b * v0[Z];
    p2 = -a * v2[X] + b * v2[Z];
    pmin = std::min(p0, p2);
    pmax = std::max(p0, p2);
    rad = fa * box_half_dx[X] + fb * box_half_dx[Z];
    if (pmin > rad || pmax < -rad) return false;

    // test 6
    a = e1[Y];
    b = e1[X];
    fa = fey;
    fb = fex;
    p0 = a * v0[X] - b * v0[Y];
    p1 = a * v1[X] - b * v1[Y];
    pmin = std::min(p0, p1);
    pmax = std::max(p0, p1);
    rad = fa * box_half_dx[X] + fb * box_half_dx[Y];
    if (pmin > rad || pmax < -rad) return false;

    fex = std::abs(e2(X));
    fey = std::abs(e2(Y));
    fez = std::abs(e2(Z));

    // test 7
    a = e2[Z];
    b = e2[Y];
    fa = fez;
    fb = fey;
    p0 = a * v0[Y] - b * v0[Z];
    p1 = a * v1[Y] - b * v1[Z];
    pmin = std::min(p0, p1);
    pmax = std::max(p0, p1);
    rad = fa * box_half_dx[Y] + fb * box_half_dx[Z];
    if (pmin > rad || pmax < -rad) return false;

    // test 8
    a = e2[Z];
    b = e2[X];
    fa = fez;
    fb = fex;
    p0 = -a * v0[X] + b * v0[Z];
    p1 = -a * v1[X] + b * v1[Z];
    pmin = std::min(p0, p1);
    pmax = std::max(p0, p1);
    rad = fa * box_half_dx[X] + fb * box_half_dx[Z];
    if (pmin > rad || pmax < -rad) return false;

    // test 9
    a = e2[Y];
    b = e2[X];
    fa = fey;
    fb = fex;
    p1 = a * v1[X] - b * v1[Y];
    p2 = a * v2[X] - b * v2[Y];
    pmin = std::min(p1, p2);
    pmax = std::max(p1, p2);
    rad = fa * box_half_dx[X] + fb * box_half_dx[Y];
    if (pmin > rad || pmax < -rad) return false;

    // Outer test 1 (3 tests total)
    // test 1
    pmin = std::min(v0[X], std::min(v1[X], v2[X]));
    pmax = std::max(v0[X], std::max(v1[X], v2[X]));
    if (pmin > box_half_dx[X] || pmax < -box_half_dx[X]) return false;

    // test 2
    pmin = std::min(v0[Y], std::min(v1[Y], v2[Y]));
    pmax = std::max(v0[Y], std::max(v1[Y], v2[Y]));
    if (pmin > box_half_dx[Y] || pmax < -box_half_dx[Y]) return false;

    // test 3
    pmin = std::min(v0[Z], std::min(v1[Z], v2[Z]));
    pmax = std::max(v0[Z], std::max(v1[Z], v2[Z]));
    if (pmin > box_half_dx[Z] || pmax < -box_half_dx[Z]) return false;

    // Outer test 2 (1 test total)
    IBTK::Vector3d vmin;
    normal = e0.cross(e1);
    for (int q = X; q < Z; ++q)
    {
        const double v = v0[q];
        if (normal[q] > 0.0)
        {
            vmin[q] = -box_half_dx[q] - v;
        }
        else
        {
            vmin[q] = box_half_dx[q] - v;
        }
    }
    if (normal.dot(vmin) > 0.0) return false;

    return true;
} // checkIntersection3D
#endif

/////////////////////////////// PRIVATE ////////////////////////////
#if (NDIM == 3)
std::pair<IBTK::Vector3d, IBTK::Vector3d>
FESurfaceDistanceEvaluator::getClosestPointandAngleWeightedNormal3D(const IBTK::Vector3d& P, libMesh::Elem* elem)
{
    const libMesh::Point& n0 = elem->point(0);
    const libMesh::Point& n1 = elem->point(1);
    const libMesh::Point& n2 = elem->point(2);
    const IBTK::Vector3d B(n0(0), n0(1), n0(2));
    const IBTK::Vector3d D = B - P;

    const IBTK::Vector3d E0 = make_edge<IBTK::Vector3d, libMesh::Point>(n0, n1);
    const IBTK::Vector3d E1 = make_edge<IBTK::Vector3d, libMesh::Point>(n0, n2);

    // Compute the required dot products.
    const double a = E0.dot(E0);
    const double b = E0.dot(E1);
    const double c = E1.dot(E1);
    const double d = E0.dot(D);
    const double e = E1.dot(D);
    const double det = a * c - b * b;
    double s = b * e - c * d;
    double t = b * d - a * e;

    // brief Enumerated type for different triangle feature.
    enum TRIANGLEFEATURE
    {
        FACE,
        EDGE_0,
        EDGE_1,
        EDGE_2,
        VERTEX_0,
        VERTEX_1,
        VERTEX_2,
        UNKNOWN_FEATURE = -1
    };

    TRIANGLEFEATURE feature = UNKNOWN_FEATURE;

    if (s + t <= det)
    {
        if (s < 0)
        {
            if (t < 0)
            {
                // region 4
                if (d < 0)
                {
                    t = 0;
                    if (-d >= a)
                    {
                        s = 1;
                        feature = VERTEX_1;
                    }
                    else
                    {
                        s = -d / a;
                        feature = EDGE_0;
                    }
                }
                else
                {
                    s = 0;
                    if (e >= 0)
                    {
                        t = 0;
                        feature = VERTEX_0;
                    }
                    else
                    {
                        if (-e >= c)
                        {
                            t = 1;
                            feature = VERTEX_2;
                        }
                        else
                        {
                            t = -e / c;
                            feature = EDGE_2;
                        }
                    }
                }
            }
            else
            {
                // region 3
                s = 0;
                if (e >= 0)
                {
                    t = 0;
                    feature = VERTEX_0;
                }
                else
                {
                    if (-e >= c)
                    {
                        t = 1;
                        feature = VERTEX_2;
                    }
                    else
                    {
                        t = -e / c;
                        feature = EDGE_2;
                    }
                }
            }
        }
        else
        {
            if (t < 0)
            {
                // region 5
                t = 0;
                if (d >= 0)
                {
                    s = 0;
                    feature = VERTEX_0;
                }
                else
                {
                    if (-d >= a)
                    {
                        s = 1;
                        feature = VERTEX_1;
                    }
                    else
                    {
                        s = -d / a;
                        feature = EDGE_0;
                    }
                }
            }
            else
            {
                // region 0
                const double invDet = 1.0 / det;
                s = s * invDet;
                t = t * invDet;
                feature = FACE;
            }
        }
    }
    else
    {
        if (s < 0)
        {
            // region 2
            const double tmp0 = b + d;
            const double tmp1 = c + e;
            if (tmp1 > tmp0) // minimum on edge s+t=1
            {
                const double numer = tmp1 - tmp0;
                const double denom = a - 2.0 * b + c;
                if (numer >= denom)
                {
                    s = 1.0;
                    t = 0.0;
                    feature = VERTEX_1;
                }
                else
                {
                    s = numer / denom;
                    t = 1.0 - s;
                    feature = EDGE_1;
                }
            }
            else
            {
                // minimum on edge s=0
                s = 0.0;
                if (tmp1 <= 0.0)
                {
                    t = 1.0;
                    feature = VERTEX_2;
                }
                else
                {
                    if (e >= 0.0)
                    {
                        t = 0.0;
                        feature = VERTEX_0;
                    }
                    else
                    {
                        t = -e / c;
                        feature = EDGE_2;
                    }
                }
            }
        }
        else
        {
            if (t < 0.0)
            {
                // region 6
                const double tmp0 = b + e;
                const double tmp1 = a + d;
                if (tmp1 > tmp0)
                {
                    const double numer = tmp1 - tmp0;
                    const double denom = a - 2.0 * b + c;
                    if (numer >= denom)
                    {
                        t = 1.0;
                        s = 0.0;
                        feature = VERTEX_2;
                    }
                    else
                    {
                        t = numer / denom;
                        s = 1.0 - t;
                        feature = EDGE_1;
                    }
                }
                else
                {
                    t = 0.0;
                    if (tmp1 <= 0.0)
                    {
                        s = 1.0;
                        feature = VERTEX_1;
                    }
                    else
                    {
                        if (d >= 0.0)
                        {
                            s = 0.0;
                            feature = VERTEX_0;
                        }
                        else
                        {
                            s = -d / a;
                            feature = EDGE_0;
                        }
                    }
                }
            }
            else
            {
                // region 1
                const double numer = c + e - b - d;
                if (numer <= 0.0)
                {
                    s = 0.0;
                    t = 1.0;
                    feature = VERTEX_2;
                }
                else
                {
                    const double denom = a - 2.0 * b + c;
                    if (numer >= denom)
                    {
                        s = 1.0;
                        t = 0.0;
                        feature = VERTEX_1;
                    }
                    else
                    {
                        s = numer / denom;
                        t = 1.0 - s;
                        feature = EDGE_1;
                    }
                }
            }
        }
    }

#if !defined(NDEBUG)
    TBOX_ASSERT(feature != UNKNOWN_FEATURE);
    TBOX_ASSERT(s >= 0.0 && s <= 1.0);
    TBOX_ASSERT(t >= 0.0 && t <= 1.0);
#endif

    std::vector<Node*> edge_nodes(2);
    std::set<Elem*> elem_set;
    std::pair<Node*, Node*> edge;
    IBTK::Vector3d nbr_e0, nbr_e1, angle_weighted_normal = IBTK::Vector3d::Zero();

    // Depending on the identified feature, angle-weighted pseudo-normal is computed.
    // If the identified feature is face, the face normal of the element is returned.
    // If the identified feature is edge, then the average of the face normal of the
    // elements share the edge is returned.
    // If the identified feature is vertex, then the angle-weighted average of the
    // face normal of the elements share the vertex is returned.
    switch (feature)
    {
    case FACE:
        angle_weighted_normal = d_elem_face_normal[elem];
        break;
    case EDGE_0:
        edge_nodes[0] = elem->node_ptr(0);
        edge_nodes[1] = elem->node_ptr(1);
        std::sort(edge_nodes.begin(), edge_nodes.end());
        edge = std::make_pair(edge_nodes[0], edge_nodes[1]);
        elem_set = d_edge_to_elem[edge];
        for (const auto& itr : elem_set)
        {
            Elem* nbr_elem = itr;
            angle_weighted_normal += d_elem_face_normal[nbr_elem];
        }
        angle_weighted_normal /= angle_weighted_normal.norm();
        break;
    case EDGE_1:
        edge_nodes[0] = elem->node_ptr(1);
        edge_nodes[1] = elem->node_ptr(2);
        std::sort(edge_nodes.begin(), edge_nodes.end());
        edge = std::make_pair(edge_nodes[0], edge_nodes[1]);
        elem_set = d_edge_to_elem[edge];
        for (const auto& itr : elem_set)
        {
            Elem* nbr_elem = itr;
            angle_weighted_normal += d_elem_face_normal[nbr_elem];
        }
        angle_weighted_normal /= angle_weighted_normal.norm();
        break;
    case EDGE_2:
        edge_nodes[0] = elem->node_ptr(2);
        edge_nodes[1] = elem->node_ptr(0);
        std::sort(edge_nodes.begin(), edge_nodes.end());
        edge = std::make_pair(edge_nodes[0], edge_nodes[1]);
        elem_set = d_edge_to_elem[edge];
        for (const auto& itr : elem_set)
        {
            Elem* nbr_elem = itr;
            angle_weighted_normal += d_elem_face_normal[nbr_elem];
        }
        angle_weighted_normal /= angle_weighted_normal.norm();
        break;
    case VERTEX_0:
        elem_set = d_node_to_elem[elem->node_ptr(0)];
        for (const auto& itr : elem_set)
        {
            Elem* nbr_elem = itr;
            const libMesh::Point& nbr_n0 = nbr_elem->point(0);
            const libMesh::Point& nbr_n1 = nbr_elem->point(1);
            const libMesh::Point& nbr_n2 = nbr_elem->point(2);
            if (nbr_elem->node_ptr(0) == elem->node_ptr(0))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n0, nbr_n1);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n0, nbr_n2);
            }
            else if (nbr_elem->node_ptr(1) == elem->node_ptr(0))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n1, nbr_n0);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n1, nbr_n2);
            }
            else if (nbr_elem->node_ptr(2) == elem->node_ptr(0))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n2, nbr_n0);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n2, nbr_n1);
            }
            const double nbr_e0_norm = nbr_e0.norm();
            const double nbr_e1_norm = nbr_e1.norm();
            const double angle = acos(nbr_e0.dot(nbr_e1) / (nbr_e0_norm * nbr_e1_norm));
            angle_weighted_normal += (angle * d_elem_face_normal[nbr_elem]);
        }
        angle_weighted_normal /= angle_weighted_normal.norm();
        break;
    case VERTEX_1:
        elem_set = d_node_to_elem[elem->node_ptr(1)];
        for (const auto& itr : elem_set)
        {
            Elem* nbr_elem = itr;
            const libMesh::Point& nbr_n0 = nbr_elem->point(0);
            const libMesh::Point& nbr_n1 = nbr_elem->point(1);
            const libMesh::Point& nbr_n2 = nbr_elem->point(2);
            if (nbr_elem->node_ptr(0) == elem->node_ptr(1))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n0, nbr_n1);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n0, nbr_n2);
            }
            else if (nbr_elem->node_ptr(1) == elem->node_ptr(1))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n1, nbr_n0);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n1, nbr_n2);
            }
            else if (nbr_elem->node_ptr(2) == elem->node_ptr(1))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n2, nbr_n0);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n2, nbr_n1);
            }
            const double nbr_e0_norm = nbr_e0.norm();
            const double nbr_e1_norm = nbr_e1.norm();
            const double angle = acos(nbr_e0.dot(nbr_e1) / (nbr_e0_norm * nbr_e1_norm));
            angle_weighted_normal += (angle * d_elem_face_normal[nbr_elem]);
        }
        angle_weighted_normal /= angle_weighted_normal.norm();
        break;
    case VERTEX_2:
        elem_set = d_node_to_elem[elem->node_ptr(2)];
        for (const auto& itr : elem_set)
        {
            Elem* nbr_elem = itr;
            const libMesh::Point& nbr_n0 = nbr_elem->point(0);
            const libMesh::Point& nbr_n1 = nbr_elem->point(1);
            const libMesh::Point& nbr_n2 = nbr_elem->point(2);
            if (nbr_elem->node_ptr(0) == elem->node_ptr(2))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n0, nbr_n1);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n0, nbr_n2);
            }
            else if (nbr_elem->node_ptr(1) == elem->node_ptr(2))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n1, nbr_n0);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n1, nbr_n2);
            }
            else if (nbr_elem->node_ptr(2) == elem->node_ptr(2))
            {
                nbr_e0 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n2, nbr_n0);
                nbr_e1 = make_edge<IBTK::Vector3d, libMesh::Point>(nbr_n2, nbr_n1);
            }
            const double nbr_e0_norm = nbr_e0.norm();
            const double nbr_e1_norm = nbr_e1.norm();
            const double angle = acos(nbr_e0.dot(nbr_e1) / (nbr_e0_norm * nbr_e1_norm));
            angle_weighted_normal += (angle * d_elem_face_normal[nbr_elem]);
        }
        angle_weighted_normal /= angle_weighted_normal.norm();
        break;
    case UNKNOWN_FEATURE:
    default:
        TBOX_ERROR("FESurfaceDistanceEvaluator::getClosestPointandAngleWeightedNormal3D()::Unknown geometric feature");
    }
    return std::make_pair((B + s * E0 + t * E1), angle_weighted_normal);
} // getClosestPointandAngleWeightedNormal3D
#endif

void
FESurfaceDistanceEvaluator::collectNeighboringPatchElements(int level_number)
{
    // Clear out the data structure.
    d_active_neighbor_patch_bdry_elem_map.clear();

    // Setup data structures used to assign elements to patches.
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(level_number);
    const int num_local_patches = level->getProcessorMapping().getNumberOfLocalIndices();
    d_active_neighbor_patch_bdry_elem_map.resize(num_local_patches);
    IntVector<NDIM> ghost_width = d_gcw;

    int local_patch_num = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++, ++local_patch_num)
    {
        d_active_neighbor_patch_bdry_elem_map[local_patch_num] = std::vector<Elem*>(0);
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
        IBTK::Point x_lower;
        for (unsigned int d = 0; d < NDIM; ++d) x_lower[d] = pgeom->getXLower()[d];
        IBTK::Point x_upper;
        for (unsigned int d = 0; d < NDIM; ++d) x_upper[d] = pgeom->getXUpper()[d];
        const double* const dx = pgeom->getDx();
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            x_lower[d] -= dx[d] * ghost_width[d];
            x_upper[d] += dx[d] * ghost_width[d];
        }

        // Loop over the active elements and see if their centroids or nodes reside
        // in the patch interior grown by the specified ghost cell width.
        MeshBase::const_element_iterator el_it = d_mesh.active_elements_begin();
        MeshBase::const_element_iterator el_end = d_mesh.active_elements_end();

        if (d_use_vol_extracted_bdry_mesh)
        {
            el_it = d_bdry_mesh.active_elements_begin();
            el_end = d_bdry_mesh.active_elements_end();
        }

        for (; el_it != el_end; ++el_it)
        {
            bool in_patch = false;
            Elem* const elem = *el_it;

            // Error checking for element type.
            if (elem->type() != d_supported_elem_type)
            {
                TBOX_ERROR(
                    "FESurfaceDistanceEvaluator presently does not support "
                    "elements of type "
                    << Utility::enum_to_string<ElemType>(elem->type()) << " for NDIM = " << NDIM
                    << ".\nSupported type is " << Utility::enum_to_string<ElemType>(d_supported_elem_type));
            }

            // First check the centroids.
            const libMesh::Point& c = elem->centroid();

            bool centroid_in_patch = true;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                in_patch = centroid_in_patch && (x_lower[d] <= c(d) && c(d) <= x_upper[d]);
            }

            // Next, check the nodes.
            const unsigned int n_nodes = elem->n_nodes();
            for (unsigned int k = 0; k < n_nodes && !in_patch; ++k)
            {
                const libMesh::Point& n = elem->point(k);

                bool node_in_patch = true;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    in_patch = node_in_patch && (x_lower[d] <= n(d) && n(d) <= x_upper[d]);
                }
            }

            // Add the element to the patch list if necessary.
            if (in_patch)
            {
                d_active_neighbor_patch_bdry_elem_map[local_patch_num].push_back(elem);
            }
        }
    }
    return;
} // collectNeighboringPatchElements
/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
