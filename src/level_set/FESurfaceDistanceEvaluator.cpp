// Filename FESurfaceDistanceEvaluator.cpp
// Created on Sep 5, 2018 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/FESurfaceDistanceEvaluator.h"
#include "CartesianGridGeometry.h"
#include "HierarchyCellDataOpsReal.h"
#include "IBAMR_config.h"
#include "ibamr/app_namespaces.h"
#include "ibtk/FEDataManager.h"
#include "ibtk/IndexUtilities.h"
#include "libmesh/equation_systems.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

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

inline double
line_equation(const IBTK::Vector3d& coord, const libMesh::Point& n0, const libMesh::Point& n1)
{
    const double x0 = n0(0);
    const double y0 = n0(1);
    const double x1 = n1(0);
    const double y1 = n1(1);

    return (y1 - y0) * coord(0) + (x0 - x1) * coord(1) + (x1 * y0 - x0 * y1);
} // line_equation

} // namespace

const double FESurfaceDistanceEvaluator::s_large_distance = 1234567.0;

/////////////////////////////// PUBLIC //////////////////////////////////////
FESurfaceDistanceEvaluator::FESurfaceDistanceEvaluator(const std::string& object_name,
                                                       Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                       Pointer<IBFEMethod> ibfe_method,
                                                       const BoundaryMesh& bdry_mesh,
                                                       const int part,
                                                       const int gcw)
    : d_object_name(object_name),
      d_patch_hierarchy(patch_hierarchy),
      d_ibfe_method(ibfe_method),
      d_bdry_mesh(bdry_mesh),
      d_part(part),
      d_gcw(gcw)
{
// The only supported element type for this class
#if (NDIM == 2)
    d_supported_elem_type = EDGE2;
#endif
#if (NDIM == 3)
    d_supported_elem_type = TRI3;
#endif

    // Note that this class is specialized to work on a boundary mesh with dim = NDIM - 1
    // derived from a volumetric mesh.
    if (d_bdry_mesh.mesh_dimension() != NDIM - 1)
    {
        TBOX_ERROR(
            "FESurfaceDistanceEvaluator presently requires a boundary mesh with dim = NDIM - 1 to be registered");
    }

    // Set up timers
    IBTK_DO_ONCE(t_collectNeighboringPatchElements = TimerManager::getManager()->getTimer(
                     "FESurfaceDistanceEvaluator::collectNeighboringPatchElements()", true);
                 t_buildIntersectionMap =
                     TimerManager::getManager()->getTimer("FESurfaceDistanceEvaluator::buildIntersectionMap()", true););

    return;
} // FESurfaceDistanceEvaluator

FESurfaceDistanceEvaluator::~FESurfaceDistanceEvaluator()
{
    // intentionally left blank
    return;
}

void
FESurfaceDistanceEvaluator::mapIntersections()
{
    // Loop over each cell and a set of triangles and map cell-triangle intersections.
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    IBTK_TIMER_START(t_collectNeighboringPatchElements);
    collectNeighboringPatchElements(finest_ln);
    IBTK_TIMER_STOP(t_collectNeighboringPatchElements);

    IBTK_TIMER_START(t_buildIntersectionMap);
    // Loop over patches on finest level, while keeping track of the local patch indexing
    IBFEMethod::CoordinateMappingFcnData mapping = d_ibfe_method->getInitialCoordinateMappingFunction(d_part);
    const bool identity_mapping = !(mapping.fcn);
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(finest_ln);
    int local_patch_num = 0;

    // Desired ghost cell width
    IntVector<NDIM> ghost_width = d_gcw;

    // Compute a bounding box for the entire structure, relying on LibMesh parallel decomposition
    MeshBase::const_element_iterator el_it = d_bdry_mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = d_bdry_mesh.active_local_elements_end();
    IBTK::Vector3d elem_bl, elem_tr;
    elem_bl << std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), std::numeric_limits<int>::max();
    elem_tr << std::numeric_limits<int>::min(), std::numeric_limits<int>::min(), std::numeric_limits<int>::min();
    for (; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        // Get the coordinates of the nodes
        const libMesh::Point& s0 = elem->point(0);
        libMesh::Point n0 = s0;
        const libMesh::Point& s1 = elem->point(1);
        libMesh::Point n1 = s1;
#if (NDIM == 3)
        const libMesh::Point& s2 = elem->point(2);
        libMesh::Point n2 = s2;
#endif

        if (!identity_mapping)
        {
            mapping.fcn(n0, s0, mapping.ctx);
            mapping.fcn(n1, s1, mapping.ctx);
#if (NDIM == 3)
            mapping.fcn(n2, s2, mapping.ctx);
#endif
        }
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
    SAMRAI_MPI::minReduction(elem_bl.data(), 3);
    SAMRAI_MPI::maxReduction(elem_tr.data(), 3);

    // Structure bounding box, taking into account ghost cell width
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = level->getGridGeometry();
    const double* const dx0 = grid_geom->getDx();
    const IntVector<NDIM>& level_ratio = level->getRatio();
    double level_dx[NDIM] = { 0.0 };
    for (int d = 0; d < NDIM; ++d) level_dx[d] = dx0[d] / level_ratio(d);

    IBTK::Vector3d struct_bl, struct_tr;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        struct_bl[d] = elem_bl[d] - (ghost_width[d] * level_dx[d]);
        struct_tr[d] = elem_tr[d] + (ghost_width[d] * level_dx[d]);
    }
    Box<NDIM> struct_box(IndexUtilities::getCellIndex(struct_bl.data(), grid_geom, level_ratio),
                         IndexUtilities::getCellIndex(struct_tr.data(), grid_geom, level_ratio));

    // Map the neighbor intersections
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
        const Index<NDIM>& patch_lower_index = patch_box.lower();
        const double* const patch_dx = patch_geom->getDx();

        // If the patch box doesn't intersect the structure box, no need to do computations
        if (!patch_box.intersects(struct_box)) continue;

        // Loop over cells
        for (Box<NDIM>::Iterator it(patch_box); it; it++)
        {
            // Get the coordinates/dimensions a box grown out of the cell center of ci
            // by the ghost cell width. This will ensure that we capture not only the elements
            // intersecting the cell, but also the elements intersecting the neighbor
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
            // checkIntersection routines
            Box<NDIM> ghost_box(IndexUtilities::getCellIndex(r_bl.data(), grid_geom, level_ratio),
                                IndexUtilities::getCellIndex(r_tr.data(), grid_geom, level_ratio));
            if (!ghost_box.intersects(struct_box)) continue;

                // Prepare the required vectors
#if (NDIM == 2)
            IBTK::Vector3d r_br, r_tl;
            r_br(1) = r_bl(1);
            r_br(0) = r_tr(0);
            r_tl(1) = r_tr(1);
            r_tl(0) = r_bl(0);
#endif
#if (NDIM == 3)
            // Get the coordinates of the cell center
            IBTK::Vector3d grown_box_center, grown_box_half_dx;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                grown_box_center[d] =
                    patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_index(d)) + 0.5);
                grown_box_half_dx[d] = (0.5 + ghost_width[d]) * patch_dx[d];
            }
#endif

            // Loop over elements in the patch
            for (std::vector<Elem*>::const_iterator eit = patch_elems.begin(); eit != patch_elems.end(); ++eit)
            {
                // Get the coordinates of the nodes
                Elem* const elem = *eit;
                const libMesh::Point& s0 = elem->point(0);
                libMesh::Point n0 = s0;
                const libMesh::Point& s1 = elem->point(1);
                libMesh::Point n1 = s1;
#if (NDIM == 3)
                const libMesh::Point& s2 = elem->point(2);
                libMesh::Point n2 = s2;
#endif

                if (!identity_mapping)
                {
                    mapping.fcn(n0, s0, mapping.ctx);
                    mapping.fcn(n1, s1, mapping.ctx);
#if (NDIM == 3)
                    mapping.fcn(n2, s2, mapping.ctx);
#endif
                }

                // Intersection detection routines
#if (NDIM == 2)
                const bool found_intersection = checkIntersection2D(r_bl, r_tr, r_br, r_tl, n0, n1);
#endif
#if (NDIM == 3)
                // In 3D, it is more convenient to pass vertices as IBTK::Vector3d
                IBTK::Vector3d vert0, vert1, vert2;
                vert0 << n0(0), n0(1), n0(2);
                vert1 << n1(0), n1(1), n1(2);
                vert2 << n2(0), n2(1), n2(2);
                const bool found_intersection =
                    checkIntersection3D(grown_box_center, grown_box_half_dx, vert0, vert1, vert2);
#endif
                if (found_intersection)
                {
                    // Store the element corresponding to the volumetric mesh
                    d_cell_elem_neighbor_map[ci].insert(elem->interior_parent());
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
    // Get the normals for the boundary mesh using a surface quadrature rule
    EquationSystems* equation_systems = d_ibfe_method->getFEDataManager()->getEquationSystems();
    System& X_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    X_system.solution->localize(*X_system.current_local_solution);
    DofMap& X_dof_map = X_system.get_dof_map();
    FEType fe_type = X_dof_map.variable_type(0);
    const int bdry_mesh_dim = d_bdry_mesh.mesh_dimension();
    libMesh::UniquePtr<FEBase> fe_bdry(FEBase::build(NDIM, fe_type));

    // Ensures only a single quadrature point is used for the normals, which is sufficient for linear elements.
    libMesh::UniquePtr<QBase> qrule_bdry = QBase::build(QGAUSS, bdry_mesh_dim, CONSTANT);
    fe_bdry->attach_quadrature_rule(qrule_bdry.get());

    // Loop over patches on finest level
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    IBFEMethod::CoordinateMappingFcnData mapping = d_ibfe_method->getInitialCoordinateMappingFunction(/*part*/ 0);
    const bool identity_mapping = !(mapping.fcn);
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Index<NDIM>& patch_lower_index = patch_box.lower();
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
                std::set<Elem*> elem_set = d_cell_elem_neighbor_map[ci];
                const int num_elements = static_cast<int>(elem_set.size());
                (*n_data)(ci) = num_elements;

                // Loop over the cutting elements and find the minimum distance
                IBTK::VectorNd P;
                for (int d = 0; d < NDIM; ++d)
                {
                    P[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_index(d)) + 0.5);
                }
                double min_dist = std::numeric_limits<double>::max();

                // Create a pair to take care of normals for cells equidistant to multiple elements.
                std::vector<std::pair<libMesh::Elem*, IBTK::VectorNd> > vec_equidistant_pair;
                for (std::set<Elem*>::const_iterator it = elem_set.begin(); it != elem_set.end(); ++it)
                {
                    Elem* elem = *it;

                    // Loop over the sides of the element. If it has no neighbors on a side,
                    // then it MUST live on the boundary of the mesh
                    for (unsigned int s = 0; s < elem->n_sides(); ++s)
                    {
                        if (elem->neighbor(s) != NULL) continue;
                        UniquePtr<Elem> side_elem = elem->build_side(s, /*proxy*/ false);
                        IBTK::VectorNd v, w, proj;
                        double dist = std::numeric_limits<double>::max();

                        // Get the nodes
                        const libMesh::Point& s0 = side_elem->point(0);
                        libMesh::Point n0 = s0;
                        const libMesh::Point& s1 = side_elem->point(1);
                        libMesh::Point n1 = s1;
#if (NDIM == 3)
                        const libMesh::Point& s2 = side_elem->point(2);
                        libMesh::Point n2 = s2;
#endif
                        if (!identity_mapping)
                        {
                            mapping.fcn(n0, s0, mapping.ctx);
                            mapping.fcn(n1, s1, mapping.ctx);
#if (NDIM == 3)
                            mapping.fcn(n2, s2, mapping.ctx);
#endif
                        }
#if (NDIM == 2)
                        v << n0(0), n0(1);
                        w << n1(0), n1(1);

                        const double L2 = (v - w).squaredNorm();
                        if (MathUtilities<double>::equalEps(L2, 0.0))
                        {
                            // Special case where line element collapses to a point.
                            // Shouldn't happen.
                            proj = v;
                            dist = (proj - P).norm();
                        }
                        else
                        {
                            // Parameterize and project
                            // Note that this will take care of the edge case where the projection
                            // does not fall on the line
                            const double t = std::max(0.0, std::min(1.0, (P - v).dot(w - v) / L2));
                            proj = v + t * (w - v);
                            dist = (P - proj).norm();
                        }
#endif
#if (NDIM == 3)
                        proj = getClosestPoint3D(P, n0, n1, n2);
                        dist = (P - proj).norm();
#endif

                        // If the distance is the same as the minimal distance, then the cell is equidistant
                        // to multiple elements, so add it to the set.
                        if (MathUtilities<double>::equalEps(dist, min_dist))
                        {
                            vec_equidistant_pair.push_back(std::make_pair(elem, proj));
                        }
                        else if (dist < min_dist)
                        {
                            // If a new minimal element is found, clear the previous vector of pairs and simply keep
                            // this one
                            min_dist = dist;
                            vec_equidistant_pair.clear();
                            vec_equidistant_pair.push_back(std::make_pair(elem, proj));
                        }
                    }
                }

                // Determine the sign based on the normal vectors at the quadrature points
                // For linear elements considered in this class, the normal is the same no matter the quadrature point
                // If the cell is equidistant to multiple elements, take the average normal from those elements
                double sgn = 0.0;
                const size_t vec_length = vec_equidistant_pair.size();
                TBOX_ASSERT(vec_length > 0);
                IBTK::VectorNd avg_unit_normal, avg_proj;
                avg_unit_normal.setZero();
                avg_proj.setZero();
                std::vector<std::pair<Elem*, IBTK::VectorNd> >::const_iterator it;
                for (it = vec_equidistant_pair.begin(); it != vec_equidistant_pair.end(); ++it)
                {
                    Elem* elem = (*it).first;
                    IBTK::VectorNd proj = (*it).second;
                    avg_proj += proj;

                    // Loop over the sides of the element. If it has no neighbors on a side,
                    // then it MUST live on the boundary of the mesh
                    for (unsigned int s = 0; s < elem->n_sides(); ++s)
                    {
                        if (elem->neighbor(s) != NULL) continue;
                        {
                            fe_bdry->reinit(elem, s);
                            const unsigned int n_qp = qrule_bdry->n_points();
                            const std::vector<libMesh::Point>& bdry_normals = fe_bdry->get_normals();
                            TBOX_ASSERT(n_qp == 1); // for now
                            unsigned int qp = 0;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                avg_unit_normal(d) += bdry_normals[qp](d);
                            }
                        }
                    }
                }

                // Take the average normal and normalize it
                avg_unit_normal /= (double)vec_length;
                avg_unit_normal /= avg_unit_normal.norm();

                // Average the proj point
                avg_proj /= (double)vec_length;

                // Compute the signed distance function.
                sgn = avg_unit_normal.dot(P - avg_proj) <= 0.0 ? -1.0 : 1.0;
                (*d_data)(ci) = sgn * min_dist;
            }
        }
    }
} // computeSignedDistance

void
FESurfaceDistanceEvaluator::updateSignAwayFromInterface(int D_idx)
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
    const int D_iter_idx =
        var_db->registerVariableAndContext(D_var, var_db->getContext(d_object_name + "::ITER"), cell_ghosts);
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(finest_ln);
    level->allocatePatchData(D_iter_idx, /*time*/ 0.0);

    // Copy d_idx to D_iter_idx
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(d_patch_hierarchy, finest_ln, finest_ln);
    hier_cc_data_ops.setToScalar(D_iter_idx, s_large_distance, /*interior_only*/ false);
    hier_cc_data_ops.copyData(D_iter_idx, D_idx);

    // Fill ghost cells
    RefineAlgorithm<NDIM> ghost_fill_alg;
    ghost_fill_alg.registerRefine(D_iter_idx, D_iter_idx, D_iter_idx, NULL);
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
            const Index<NDIM>& patch_lower_index = patch_box.lower();
            const Index<NDIM>& patch_upper_index = patch_box.upper();

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
                          s_large_distance,
                          n_local_updates);
        }
        n_global_updates = SAMRAI_MPI::sumReduction(n_local_updates);
    }

    // Copy D_iter_idx to D_idx and deallocate D_iter_idx
    hier_cc_data_ops.copyData(D_idx, D_iter_idx);
    level->deallocatePatchData(D_iter_idx);
    var_db->removePatchDataIndex(D_iter_idx);

    return;
} // updateSignAwayFromInterface

/////////////////////////////// STATIC //////////////////////////////////////

bool
FESurfaceDistanceEvaluator::checkIntersection2D(const IBTK::Vector3d& box_bl,
                                                const IBTK::Vector3d& box_tr,
                                                const IBTK::Vector3d& box_br,
                                                const IBTK::Vector3d& box_tl,
                                                const libMesh::Point& n0,
                                                const libMesh::Point& n1)
{
    // If the line element is entirely contained within the box, then count as "intersected"
    if ((box_bl(0) <= n0(0) && n0(0) <= box_tr(0)) && (box_bl(1) <= n0(1) && n0(1) <= box_tr(1)) &&
        (box_bl(0) <= n1(0) && n1(0) <= box_tr(0)) && (box_bl(1) <= n1(1) && n1(1) <= box_tr(1)))
    {
        return true;
    }

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

bool
FESurfaceDistanceEvaluator::checkIntersection3D(const IBTK::Vector3d& box_center,
                                                const IBTK::Vector3d& box_half_dx,
                                                const IBTK::Vector3d& vert0,
                                                const IBTK::Vector3d& vert1,
                                                const IBTK::Vector3d& vert2)
{
    // If the tri element is entirely contained within the box, then count as "intersected"
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

    // Use the separating axis theorem to test the overlap between the tri element and the box
    double pmin, pmax, p0, p1, p2, rad, fex, fey, fez, a, b, fa, fb;
    IBTK::Vector3d v0, v1, v2, normal, e0, e1, e2;
    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    // Shift everything so that box_center is the origin
    v0 = vert0 - box_center;
    v1 = vert1 - box_center;
    v2 = vert2 - box_center;

    // Compute the triangle edges
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

IBTK::Vector3d
FESurfaceDistanceEvaluator::getClosestPoint3D(const IBTK::Vector3d& P,
                                              const libMesh::Point& vert0,
                                              const libMesh::Point& vert1,
                                              const libMesh::Point& vert2)
{
    IBTK::Vector3d B, E0, E1, D;
    B << vert0(0), vert0(1), vert0(2);
    E0 << vert1(0) - B(0), vert1(1) - B(1), vert1(2) - B(2);
    E1 << vert2(0) - B(0), vert2(1) - B(1), vert2(2) - B(2);
    D = B - P;

    // Compute the required dot products
    const double a = E0.dot(E0);
    const double b = E0.dot(E1);
    const double c = E1.dot(E1);
    const double d = E0.dot(D);
    const double e = E1.dot(D);
    const double det = a * c - b * b;
    double s = b * e - c * d;
    double t = b * d - a * e;

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
                        s = 1;
                    else
                        s = -d / a;
                }
                else
                {
                    s = 0;
                    if (e >= 0)
                    {
                        t = 0;
                    }
                    else
                    {
                        if (-e >= c)
                            t = 1;
                        else
                            t = -e / c;
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
                }
                else
                {
                    if (-e >= c)
                        t = 1;
                    else
                        t = -e / c;
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
                }
                else
                {
                    if (-d >= a)
                        s = 1;
                    else
                        s = -d / a;
                }
            }
            else
            {
                // region 0
                const double invDet = 1.0 / det;
                s = s * invDet;
                t = t * invDet;
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
                }
                else
                {
                    s = numer / denom;
                    t = 1.0 - s;
                }
            }
            else
            {
                // minimum on edge s=0
                s = 0.0;
                if (tmp1 <= 0.0)
                {
                    t = 1.0;
                }
                else
                {
                    if (e >= 0.0)
                        t = 0.0;
                    else
                        t = -e / c;
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
                    }
                    else
                    {
                        t = numer / denom;
                        s = 1.0 - t;
                    }
                }
                else
                {
                    t = 0.0;
                    if (tmp1 <= 0.0)
                    {
                        s = 1.0;
                    }
                    else
                    {
                        if (d >= 0.0)
                            s = 0.0;
                        else
                            s = -d / a;
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
                }
                else
                {
                    const double denom = a - 2.0 * b + c;
                    if (numer >= denom)
                    {
                        s = 1.0;
                        t = 0.0;
                    }
                    else
                    {
                        s = numer / denom;
                        t = 1.0 - s;
                    }
                }
            }
        }
    }

    return B + s * E0 + t * E1;
} // getClosestPoint3D

/////////////////////////////// PRIVATE //////////////////////////////////////

void
FESurfaceDistanceEvaluator::collectNeighboringPatchElements(int level_number)
{
    IBFEMethod::CoordinateMappingFcnData mapping = d_ibfe_method->getInitialCoordinateMappingFunction(d_part);
    const bool identity_mapping = !(mapping.fcn);

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
        // in the patch interior grown by the specified ghost cell width
        MeshBase::const_element_iterator el_it = d_bdry_mesh.active_elements_begin();
        const MeshBase::const_element_iterator el_end = d_bdry_mesh.active_elements_end();
        for (; el_it != el_end; ++el_it)
        {
            bool in_patch = false;
            Elem* const elem = *el_it;

            // Error checking for element type
            if (elem->type() != d_supported_elem_type)
            {
                TBOX_ERROR("FESurfaceDistanceEvaluator presently does not support elements of type "
                           << Utility::enum_to_string<ElemType>(elem->type()) << " for NDIM = " << NDIM
                           << ".\nSupported type is " << Utility::enum_to_string<ElemType>(d_supported_elem_type));
            }

            // First check the centroids
            const libMesh::Point& c0 = elem->centroid();
            libMesh::Point c = c0;
            if (!identity_mapping) mapping.fcn(c, c0, mapping.ctx);

            bool centroid_in_patch = true;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                in_patch = centroid_in_patch && (x_lower[d] <= c(d) && c(d) <= x_upper[d]);
            }

            // Next, check the nodes
            const unsigned int n_nodes = elem->n_nodes();
            for (unsigned int k = 0; k < n_nodes && !in_patch; ++k)
            {
                const libMesh::Point& n0 = elem->point(k);
                libMesh::Point n = n0;
                if (!identity_mapping) mapping.fcn(n, n0, mapping.ctx);

                bool node_in_patch = true;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    in_patch = node_in_patch && (x_lower[d] <= c(d) && c(d) <= x_upper[d]);
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
