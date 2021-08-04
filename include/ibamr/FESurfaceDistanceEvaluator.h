// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#ifndef included_FESurfaceDistanceEvaluator
#define included_FESurfaceDistanceEvaluator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#ifdef IBAMR_HAVE_LIBMESH

#include "ibamr/IBFEMethod.h"

#include "ibtk/IndexUtilities.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/libmesh_utilities.h"

#include "PatchHierarchy.h"
#include "tbox/Pointer.h"

#include "libmesh/boundary_mesh.h"

#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief class FESurfaceDistanceEvaluator is a utility class which is used to identify
 * which line (triangle) in 2D (3D) elements are intersecting which grid cells.
 */
class FESurfaceDistanceEvaluator
{
public:
    static const double s_large_distance;

    /*!
     * \brief The only constructor of this class.
     */
    FESurfaceDistanceEvaluator(std::string object_name,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                               const libMesh::Mesh& mesh,
                               const libMesh::BoundaryMesh& bdry_mesh,
                               const int gcw = 1,
                               bool use_extracted_bdry_mesh = true);

    /*!
     * \brief Default dstructor.
     */
    ~FESurfaceDistanceEvaluator() = default;

    /*!
     * \brief Map the triangles intersecting a particular grid cell.
     */
    void mapIntersections();

    /*!
     * \brief Compute the face normal of the surface elements.
     */
    void calculateSurfaceNormals();

    /*!
     * \brief Get the map maintaining triangle-cell intersection and neighbors.
     */
    const std::map<SAMRAI::pdat::CellIndex<NDIM>, std::set<libMesh::Elem*>, IBTK::CellIndexFortranOrder>&
    getNeighborIntersectionsMap();

    /*!
     * \brief Compute the signed distance in the viscinity of the finite element mesh.
     */
    void computeSignedDistance(int n_idx, int d_idx);

    /*!
     * \brief Update the sign of the distance function away from the finite element mesh.
     */
    void updateSignAwayFromInterface(int d_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                                     double large_distance = s_large_distance);

    /*!
     * \brief Check whether the grown box and line element intersect in 2D.
     *
     *     (box_tl)   *--------------*  (box_tr)
     *                |              |
     *                |              |
     *                |              |
     *                |              |
     *                |              |
     *     (box_bl)   *--------------*  (box_br)
     */
    static bool checkIntersection2D(const IBTK::Vector3d& box_bl,
                                    const IBTK::Vector3d& box_tr,
                                    const IBTK::Vector3d& box_br,
                                    const IBTK::Vector3d& box_tl,
                                    const libMesh::Point& n0,
                                    const libMesh::Point& n1);

    /*!
     * \brief Check whether the grown box and triangle element intersect in 3D.
     *
     *                  ^   --------------
     *                  |  |              |                  vert0
     *  box_half_dx[1]  |  |              |                    *
     *                  |  |              |                   / \
     *                  ^  |      *       |                  /   \
     *                     |   box_center |                 /     \
     *                     |              |                /       \
     *                     |              |               /         \
     *                      --------------               *-----------*
     *                           <------->
     *                            box_half_dx[0]      vert1         vert2
     */
    static bool checkIntersection3D(const IBTK::Vector3d& box_center,
                                    const IBTK::Vector3d& box_half_dx,
                                    const IBTK::Vector3d& vert0,
                                    const IBTK::Vector3d& vert1,
                                    const IBTK::Vector3d& vert2);

private:
    /*!
     * \brief Default constructor is not implemented and should not be used.
     */
    FESurfaceDistanceEvaluator() = delete;

    /*!
     * \brief Default assignment operator is not implemented and should not be used.
     */
    FESurfaceDistanceEvaluator& operator=(const FESurfaceDistanceEvaluator& that) = delete;

    /*!
     * \brief Default copy constructor is not implemented and should not be used.
     */
    FESurfaceDistanceEvaluator(const FESurfaceDistanceEvaluator& from) = delete;

    /*!
     * \brief Collect all of the neighboring elements which are located within a local
     * Cartesian grid patch grown by the specified ghost cell width.
     *
     * In this method, the determination as to whether an element is local or
     * not is based on the physical location of its nodes and centroid.
     */
    void collectNeighboringPatchElements(int level_number);

    /*!
     * \brief Compute the closest point on a triangle for a given Eulerian point P and obtain
     * its corresponding angle-weighted pseudo-normal from the surface mesh elements.
     */
    std::pair<IBTK::Vector3d, IBTK::Vector3d> getClosestPointandAngleWeightedNormal3D(const IBTK::Vector3d& P,
                                                                                      libMesh::Elem* elem);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to Patch Hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_patch_hierarchy;

    /*!
     * Volume mesh object
     */
    const libMesh::Mesh& d_mesh;

    /*!
     * Boundary mesh object.
     */
    const libMesh::BoundaryMesh& d_bdry_mesh;

    /*!
     * The desired ghost cell width.
     */
    int d_gcw;

    /*!
     * Check if we are using volume extracted boundary mesh.
     */
    bool d_use_vol_extracted_bdry_mesh;

    /*!
     * The supported element type for this class.
     */
    libMesh::ElemType d_supported_elem_type;

    /*!
     * Data to manage mapping between boundary mesh elements and grid patches.
     */
    std::vector<std::vector<libMesh::Elem*> > d_active_neighbor_patch_bdry_elem_map;

    /*!
     * Map object keeping track of element-cell intersections as well as elements intersecting that cell
     * and its neighboring cells within the ghost cell width. Note that the elements contained in thie
     * map belong to the original solid mesh.
     */
    std::map<SAMRAI::pdat::CellIndex<NDIM>, std::set<libMesh::Elem*>, IBTK::CellIndexFortranOrder>
        d_cell_elem_neighbor_map;

    /*!
     * Map the node and the set of elements sharing this node.
     */
    std::map<libMesh::Node*, std::set<libMesh::Elem*> > d_node_to_elem;
    /*!
     * Map the edge and the set of elements sharing this edge.
     */
    std::map<std::pair<libMesh::Node*, libMesh::Node*>, std::set<libMesh::Elem*> > d_edge_to_elem;

    /*!
     * Map the element and its face normal.
     */
    std::map<libMesh::Elem*, IBTK::VectorNd> d_elem_face_normal;

    /*!
     * Object to create a bounding box for sign update sweeping algorithm.
     */
    SAMRAI::hier::Box<NDIM> d_large_struct_box;
};

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifdef IBAMR_HAVE_LIBMESH
#endif // #ifndef included_FESurfaceDistanceEvaluator
