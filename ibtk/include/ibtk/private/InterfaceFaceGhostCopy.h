// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_InterfaceFaceGhostCopy
#define included_IBTK_InterfaceFaceGhostCopy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/ArrayDataTransaction.h>
#include <ibtk/private/HierarchyLevelsIdentity.h>

#include <tbox/Pointer.h>
#include <tbox/Schedule.h>

#include <Box.h>
#include <BoxList.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <VariableFillPattern.h>

#include <map>
#include <memory>
#include <vector>

namespace IBTK
{
namespace detail
{
/*!
 * \brief Class InterfaceFaceGhostCopy sets the faces of side-centered double precision data on the patch boundaries
 * along the coarse-fine interface of a fine level. It does two things:
 *
 * - It copies each such face to the ghost cells of every patch whose ghost cells contain it, on any process and across
 * periodic boundaries, within the fill pattern of the component. The same-level copy of a refine schedule reaches a
 * ghost face only if the ghost cells overlap the interior cells of the patch that owns the face, and it misses the
 * faces of patches that are separated by a coarse cell or touch only at a corner. Without this copy, such a ghost face
 * has its value from before the fill, or the interpolated coarse value, while the owner's face has the fine value. This
 * applies to every refine operator, including none.
 * - If protect_faces is true, it protects the faces against the refine operator of an in-place schedule, which
 * overwrites them with refined coarse values (IBSAMRAI2#67), although they are fine-level data that the fill must not
 * change. save() and restore() are a workaround for IBSAMRAI2#67; remove them once IBAMR requires an IBSAMRAI2 release
 * that fixes #67.
 *
 * Call save() before the fill and restore() for each patch right after the refine operator fills it. The copy is a
 * set of transactions that appendTransactions() adds to a schedule, which runs after the fill. save() and restore() do
 * nothing unless protect_faces is true.
 */
class InterfaceFaceGhostCopy
{
public:
    /*!
     * \brief Constructor. Finds the faces on the levels coarsest_ln to finest_ln of the hierarchy (level 0 has no
     * coarse-fine interface). The fill copies the values of the data src_data_idx to the data dst_data_idx, within
     * fill_pattern (null for the default pattern). Objects on the same levels with the same ghost width, fill pattern
     * shape, and protect_faces share the faces and copies while any of them exists.
     */
    InterfaceFaceGhostCopy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                           int coarsest_ln,
                           int finest_ln,
                           int dst_data_idx,
                           int src_data_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                           bool protect_faces);

    /*!
     * \brief Use other patch data indices and fill pattern. The ghost width of the data and the shape of the fill
     * pattern must not change; a change stops with an error.
     */
    void resetPatchDataIndices(int dst_data_idx,
                               int src_data_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern);

    /*!
     * \brief Whether the constructor was asked to protect the faces.
     */
    bool protectsFaces() const;

    /*!
     * \brief Record the values of the faces on level ln.
     */
    void save(int ln);

    /*!
     * \brief Set the faces of a patch to the recorded values. Call it after the refine operator has filled the patch,
     * and before anything else reads the faces.
     */
    void restore(const SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Append to a schedule the transactions that copy the faces of level ln to the ghost cells of the patches
     * whose ghost cells contain them, on any process. Every process must append the transactions of the same objects
     * in the same order.
     */
    void appendTransactions(int ln, SAMRAI::tbox::Schedule& schedule);

private:
    InterfaceFaceGhostCopy() = delete;
    InterfaceFaceGhostCopy(const InterfaceFaceGhostCopy& from) = delete;
    InterfaceFaceGhostCopy& operator=(const InterfaceFaceGhostCopy& that) = delete;

    /*!
     * \brief The faces and copies that the object finds: boxes, axes, patch numbers, and periodic offsets. The
     * geometry contains no data index and no values, so that objects with the same key share one geometry. The key is
     * the hierarchy levels, the ghost width, the shape of the fill pattern, and protect_faces.
     */
    struct Geometry
    {
        /*!
         * \brief A copy of faces of the patch src_patch to the ghost cells of the patch dst_patch, shifted by offset.
         */
        struct Copy
        {
            //! The patch numbers of the destination and the source.
            int dst_patch = -1;
            int src_patch = -1;
            //! The axis of the faces.
            int axis = 0;
            //! The periodic shift of the source.
            SAMRAI::hier::IntVector<NDIM> offset = SAMRAI::hier::IntVector<NDIM>(0);
            //! The cells of the destination that get the faces.
            SAMRAI::hier::Box<NDIM> box;
        };

        //! The geometry of one level.
        struct Level
        {
            //! For each local patch, the faces of each axis.
            std::map<int, std::vector<SAMRAI::hier::BoxList<NDIM>>> faces;
            //! The copies that involve a local patch, in the same order on every process.
            std::vector<Copy> copies;
        };

        //! The levels max(1, coarsest_ln) through finest_ln, recorded when the geometry was built.
        HierarchyLevelsIdentity identity;
        int first_ln = 1;
        int finest_ln = 0;
        SAMRAI::hier::IntVector<NDIM> ghost_width = SAMRAI::hier::IntVector<NDIM>(0);
        detail::FillPatternShape shape;
        bool protect_faces = false;
        //! The levels, indexed by level number.
        std::vector<Level> levels;
    };

    /*!
     * \brief The values of one level.
     */
    struct LevelState
    {
        //! For each local patch, the recorded values of the faces of each axis.
        std::map<int, std::vector<std::vector<double>>> saved_values;
    };

    /*!
     * \brief Return the geometry with this object's key: the registered one if it exists, else a new one, which is
     * registered unless the fill pattern has the kind OTHER.
     */
    std::shared_ptr<const Geometry>
    findGeometry(SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern,
                 int coarsest_ln,
                 int finest_ln) const;

    //! The patch hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;

    //! The patch data indices of the destination and of the source.
    int d_data_idx = invalid_index;
    int d_src_idx = invalid_index;

    //! The ghost width of the data and the shape of the fill pattern that the copies were found with. The shape stores
    //! a pattern of the kind OTHER.
    SAMRAI::hier::IntVector<NDIM> d_ghost_width = SAMRAI::hier::IntVector<NDIM>(0);
    detail::FillPatternShape d_shape;

    //! Whether to protect the faces against the refine operator.
    bool d_protect = false;

    //! The faces and copies, which objects with the same key share.
    std::shared_ptr<const Geometry> d_geometry;

    //! The values of each level.
    std::vector<LevelState> d_levels;
};
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_InterfaceFaceGhostCopy
