// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_IndexUtilities
#define included_IBTK_IndexUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"

#include <functional>
#include <vector>

namespace SAMRAI
{
namespace pdat
{
template <int DIM>
class CellIndex;
template <int DIM>
class SideIndex;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
template <typename T>
struct IndexOrder
{
    inline bool operator()(const T& lhs, const T& rhs) const
    {
        return (lhs(0) < rhs(0)
#if (NDIM > 1)
                || (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                || (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
#endif
#endif
        );
    }
};

using IndexFortranOrder = struct IndexOrder<SAMRAI::hier::IndexNd>;
using CellIndexFortranOrder = struct IndexOrder<SAMRAI::pdat::CellIndexNd>;

/*!
 * \brief Class IndexUtilities is a utility class that defines simple functions
 * such as conversion routines between physical coordinates and Cartesian index
 * space.
 */
class IndexUtilities
{
public:
    /*
     * \return The coarsened version of a cell-centered index.
     */
    static SAMRAI::hier::IndexNd coarsen(const SAMRAI::hier::IndexNd& i_fine, const SAMRAI::hier::IndexNd& ratio);

    /*
     * \return The refined version of a cell-centered index.
     */
    static SAMRAI::hier::IndexNd refine(const SAMRAI::hier::IndexNd& i_coarsen, const SAMRAI::hier::IndexNd& ratio);

    /*!
     * \return The cell index corresponding to location \p X relative
     * to \p x_lower and \p x_upper for the specified Cartesian grid
     * spacings \p dx and box extents \p ilower and \p iupper.
     *
     * \note Because of round-off error in floating point, this routine
     * cannot guarantee that the same spatial location X is assigned the
     * same cell index for different patches.  To obtain a unique index,
     * use the globalized version.
     *
     * \see SAMRAI::geom::CartesianPatchGeometry
     * \see SAMRAI::geom::CartesianPatchGeometry
     */
    template <class DoubleArray>
    static SAMRAI::hier::IndexNd getCellIndex(const DoubleArray& X,
                                              const double* x_lower,
                                              const double* x_upper,
                                              const double* dx,
                                              const SAMRAI::hier::IndexNd& ilower,
                                              const SAMRAI::hier::IndexNd& iupper);

    /*!
     * \return The cell index corresponding to location \p X relative
     * to the extents of the supplied Cartesian grid patch geometry and
     * patch box.
     *
     * \note Because of round-off error in floating point, this routine cannot
     * guarantee that the same spatial location X is assigned the same cell
     * index for different patches. To obtain a unique index, use the globalized
     * version (i.e., the one that uses a SAMRAI::geom::CartesianGridGeometry)
     * instead.
     *
     * \see SAMRAI::geom::CartesianPatchGeometry
     */
    template <class DoubleArray>
    static SAMRAI::hier::IndexNd
    getCellIndex(const DoubleArray& X,
                 const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianPatchGeometryNd>& patch_geom,
                 const SAMRAI::hier::BoxNd& patch_box);

    /*!
     * \return The cell index corresponding to location \p X relative
     * to the corner of the computational domain specified by the grid
     * geometry object.
     *
     * \see SAMRAI::geom::CartesianGridGeometry
     */
    template <class DoubleArray>
    static SAMRAI::hier::IndexNd
    getCellIndex(const DoubleArray& X,
                 const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd>& grid_geom,
                 const SAMRAI::hier::IntVectorNd& ratio);

    /*!
     * \return The cell index corresponding to location \p X relative to the
     * corner of the computational domain specified by the grid geometry object.
     * Unlike getCellIndex(), this function assigns points on upper boundaries
     * to the cells inside the domain (instead of in the ghost region) which
     * makes it more suitable for partitioning purposes.
     *
     * This function does not shift points on periodic domains. It is up to the
     * caller to correctly shift @p X with respect to periodicity before calling
     * this function.
     *
     * \see SAMRAI::geom::CartesianGridGeometry
     */
    template <class DoubleArray>
    static SAMRAI::hier::IndexNd
    getAssignedCellIndex(const DoubleArray& X,
                         const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd>& grid_geom,
                         const SAMRAI::hier::IntVectorNd& ratio);

    /*!
     * \return The spatial coordinate of the given cell center.
     *
     * @param patch The patch on which the cell lives.
     *
     * @param cell_idx The CellIndex describing the current cell.
     */
    template <typename Vector>
    static Vector getCellCenter(const SAMRAI::hier::PatchNd& patch, const SAMRAI::pdat::CellIndexNd& cell_idx);

    /*!
     * \return The spatial coordinate of the given cell center.
     *
     * @param patch The patch on which the cell lives.
     *
     * @param cell_idx The CellIndex describing the current cell.
     */
    static IBTK::VectorNd getCellCenter(const SAMRAI::hier::PatchNd& patch, const SAMRAI::pdat::CellIndexNd& cell_idx);

    /*!
     * \return The spatial coordinate of the given cell center.
     *
     * @param grid_geom The grid geometry provides the extents of the computational domain.
     *
     * @param ratio Refinement ratio.
     *
     * @param cell_idx The CellIndex describing the current cell.
     */
    template <typename Vector>
    static Vector getCellCenter(const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd>& grid_geom,
                                const SAMRAI::hier::IntVectorNd& ratio,
                                const SAMRAI::pdat::CellIndexNd& cell_idx);

    /*!
     * \return The spatial coordinate of the given side center.
     *
     * @param patch The patch on which the cell lives.
     *
     * @param side_idx The SideIndex describing the current side.
     */
    template <typename Vector>
    static Vector getSideCenter(const SAMRAI::hier::PatchNd& patch, const SAMRAI::pdat::SideIndexNd& side_idx);

    /*!
     * \return The spatial coordinate of the given side center.
     *
     * @param patch The patch on which the cell lives.
     *
     * @param side_idx The SideIndex describing the current side.
     */
    static IBTK::VectorNd getSideCenter(const SAMRAI::hier::PatchNd& patch, const SAMRAI::pdat::SideIndexNd& side_idx);

    /*!
     * \return The spatial coordinate of the given side center.
     *
     * @param grid_geom The grid geometry provides the extents of the computational domain.
     *
     * @param ratio Refinement ratio.
     *
     * @param side_idx The SideIndex describing the current side.
     */
    template <typename Vector>
    static Vector getSideCenter(const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd>& grid_geom,
                                const SAMRAI::hier::IntVectorNd& ratio,
                                const SAMRAI::pdat::SideIndexNd& side_idx);

    /*!
     * \return The spatial coordinate of the given side center.
     *
     * @param grid_geom The grid geometry provides the extents of the computational domain.
     *
     * @param ratio Refinement ratio.
     *
     * @param side_idx The SideIndex describing the current side.
     */
    static IBTK::VectorNd getSideCenter(const SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometryNd>& grid_geom,
                                        const SAMRAI::hier::IntVectorNd& ratio,
                                        const SAMRAI::pdat::SideIndexNd& side_idx);

    /*!
     * \brief Map (i,j,k,d) index for a DOF defined for a SAMRAI variable
     * on a particular patch level to a positive integer. Such a mapping can
     * be useful for creating an application ordering (AO) between SAMRAI and
     * PETSc data structures.
     *
     * \param i AMR index representing the (i,j,k) array data index for a
     * variable on particular patch level.
     *
     * \param domain_lower Lower index of the domain for that patch level,
     * assuming that the patch level covers the entire domain.
     *
     * \param num_cells Number of data array cells for a patch level, which is
     * assumed to cover the entire domain. It can be thought of size of the
     * rectangular array that can store the variable data for the patch
     * level that covers the entire domain. For a cc-variable the number of data
     * array cells are same as patch level cells. For a sc-variable, the number
     * of cells for the normal component exceeds the patch level cells by 1 in
     * the normal direction.
     *
     * \param depth Data depth.
     *
     * \param offset Component offset. This is useful for getting unique values
     * for different components of the variable, e.g., a sc-variable. Different
     * components can have different depth.
     *
     * \param periodic_shift Periodic shift in each direction.
     *
     * \return The linear mapping of an AMR index to a continuous non-negative
     * integer space.
     */
    static int mapIndexToInteger(const SAMRAI::hier::IndexNd& i,
                                 const SAMRAI::hier::IndexNd& domain_lower,
                                 const SAMRAI::hier::IndexNd& num_cells,
                                 const int depth,
                                 const int offset = 0,
                                 const SAMRAI::hier::IntVectorNd& periodic_shift = SAMRAI::hier::IntVectorNd(0));

    /*!
     * \brief Partition a patch box into subdomains of size \em box_size
     * and into equal number of overlapping subdomains whose overlap region
     * is defined by \em overlap_size.
     *
     * \return Total number of subdomains the patch box got
     * partitioned intoin various dimensions.
     *
     * \note The overlap boxes are obtained from nonoverlap_boxes by growing
     * them suitably.
     */
    static SAMRAI::hier::IntVectorNd partitionPatchBox(std::vector<SAMRAI::hier::BoxNd>& overlap_boxes,
                                                       std::vector<SAMRAI::hier::BoxNd>& nonoverlap_boxes,
                                                       const SAMRAI::hier::BoxNd& patch_box,
                                                       const SAMRAI::hier::IntVectorNd& box_size,
                                                       const SAMRAI::hier::IntVectorNd& overlap_size);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    IndexUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    IndexUtilities(const IndexUtilities& from) = delete;

    /*!
     * \brief Unimplemented destructor.
     */
    ~IndexUtilities() = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IndexUtilities& operator=(const IndexUtilities& that) = delete;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/IndexUtilities-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_IndexUtilities
