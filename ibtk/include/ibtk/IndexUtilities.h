// Filename: IndexUtilities.h
// Created on 06 Mar 2004 by Boyce Griffith
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

#ifndef included_IndexUtilities
#define included_IndexUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <functional>
#include <vector>

#include "Box.h"
#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace geom
{
template <int DIM>
class CartesianGridGeometry;
}
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
struct CellIndexFortranOrder : std::binary_function<SAMRAI::pdat::CellIndex<NDIM>, SAMRAI::pdat::CellIndex<NDIM>, bool>
{
    inline bool operator()(const SAMRAI::pdat::CellIndex<NDIM>& lhs, const SAMRAI::pdat::CellIndex<NDIM>& rhs) const
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

/*!
 * \brief Class IndexUtilities is a utility class that defines simple functions
 * such as conversion routines between physical coordinates and Cartesian index
 * space.
 */
class IndexUtilities
{
public:
    /*
     * \brief Initialize global values for IndexUtilities class.
     */
    static void init(SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geometry);

    /*
     * \return The coarsened version of a cell-centered index.
     */
    static SAMRAI::hier::Index<NDIM> coarsen(const SAMRAI::hier::Index<NDIM>& i_fine,
                                             const SAMRAI::hier::Index<NDIM>& ratio);

    /*
     * \return The refined version of a cell-centered index.
     */
    static SAMRAI::hier::Index<NDIM> refine(const SAMRAI::hier::Index<NDIM>& i_coarsen,
                                            const SAMRAI::hier::Index<NDIM>& ratio);

    /*!
     * \return The cell index corresponding to location \p X relative
     * to \p XLower and \p XUpper for the specified Cartesian grid
     * spacings \p dx and box extents \p ilower and \p iupper.
     *
     * \see SAMRAI::geom::CartesianPatchGeometry
     */
    template <class DoubleArray>
    static SAMRAI::hier::Index<NDIM> getCellIndex(const DoubleArray& X,
                                                  const double* x_lower,
                                                  const double* x_upper,
                                                  const double* dx,
                                                  const SAMRAI::hier::Index<NDIM>& ilower,
                                                  const SAMRAI::hier::Index<NDIM>& iupper);

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
     * \return The linear mapping of an AMR index to a continuous non-negative
     * integer space.
     */
    static int mapIndexToInteger(const SAMRAI::hier::Index<NDIM>& i,
                                 const SAMRAI::hier::Index<NDIM>& domain_lower,
                                 const SAMRAI::hier::Index<NDIM>& num_cells,
                                 const int depth,
                                 const int offset = 0);

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
    static SAMRAI::hier::IntVector<NDIM> partitionPatchBox(std::vector<SAMRAI::hier::Box<NDIM> >& overlap_boxes,
                                                           std::vector<SAMRAI::hier::Box<NDIM> >& nonoverlap_boxes,
                                                           const SAMRAI::hier::Box<NDIM>& patch_box,
                                                           const SAMRAI::hier::IntVector<NDIM>& box_size,
                                                           const SAMRAI::hier::IntVector<NDIM>& overlap_size);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    IndexUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    IndexUtilities(const IndexUtilities& from);

    /*!
     * \brief Unimplemented destructor.
     */
    ~IndexUtilities();

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IndexUtilities& operator=(const IndexUtilities& that);

    /*!
     * \brief Pointer to Cartesian grid geometry.
     */
    static SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > s_grid_geometry;

    /*!
     * \brief Grid spacing at the coarsest level.
     */
    static const double* s_dx0;

    /*!
     * \brief Domain lower and upper coordinates.
     */
    static const double *s_x_lower, *s_x_upper;

    /*!
     * \brief Domain lower and upper index.
     */
    static SAMRAI::hier::IntVector<NDIM> s_ilower, s_iupper;

    /*!
     * \brief Indicate if the IndexUtilities class has been initialized.
     */
    static bool s_is_initialized;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/IndexUtilities-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IndexUtilities
