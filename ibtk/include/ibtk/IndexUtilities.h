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

#include "CellIndex.h"
#include "Index.h"
#include "IntVector.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
struct CellIndexFortranOrder : std::binary_function<SAMRAI::pdat::CellIndex<NDIM>, SAMRAI::pdat::CellIndex<NDIM>, bool>
{
    inline bool operator()(const SAMRAI::pdat::CellIndex<NDIM>& lhs, const SAMRAI::pdat::CellIndex<NDIM>& rhs) const
    {
        return (lhs(0) < rhs(0)
#if (NDIM > 1)
                ||
                (lhs(0) == rhs(0) && lhs(1) < rhs(1))
#if (NDIM > 2)
                ||
                (lhs(0) == rhs(0) && lhs(1) == rhs(1) && lhs(2) < rhs(2))
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
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/IndexUtilities-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IndexUtilities
