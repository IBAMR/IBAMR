// Filename: IndexUtilities-inl.h
// Created on 18 Jun 2005 by Boyce Griffith
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

#ifndef included_IndexUtilities_inl_h
#define included_IndexUtilities_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cmath>

#include "ibtk/IndexUtilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline SAMRAI::hier::Index<NDIM> IndexUtilities::coarsen(const SAMRAI::hier::Index<NDIM>& i_fine,
                                                         const SAMRAI::hier::Index<NDIM>& ratio)
{
    SAMRAI::hier::Index<NDIM> i_coarse;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        i_coarse(d) = i_fine(d) < 0 ? (i_fine(d) + 1) / ratio(d) - 1 : i_fine(d) / ratio(d);
    }
    return i_coarse;
} // coarsen

inline SAMRAI::hier::Index<NDIM> IndexUtilities::refine(const SAMRAI::hier::Index<NDIM>& i_coarse,
                                                        const SAMRAI::hier::Index<NDIM>& ratio)
{
    return i_coarse * ratio;
} // refine

template <class DoubleArray>
inline SAMRAI::hier::Index<NDIM> IndexUtilities::getCellIndex(const DoubleArray& X,
                                                              const double* const x_lower,
                                                              const double* const x_upper,
                                                              const double* const dx,
                                                              const SAMRAI::hier::Index<NDIM>& ilower,
                                                              const SAMRAI::hier::Index<NDIM>& iupper)
{
    // TODO: This expression guarantees consistency between neighboring patches, but it is still possible to get
    // inconsitent mappings on disjoint patches.
    SAMRAI::hier::Index<NDIM> idx;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        double dX_lower = X[d] - x_lower[d], dX_upper = X[d] - x_upper[d];
        if (abs(dX_lower) < abs(dX_upper))
        {
            idx(d) = ilower(d) + static_cast<int>(floor((dX_lower) / dx[d]));
        }
        else
        {
            idx(d) = iupper(d) + static_cast<int>(floor((dX_upper) / dx[d])) + 1;
        }
    }
    return idx;
} // getCellIndex

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IndexUtilities_inl_h
