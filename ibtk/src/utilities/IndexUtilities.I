// Filename: IndexUtilities.I
// Created on 18 Jun 2005 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

// C++ STDLIB INCLUDES
#include <cmath>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getCellIndex(
    const double* const X,
    const double* const XLower,
    const double* const XUpper,
    const double* const dx,
    const SAMRAI::hier::Index<NDIM>& ilower,
    const SAMRAI::hier::Index<NDIM>& iupper)
{
    SAMRAI::hier::Index<NDIM> idx(
        int( floor((X[0]-XLower[0])/dx[0]))+ilower(0)
#if (NDIM > 1)
        ,int(floor((X[1]-XLower[1])/dx[1]))+ilower(1)
#if (NDIM > 2)
        ,int(floor((X[2]-XLower[2])/dx[2]))+ilower(2)
#endif
#endif
                                  );
    return idx;
}// getCellIndex

inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getCellIndex(
    const std::vector<double>& X,
    const double* const XLower,
    const double* const XUpper,
    const double* const dx,
    const SAMRAI::hier::Index<NDIM>& ilower,
    const SAMRAI::hier::Index<NDIM>& iupper)
{
    return getCellIndex(&X[0],XLower,XUpper,dx,ilower,iupper);
}// getCellIndex

inline SAMRAI::hier::Index<NDIM>
IndexUtilities::getCellIndex(
    const blitz::TinyVector<double,NDIM>& X,
    const double* const XLower,
    const double* const XUpper,
    const double* const dx,
    const SAMRAI::hier::Index<NDIM>& ilower,
    const SAMRAI::hier::Index<NDIM>& iupper)
{
    return getCellIndex(X.data(),XLower,XUpper,dx,ilower,iupper);
}// getCellIndex

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
