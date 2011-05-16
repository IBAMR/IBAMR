// Filename: Interpolator.I
// Created on 16 Jun 2010 by Boyce Griffith
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

// IBTK INCLUDES
#include <ibtk/compiler_hints.h>

// C++ STDLIB INCLUDES
#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

template<typename XInputIterator,typename YInputIterator>
inline
Interpolator::Interpolator(
    XInputIterator x_first,
    XInputIterator x_last,
    YInputIterator y_first,
    YInputIterator y_last)
    : d_x(x_first,x_last),
      d_y(y_first,y_last)
{
    // intentionally blank
    return;
}// Interpolator

template<typename InputIterator>
inline
Interpolator::Interpolator(
    double x_first,
    double x_step,
    InputIterator y_first,
    InputIterator y_last)
    : d_y(y_first,y_last)
{
    d_x.reserve(d_y.size());
    d_x.push_back(x_first);
    for (unsigned int k = 1; k < d_y.size(); ++k)
    {
        d_x.push_back(d_x.back()+x_step);
    }
    return;
}// Interpolator

inline
Interpolator::~Interpolator()
{
    // intentionally blank
    return;
}// ~Interpolator

inline double
Interpolator::operator()(
    const double& x) const
{
    std::vector<double>::const_iterator x_it = std::lower_bound(d_x.begin(), d_x.end(), x);
    if (UNLIKELY(x_it == d_x.begin())) return d_x.front();
    if (UNLIKELY(x_it == d_x.end())) return d_x.back();

    const double x_up = *x_it;
    const double y_up = *(d_y.begin() + static_cast<int>(x_it-d_x.begin()));
    --x_it;
    const double x_lo = *x_it;
    const double y_lo = *(d_y.begin() + static_cast<int>(x_it-d_x.begin()));

    return y_lo + (y_up-y_lo)*(x-x_lo)/(x_up-x_lo);
}// operator()

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////
