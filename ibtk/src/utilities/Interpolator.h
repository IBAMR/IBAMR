// Filename: Interpolator.h
// Created on 30 Apr 2010 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

#ifndef included_Interpolator
#define included_Interpolator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class Interpolator interpolates a stream of values using linear
 * interpolation.
 *
 * \note Currently, only piecewise-linear interpolation is implemented.
 */
class Interpolator
{
public:
    /*!
     * \brief Constructor.
     *
     * \note The range of x values is required to be ascending with no
     * duplicates.
     */
    template<typename XInputIterator,typename YInputIterator>
    Interpolator(
        XInputIterator x_first,
        XInputIterator x_last,
        YInputIterator y_first,
        YInputIterator y_last);

    /*!
     * \brief Constructor.
     */
    template<typename InputIterator>
    Interpolator(
        double x_first,
        double x_step,
        InputIterator y_first,
        InputIterator y_last);

    /*!
     * \brief Destructor.
     */
    ~Interpolator();

    /*!
     * \brief Interpolate the data to the specified location x.
     *
     * \note If x < x_first, then the value y_first is returned.  If x > x_last,
     * then the value y_last is returned.
     */
    double
    operator()(
        double x) const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    Interpolator();

    /*!
     * \brief Copy constructor.
     *
     * \param from The value to copy to this object.
     *
     * \note This constructor is not implemented and should not be used.
     */
    Interpolator(
        const Interpolator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    Interpolator&
    operator=(
        const Interpolator& that);

    /*
     * The data to interpolate.
     */
    std::vector<double> d_x, d_y;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/Interpolator-inl.h"  // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_Interpolator
