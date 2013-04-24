// Filename: Vector.h
// Created on 27 Jan 2011 by Boyce Griffith
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

#ifndef included_IBTK_Vector
#define included_IBTK_Vector

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "boost/array.hpp"     // IBTK pragma: export
#include "boost/operators.hpp" // IBTK pragma: export

//////////////////////////////////////////////////////////////////////////////

namespace IBTK
{
/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * Class Vector is a simple fixed-sized vector class similar to class
 * blitz::TinyVector<T,N>.
 */
template <class T,size_t N>
class Vector
    : public boost::array<T,N>
    , boost::addable< Vector<T,N>         // Vector + Vector
                      , boost::subtractable< Vector<T,N>    // Vector - Vector
                                             , boost::dividable< Vector<T,N>, T    // Vector / T
                                                                 , boost::multipliable< Vector<T,N>, T // Vector * T, T * Vector
                                                                                        > > > >
{
public:
    inline Vector() { }

    inline Vector(const T& v) {
        this->assign(v);
    }

    inline Vector(const T& v0,const T& v1) {
        (*this)[0] = v0;
        (*this)[1] = v1;
    }

    inline Vector(const T& v0,const T& v1,const T& v2) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3,const T& v4) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
        (*this)[4] = v4;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3,const T& v4,const T& v5) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
        (*this)[4] = v4;
        (*this)[5] = v5;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
        (*this)[4] = v4;
        (*this)[5] = v5;
        (*this)[6] = v6;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
        (*this)[4] = v4;
        (*this)[5] = v5;
        (*this)[6] = v6;
        (*this)[7] = v7;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7,const T& v8) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
        (*this)[4] = v4;
        (*this)[5] = v5;
        (*this)[6] = v6;
        (*this)[7] = v7;
        (*this)[8] = v8;
    }

    inline Vector(const T& v0,const T& v1,const T& v2,const T& v3,const T& v4,const T& v5,const T& v6,const T& v7,const T& v8,const T& v9) {
        (*this)[0] = v0;
        (*this)[1] = v1;
        (*this)[2] = v2;
        (*this)[3] = v3;
        (*this)[4] = v4;
        (*this)[5] = v5;
        (*this)[6] = v6;
        (*this)[7] = v7;
        (*this)[8] = v8;
        (*this)[9] = v9;
    }

    inline Vector<T,N>& operator+=(const Vector<T,N>& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] += v[k];
        }
        return *this;
    }

    inline Vector<T,N>& operator+=(const T& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] += v;
        }
        return *this;
    }

    inline Vector<T,N>& operator-=(const Vector<T,N>& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] -= v[k];
        }
        return *this;
    }

    inline Vector<T,N>& operator-=(const T& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] -= v;
        }
        return *this;
    }

    inline Vector<T,N>& operator*=(const Vector<T,N>& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] *= v[k];
        }
        return *this;
    }

    inline Vector<T,N>& operator*=(const T& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] *= v;
        }
        return *this;
    }

    inline Vector<T,N>& operator/=(const Vector<T,N>& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] /= v[k];
        }
        return *this;
    }

    inline Vector<T,N>& operator/=(const T& v) {
        for (size_t k = 0; k < N; ++k)
        {
            (*this)[k] /= v;
        }
        return *this;
    }

};

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

template <class T>
inline Vector<T,3>
cross(const Vector<T,3>& x, const Vector<T,3>& y)
{
    Vector<T,3> x_cross_y;
    x_cross_y[0] = x[1]*y[2] - x[2]*y[1];
    x_cross_y[1] = x[2]*y[0] - x[0]*y[2];
    x_cross_y[2] = x[0]*y[1] - x[1]*y[0];
    return x_cross_y;
}// cross

template <class T,size_t N>
inline T
dot(const Vector<T,N>& x, const Vector<T,N>& y)
{
    T x_dot_y = 0;
    for (size_t k = 0; k < N; ++k)
    {
        x_dot_y += x[k] * y[k];
    }
    return x_dot_y;
}// norm

template <size_t N>
inline double
norm(const Vector<double,N> x)
{
    return sqrt(dot(x,x));
}// norm

template <size_t N>
inline float
norm(const Vector<float,N> x)
{
    return sqrt(dot(x,x));
}// norm

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_Vector
