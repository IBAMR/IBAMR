// Filename: TGACoefs.h
// Created on 26 Aug 2007 by Boyce Griffith
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

#ifndef included_TGACoefs
#define included_TGACoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class TGACoefs is a lightweight utility class which is used to specify
 * the coefficients employed by the TGA implicit time integration algorithm.
 *
 * The TGA discretization is:
 *
 *     (I-nu2*dt*mu*L(t_int)) (I-nu1*dt*mu*L(t_new)) Q(n+1) = [(I+nu3*dt*mu*L(t_old)) Q(n) + (I+nu4*dt*mu*L) F(t_avg) dt]
 *
 * where
 *
 *    t_old = n dt
 *    t_new = (n+1) dt
 *    t_int = t_new - nu1*dt = t_old + (nu2+nu3)*dt
 *    t_avg = (t_new+t_old)/2 = t_old + (nu1+nu2+nu4)*dt
 *
 * Following McCorquodale et al., the coefficients for the TGA discretization
 * are:
 *
 *     nu1 = (a - sqrt(a^2-4*a+2))/2
 *     nu2 = (a + sqrt(a^2-4*a+2))/2
 *     nu3 = (1-a)
 *     nu4 = (0.5-a)
 *
 * Note that by choosing a = 2 - sqrt(2), nu1 == nu2.
 *
 * Ref: McCorquodale, Colella, Johansen.  "A Cartesian grid embedded boundary
 * method for the heat equation on irregular domains." JCP 173, pp. 620-635
 * (2001)
 */
class TGACoefs
{
public:
    static const double a;
    static const double nu1;
    static const double nu2;
    static const double nu3;
    static const double nu4;

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    TGACoefs();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    TGACoefs(
        const TGACoefs& from);

    /*!
     * \brief Destructor.
     */
    ~TGACoefs();

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    TGACoefs& operator=(
        const TGACoefs& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "TGACoefs.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_TGACoefs
