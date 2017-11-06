// Filename: CartGridBodyForce.h
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

// INCLUDE GAURDS
#ifndef included_CartGridBodyForce
#define included_CartGridBodyForce

//////////////////////////////////INCLUDES///////////////////////////////////////////////

#include "Patch.h"
#include "PatchLevel.h"
#include "Variable.h"
#include "ibtk/CartGridFunction.h"
#include "tbox/Pointer.h"

namespace IBTK
{
/*!
        * \brief CartGridBodyForce class.
        *
        * This class is a concrete implementation of IBTK::CartGridFunction class and is used to set the
        * body force on patches for the fluid solver.
        */

class CartGridBodyForce : public IBTK::CartGridFunction
{
private:
    ///////////////////////////// PRIVATE DATA MEMBERS /////////////////////
    /*!
     * Patch index for body force var.
     */
    int d_body_force_idx;

    ////////////////////////////// PRIVATE MEMBER FUNCTIONS //////////////////
    /*!
     * \note The default ctor is not implemented and should not be used.
     */
    CartGridBodyForce();

    /*!
     * \note The default copy ctor is not implemented and should not be used.
     */
    CartGridBodyForce(const CartGridBodyForce& from);

    /*!
     * \note The default assignment operator is not implemented and should not be used.
     */
    CartGridBodyForce& operator=(const CartGridBodyForce& that);

    //////////////////////////// PUBLIC MEMBER FUNCTIONS ////////////////////////
public:
    /*!
     * \brief Constructor
     *
     * \param body_force_idx is the array index of the patch where an external body force for INS is calculated.
     */
    CartGridBodyForce(const int body_force_idx);

    /*!
     * Returns if the body force cartesian grid function is time dependent or not.
     */
    virtual bool isTimeDependent() const;

    /*!
     * This method copies the body force calculated a \textit{piori} on the  patch interior passed to this method.
     */
    virtual void setDataOnPatch(const int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                const double data_time,
                                const bool initial_time = false,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level =
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

}; // CartGridBodyForce

} // namespace IBTK
#endif //#ifndef included_CartGridBodyForce
