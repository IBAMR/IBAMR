// Filename: CartGridFunction.h
// Created on 15 Mar 2004 by Boyce Griffith
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

#ifndef included_CartGridFunction
#define included_CartGridFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>

#include "PatchLevel.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CartGridFunction provides an abstract interface for objects for
 * evaluating functions to set values in SAMRAI::hier::PatchData objects.
 */
class CartGridFunction : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief The default constructor sets the name of the strategy object.
     */
    CartGridFunction(const std::string& object_name = "");

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~CartGridFunction();

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    virtual bool isTimeDependent() const = 0;

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * levels of the patch hierarchy using the virtual function
     * setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    virtual void setDataOnPatchHierarchy(int data_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                         double data_time,
                                         bool initial_time = false,
                                         int coarsest_ln = -1,
                                         int finest_ln = -1);

    /*!
     * \brief Evaluate the function on the patch interiors on the specified
     * level of the patch hierarchy using the virtual function setDataOnPatch().
     *
     * \see setDataOnPatch
     */
    virtual void setDataOnPatchLevel(int data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                     double data_time,
                                     bool initial_time = false);

    /*!
     * \brief Pure virtual function to evaluate the function on the patch
     * interior.
     */
    virtual void setDataOnPatch(int data_idx,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                double data_time,
                                bool initial_time = false,
                                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >
                                    patch_level = SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) = 0;

    //\}

protected:
    /*
     * The object name is used for error/warning reporting.
     */
    std::string d_object_name;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CartGridFunction(const CartGridFunction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CartGridFunction& operator=(const CartGridFunction& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CartGridFunction
