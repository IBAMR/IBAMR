// Filename: PatchVecCellDataBasicOps.h
// Created on 09 Apr 2010 by Boyce Griffith
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

#ifndef included_PatchVecCellDataBasicOps
#define included_PatchVecCellDataBasicOps

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/VecCellData.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class PatchVecCellDataBasicOps provides access to a collection of basic
 * numerical operations that may be applied to numerical cell-centered patch
 * data.  These operations include simple arithmetic operations as well as min
 * and max, etc.  The primary intent of this class is to provide the interface
 * to these standard operations for an PatchVecCellDataOps object which provides
 * access to a complete set of operations that may be used to manipulate
 * cell-centered patch data objects.  Each member function accepts a box
 * argument indicating the region of index space on which the operation should
 * be performed.  The operation will be performed on the intersection of this
 * box and those boxes corresponding to the patch data objects involved.
 *
 * These operations typically apply only to the numerical standard built-in
 * types, such as double, float, and int, and the complex type (which may or may
 * not be a built-in type depending on the C++ compiler).  Thus, this templated
 * class should only be used to instantiate objects with those types as the
 * template parameter.  None of the operations are implemented for any other
 * type.
 */
template<class TYPE>
class PatchVecCellDataBasicOps
{
public:
    /*!
     * Empty constructor and destructor.
     */
    PatchVecCellDataBasicOps();

    ~PatchVecCellDataBasicOps<TYPE>();

    /*!
     * Set dst = alpha * src, elementwise.
     */
    void
    scale(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const TYPE& alpha,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = src + alpha, elementwise.
     */
    void
    addScalar(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src,
        const TYPE& alpha,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = src1 + src2, elementwise.
     */
    void
    add(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = src1 - src2, elementwise.
     */
    void
    subtract(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = src1 * src2, elementwise.
     */
    void
    multiply(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = src1 / src2, elementwise.  No check for division by zero.
     */
    void
    divide(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = 1 / src, elementwise.  No check for division by zero.
     */
    void
    reciprocal(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = alpha * src1 + beta * src2, elementwise.
     */
    void
    linearSum(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const TYPE& alpha,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const TYPE& beta,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = alpha * src1 + src2, elementwise.
     */
    void
    axpy(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const TYPE& alpha,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Set dst = alpha * src1 - src2, elementwise.
     */
    void
    axmy(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const TYPE& alpha,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src1,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src2,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Return the minimum patch data component entry  When the data is
     * complex, the result is the data element with the smallest norm.
     */
    TYPE
    min(
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& data,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Return the maximum patch data component entry  When the data is
     * complex, the result is the data element with the largest norm.
     */
    TYPE
    max(
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& data,
        const SAMRAI::hier::Box<NDIM>& box) const;

private:
    // The following are not implemented:
    PatchVecCellDataBasicOps(const PatchVecCellDataBasicOps<TYPE>&);
    void operator=(const PatchVecCellDataBasicOps<TYPE>&);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/PatchVecCellDataBasicOps.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PatchVecCellDataBasicOps
