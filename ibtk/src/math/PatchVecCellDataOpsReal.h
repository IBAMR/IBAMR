// Filename: PatchVecCellDataOpsReal.h
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

#ifndef included_PatchVecCellDataOpsReal
#define included_PatchVecCellDataOpsReal

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/PatchVecCellDataBasicOps.h>

// SAMRAI INCLUDES
#include <Patch.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * Class PatchVecCellDataOpsReal provides a collection of operations to
 * manipulate float and double numerical cell-centered patch data.
 *
 * Note that this templated class should only be used to instantiate objects
 * with double or float as the template parameter.
 *
 * \see PatchVecCellDataBasicOps
 */
template<class TYPE>
class PatchVecCellDataOpsReal
    : public SAMRAI::tbox::DescribedClass,
      public PatchVecCellDataBasicOps<TYPE>
{
public:
    /*!
     * Empty constructor and destructor.
     */
    PatchVecCellDataOpsReal();

    ~PatchVecCellDataOpsReal<TYPE>();

    /*!
     * Copy dst data to src data over given box.
     */
    void
    copyData(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& src,
        const SAMRAI::hier::Box<NDIM>& box) const;

    /*!
     * Swap pointers for patch data objects.  Objects are checked for
     * consistency of depth, box, and ghost box.
     */
    void
    swapData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        int data1_id,
        int data2_id) const;

    /*!
     * Print data entries over given box to given output stream.
     */
    void
    printData(
        const SAMRAI::tbox::Pointer<VecCellData<TYPE> >& data,
        const SAMRAI::hier::Box<NDIM>& box,
        std::ostream& s=SAMRAI::tbox::plog) const;

    /*!
     * Initialize data to given scalar over given box.
     */
    void
    setToScalar(
        SAMRAI::tbox::Pointer<VecCellData<TYPE> >& dst,
        const TYPE& alpha,
        const SAMRAI::hier::Box<NDIM>& box) const;

private:
    // The following are not implemented:
    PatchVecCellDataOpsReal(const PatchVecCellDataOpsReal<TYPE>&);
    void operator=(const PatchVecCellDataOpsReal<TYPE>&);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/PatchVecCellDataOpsReal.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PatchVecCellDataOpsReal
