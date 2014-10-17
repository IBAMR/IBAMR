// Filename: CopyToRootSchedule.h
// Created on 04 May 2011 by Boyce Griffith
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

#ifndef included_CopyToRootSchedule
#define included_CopyToRootSchedule

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "IntVector.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Schedule.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchData;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CopyToRootSchedule is used to communicate distributed patch data
 * to a unified patch data object on a root MPI process.
 *
 * \note This class is designed to be used with uniform grid data only.
 */
class CopyToRootSchedule
{
public:
    /*!
     * \brief Constructor
     */
    CopyToRootSchedule(int root_proc,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                       int src_patch_data_idx);

    /*!
     * \brief Constructor
     */
    CopyToRootSchedule(int root_proc,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                       const std::vector<int>& src_patch_data_idxs);

    /*!
     * \brief Destructor
     */
    ~CopyToRootSchedule();

    /*!
     * \brief Communicate data.
     */
    void communicate();

    /*!
     * \brief Get unified patch data.
     *
     * \note Patch data objects are allocated only on the root MPI process.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > >& getRootPatchData() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CopyToRootSchedule();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CopyToRootSchedule(const CopyToRootSchedule& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CopyToRootSchedule& operator=(const CopyToRootSchedule& that);

    void commonClassCtor();

    const int d_root_proc;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_patch_level;
    const std::vector<int> d_src_patch_data_idxs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > > d_root_patch_data;
    SAMRAI::tbox::Schedule d_schedule;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CopyToRootSchedule
