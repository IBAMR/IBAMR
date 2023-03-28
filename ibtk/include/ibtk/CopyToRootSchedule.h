// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_CopyToRootSchedule
#define included_IBTK_CopyToRootSchedule

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "IntVector.h"
#include "PatchData.h"
#include "PatchLevel.h"
#include "tbox/Pointer.h"
#include "tbox/Schedule.h"

#include <string>
#include <vector>

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
                       std::vector<int> src_patch_data_idxs);

    /*!
     * \brief Destructor
     */
    ~CopyToRootSchedule() = default;

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
    CopyToRootSchedule() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CopyToRootSchedule(const CopyToRootSchedule& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CopyToRootSchedule& operator=(const CopyToRootSchedule& that) = delete;

    void commonClassCtor();

    const int d_root_proc;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > d_patch_level;
    const std::vector<int> d_src_patch_data_idxs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::hier::PatchData<NDIM> > > d_root_patch_data;
    SAMRAI::tbox::Schedule d_schedule;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CopyToRootSchedule
