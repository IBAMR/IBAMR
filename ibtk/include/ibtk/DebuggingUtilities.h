// Filename: DebuggingUtilities.h
// Created on 12 Dec 2008 by Boyce Griffith
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

#ifndef included_DebuggingUtilities
#define included_DebuggingUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "tbox/Pointer.h"

namespace IBTK
{
class LData;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class DebuggingUtilities provides debugging functionality.
 */
class DebuggingUtilities
{
public:
    /*!
     * \brief Check a cell-centered variable for NaN or unusually large values.
     */
    static bool checkCellDataForNaNs(int patch_data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     bool interior_only = true,
                                     int coarsest_ln = -1,
                                     int finest_ln = -1);

    /*!
     * \brief Check a face-centered variable for NaN or unusually large values.
     */
    static bool checkFaceDataForNaNs(int patch_data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     bool interior_only = true,
                                     int coarsest_ln = -1,
                                     int finest_ln = -1);

    /*!
     * \brief Check a node-centered variable for NaN or unusually large values.
     */
    static bool checkNodeDataForNaNs(int patch_data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     bool interior_only = true,
                                     int coarsest_ln = -1,
                                     int finest_ln = -1);

    /*!
     * \brief Check a side-centered variable for NaN or unusually large values.
     */
    static bool checkSideDataForNaNs(int patch_data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                     bool interior_only = true,
                                     int coarsest_ln = -1,
                                     int finest_ln = -1);

    /*!
     * \brief Save the local portion of a cell-centered variable to disk.
     */
    static void saveCellData(int patch_data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const std::string& filename,
                             const std::string& dirname);

    /*!
     * \brief Save the local portion of a face-centered variable to disk.
     */
    static void saveFaceData(int patch_data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const std::string& filename,
                             const std::string& dirname);

    /*!
     * \brief Save the local portion of a node-centered variable to disk.
     */
    static void saveNodeData(int patch_data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const std::string& filename,
                             const std::string& dirname);

    /*!
     * \brief Save the local portion of a side-centered variable to disk.
     */
    static void saveSideData(int patch_data_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                             const std::string& filename,
                             const std::string& dirname);

    /*!
     * \brief Save the local portion of a Lagrangian variable to disk.
     */
    static void saveLagrangianData(SAMRAI::tbox::Pointer<LData> lag_data,
                                   bool save_ghost_nodes,
                                   const std::string& filename,
                                   const std::string& dirname);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    DebuggingUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    DebuggingUtilities(const DebuggingUtilities& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    DebuggingUtilities& operator=(const DebuggingUtilities& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_DebuggingUtilities
