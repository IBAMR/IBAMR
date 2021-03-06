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

#ifndef included_IBTK_DebuggingUtilities
#define included_IBTK_DebuggingUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/Pointer.h"

#include <string>

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
    DebuggingUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    DebuggingUtilities(const DebuggingUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    DebuggingUtilities& operator=(const DebuggingUtilities& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_DebuggingUtilities
