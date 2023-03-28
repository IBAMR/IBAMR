// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_solver_utilities
#define included_IBTK_solver_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "petscksp.h"
#include "petscsnes.h"

IBTK_DISABLE_EXTRA_WARNINGS
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_sstruct_mv.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_struct_mv.h"
IBTK_ENABLE_EXTRA_WARNINGS

#include "CellData.h"
#include "SideData.h"

#include <array>
#include <vector>

namespace IBTK
{
/*!
 * \brief Report the KSPConvergedReason.
 */
void reportPETScKSPConvergedReason(const std::string& object_name, const KSPConvergedReason& reason, std::ostream& os);

/*!
 * \brief Report the SNESConvergedReason.
 */
void
reportPETScSNESConvergedReason(const std::string& object_name, const SNESConvergedReason& reason, std::ostream& os);

/*!
 * \brief Helper function to convert SAMRAI indices to Hypre integers.
 *
 * \note Hypre can use 64 bit indices, but SAMRAI IntVectors are always 32.
 */
std::array<HYPRE_Int, NDIM> hypre_array(const SAMRAI::hier::Index<NDIM>& index);

/*!
 * \brief Copy data from a vector of Hypre vectors to SAMRAI cell centered data with depth equal to number of Hypre
 * vectors.
 *
 * \param[out] dst_data Reference to destination for data to be copied.
 * \param[in] vectors Vector of Hypre data to copy.
 * \param[in] box Box over which to copy.
 */
void copyFromHypre(SAMRAI::pdat::CellData<NDIM, double>& dst_data,
                   const std::vector<HYPRE_StructVector>& vectors,
                   const SAMRAI::hier::Box<NDIM>& box);

/*!
 * \brief Copy data from a Hypre vector to SAMRAI side centered data.
 *
 * \note This function is specialized for cases when the Hypre vector has one part with number of variables equal to the
 * spatial dimension.
 *
 * \param[out] dst_data Reference to destination for data to be copied.
 * \param[in] vector Vector of Hypre data to copy
 * \param[in] box Box over which to copy.
 */
void copyFromHypre(SAMRAI::pdat::SideData<NDIM, double>& dst_data,
                   HYPRE_SStructVector vector,
                   const SAMRAI::hier::Box<NDIM>& box);

/*!
 * \brief Copy data from SAMRAI cell centered data to Hypre vectors.
 *
 * \param[out] vectors Reference to vector of Hypre vectors to be copied to.
 * \param[in] src_data Reference to cell centered data to be copied.
 * \param[in] box Box over which to copy.
 */
void copyToHypre(const std::vector<HYPRE_StructVector>& vectors,
                 SAMRAI::pdat::CellData<NDIM, double>& src_data,
                 const SAMRAI::hier::Box<NDIM>& box);

/*!
 * \brief Copy data from SAMRAI side centered data to a Hypre vector.
 *
 * \note This function is specialized for cases when the Hypre vector has one part with number of variables equal to the
 * spatial dimension.
 *
 * \param[out] vectors Reference to Hypre vector to be copied to.
 * \param[in] src_data Reference to side centered data to be copied.
 * \param[in] box Box over which to copy.
 */
void copyToHypre(HYPRE_SStructVector& vector,
                 SAMRAI::pdat::SideData<NDIM, double>& src_data,
                 const SAMRAI::hier::Box<NDIM>& box);
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_solver_utilities
