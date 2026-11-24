// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2026 by the IBAMR developers
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

#include <petscksp.h>
#include <petscsnes.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <HYPRE_sstruct_ls.h>
#include <HYPRE_sstruct_mv.h>
#include <HYPRE_struct_ls.h>
#include <HYPRE_struct_mv.h>
IBTK_ENABLE_EXTRA_WARNINGS

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <CellData.h>
#include <SideData.h>

#include <array>
#include <vector>

namespace IBTK
{
class GeneralSolver;

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
 * \brief Make \a solver apply a fixed number of Richardson iterations with no convergence test if it is a
 * PETScLevelSolver and \a solver_db does not select a KSP type.
 *
 * Call this function before the solver state is initialized.
 */
void set_fixed_iteration_ksp_defaults(GeneralSolver* solver,
                                      const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& solver_db);

/*!
 * \brief Initialize Hypre if it is not already initialized, and then finalize it when PETSc is finalized.
 *
 * Call this function before creating Hypre objects. This function calls HYPRE_Initialize() with Hypre 2.29 and newer
 * and HYPRE_Init() with Hypre 2.21 through 2.28, and does nothing with older versions.
 */
void initialize_hypre();

/*!
 * \brief Helper function to convert SAMRAI indices to Hypre integers.
 *
 * \note Hypre can use 64 bit indices, but SAMRAI IntVectors are always 32.
 */
std::array<HYPRE_Int, NDIM> hypre_array(const SAMRAI::hier::Index<NDIM>& index);

/*!
 * \brief Return the periodic shift to give hypre in direction \a d for the grid of \a level: the periodic shift of the
 * domain if the level has cells at both ends of the domain in that direction, and zero otherwise. If
 * \a require_periodic_images, throw an error naming \a object_name unless those cells are periodic images of each
 * other, as the side-centered solvers require: hypre does not store the faces at the upper end of a periodic grid.
 */
int get_hypre_periodic_shift(const SAMRAI::hier::PatchLevel<NDIM>& level,
                             const std::string& object_name,
                             unsigned int d,
                             bool require_periodic_images);

/*!
 * \brief Copy data from a vector of Hypre vectors to SAMRAI cell centered data with depth equal to number of Hypre
 * vectors.
 *
 * \param[out] dst_data Reference to destination for data to be copied.
 * \param[in] vectors Vector of Hypre data to copy.
 * \param[in] box Box of cells to copy. Cells that \p dst_data does not store are skipped.
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
 * \param[in] box Box of cells to copy. All sides of these cells are copied, including the sides on the boundary
 * of the box. Sides that \p dst_data does not store are skipped.
 */
void copyFromHypre(SAMRAI::pdat::SideData<NDIM, double>& dst_data,
                   HYPRE_SStructVector vector,
                   const SAMRAI::hier::Box<NDIM>& box);

/*!
 * \brief Copy data from SAMRAI cell centered data to Hypre vectors.
 *
 * \param[out] vectors Reference to vector of Hypre vectors to be copied to.
 * \param[in] src_data Reference to cell centered data to be copied.
 * \param[in] box Box of cells to copy. Cells that \p src_data does not store are skipped.
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
 * \param[in] box Box of cells to copy. All sides of these cells are copied, including the sides on the boundary
 * of the box. Sides that \p src_data does not store are skipped.
 */
void copyToHypre(HYPRE_SStructVector& vector,
                 SAMRAI::pdat::SideData<NDIM, double>& src_data,
                 const SAMRAI::hier::Box<NDIM>& box);
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_solver_utilities
