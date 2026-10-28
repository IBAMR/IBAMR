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

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

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
 * Call this function before creating Hypre objects. It calls HYPRE_Initialize() with Hypre 2.29 and newer and
 * HYPRE_Init() with Hypre 2.21 through 2.28. It does nothing with older versions of Hypre.
 */
void initialize_hypre();

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

/*!
 * \brief Set to zero the entries of a cell-centered matrix that couple a cell of a patch to a cell of the physical
 * domain that is not a cell of the patch level, which is a cell on the coarse side of a coarse-fine boundary.
 *
 * hypre's multigrid solvers form the matrices of the coarser grids from such entries and converge only if they are
 * zero. Cells of other patches of the level, including across a periodic boundary, are cells of the level. Entries
 * that couple to cells outside the physical domain are left to the boundary condition treatment. The diagonal entry is
 * unchanged, and PoissonUtilities::adjustRHSAtCoarseFineBoundary() moves the dropped terms to the right-hand side.
 * This function does nothing on the coarsest level.
 *
 * \param[in,out] matrix_coefficients Matrix coefficients on the box of a patch of \p level with one depth per entry of
 * \p stencil and no ghost cells.
 * \param[in] level Patch level that contains the patch.
 * \param[in] stencil Offsets of the cells that the entries couple to.
 */
void clearOffLevelMatrixEntries(SAMRAI::pdat::CellData<NDIM, double>& matrix_coefficients,
                                const SAMRAI::hier::PatchLevel<NDIM>& level,
                                const std::vector<SAMRAI::hier::Index<NDIM>>& stencil);

/*!
 * \brief Set to zero the entries of a side-centered matrix that couple a side of a patch to a side of the same
 * component that lies in the physical domain but is not a side of the patch level.
 *
 * The sides of the level are the sides of the cells of its patches, including the sides that patches share. See the
 * cell-centered version for the reason.
 *
 * \param[in,out] matrix_coefficients Matrix coefficients on the sides of the box of a patch of \p level with one depth
 * per entry of \p stencil and no ghost sides.
 * \param[in] level Patch level that contains the patch.
 * \param[in] stencil Offsets of the sides that the entries couple to, along the axes.
 */
void clearOffLevelMatrixEntries(SAMRAI::pdat::SideData<NDIM, double>& matrix_coefficients,
                                const SAMRAI::hier::PatchLevel<NDIM>& level,
                                const std::vector<SAMRAI::hier::Index<NDIM>>& stencil);
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_solver_utilities
