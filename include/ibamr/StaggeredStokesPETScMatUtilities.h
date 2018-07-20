// Filename: StaggeredStokesPETScMatUtilities.h
// Created on 03 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBAMR_StaggeredStokesPETScMatUtilities
#define included_IBAMR_StaggeredStokesPETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "PoissonSpecifications.h"
#include "petscao.h"
#include "petscmat.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
template <int DIM>
class CoarseFineBoundary;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesPETScMatUtilities provides utility functions for
 * <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Mat objects.
 */
class StaggeredStokesPETScMatUtilities
{
public:
    /*!
     * \name Methods acting on SAMRAI::hier::PatchLevel and
     * SAMRAI::hier::Variable objects.
     */
    //\{

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to a MAC
     * discretization of the time-dependent incompressible Stokes equations on a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelMACStokesOp(Mat& mat,
                                               const SAMRAI::solv::PoissonSpecifications& u_problem_coefs,
                                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                               double data_time,
                                               const std::vector<int>& num_dofs_per_proc,
                                               int u_dof_index_idx,
                                               int p_dof_index_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * additive Schwarz method.
     */
    static void
    constructPatchLevelASMSubdomains(std::vector<std::set<int> >& is_overlap,
                                     std::vector<std::set<int> >& is_nonoverlap,
                                     const SAMRAI::hier::IntVector<NDIM>& box_size,
                                     const SAMRAI::hier::IntVector<NDIM>& overlap_size,
                                     const std::vector<int>& num_dofs_per_proc,
                                     int u_dof_index_idx,
                                     int p_dof_index_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM> > cf_boundary);

    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * PCFieldSplit preconditioner.
     */
    static void constructPatchLevelFields(std::vector<std::set<int> >& is_field,
                                          std::vector<std::string>& is_field_name,
                                          const std::vector<int>& num_dofs_per_proc,
                                          int u_dof_index_idx,
                                          int p_dof_index_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to data
     * prolongation from a coarser level to a finer level.
     */
    static void constructProlongationOp(Mat& mat,
                                        const std::string& u_op_type,
                                        const std::string& p_op_type,
                                        int u_dof_index_idx,
                                        int p_dof_index_idx,
                                        const std::vector<int>& num_fine_dofs_per_proc,
                                        const std::vector<int>& num_coarse_dofs_per_proc,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > fine_patch_level,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarse_patch_level,
                                        const AO& coarse_level_ao,
                                        const int u_coarse_ao_offset,
                                        const int p_coarse_ao_offset);

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesPETScMatUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesPETScMatUtilities(const StaggeredStokesPETScMatUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesPETScMatUtilities& operator=(const StaggeredStokesPETScMatUtilities& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StaggeredStokesPETScMatUtilities
