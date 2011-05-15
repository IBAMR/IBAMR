// Filename: PETScMatUtilities.h
// Created on 24 Aug 2010 by Boyce Griffith
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

#ifndef included_PETScMatUtilities
#define included_PETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscmat.h>

// SAMRAI INCLUDES
#include <CellData.h>
#include <CellVariable.h>
#include <Patch.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideVariable.h>

// BLITZ++ INCLUDES
#include <blitz/tinyvec.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScMatUtilities provides utility functions for <A
 * HREF="http://www-unix.mcs.anl.gov/petsc">PETSc</A> Mat objects.
 */
class PETScMatUtilities
{
public:
    /*!
     * \name Methods acting on SAMRAI::hier::Patch and SAMRAI::hier::PatchData
     * objects.
     */
    //\{

    /*!
     * \brief Construct a sequential PETSc Mat object corresponding to the
     * cell-centered Laplacian of a cell-centered variable restricted to a
     * single SAMRAI::hier::Patch.
     */
    static void
    constructPatchLaplaceOp(
        Mat& mat,
        const double C,
        const double D,
        SAMRAI::pdat::CellData<NDIM,double>& src_data,
        SAMRAI::pdat::CellData<NDIM,double>& dst_data,
        SAMRAI::hier::Patch<NDIM>& patch);

    /*!
     * \brief Construct sequential PETSc Mat objects corresponding to the
     * side-centered Laplacian of a side-centered variable restricted to a
     * single SAMRAI::hier::Patch.
     */
    static void
    constructPatchLaplaceOps(
        blitz::TinyVector<Mat,NDIM>& mats,
        const double C,
        const double D,
        SAMRAI::pdat::SideData<NDIM,double>& src_data,
        SAMRAI::pdat::SideData<NDIM,double>& dst_data,
        SAMRAI::hier::Patch<NDIM>& patch);

    //\}

    /*!
     * \name Methods acting on SAMRAI::hier::PatchLevel and
     * SAMRAI::hier::Variable objects.
     */
    //\{

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * cell-centered Laplacian of a cell-centered variable restricted to a
     * single SAMRAI::hier::PatchLevel.
     *
     * \note Rows of this operator corresponding to slaved ghost DOFs are
     * constrained to equal their master value.
     */
    static void
    constructPatchLevelLaplaceOp(
        Mat& mat,
        const double C,
        const double D,
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > data_var,
        const int dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,int> > dof_index_var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > dof_index_fill=NULL);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * side-centered Laplacian of a side-centered variable restricted to a
     * single SAMRAI::hier::PatchLevel.
     *
     * \note Rows of this operator corresponding to slaved ghost DOFs are
     * constrained to equal their master value.
     */
    static void
    constructPatchLevelLaplaceOp(
        Mat& mat,
        const double C,
        const double D,
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > data_var,
        const int dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,int> > dof_index_var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > dof_index_fill=NULL);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the IB
     * interpolation operator for side-centered variable restricted to a single
     * SAMRAI::hier::PatchLevel.
     *
     * \note Rows of this operator corresponding to slaved ghost DOFs are left
     * empty.
     */
    static void
    constructPatchLevelInterpOp(
        Mat& mat,
        Vec& X_vec,
        const int data_idx,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > data_var,
        const int dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,int> > dof_index_var,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > dof_index_fill=NULL);

    //\}

protected:

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScMatUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScMatUtilities(
        const PETScMatUtilities& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScMatUtilities&
    operator=(
        const PETScMatUtilities& that);

    /*!
     * \brief Construct a sequential PETSc Mat object corresponding to the
     * Laplacian on a box.
     */
    static void
    constructBoxLaplaceOp(
        Mat& mat,
        const double C,
        const double D,
        const SAMRAI::hier::Box<NDIM>& src_ghost_box,
        const SAMRAI::hier::Box<NDIM>& dst_ghost_box,
        const SAMRAI::hier::Box<NDIM>& interior_box,
        const int data_depth,
        const double* const dx);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include "PETScMatUtilities.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScMatUtilities
