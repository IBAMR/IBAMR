// Filename: HierarchyMathOps.h
// Created on 11 Jun 2003 by Boyce Griffith
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

#ifndef included_HierarchyMathOps
#define included_HierarchyMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/PatchMathOps.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenOperator.h>
#include <CoarsenSchedule.h>
#include <FaceVariable.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchyFaceDataOpsReal.h>
#include <HierarchyNodeDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <NodeVariable.h>
#include <OuterfaceVariable.h>
#include <OutersideVariable.h>
#include <PatchHierarchy.h>
#include <PoissonSpecifications.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <tbox/Pointer.h>
#include <tbox/DescribedClass.h>

// C++ STDLIB INCLUDES
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class HierarchyMathOps provides functionality to perform
 * "composite-grid" mathematical operations on a range of levels in a
 * SAMRAI::hier::PatchHierarchy object.
 *
 * \note All specified variable descriptor indices must refer to
 * SAMRAI::hier::Variable / SAMRAI::hier::VariableContext pairs that have been
 * registered with the SAMRAI::hier::VariableDatabase.
 */
class HierarchyMathOps
    : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     *
     * Does nothing interesting.
     */
    HierarchyMathOps(
        const std::string& name,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_ln=-1,
        const int finest_ln=-1,
        const std::string& coarsen_op_name="CONSERVATIVE_COARSEN");

    /*!
     * \brief Destructor.
     */
    ~HierarchyMathOps();

    /*!
     * \name Methods to set the hierarchy and range of levels.
     */
    //\{

    /*!
     * \brief Reset the patch hierarchy over which operations occur.
     */
    void
    setPatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * \brief Reset range of patch levels over which operations occur.
     *
     * The specified levels must exist in the hierarchy or an assertion will
     * result.
     */
    void
    resetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //\}

    /*!
     * \name Methods to get cell weight data.
     */
    //\{

    /*!
     * \brief Access the SAMRAI::pdat::CellVariable that is used to store cell
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a cell is not refined in the hierarchy, its weight is set to its
     * volume.  If a cell is refined, its weight is set to zero.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
    getCellWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::CellVariable that is used to store cell weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a cell is not refined in the hierarchy, its weight is set to its
     * volume.  If a cell is refined, its weight is set to zero.
     */
    int
    getCellWeightPatchDescriptorIndex() const;

    /*!
     * \brief Access the SAMRAI::pdat::FaceVariable that is used to store face
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a face is not refined in the hierarchy, its weight is set to its
     * volume.  If a face is refined, its weight is set to zero.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >
    getFaceWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::FaceVariable that is used to store face weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a face is not refined in the hierarchy, its weight is set to its
     * volume.  If a face is refined, its weight is set to zero.
     */
    int
    getFaceWeightPatchDescriptorIndex() const;

    /*!
     * \brief Access the SAMRAI::pdat::SideVariable that is used to store side
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a side is not refined in the hierarchy, its weight is set to its
     * volume.  If a side is refined, its weight is set to zero.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
    getSideWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::SideVariable that is used to store side weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a side is not refined in the hierarchy, its weight is set to its
     * volume.  If a side is refined, its weight is set to zero.
     */
    int
    getSideWeightPatchDescriptorIndex() const;

    /*!
     * \brief Returns the volume of the physical domain.
     */
    double
    getVolumeOfPhysicalDomain() const;

    //\}

    //\{ \name Mathematical operations.

    /*!
     * \brief Set the name of the coarsen operator used for synchronizing the
     * coarse-fine boundary.
     */
    void
    setCoarsenOperatorName(
        const std::string& coarsen_op_name);

    /*!
     * \brief Compute the cell-centered curl of a cell-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    curl(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Compute the cell-centered curl of a face-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    curl(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Compute the cell-centered curl of a face-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    curl(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Compute the cell-centered curl of a side-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    curl(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Compute the cell-centered curl of a side-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    curl(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Compute the cell-centered divergence of a cell-centered vector
     * field using centered differences.
     *
     * Sets dst = alpha div src1 + beta src2.
     *
     * Compute the divergence of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    div(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute the cell-centered divergence of a face-centered normal
     * vector field using centered differences.
     *
     * Sets dst = alpha div src1 + beta src2.
     *
     * Compute the divergence of a vector field using centered differences.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized prior to performing the differencing.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    div(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const bool src1_cf_bdry_synch,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute the cell-centered divergence of a side-centered normal
     * vector field using centered differences.
     *
     * Sets dst = alpha div src1 + beta src2.
     *
     * Compute the divergence of a vector field using centered differences.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized prior to performing the differencing.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    div(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const bool src1_cf_bdry_synch,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized after performing the differencing.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const bool dst_cf_bdry_synch,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized after performing the differencing.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const bool dst_cf_bdry_synch,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized after performing the differencing.
     *
     * \note Homogeneous Neumann boundary conditions are explicitly enforced in
     * the case of non-grid aligned anisotropy.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const bool dst_cf_bdry_synch,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Compute the gradient of a scalar quantity using centered
     * differences.
     *
     * Sets dst = alpha grad src1 + beta src2.
     *
     * Compute the gradient of a scalar quantity using centered differences.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized after performing the differencing.
     *
     * \note Homogeneous Neumann boundary conditions are explicitly enforced in
     * the case of non-grid aligned anisotropy.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    grad(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const bool dst_cf_bdry_synch,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL,
        const int src1_depth=0);

    /*!
     * \brief Interpolate to a cell-centered vector/tensor field from a
     * face-centered normal vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When specified, coarse values
     * on each coarse-fine interface are synchronized prior to performing the
     * interpolation.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    interp(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time,
        const bool src_cf_bdry_synch);

    /*!
     * \brief Interpolate to a cell-centered vector/tensor field from a
     * side-centered normal vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When specified, coarse values
     * on each coarse-fine interface are synchronized prior to performing the
     * interpolation.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    interp(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time,
        const bool src_cf_bdry_synch);

    /*!
     * \brief Interpolate to a face-centered normal vector/tensor field from a
     * cell-centered vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When the interpolation occurs
     * over over multiple levels of the hierarchy, second order interpolation is
     * used at the coarse-fine interface to correct fine values along the
     * interface.  When specified, coarse values on each coarse-fine interface
     * are synchronized after performing the interpolation.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    interp(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const bool dst_cf_bdry_synch,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Interpolate to a side-centered normal vector/tensor field from a
     * cell-centered vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When the interpolation occurs
     * over multiple levels of the hierarchy, second order interpolation is used
     * at the coarse-fine interface to correct fine values along the interface.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized after performing the interpolation.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    interp(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const bool dst_cf_bdry_synch,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
        const double& src_ghost_fill_time);

    /*!
     * \brief Compute the Laplacian of a scalar quantity using centered
     * differences.
     *
     * Sets dst = C src1 + div D grad src1 + gamma src2.
     *
     * Compute the Laplacian of a scalar quantity using centered differences,
     * where C and D are specified by the SAMRAI::solv::PoissonSpecifications
     * object poisson_spec.  When the operation is being performed on multiple
     * levels of the hierarchy, the appropriately synchronized gradient and
     * divergence operators are employed to obtain a consistent discretization
     * of the Laplace operator.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    laplace(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& gamma=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute the Laplacian of a scalar quantity using centered
     * differences.
     *
     * Sets dst = C src1 + div D grad src1 + gamma src2.
     *
     * Compute the Laplacian of a scalar quantity using centered differences,
     * where C and D are specified by the SAMRAI::solv::PoissonSpecifications
     * object poisson_spec.
     *
     * \note The present implementation of this operator \em requires that
     * damping factor C and diffusivity D be spatially constant and
     * scalar-valued.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    laplace(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& gamma=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL);

    /*!
     * \brief Compute dst = alpha div coef ((grad src1) + (grad src1)^T) + beta
     * src1 + gamma src2, the variable coefficient generalized Laplacian of
     * src1.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    vc_laplace(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const double& beta,
        const int coef_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& coef_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src1_var,
        SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
        const double& src1_ghost_fill_time,
        const double& gamma=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src1_var,
        const int beta_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& beta_var=NULL,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0,
        const int beta_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src1_var,
        const int beta_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& beta_var=NULL,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0,
        const int beta_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src1_var,
        const int beta_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& beta_var=NULL,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0,
        const int beta_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const double& alpha,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src1_var,
        const double& beta=0.0,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMultiply(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& dst_var,
        const int alpha_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& alpha_var,
        const int src1_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src1_var,
        const int beta_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& beta_var=NULL,
        const int src2_idx=-1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& src2_var=NULL,
        const int dst_depth=0,
        const int src1_depth=0,
        const int src2_depth=0,
        const int alpha_depth=0,
        const int beta_depth=0);

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseL1Norm(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src_var);

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseL2Norm(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src_var);

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMaxNorm(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& src_var);

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseL1Norm(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src_var);

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseL2Norm(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src_var);

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void
    pointwiseMaxNorm(
        const int dst_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& dst_var,
        const int src_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM,double> >& src_var);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyMathOps();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyMathOps(
        const HierarchyMathOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyMathOps&
    operator=(
        const HierarchyMathOps& that);

    /*!
     * \brief Reset the coarsen operators.
     */
    void
    resetCoarsenOperators();

    /*!
     * \brief Reset the refine operators.
     */
    void
    resetRefineOperators();

    /*!
     * \brief Execute schedule for restricting Outerface data to the specified
     * level from the next finer level.
     */
    void
    xeqScheduleOuterfaceRestriction(
        const int dst_idx,
        const int src_idx,
        const int dst_ln);

    /*!
     * \brief Execute schedule for restricting Outerside data to the specified
     * level from the next finer level.
     */
    void
    xeqScheduleOutersideRestriction(
        const int dst_idx,
        const int src_idx,
        const int dst_ln);

    // Housekeeping.
    std::string d_object_name;

    // Patch hierarchy information.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln, d_finest_ln;

    // Scratch Variables.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > d_fc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::OuterfaceVariable<NDIM,double> > d_of_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::OutersideVariable<NDIM,double> > d_os_var;
    int d_fc_idx, d_sc_idx, d_of_idx, d_os_idx;

    // Communications operators, algorithms, and schedules.
    std::string d_coarsen_op_name;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_of_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_os_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_of_coarsen_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_os_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_of_coarsen_scheds;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_os_coarsen_scheds;

    // Hierarchy data operations.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM,double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;

    // Patch math operations.
    PatchMathOps d_patch_math_ops;

    // The cell weights are used to compute norms of data defined on the patch
    // hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > d_wgt_fc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_wgt_sc_var;
    int d_wgt_cc_idx, d_wgt_fc_idx, d_wgt_sc_idx;
    double d_volume;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/HierarchyMathOps.I>

//////////////////////////////////////////////////////////////////////////////

#endif  //#ifndef included_HierarchyMathOps
