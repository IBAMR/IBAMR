// Filename: HierarchyMathOps.h
// Created on 11 Jun 2003 by Boyce Griffith
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

#ifndef included_HierarchyMathOps
#define included_HierarchyMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyFaceDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "ibtk/PatchMathOps.h"

namespace SAMRAI
{
namespace pdat
{
template <class TYPE>
class EdgeVariable;
} // namespace pdat
} // namespace SAMRAI

namespace IBTK
{
class HierarchyGhostCellInterpolation;
} // namespace IBTK
namespace SAMRAI
{
namespace xfer
{

class CoarsenSchedule;
} // namespace xfer
} // namespace SAMRAI

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
{
public:
    /*!
     * \brief Constructor.
     *
     * Does nothing interesting.
     */
    HierarchyMathOps(const std::string& name,
                     boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
                     int coarsest_ln = -1,
                     int finest_ln = -1,
                     const std::string& coarsen_op_name = "CONSERVATIVE_COARSEN");

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
    void setPatchHierarchy(boost::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy);

    /*!
     * \brief Reset range of patch levels over which operations occur.
     *
     * The specified levels must exist in the hierarchy or an assertion will
     * result.
     */
    void resetLevels(int coarsest_ln, int finest_ln);

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
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > getCellWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::CellVariable that is used to store cell weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a cell is not refined in the hierarchy, its weight is set to its
     * volume.  If a cell is refined, its weight is set to zero.
     */
    int getCellWeightPatchDescriptorIndex() const;

    /*!
     * \brief Access the SAMRAI::pdat::FaceVariable that is used to store face
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a face is not refined in the hierarchy, its weight is set to its
     * volume.  If a face is refined, its weight is set to zero.
     */
    boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > getFaceWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::FaceVariable that is used to store face weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a face is not refined in the hierarchy, its weight is set to its
     * volume.  If a face is refined, its weight is set to zero.
     */
    int getFaceWeightPatchDescriptorIndex() const;

    /*!
     * \brief Access the SAMRAI::pdat::SideVariable that is used to store side
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a side is not refined in the hierarchy, its weight is set to its
     * volume.  If a side is refined, its weight is set to zero.
     */
    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > getSideWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::SideVariable that is used to store side weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a side is not refined in the hierarchy, its weight is set to its
     * volume.  If a side is refined, its weight is set to zero.
     */
    int getSideWeightPatchDescriptorIndex() const;

    /*!
     * \brief Returns the volume of the physical domain.
     */
    double getVolumeOfPhysicalDomain() const;

    //\}

    //\{ \name Mathematical operations.

    /*!
     * \brief Set the name of the coarsen operator used for synchronizing the
     * coarse-fine boundary.
     */
    void setCoarsenOperatorName(const std::string& coarsen_op_name);

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
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

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
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

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
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

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
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

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
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

    /*!
     * \brief Compute the node-centered curl of a side-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

    /*!
     * \brief Compute the edge-centered curl of a side-centered vector field
     * using centered differences.
     *
     * Sets dst = curl src.
     *
     * Compute the curl of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void curl(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::EdgeVariable<double> > dst_var,
              int src_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
              double src_ghost_fill_time);

    /*!
     * \brief Compute the side-centered rot of a node-centered scalar field
     * using centered differences.
     *
     * Sets dst = rot src.
     *
     * Compute the rot of a 2d scalar field defined on nodes using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void rot(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
             int src_idx,
             boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time);

    /*!
     * \brief Compute the side-centered rot of a node-centered scalar field
     * using centered differences.
     *
     * Sets dst = rot src.
     *
     * Compute the rot of a 2d scalar field defined on nodes using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void rot(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
             int src_idx,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time);

    /*!
     * \brief Compute the side-centered rot of a edge-centered vector field
     * using centered differences.
     *
     * Sets dst = rot src.
     *
     * Compute the rot of a 3d vector field defined on edges using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void rot(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
             int src_idx,
             boost::shared_ptr<SAMRAI::pdat::EdgeVariable<double> > src_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time);

    /*!
     * \brief Compute the side-centered rot of a side-centered vector field
     * using centered differences.
     *
     * Sets dst = rot src.
     *
     * Compute the rot of a 3d vector field defined on sides using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void rot(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
             int src_idx,
             boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time);

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
    void div(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
             double alpha,
             int src1_idx,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
             double src1_ghost_fill_time,
             double beta = 0.0,
             int src2_idx = -1,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
             int dst_depth = 0,
             int src2_depth = 0);

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
    void div(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
             double alpha,
             int src1_idx,
             boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src1_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
             double src1_ghost_fill_time,
             bool src1_cf_bdry_synch,
             double beta = 0.0,
             int src2_idx = -1,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
             int dst_depth = 0,
             int src2_depth = 0);

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
    void div(int dst_idx,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
             double alpha,
             int src1_idx,
             boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src1_var,
             boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
             double src1_ghost_fill_time,
             bool src1_cf_bdry_synch,
             double beta = 0.0,
             int src2_idx = -1,
             boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
             int dst_depth = 0,
             int src2_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
              double alpha,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
              bool dst_cf_bdry_synch,
              double alpha,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
              bool dst_cf_bdry_synch,
              double alpha,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
              int alpha_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > alpha_var,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
              int alpha_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > alpha_var,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
              bool dst_cf_bdry_synch,
              int alpha_idx,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > alpha_var,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void grad(int dst_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
              bool dst_cf_bdry_synch,
              int alpha_idx,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > alpha_var,
              int src1_idx,
              boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
              boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL,
              int src1_depth = 0);

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
    void interp(int dst_idx,
                boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                int src_idx,
                boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src_var,
                boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time,
                bool src_cf_bdry_synch);

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
    void interp(int dst_idx,
                boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                int src_idx,
                boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src_var,
                boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time,
                bool src_cf_bdry_synch);

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
    void interp(int dst_idx,
                boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
                bool dst_cf_bdry_synch,
                int src_idx,
                boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var,
                boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

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
    void interp(int dst_idx,
                boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
                bool dst_cf_bdry_synch,
                int src_idx,
                boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var,
                boost::shared_ptr<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

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
    void laplace(int dst_idx,
                 boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                 const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                 int src1_idx,
                 boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
                 boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
                 double src1_ghost_fill_time,
                 double gamma = 0.0,
                 int src2_idx = -1,
                 boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
                 int dst_depth = 0,
                 int src1_depth = 0,
                 int src2_depth = 0);

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
    void laplace(int dst_idx,
                 boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
                 const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                 int src1_idx,
                 boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src1_var,
                 boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
                 double src1_ghost_fill_time,
                 double gamma = 0.0,
                 int src2_idx = -1,
                 boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL);

    /*!
     * \brief Compute dst = alpha div coef ((grad src1) + (grad src1)^T) + beta
     * src1 + gamma src2, the variable coefficient generalized Laplacian of
     * src1.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void vc_laplace(int dst_idx,
                    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
                    double alpha,
                    double beta,
                    int coef_idx,
                    boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > coef_var,
                    int src1_idx,
                    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src1_var,
                    boost::shared_ptr<HierarchyGhostCellInterpolation> src1_ghost_fill,
                    double src1_ghost_fill_time,
                    double gamma = 0.0,
                    int src2_idx = -1,
                    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                           double alpha,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src1_var,
                           int beta_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > beta_var = NULL,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0,
                           int beta_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
                           double alpha,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src1_var,
                           int beta_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > beta_var = NULL,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0,
                           int beta_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
                           double alpha,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src1_var,
                           int beta_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > beta_var = NULL,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0,
                           int beta_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
                           double alpha,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > dst_var,
                           int alpha_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > alpha_var,
                           int src1_idx,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src1_var,
                           int beta_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > beta_var = NULL,
                           int src2_idx = -1,
                           boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > src2_var = NULL,
                           int dst_depth = 0,
                           int src1_depth = 0,
                           int src2_depth = 0,
                           int alpha_depth = 0,
                           int beta_depth = 0);

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL1Norm(int dst_idx,
                         boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                         int src_idx,
                         boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var);

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL2Norm(int dst_idx,
                         boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                         int src_idx,
                         boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var);

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMaxNorm(int dst_idx,
                          boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > dst_var,
                          int src_idx,
                          boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > src_var);

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL1Norm(int dst_idx,
                         boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
                         int src_idx,
                         boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src_var);

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL2Norm(int dst_idx,
                         boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
                         int src_idx,
                         boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src_var);

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMaxNorm(int dst_idx,
                          boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > dst_var,
                          int src_idx,
                          boost::shared_ptr<SAMRAI::pdat::NodeVariable<double> > src_var);

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
    HierarchyMathOps(const HierarchyMathOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyMathOps& operator=(const HierarchyMathOps& that);

    /*!
     * \brief Reset the coarsen operators.
     */
    void resetCoarsenOperators();

    /*!
     * \brief Reset the refine operators.
     */
    void resetRefineOperators();

    /*!
     * \brief Execute schedule for restricting Outerface data to the specified
     * level from the next finer level.
     */
    void xeqScheduleOuterfaceRestriction(int dst_idx, int src_idx, int dst_ln);

    /*!
     * \brief Execute schedule for restricting Outerside data to the specified
     * level from the next finer level.
     */
    void xeqScheduleOutersideRestriction(int dst_idx, int src_idx, int dst_ln);

    // Housekeeping.
    std::string d_object_name;

    // Patch hierarchy information.
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy> d_hierarchy;
    boost::shared_ptr<SAMRAI::geom::CartesianGridGeometry> d_grid_geom;
    int d_coarsest_ln, d_finest_ln;

    // Scratch Variables.
    boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > d_fc_var;
    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > d_sc_var;
    boost::shared_ptr<SAMRAI::pdat::OuterfaceVariable<double> > d_of_var;
    boost::shared_ptr<SAMRAI::pdat::OutersideVariable<double> > d_os_var;
    int d_fc_idx, d_sc_idx, d_of_idx, d_os_idx;

    // Communications operators, algorithms, and schedules.
    std::string d_coarsen_op_name;
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator> d_of_coarsen_op;
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator> d_os_coarsen_op;
    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> d_of_coarsen_alg;
    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> d_os_coarsen_alg;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> > d_of_coarsen_scheds;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule> > d_os_coarsen_scheds;

    // Hierarchy data operations.
    boost::shared_ptr<SAMRAI::math::HierarchyCellDataOpsReal<double> > d_hier_cc_data_ops;
    boost::shared_ptr<SAMRAI::math::HierarchyFaceDataOpsReal<double> > d_hier_fc_data_ops;
    boost::shared_ptr<SAMRAI::math::HierarchySideDataOpsReal<double> > d_hier_sc_data_ops;

    // Patch math operations.
    PatchMathOps d_patch_math_ops;

    // The cell weights are used to compute norms of data defined on the patch
    // hierarchy.
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > d_wgt_cc_var;
    boost::shared_ptr<SAMRAI::pdat::FaceVariable<double> > d_wgt_fc_var;
    boost::shared_ptr<SAMRAI::pdat::SideVariable<double> > d_wgt_sc_var;
    int d_wgt_cc_idx, d_wgt_fc_idx, d_wgt_sc_idx;
    double d_volume;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HierarchyMathOps
