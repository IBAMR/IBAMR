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

#ifndef included_IBTK_HierarchyMathOps
#define included_IBTK_HierarchyMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/PatchMathOps.h"
#include "ibtk/SAMRAIDataCache.h"
#include "ibtk/ibtk_enums.h"
#include "ibtk/ibtk_utilities.h"

#include "CartesianGridGeometry.h"
#include "CellVariable.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "FaceVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "NodeVariable.h"
#include "OuterfaceVariable.h"
#include "OutersideVariable.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SideVariable.h"
#include "VariableContext.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace pdat
{
template <int DIM, class TYPE>
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
template <int DIM>
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
class HierarchyMathOps : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     *
     * Does nothing interesting.
     */
    HierarchyMathOps(std::string name,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                     int coarsest_ln = -1,
                     int finest_ln = -1,
                     std::string coarsen_op_name = "CONSERVATIVE_COARSEN");

    /*!
     * \brief Destructor.
     */
    virtual ~HierarchyMathOps() = default;

    /*!
     * \name Methods to set the hierarchy and range of levels.
     */
    //\{

    /*!
     * \brief Reset the patch hierarchy over which operations occur.
     */
    void setPatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

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
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > getCellWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::CellVariable that is used to store cell weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a cell is not refined in the hierarchy, its weight is set to its
     * volume.  If a cell is refined, its weight is set to zero.
     */
    int getCellWeightPatchDescriptorIndex();

    /*!
     * \brief Access the SAMRAI::pdat::FaceVariable that is used to store face
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a face is not refined in the hierarchy, its weight is set to its
     * volume.  If a face is refined, its weight is set to zero.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > getFaceWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::FaceVariable that is used to store face weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a face is not refined in the hierarchy, its weight is set to its
     * volume.  If a face is refined, its weight is set to zero.
     */
    int getFaceWeightPatchDescriptorIndex();

    /*!
     * \brief Access the SAMRAI::pdat::SideVariable that is used to store side
     * weights.  These weights are used in computing discrete vector norms.
     *
     * If a side is not refined in the hierarchy, its weight is set to its
     * volume.  If a side is refined, its weight is set to zero.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > getSideWeightVariable() const;

    /*!
     * \brief Access the patch descriptor index that corresponds to the
     * SAMRAI::pdat::SideVariable that is used to store side weights.  These
     * weights are used in computing discrete vector norms.
     *
     * If a side is not refined in the hierarchy, its weight is set to its
     * volume.  If a side is refined, its weight is set to zero.
     */
    int getSideWeightPatchDescriptorIndex();

    /*!
     * \brief Returns the volume of the physical domain.
     */
    double getVolumeOfPhysicalDomain() const;

    //\}

    /*!
     * \brief Get the patch hierarchy object associated with the class.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > getPatchHierarchy() const;

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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > dst_var,
              int src_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
             int src_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time,
             const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs =
                 std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>());

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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
             int src_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time,
             const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs =
                 std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>());

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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
             int src_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > src_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time,
             const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs =
                 std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>());

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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
             int src_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
             double src_ghost_fill_time,
             const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs =
                 std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>());

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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
             double alpha,
             int src1_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
             double src1_ghost_fill_time,
             double beta = 0.0,
             int src2_idx = -1,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
             double alpha,
             int src1_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src1_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
             double src1_ghost_fill_time,
             bool src1_cf_bdry_synch,
             double beta = 0.0,
             int src2_idx = -1,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
             double alpha,
             int src1_idx,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
             SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
             double src1_ghost_fill_time,
             bool src1_cf_bdry_synch,
             double beta = 0.0,
             int src2_idx = -1,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
              double alpha,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
              bool dst_cf_bdry_synch,
              double alpha,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
              bool dst_cf_bdry_synch,
              double alpha,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
              int alpha_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > alpha_var,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
              int alpha_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > alpha_var,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
              bool dst_cf_bdry_synch,
              int alpha_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > alpha_var,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src2_var = NULL,
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
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
              bool dst_cf_bdry_synch,
              int alpha_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > alpha_var,
              int src1_idx,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
              SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
              double src1_ghost_fill_time,
              double beta = 0.0,
              int src2_idx = -1,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL,
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
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
                bool dst_cf_bdry_synch,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                bool dst_cf_bdry_synch,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

    /*!
     * \brief Interpolate to a cell-centered vector/tensor field from a
     * node-centered vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void interp(int dst_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

    /*!
     * \brief Interpolate to a cell-centered vector/tensor field from an
     * edge-centered vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void interp(int dst_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

    /*!
     * \brief Interpolate to a node-centered vector/tensor field from a
     * cell-centered vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When specified, the ghost cells
     * of the node centered variable are computed as averages of the cell centered
     * variable
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void interp(int dst_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                bool dst_ghost_interp,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

    /*!
     * \brief Interpolate to a edge-centered vector/tensor field from a
     * cell-centered vector/tensor field.
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When specified, the ghost cells
     * of the edge centered variable are computed as averages of the cell centered
     * variable
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void interp(int dst_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > dst_var,
                bool dst_ghost_interp,
                int src_idx,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                double src_ghost_fill_time);

    /*!
     * \brief Harmonically interpolate to a side-centered normal vector/tensor field from a
     * cell-centered vector/tensor field. This should be used when the gradients in the src
     * quantity are large
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) harmonic averaging.  When the interpolation occurs
     * over multiple levels of the hierarchy, second order interpolation is used
     * at the coarse-fine interface to correct fine values along the interface.
     * When specified, coarse values on each coarse-fine interface are
     * synchronized after performing the interpolation.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void harmonic_interp(int dst_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                         bool dst_cf_bdry_synch,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                         SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                         double src_ghost_fill_time);

    /*!
     * \brief Harmonically interpolate to a node-centered vector/tensor field from a
     * cell-centered vector/tensor field. This should be used when the gradients in the src
     * quantity are large
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When specified, the ghost cells
     * of the node centered variable are computed as averages of the cell centered
     * variable
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void harmonic_interp(int dst_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                         bool dst_ghost_interp,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                         SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                         double src_ghost_fill_time);

    /*!
     * \brief Harmonically interpolate to a edge-centered vector/tensor field from a
     * cell-centered vector/tensor field. This should be used when the gradients in the src
     * quantity are large
     *
     * Interpolate a vector or tensor field from one variable type to another
     * using (second-order accurate) averaging.  When specified, the ghost cells
     * of the edge centered variable are computed as averages of the cell centered
     * variable
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void harmonic_interp(int dst_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > dst_var,
                         bool dst_ghost_interp,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var,
                         SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
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
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                 const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                 int src1_idx,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
                 SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                 double src1_ghost_fill_time,
                 double gamma = 0.0,
                 int src2_idx = -1,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                 const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                 int src1_idx,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
                 SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                 double src1_ghost_fill_time,
                 double gamma = 0.0,
                 int src2_idx = -1,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL);

    /*!
     * \brief Compute dst = alpha div coef1 ((grad src1) + (grad src1)^T) + beta coef2
     * src1 + gamma src2, the variable coefficient generalized Laplacian of
     * src1.
     *
     * \note This routine is only appropriate for NDIM == 2 and uses node centered coef1.
     * Interally, coef1 will be interpolated onto cell centers when needed. Harmonic interpolation
     * is used by default, but coef1_interp_type = VC_AVERAGE_INTERP will apply arithmetic averaging
     * from nodes to cells
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void vc_laplace(int dst_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                    double alpha,
                    double beta,
                    int coef1_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > coef1_var,
                    int src1_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
                    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                    double src1_ghost_fill_time,
                    IBTK::VCInterpType coef1_interp_type = VC_HARMONIC_INTERP,
                    int coef2_idx = -1,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > coef2_var = NULL,
                    double gamma = 0.0,
                    int src2_idx = -1,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL);

    /*!
     * \brief Compute dst = alpha div coef1 ((grad src1) + (grad src1)^T) + beta coef2
     * src1 + gamma src2, the variable coefficient generalized Laplacian of
     * src1.
     *
     * \note This routine is only appropriate for NDIM == 3 and uses edge centered coef1.
     * Interally, coef1 will be interpolated onto cell centers when needed. Harmonic interpolation
     * is used by default, but coef1_interp_type = VC_AVERAGE_INTERP will apply arithmetic averaging
     * from edges to cells
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void vc_laplace(int dst_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                    double alpha,
                    double beta,
                    int coef1_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > coef1_var,
                    int src1_idx,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
                    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src1_ghost_fill,
                    double src1_ghost_fill_time,
                    IBTK::VCInterpType coef1_interp_type = VC_HARMONIC_INTERP,
                    int coef2_idx = -1,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > coef2_var = NULL,
                    double gamma = 0.0,
                    int src2_idx = -1,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL);

    /*!
     * \brief Compute dst = alpha src1 + beta src2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMultiply(int dst_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                           double alpha,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src1_var,
                           int beta_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > beta_var = NULL,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
                           double alpha,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src1_var,
                           int beta_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > beta_var = NULL,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                           double alpha,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src1_var,
                           int beta_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > beta_var = NULL,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                           double alpha,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
                           double beta = 0.0,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL,
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
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > dst_var,
                           int alpha_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > alpha_var,
                           int src1_idx,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src1_var,
                           int beta_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > beta_var = NULL,
                           int src2_idx = -1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src2_var = NULL,
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
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var);

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL2Norm(int dst_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var);

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMaxNorm(int dst_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                          int src_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > src_var);

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL1Norm(int dst_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src_var);

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseL2Norm(int dst_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                         int src_idx,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src_var);

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void pointwiseMaxNorm(int dst_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > dst_var,
                          int src_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > src_var);

    /*!
     * \brief Compute the cell-centered strain rate of a side-centered vector field
     * using centered differences.
     *
     * Sets dst1 = strain src (diagonal components).
     * Sets dst2 = strain src (off diagonal components).
     *
     * Compute the strain of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void strain_rate(int dst1_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst1_var,
                     int dst2_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst2_var,
                     int src_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
                     SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                     double src_ghost_fill_time);

    /*!
     * \brief Compute the cell-centered strain rate of a side-centered vector field
     * using centered differences.

     * Sets dst = strain src.
     *
     * Compute the strain of a vector field using centered differences.
     *
     * \see setPatchHierarchy
     * \see resetLevels
     */
    void strain_rate(int dst_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > dst_var,
                     int src_idx,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > src_var,
                     SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> src_ghost_fill,
                     double src_ghost_fill_time);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyMathOps() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyMathOps(const HierarchyMathOps& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyMathOps& operator=(const HierarchyMathOps& that) = delete;

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

    /*!
     * \brief Reset cell-centered weights, allocating patch data if needed.
     */
    void resetCellWeights(int coarsest_ln, int finest_ln);

    /*!
     * \brief Reset face-centered weights, allocating patch data if needed.
     */
    void resetFaceWeights(int coarsest_ln, int finest_ln);

    /*!
     * \brief Reset side-centered weights, allocating patch data if needed.
     */
    void resetSideWeights(int coarsest_ln, int finest_ln);

    // Housekeeping.
    std::string d_object_name;

    // Patch hierarchy information.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;
    int d_coarsest_ln, d_finest_ln;

    // Scratch Variables.
    SAMRAIDataCache d_cached_eulerian_data;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_fc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::OuterfaceVariable<NDIM, double> > d_of_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::OutersideVariable<NDIM, double> > d_os_var;
    int d_fc_idx = IBTK::invalid_index, d_sc_idx = IBTK::invalid_index, d_of_idx = IBTK::invalid_index,
        d_os_idx = IBTK::invalid_index;

    // Communications operators, algorithms, and schedules.
    std::string d_coarsen_op_name;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_of_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_os_coarsen_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_of_coarsen_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_os_coarsen_alg;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_of_coarsen_scheds;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_os_coarsen_scheds;

    // Hierarchy data operations.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;

    // Patch math operations.
    PatchMathOps d_patch_math_ops;

    // The cell weights are used to compute norms of data defined on the patch
    // hierarchy.
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_wgt_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_wgt_fc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_wgt_sc_var;
    int d_wgt_cc_idx = IBTK::invalid_index, d_wgt_fc_idx = IBTK::invalid_index, d_wgt_sc_idx = IBTK::invalid_index;
    bool d_using_wgt_cc = false, d_using_wgt_fc = false, d_using_wgt_sc = false;
    double d_volume = 0.0;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_HierarchyMathOps
