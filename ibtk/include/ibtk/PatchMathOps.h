// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBTK_PatchMathOps
#define included_IBTK_PatchMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAICellData.h"
#include "SAMRAIEdgeData.h"
#include "SAMRAIFaceData.h"
#include "SAMRAINodeData.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideData.h"
#include "tbox/DescribedClass.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class CellData;
template <int DIM, class TYPE>
class FaceData;
template <int DIM, class TYPE>
class NodeData;
template <int DIM, class TYPE>
class EdgeData;
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PatchMathOps provides functionality to perform mathematical
 * operations on \em individual SAMRAI::hier::Patch objects.
 *
 * \note Coarse-fine interface discretizations are handled in an implicit manner
 * via ghost cells.
 */
class PatchMathOps : public SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Default constructor.
     */
    PatchMathOps() = default;

    /*!
     * \brief Destructor.
     */
    virtual ~PatchMathOps() = default;

    /*!
     * \name Mathematical operations.
     */
    //\{

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAICellData<double>> dst,
              SAMRAIPointer<SAMRAICellData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAICellData<double>> dst,
              SAMRAIPointer<SAMRAIFaceData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAIFaceData<double>> dst,
              SAMRAIPointer<SAMRAIFaceData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAICellData<double>> dst,
              SAMRAIPointer<SAMRAISideData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAISideData<double>> dst,
              SAMRAIPointer<SAMRAISideData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAINodeData<double>> dst,
              SAMRAIPointer<SAMRAISideData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAIEdgeData<double>> dst,
              SAMRAIPointer<SAMRAISideData<double>> src,
              SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAISideData<double>> dst,
             SAMRAIPointer<SAMRAINodeData<double>> src,
             SAMRAIPointer<SAMRAIPatch> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAISideData<double>> dst,
             SAMRAIPointer<SAMRAICellData<double>> src,
             SAMRAIPointer<SAMRAIPatch> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAISideData<double>> dst,
             SAMRAIPointer<SAMRAIEdgeData<double>> src,
             SAMRAIPointer<SAMRAIPatch> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAISideData<double>> dst,
             SAMRAIPointer<SAMRAISideData<double>> src,
             SAMRAIPointer<SAMRAIPatch> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAIPointer<SAMRAICellData<double>> dst,
             double alpha,
             SAMRAIPointer<SAMRAICellData<double>> src1,
             double beta,
             SAMRAIPointer<SAMRAICellData<double>> src2,
             SAMRAIPointer<SAMRAIPatch> patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAIPointer<SAMRAICellData<double>> dst,
             double alpha,
             SAMRAIPointer<SAMRAIFaceData<double>> src1,
             double beta,
             SAMRAIPointer<SAMRAICellData<double>> src2,
             SAMRAIPointer<SAMRAIPatch> patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAIPointer<SAMRAICellData<double>> dst,
             double alpha,
             SAMRAIPointer<SAMRAISideData<double>> src1,
             double beta,
             SAMRAIPointer<SAMRAICellData<double>> src2,
             SAMRAIPointer<SAMRAIPatch> patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAICellData<double>> dst,
              double alpha,
              SAMRAIPointer<SAMRAICellData<double>> src1,
              double beta,
              SAMRAIPointer<SAMRAICellData<double>> src2,
              SAMRAIPointer<SAMRAIPatch> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAIFaceData<double>> dst,
              double alpha,
              SAMRAIPointer<SAMRAICellData<double>> src1,
              double beta,
              SAMRAIPointer<SAMRAIFaceData<double>> src2,
              SAMRAIPointer<SAMRAIPatch> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAISideData<double>> dst,
              double alpha,
              SAMRAIPointer<SAMRAICellData<double>> src1,
              double beta,
              SAMRAIPointer<SAMRAISideData<double>> src2,
              SAMRAIPointer<SAMRAIPatch> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAIFaceData<double>> dst,
              SAMRAIPointer<SAMRAIFaceData<double>> alpha,
              SAMRAIPointer<SAMRAICellData<double>> src1,
              double beta,
              SAMRAIPointer<SAMRAIFaceData<double>> src2,
              SAMRAIPointer<SAMRAIPatch> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAISideData<double>> dst,
              SAMRAIPointer<SAMRAISideData<double>> alpha,
              SAMRAIPointer<SAMRAICellData<double>> src1,
              double beta,
              SAMRAIPointer<SAMRAISideData<double>> src2,
              SAMRAIPointer<SAMRAIPatch> patch,
              int l = 0) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the face-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAICellData<double>> dst,
                SAMRAIPointer<SAMRAIFaceData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the side-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAICellData<double>> dst,
                SAMRAIPointer<SAMRAISideData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the face-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAIFaceData<double>> dst,
                SAMRAIPointer<SAMRAICellData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAISideData<double>> dst,
                SAMRAIPointer<SAMRAICellData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the node-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAICellData<double>> dst,
                SAMRAIPointer<SAMRAINodeData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the edge-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAICellData<double>> dst,
                SAMRAIPointer<SAMRAIEdgeData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAINodeData<double>> dst,
                SAMRAIPointer<SAMRAICellData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the node-centered vector field dst from the face-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAINodeData<double>> dst,
                SAMRAIPointer<SAMRAIFaceData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the side-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAINodeData<double>> dst,
                SAMRAIPointer<SAMRAISideData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAIEdgeData<double>> dst,
                SAMRAIPointer<SAMRAICellData<double>> src,
                SAMRAIPointer<SAMRAIPatch> patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAIPointer<SAMRAISideData<double>> dst,
                         SAMRAIPointer<SAMRAICellData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAIPointer<SAMRAINodeData<double>> dst,
                         SAMRAIPointer<SAMRAICellData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAIPointer<SAMRAIEdgeData<double>> dst,
                         SAMRAIPointer<SAMRAICellData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAICellData<double>> dst,
                 double alpha,
                 double beta,
                 SAMRAIPointer<SAMRAICellData<double>> src1,
                 double gamma,
                 SAMRAIPointer<SAMRAICellData<double>> src2,
                 SAMRAIPointer<SAMRAIPatch> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAISideData<double>> dst,
                 double alpha,
                 double beta,
                 SAMRAIPointer<SAMRAISideData<double>> src1,
                 double gamma,
                 SAMRAIPointer<SAMRAISideData<double>> src2,
                 SAMRAIPointer<SAMRAIPatch> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAICellData<double>> dst,
                 SAMRAIPointer<SAMRAIFaceData<double>> alpha,
                 double beta,
                 SAMRAIPointer<SAMRAICellData<double>> src1,
                 double gamma,
                 SAMRAIPointer<SAMRAICellData<double>> src2,
                 SAMRAIPointer<SAMRAIPatch> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAICellData<double>> dst,
                 SAMRAIPointer<SAMRAISideData<double>> alpha,
                 double beta,
                 SAMRAIPointer<SAMRAICellData<double>> src1,
                 double gamma,
                 SAMRAIPointer<SAMRAICellData<double>> src2,
                 SAMRAIPointer<SAMRAIPatch> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAIPointer<SAMRAISideData<double>> dst,
                    double alpha,
                    double beta,
                    SAMRAIPointer<SAMRAINodeData<double>> coef1,
                    SAMRAIPointer<SAMRAISideData<double>> coef2,
                    SAMRAIPointer<SAMRAISideData<double>> src1,
                    double gamma,
                    SAMRAIPointer<SAMRAISideData<double>> src2,
                    SAMRAIPointer<SAMRAIPatch> patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAIPointer<SAMRAISideData<double>> dst,
                    double alpha,
                    double beta,
                    SAMRAIPointer<SAMRAIEdgeData<double>> coef1,
                    SAMRAIPointer<SAMRAISideData<double>> coef2,
                    SAMRAIPointer<SAMRAISideData<double>> src1,
                    double gamma,
                    SAMRAIPointer<SAMRAISideData<double>> src2,
                    SAMRAIPointer<SAMRAIPatch> patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAICellData<double>> dst,
                           double alpha,
                           SAMRAIPointer<SAMRAICellData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAICellData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAICellData<double>> dst,
                           SAMRAIPointer<SAMRAICellData<double>> alpha,
                           SAMRAIPointer<SAMRAICellData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAICellData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAICellData<double>> dst,
                           SAMRAIPointer<SAMRAICellData<double>> alpha,
                           SAMRAIPointer<SAMRAICellData<double>> src1,
                           SAMRAIPointer<SAMRAICellData<double>> beta,
                           SAMRAIPointer<SAMRAICellData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAIFaceData<double>> dst,
                           double alpha,
                           SAMRAIPointer<SAMRAIFaceData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAIFaceData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAIFaceData<double>> dst,
                           SAMRAIPointer<SAMRAIFaceData<double>> alpha,
                           SAMRAIPointer<SAMRAIFaceData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAIFaceData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAIFaceData<double>> dst,
                           SAMRAIPointer<SAMRAIFaceData<double>> alpha,
                           SAMRAIPointer<SAMRAIFaceData<double>> src1,
                           SAMRAIPointer<SAMRAIFaceData<double>> beta,
                           SAMRAIPointer<SAMRAIFaceData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAINodeData<double>> dst,
                           double alpha,
                           SAMRAIPointer<SAMRAINodeData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAINodeData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAINodeData<double>> dst,
                           SAMRAIPointer<SAMRAINodeData<double>> alpha,
                           SAMRAIPointer<SAMRAINodeData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAINodeData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAINodeData<double>> dst,
                           SAMRAIPointer<SAMRAINodeData<double>> alpha,
                           SAMRAIPointer<SAMRAINodeData<double>> src1,
                           SAMRAIPointer<SAMRAINodeData<double>> beta,
                           SAMRAIPointer<SAMRAINodeData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAISideData<double>> dst,
                           double alpha,
                           SAMRAIPointer<SAMRAISideData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAISideData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAISideData<double>> dst,
                           SAMRAIPointer<SAMRAISideData<double>> alpha,
                           SAMRAIPointer<SAMRAISideData<double>> src1,
                           double beta,
                           SAMRAIPointer<SAMRAISideData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAISideData<double>> dst,
                           SAMRAIPointer<SAMRAISideData<double>> alpha,
                           SAMRAIPointer<SAMRAISideData<double>> src1,
                           SAMRAIPointer<SAMRAISideData<double>> beta,
                           SAMRAIPointer<SAMRAISideData<double>> src2,
                           SAMRAIPointer<SAMRAIPatch> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL1Norm(SAMRAIPointer<SAMRAICellData<double>> dst,
                         SAMRAIPointer<SAMRAICellData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(SAMRAIPointer<SAMRAICellData<double>> dst,
                         SAMRAIPointer<SAMRAICellData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(SAMRAIPointer<SAMRAICellData<double>> dst,
                          SAMRAIPointer<SAMRAICellData<double>> src,
                          SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL1Norm(SAMRAIPointer<SAMRAINodeData<double>> dst,
                         SAMRAIPointer<SAMRAINodeData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(SAMRAIPointer<SAMRAINodeData<double>> dst,
                         SAMRAIPointer<SAMRAINodeData<double>> src,
                         SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(SAMRAIPointer<SAMRAINodeData<double>> dst,
                          SAMRAIPointer<SAMRAINodeData<double>> src,
                          SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst1 = strain src (diagonal), and dst2 = strain src (off diagonal).
     *
     * Uses centered differences.
     */
    void strain_rate(SAMRAIPointer<SAMRAICellData<double>> dst1,
                     SAMRAIPointer<SAMRAICellData<double>> dst2,
                     SAMRAIPointer<SAMRAISideData<double>> src,
                     SAMRAIPointer<SAMRAIPatch> patch) const;

    /*!
     * \brief Computes dst = strain src.
     *
     * Uses centered differences.
     */
    void strain_rate(SAMRAIPointer<SAMRAICellData<double>> dst,
                     SAMRAIPointer<SAMRAISideData<double>> src,
                     SAMRAIPointer<SAMRAIPatch> patch) const;

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PatchMathOps(const PatchMathOps& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PatchMathOps& operator=(const PatchMathOps& that) = delete;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PatchMathOps
