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

#ifndef included_IBTK_PatchMathOps
#define included_IBTK_PatchMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

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
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
             double alpha,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
             double beta,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
             double alpha,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src1,
             double beta,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
             double alpha,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
             double beta,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
              double alpha,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
              double beta,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
              double alpha,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
              double beta,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src2,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
              double alpha,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
              double beta,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > alpha,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
              double beta,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src2,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > alpha,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
              double beta,
              SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
              SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
              int l = 0) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the face-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the side-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the face-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the node-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the edge-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > dst,
                SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                 double alpha,
                 double beta,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
                 double gamma,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
                 SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                 double alpha,
                 double beta,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                 double gamma,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                 SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > alpha,
                 double beta,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
                 double gamma,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
                 SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > alpha,
                 double beta,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
                 double gamma,
                 SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
                 SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                    double alpha,
                    double beta,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > coef1,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > coef2,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                    double gamma,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                    double alpha,
                    double beta,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > coef1,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > coef2,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                    double gamma,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                           double alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
                           double alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                           double alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                           double alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                           double beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > alpha,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > beta,
                           SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
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
    void pointwiseL1Norm(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL1Norm(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                         SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > dst,
                          SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > src,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst1 = strain src (diagonal), and dst2 = strain src (off diagonal).
     *
     * Uses centered differences.
     */
    void strain_rate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst1,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst2,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * \brief Computes dst = strain src.
     *
     * Uses centered differences.
     */
    void strain_rate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > dst,
                     SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
                     SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

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

#endif //#ifndef included_IBTK_PatchMathOps
