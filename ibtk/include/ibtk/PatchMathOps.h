// Filename: PatchMathOps.h
// Created on 23 Jul 2002 by Boyce Griffith
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

#ifndef included_PatchMathOps
#define included_PatchMathOps

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"

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
    PatchMathOps();

    /*!
     * \brief Destructor.
     */
    ~PatchMathOps();

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
             CartSideRobinPhysBdryOp* bc_op = NULL,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = NULL,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = NULL,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
             SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src,
             SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
             CartSideRobinPhysBdryOp* bc_op = NULL,
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
     * \brief Computes dst_l = alpha div coef ((grad src1_m) + (grad src1_m)^T)
     * + beta src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > dst,
                    double alpha,
                    double beta,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > coef,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src1,
                    double gamma,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > src2,
                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
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

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PatchMathOps(const PatchMathOps& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PatchMathOps& operator=(const PatchMathOps& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PatchMathOps
