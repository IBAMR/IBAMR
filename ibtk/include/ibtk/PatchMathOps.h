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
    void curl(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = curl src.
     *
     * Uses centered differences.
     */
    void curl(SAMRAIPointer<SAMRAI::pdat::EdgeDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
             SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
             SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
             SAMRAIPointer<SAMRAI::pdat::EdgeDataNd<double> > src,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst = rot src.
     *
     * Uses centered differences.
     */
    void rot(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
             SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             CartSideRobinPhysBdryOp* bc_op = nullptr,
             double fill_time = 0.0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
             double alpha,
             SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
             double beta,
             SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
             double alpha,
             SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src1,
             double beta,
             SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst_l = alpha div src1 + beta src2_m.
     *
     * Uses centered differences.
     */
    void div(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
             double alpha,
             SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
             double beta,
             SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
             SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
             int l = 0,
             int m = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
              double alpha,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
              double beta,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
              double alpha,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
              double beta,
              SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src2,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
              double alpha,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
              double beta,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > alpha,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
              double beta,
              SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src2,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
              int l = 0) const;

    /*!
     * \brief Computes dst = alpha grad src1_l + beta src2.
     *
     * Uses centered differences.
     */
    void grad(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > alpha,
              SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
              double beta,
              SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
              SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
              int l = 0) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the face-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the side-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the face-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the node-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the cell-centered vector field dst from the edge-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::EdgeDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the node-centered vector field dst from the face-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the side-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial averaging.
     */
    void interp(SAMRAIPointer<SAMRAI::pdat::EdgeDataNd<double> > dst,
                SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                bool dst_ghost_interp) const;

    /*!
     * \brief Computes the side-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes the node-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes the edge-centered vector field dst from the cell-centered
     * vector field src by spatial harmonic averaging.
     */
    void harmonic_interp(SAMRAIPointer<SAMRAI::pdat::EdgeDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                         bool dst_ghost_interp) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                 double alpha,
                 double beta,
                 SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
                 double gamma,
                 SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
                 SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha L src1_m + beta src1_m + gamma src2_n.
     *
     * Uses the standard 5 point stencil in 2D (7 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                 double alpha,
                 double beta,
                 SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
                 double gamma,
                 SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
                 SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                 SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > alpha,
                 double beta,
                 SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
                 double gamma,
                 SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
                 SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = div alpha grad src1_m + beta src1_m + gamma
     * src2_n.
     *
     * Uses a 9 point stencil in 2D (19 point stencil in 3D).
     */
    void laplace(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                 SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > alpha,
                 double beta,
                 SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
                 double gamma,
                 SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
                 SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                 int l = 0,
                 int m = 0,
                 int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                    double alpha,
                    double beta,
                    SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > coef1,
                    SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > coef2,
                    SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
                    double gamma,
                    SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
                    SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Computes dst_l = alpha div coef1 ((grad src1_m) + (grad src1_m)^T)
     * + beta coef2 src1_m + gamma src2_n.
     */
    void vc_laplace(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                    double alpha,
                    double beta,
                    SAMRAIPointer<SAMRAI::pdat::EdgeDataNd<double> > coef1,
                    SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > coef2,
                    SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
                    double gamma,
                    SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
                    SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                    bool use_harmonic_interp,
                    int l = 0,
                    int m = 0,
                    int n = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                           double alpha,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src1,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > beta,
                           SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
                           double alpha,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src1,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > beta,
                           SAMRAIPointer<SAMRAI::pdat::FaceDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                           double alpha,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src1,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > beta,
                           SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0,
                           int m = 0) const;

    /*!
     * \brief Compute dst_i = alpha src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                           double alpha,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
                           double beta,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
                           int i = 0,
                           int j = 0,
                           int k = 0,
                           int l = 0) const;

    /*!
     * \brief Compute dst_i = alpha_l src1_j + beta_m src2_k, pointwise.
     */
    void pointwiseMultiply(SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > dst,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > alpha,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src1,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > beta,
                           SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src2,
                           SAMRAIPointer<SAMRAI::hier::PatchNd> patch,
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
    void pointwiseL1Norm(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                          SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > src,
                          SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Compute dst = |src|_1, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL1Norm(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Compute dst = |src|_2, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseL2Norm(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                         SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src,
                         SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Compute dst = |src|_oo, pointwise.
     *
     * \see resetLevels
     */
    void pointwiseMaxNorm(SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > dst,
                          SAMRAIPointer<SAMRAI::pdat::NodeDataNd<double> > src,
                          SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst1 = strain src (diagonal), and dst2 = strain src (off diagonal).
     *
     * Uses centered differences.
     */
    void strain_rate(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst1,
                     SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst2,
                     SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
                     SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

    /*!
     * \brief Computes dst = strain src.
     *
     * Uses centered differences.
     */
    void strain_rate(SAMRAIPointer<SAMRAI::pdat::CellDataNd<double> > dst,
                     SAMRAIPointer<SAMRAI::pdat::SideDataNd<double> > src,
                     SAMRAIPointer<SAMRAI::hier::PatchNd> patch) const;

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
