// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_LEInteractor
#define included_IBTK_LEInteractor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "Box.h"
#include "IntVector.h"
#include "tbox/Pointer.h"

#include <array>
#include <string>
#include <vector>

namespace IBTK
{
class LData;
template <class T>
class LIndexSetData;
} // namespace IBTK
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
class EdgeData;
template <int DIM, class TYPE>
class NodeData;
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LEInteractor is a utility class that defines several functions
 * to interpolate data from Eulerian grid patches onto Lagrangian meshes and to
 * spread values (\em not densities) from Lagrangian meshes to Eulerian grid
 * patches.
 */
class LEInteractor
{
public:
    /*!
     * \brief Function pointer to user-defined kernel function along with
     * corresponding stencil size and quadratic constant C.
     */
    static double (*s_kernel_fcn)(double r);
    static int s_kernel_fcn_stencil_size;

    /*!
     * \brief Set configuration options from a user-supplied database.
     */
    static void setFromDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Output class configuration.
     */
    static void printClassData(std::ostream& os);

    /*!
     * \brief Returns the interpolation/spreading stencil corresponding to the
     * specified kernel function.
     */
    static int getStencilSize(const std::string& kernel_fcn);

    /*!
     * \brief Return whether or not the provided string corresponds to a known
     * kernel function.
     */
    static bool isKnownKernel(const std::string& kernel_fcn);

    /*!
     * \brief Returns the minimum ghost width size corresponding to the
     * specified kernel function.
     *
     * The minimum ghost cell width is appropriate for simulations in which IB
     * points are allowed to move no more than one cell width between
     * regridding/redistribution operations.  Simulations in which IB points are
     * allowed to move further between regridding/redistribution operations
     * require correspondingly larger ghost cell widths.
     */
    static int getMinimumGhostWidth(const std::string& kernel_fcn);

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                            SAMRAI::tbox::Pointer<LData> X_data,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                            SAMRAI::tbox::Pointer<LData> X_data,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                            SAMRAI::tbox::Pointer<LData> X_data,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(SAMRAI::tbox::Pointer<LData> Q_data,
                            SAMRAI::tbox::Pointer<LData> X_data,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(double* Q_data,
                            int Q_depth,
                            const double* X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(double* Q_data,
                            int Q_depth,
                            const double* X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(double* Q_data,
                            int Q_depth,
                            const double* X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     */
    template <class T>
    static void interpolate(double* Q_data,
                            int Q_depth,
                            const double* X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(std::vector<double>& Q_data,
                            int Q_depth,
                            const std::vector<double>& X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");
    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(std::vector<double>& Q_data,
                            int Q_depth,
                            const std::vector<double>& X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > mask_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(std::vector<double>& Q_data,
                            int Q_depth,
                            const std::vector<double>& X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(std::vector<double>& Q_data,
                            int Q_depth,
                            const std::vector<double>& X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(std::vector<double>& Q_data,
                            int Q_depth,
                            const std::vector<double>& X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > mask_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(std::vector<double>& Q_data,
                            int Q_depth,
                            const std::vector<double>& X_data,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(double* Q_data,
                            int Q_size,
                            int Q_depth,
                            const double* X_data,
                            int X_size,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(double* Q_data,
                            int Q_size,
                            int Q_depth,
                            const double* X_data,
                            int X_size,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > mask_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(double* Q_data,
                            int Q_size,
                            int Q_depth,
                            const double* X_data,
                            int X_size,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(double* Q_data,
                            int Q_size,
                            int Q_depth,
                            const double* X_data,
                            int X_size,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(double* Q_data,
                            int Q_size,
                            int Q_depth,
                            const double* X_data,
                            int X_size,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > mask_data,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation operation.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void interpolate(double* Q_data,
                            int Q_size,
                            int Q_depth,
                            const double* X_data,
                            int X_size,
                            int X_depth,
                            SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                            SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                            const SAMRAI::hier::Box<NDIM>& interp_box,
                            const std::string& interp_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                       SAMRAI::tbox::Pointer<LData> Q_data,
                       SAMRAI::tbox::Pointer<LData> X_data,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                       SAMRAI::tbox::Pointer<LData> Q_data,
                       SAMRAI::tbox::Pointer<LData> X_data,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                       SAMRAI::tbox::Pointer<LData> Q_data,
                       SAMRAI::tbox::Pointer<LData> X_data,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                       SAMRAI::tbox::Pointer<LData> Q_data,
                       SAMRAI::tbox::Pointer<LData> X_data,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the canonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    template <class T>
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > mask_data,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > mask_data,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > mask_data,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > mask_data,
                       SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note X_data must provide the canonical location of the node---i.e.,
     * each node location must lie within the extents of the physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = q(i,j,k) + Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     *
     * \warning This method does \em not support periodic offsets for positions.
     */
    static void spread(SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeData<NDIM, double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                       const SAMRAI::hier::Box<NDIM>& spread_box,
                       const std::string& spread_fcn = "IB_4");

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LEInteractor() = delete;

    /*!
     * \brief Default destructor constructor.
     *
     * \note This destructor is not implemented and should not be used.
     */
    ~LEInteractor() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LEInteractor(const LEInteractor& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LEInteractor& operator=(const LEInteractor& that) = delete;

    /*!
     * Implementation of the IB interpolation operation.
     */
    static void interpolate(double* Q_data,
                            int Q_depth,
                            const double* X_data,
                            const double* q_data,
                            const SAMRAI::hier::Box<NDIM>& q_data_box,
                            const SAMRAI::hier::IntVector<NDIM>& q_gcw,
                            int q_depth,
                            const double* x_lower,
                            const double* x_upper,
                            const double* dx,
                            const std::array<int, NDIM>& patch_touches_lower_physical_bdry,
                            const std::array<int, NDIM>& patch_touches_upper_physical_bdry,
                            const std::vector<int>& local_indices,
                            const std::vector<double>& periodic_shifts,
                            const std::string& interp_fcn,
                            int axis = 0);

    /*!
     * Implementation of the IB spreading operation.
     */
    static void spread(double* q_data,
                       const SAMRAI::hier::Box<NDIM>& q_data_box,
                       const SAMRAI::hier::IntVector<NDIM>& q_gcw,
                       int q_depth,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       const double* x_lower,
                       const double* x_upper,
                       const double* dx,
                       const std::array<int, NDIM>& patch_touches_lower_physical_bdry,
                       const std::array<int, NDIM>& patch_touches_upper_physical_bdry,
                       const std::vector<int>& local_indices,
                       const std::vector<double>& periodic_shifts,
                       const std::string& spread_fcn,
                       int axis = 0);

    /*!
     * \brief Compute the local PETSc indices located within the provided box
     * based on the LNodeIndexSetData values.
     */
    template <class T>
    static void buildLocalIndices(std::vector<int>& local_indices,
                                  std::vector<double>& periodic_shifts,
                                  const SAMRAI::hier::Box<NDIM>& box,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                  const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
                                  SAMRAI::tbox::Pointer<LIndexSetData<T> > idx_data);

    /*!
     * \brief Compute the local PETSc indices located within the provided box
     * based on the positions of the Lagrangian mesh nodes.
     */
    static void buildLocalIndices(std::vector<int>& local_indices,
                                  const SAMRAI::hier::Box<NDIM>& box,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                                  const double* X_data,
                                  int X_size,
                                  int X_depth);

    /*!
     * Implementation of the IB interpolation operation for a user-defined
     * kernel.
     */
    static void userDefinedInterpolate(double* Q,
                                       int Q_depth,
                                       const double* X,
                                       const double* q,
                                       const SAMRAI::hier::Box<NDIM>& q_data_box,
                                       const int* q_gcw,
                                       int q_depth,
                                       const double* x_lower,
                                       const double* x_upper,
                                       const double* dx,
                                       const int* local_indices,
                                       const double* X_shift,
                                       int num_local_indices);

    /*!
     * Implementation of the IB spreading operation for a user-defined kernel.
     */
    static void userDefinedSpread(double* q,
                                  const SAMRAI::hier::Box<NDIM>& q_data_box,
                                  const int* q_gcw,
                                  int q_depth,
                                  const double* x_lower,
                                  const double* x_upper,
                                  const double* dx,
                                  const double* Q,
                                  int Q_depth,
                                  const double* X,
                                  const int* local_indices,
                                  const double* X_shift,
                                  int num_local_indices);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LEInteractor
