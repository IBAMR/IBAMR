// Filename: LEInteractor.h
// Created on 14 Jul 2004 by Boyce Griffith
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

#ifndef included_LEInteractor
#define included_LEInteractor

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <iosfwd>
#include <string>
#include <vector>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"

namespace boost
{
template <class T, std::size_t N>
class array;
} // namespace boost

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

class Patch;
} // namespace hier
namespace pdat
{
template <class TYPE>
class CellData;
template <class TYPE>
class EdgeData;
template <class TYPE>
class NodeData;
template <class TYPE>
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
    static void setFromDatabase(boost::shared_ptr<SAMRAI::tbox::Database> db);

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
    static void interpolate(boost::shared_ptr<LData> Q_data,
                            boost::shared_ptr<LData> X_data,
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
    static void interpolate(boost::shared_ptr<LData> Q_data,
                            boost::shared_ptr<LData> X_data,
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
    static void interpolate(boost::shared_ptr<LData> Q_data,
                            boost::shared_ptr<LData> X_data,
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
    static void interpolate(boost::shared_ptr<LData> Q_data,
                            boost::shared_ptr<LData> X_data,
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
                            boost::shared_ptr<LIndexSetData<T> > idx_data,
                            boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
                            const SAMRAI::hier::IntVector& periodic_shift,
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
                            boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
                            boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                            boost::shared_ptr<SAMRAI::hier::Patch> patch,
                            const SAMRAI::hier::Box& interp_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                       boost::shared_ptr<LData> Q_data,
                       boost::shared_ptr<LData> X_data,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                       boost::shared_ptr<LData> Q_data,
                       boost::shared_ptr<LData> X_data,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                       boost::shared_ptr<LData> Q_data,
                       boost::shared_ptr<LData> X_data,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                       boost::shared_ptr<LData> Q_data,
                       boost::shared_ptr<LData> X_data,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       int X_depth,
                       boost::shared_ptr<LIndexSetData<T> > idx_data,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const SAMRAI::hier::IntVector& periodic_shift,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                       const std::vector<double>& Q_data,
                       int Q_depth,
                       const std::vector<double>& X_data,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::CellData<double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::NodeData<double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::SideData<double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
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
    static void spread(boost::shared_ptr<SAMRAI::pdat::EdgeData<double> > q_data,
                       const double* Q_data,
                       int Q_size,
                       int Q_depth,
                       const double* X_data,
                       int X_size,
                       int X_depth,
                       boost::shared_ptr<SAMRAI::hier::Patch> patch,
                       const SAMRAI::hier::Box& spread_box,
                       const std::string& spread_fcn = "IB_4");

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    LEInteractor();

    /*!
     * \brief Default destructor constructor.
     *
     * \note This destructor is not implemented and should not be used.
     */
    ~LEInteractor();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LEInteractor(const LEInteractor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LEInteractor& operator=(const LEInteractor& that);

    /*!
     * Implementation of the IB interpolation operation.
     */
    static void interpolate(double* Q_data,
                            int Q_depth,
                            const double* X_data,
                            const double* q_data,
                            const SAMRAI::hier::Box& q_data_box,
                            const SAMRAI::hier::IntVector& q_gcw,
                            int q_depth,
                            const double* x_lower,
                            const double* x_upper,
                            const double* dx,
                            const boost::array<int, NDIM>& patch_touches_lower_physical_bdry,
                            const boost::array<int, NDIM>& patch_touches_upper_physical_bdry,
                            const std::vector<int>& local_indices,
                            const std::vector<double>& periodic_shifts,
                            const std::string& interp_fcn,
                            int axis = 0);

    /*!
     * Implementation of the IB spreading operation.
     */
    static void spread(double* q_data,
                       const SAMRAI::hier::Box& q_data_box,
                       const SAMRAI::hier::IntVector& q_gcw,
                       int q_depth,
                       const double* Q_data,
                       int Q_depth,
                       const double* X_data,
                       const double* x_lower,
                       const double* x_upper,
                       const double* dx,
                       const boost::array<int, NDIM>& patch_touches_lower_physical_bdry,
                       const boost::array<int, NDIM>& patch_touches_upper_physical_bdry,
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
                                  const SAMRAI::hier::Box& box,
                                  boost::shared_ptr<SAMRAI::hier::Patch> patch,
                                  const SAMRAI::hier::IntVector& periodic_shift,
                                  boost::shared_ptr<LIndexSetData<T> > idx_data);

    /*!
     * \brief Compute the local PETSc indices located within the provided box
     * based on the positions of the Lagrangian mesh nodes.
     */
    static void buildLocalIndices(std::vector<int>& local_indices,
                                  const SAMRAI::hier::Box& box,
                                  boost::shared_ptr<SAMRAI::hier::Patch> patch,
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
                                       const SAMRAI::hier::Box& q_data_box,
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
                                  const SAMRAI::hier::Box& q_data_box,
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

#endif //#ifndef included_LEInteractor
