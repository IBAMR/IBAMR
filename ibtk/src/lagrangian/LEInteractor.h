// Filename: LEInteractor.h
// Created on 14 Jul 2004 by Boyce Griffith
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

#ifndef included_LEInteractor
#define included_LEInteractor

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LNodeIndexSetData.h>
#include <ibtk/LData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellData.h>
#include <IntVector.h>
#include <Patch.h>
#include <SideData.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

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
     * \brief Function pointer to user-defined delta function kernel along with
     * corresponding stencil size and quadratic constant C.
     */
    static double (*s_delta_fcn)(double r);
    static int s_delta_fcn_stencil_size;
    static double s_delta_fcn_C;

    /*!
     * \brief Sort modes used when interpolating and spreading values.
     *
     * \note Default is: NO_SORT.
     */
    enum SortMode {NO_SORT=0, SORT_INCREASING_LAG_IDX=1, SORT_DECREASING_LAG_IDX=2};
    static SortMode s_sort_mode;

    /*!
     * \brief Floating point precision modes for spreading routines.
     *
     * \note Default is: DOUBLE.
     */
    enum PrecisionMode {DOUBLE=0, DOUBLE_DOUBLE=1};
    static PrecisionMode s_precision_mode;

    /*!
     * \brief Set configuration options from a user-supplied database.
     */
    static void
    setFromDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * \brief Output class configuration.
     */
    static void
    printClassData(
        std::ostream& os);

    /*!
     * \brief Initialize the Timer objects employed by the LEInteractor class.
     *
     * \note It is necessary to initialize the Timer objects prior to using any
     * of the functionality provided by this class.
     */
    static void
    initializeTimers();

    /*!
     * \brief Returns the interpolation/spreading stencil corresponding to the
     * specified weighting function.
     *
     * The return value is -1 for any unknown weighting function type.
     */
    static int
    getStencilSize(
        const std::string& weighting_fcn);

    /*!
     * \brief Returns the constant C associated with a particular IB delta
     * function.
     *
     * The return value is -1 for any unknown weighting function type.
     */
    static double
    getC(
        const std::string& weighting_fcn);

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
    static void
    interpolate(
        SAMRAI::tbox::Pointer<LData>& Q_data,
        const SAMRAI::tbox::Pointer<LData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        SAMRAI::tbox::Pointer<LData>& Q_data,
        const SAMRAI::tbox::Pointer<LData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        std::vector<double>& Q_data,
        const int Q_depth,
        const std::vector<double>& X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        double* const Q_data,
        const int Q_size,
        const int Q_depth,
        const double* const X_data,
        const int X_size,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        std::vector<double>& Q_data,
        const int Q_depth,
        const std::vector<double>& X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const std::string& interp_fcn="IB_4");

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
    static void
    interpolate(
        double* const Q_data,
        const int Q_size,
        const int Q_depth,
        const double* const X_data,
        const int X_size,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& interp_box,
        const std::string& interp_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<LData>& Q_data,
        const SAMRAI::tbox::Pointer<LData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<LData>& Q_data,
        const SAMRAI::tbox::Pointer<LData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const std::vector<double>& Q_data,
        const int Q_depth,
        const std::vector<double>& X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const std::vector<double>& Q_data,
        const int Q_depth,
        const std::vector<double>& X_data,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const double* const Q_data,
        const int Q_size,
        const int Q_depth,
        const double* const X_data,
        const int X_size,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const std::string& spread_fcn="IB_4");

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
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM,double> > q_data,
        const double* const Q_data,
        const int Q_size,
        const int Q_depth,
        const double* const X_data,
        const int X_size,
        const int X_depth,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& spread_box,
        const std::string& spread_fcn="IB_4");

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
    LEInteractor(
        const LEInteractor& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LEInteractor&
    operator=(
        const LEInteractor& that);

    /*!
     * Implementation of the IB interpolation operation.
     */
    static void
    interpolate(
        double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const double* const q_data,
        const SAMRAI::hier::Box<NDIM>& q_data_box,
        const SAMRAI::hier::IntVector<NDIM>& q_gcw,
        const int q_depth,
        const double* const x_lower,
        const double* const x_upper,
        const double* const dx,
        const std::vector<int>& patch_touches_lower_physical_bdry,
        const std::vector<int>& patch_touches_upper_physical_bdry,
        const std::vector<int>& use_alt_one_sided_delta,
        const std::vector<int>& local_indices,
        const std::vector<double>& periodic_offsets,
        const std::string& interp_fcn);

    /*!
     * Implementation of the IB spreading operation.
     */
    static void
    spread(
        double* const q_data,
        const SAMRAI::hier::Box<NDIM>& q_data_box,
        const SAMRAI::hier::IntVector<NDIM>& q_gcw,
        const int q_depth,
        const double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const double* const x_lower,
        const double* const x_upper,
        const double* const dx,
        const std::vector<int>& patch_touches_lower_physical_bdry,
        const std::vector<int>& patch_touches_upper_physical_bdry,
        const std::vector<int>& use_alt_one_sided_delta,
        const std::vector<int>& local_indices,
        const std::vector<double>& periodic_offsets,
        const std::string& spread_fcn);

    /*!
     * \brief Compute the local PETSc indices located within the provided box
     * based on the LNodeIndexSetData values.
     */
    static void
    buildLocalIndices(
        std::vector<int>& local_indices,
        std::vector<double>& periodic_offsets,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const SAMRAI::tbox::Pointer<LNodeIndexSetData>& idx_data);

    /*!
     * \brief Compute the local PETSc indices located within the provided box
     * based on the positions of the Lagrangian mesh nodes.
     */
    static void
    buildLocalIndices(
        std::vector<int>& local_indices,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const double* const X_data,
        const int X_size,
        const int X_depth);

    /*!
     * Implementation of the IB interpolation operation for a user-defined
     * kernel.
     */
    static void
    userDefinedInterpolate(
        double* Q,
        const int Q_depth,
        const double* const X,
        const double* const q,
        const SAMRAI::hier::Box<NDIM>& q_data_box,
        const int* const q_gcw,
        const int q_depth,
        const double* const x_lower,
        const double* const x_upper,
        const double* const dx,
        const int* const local_indices,
        const double* const X_shift,
        const int num_local_indices);

    /*!
     * Implementation of the IB spreading operation for a user-defined kernel.
     */
    static void
    userDefinedSpread(
        double* q,
        const SAMRAI::hier::Box<NDIM>& q_data_box,
        const int* const q_gcw,
        const int q_depth,
        const double* const x_lower,
        const double* const x_upper,
        const double* const dx,
        const double* const Q,
        const int Q_depth,
        const double* const X,
        const int* const local_indices,
        const double* const X_shift,
        const int num_local_indices);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/LEInteractor.I>

//////////////////////////////////////////////////////////////////////////////


#endif //#ifndef included_LEInteractor
