#ifndef included_LEInteractor
#define included_LEInteractor

// Filename: LEInteractor.h
// Last modified: <04.Jul.2007 13:36:40 boyce@bigboy.nyconnect.com>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/LNodeIndexData2.h>
#include <ibamr/LNodeLevelData.h>

// SAMRAI INCLUDES
#include <Box.h>
#include <CellData.h>
#include <IntVector.h>
#include <Patch.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class LEInteractor is a utility class that defines several functions
 * to interpolate data from Eulerian grid patches onto Lagrangian meshes and to
 * spread values (\em not densities!) from Lagrangian meshes to Eulerian grid
 * patches.
 */
class LEInteractor
{
public:
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
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the cannonical location of the
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
        SAMRAI::tbox::Pointer<LNodeLevelData>& Q_data,
        const SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& interp_fcn="IB_4");

    /*!
     * \brief Interpolate data from an Eulerian grid to a Lagrangian mesh.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the cannonical location of the
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
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& interp_fcn="IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the cannonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading operation, the
     * implemented operations spreads values, NOT densities.
     */
    static void
    spread(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > q_data,
        const SAMRAI::tbox::Pointer<LNodeLevelData>& Q_data,
        const SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& spread_fcn="IB_4");

    /*!
     * \brief Spread data from a Lagrangian mesh to an Eulerian grid.  The
     * positions of the nodes of the Lagrangian mesh are specified by X_data.
     *
     * \note This method employs periodic boundary conditions where appropriate
     * and when requested.  X_data must provide the cannonical location of the
     * node---i.e., each node location must lie within the extents of the
     * physical domain.
     *
     * \note The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).  The
     * spreading formula is
     *
     *     q(i,j,k) = Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
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
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const std::string& spread_fcn="IB_4");

    /*!
     * \brief A version of the spreading functionality that has been modified to
     * spread forces reflected across physical boundaries with reflection
     * boundary conditions.
     *
     * \warning This routine \em assumes that all physical boundaries are
     * reflection boundaries.
     *
     * \warning Presently, this function \em does \em not work correctly when
     * curvilinear mesh nodes are in the vicinity of edges or corners in the
     * physical domain.  NO ERROR IS ISSUED IF SUCH A CASE IS ATTEMPTED!
     */
    static void
    spreadReflectedForces(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > f_data,
        const SAMRAI::tbox::Pointer<LNodeLevelData>& F_data,
        const SAMRAI::tbox::Pointer<LNodeLevelData>& X_data,
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
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
     * \brief Compute the local PETSc indices located within the provided box.
     */
    static void
    getLocalIndices(
        std::vector<int>& local_indices,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data);

    /*!
     * \brief Compute the periodic shifts for the indicies located withing the
     * provided box.
     */
    static void
    getPeriodicOffsets(
        std::vector<double>& local_indices,
        const SAMRAI::hier::Box<NDIM>& box,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const SAMRAI::hier::IntVector<NDIM>& periodic_shift,
        const SAMRAI::tbox::Pointer<LNodeIndexData2>& idx_data);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/LEInteractor.I>

//////////////////////////////////////////////////////////////////////////////


#endif //#ifndef included_LEInteractor
