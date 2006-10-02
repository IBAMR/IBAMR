//
// LEInteractor.h
//
// Created on 14 Jul 2004
//         by Boyce Griffith (boyce@trasnaform.speakeasy.net).
//
// Last modified: <29.Jun.2005 18:47:16 boyce@mstu1.cims.nyu.edu>
//

#ifndef included_LEInteractor
#define included_LEInteractor

// SAMRAI-tools INCLUDES
//
#include "LNodeIndexData.h"
#include "LNodeLevelData.h"

// SAMRAI INCLUDES
//
#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

#include "Box.h"
#include "CellData.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"

using namespace SAMRAI;
using namespace std;

// CLASS DEFINITION
//

/*!
 * @brief Class LEInteractor is a Lagrangian--Eulerian utilities class
 * that defines several functions to interpolate data from Eulerian
 * grid patches onto Lagrangian meshes and to spread densities from
 * Lagrangian meshes to Eulerian grid patches.
 */
class LEInteractor
{
public:
    /*!
     * @brief Initialize the Timer objects employed by the
     * LEInteractor class.
     *
     * NOTE: It is necessary to initialize the Timer objects prior to
     * using any of the functionality provided by this class.
     */
    static void initializeTimers();
    
    /*!
     * @brief Returns the interpolation/spreading stencil
     * corresponding to the specified weighting function.
     *
     * The return value is -1 for any unknown weighting function type.
     */
    static int getStencilSize(
        const string& weighting_fcn);
    
    /*!
     * @brief Interpolate data from an Eulerian grid to a Lagrangian
     * mesh.  The positions of the nodes of the Lagrangian mesh are
     * specified by X_data.
     *
     * NOTE: This method employs periodic boundary conditions where
     * appropriate and when requested.  X_data must provide the
     * cannonical location of the node---i.e., each node location must
     * lie within the extents of the physical domain.
     * 
     * NOTE: The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation
     * operation.
     */
    static void interpolate(
        tbox::Pointer<LNodeLevelData>& Q_data,
        const tbox::Pointer<LNodeLevelData>& X_data,
        const tbox::Pointer<LNodeIndexData>& idx_data,
        const tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
        const tbox::Pointer<hier::Patch<NDIM> >& patch,
        const hier::Box<NDIM>& box,
        const string& interp_fcn="PIECEWISE_CUBIC",
        const bool enforce_periodic_bcs=false);
    
    /*!
     * @brief Interpolate data from an Eulerian grid to a Lagrangian
     * mesh.  The positions of the nodes of the Lagrangian mesh are
     * specified by X_data.
     *
     * NOTE: This method employs periodic boundary conditions where
     * appropriate and when requested.  X_data must provide the
     * cannonical location of the node---i.e., each node location must
     * lie within the extents of the physical domain.
     * 
     * NOTE: The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation
     * operation.
     */
    static void interpolate(
        double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const tbox::Pointer<LNodeIndexData>& idx_data,
        const tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
        const tbox::Pointer<hier::Patch<NDIM> >& patch,
        const hier::Box<NDIM>& box,
        const string& interp_fcn="PIECEWISE_CUBIC",
        const bool enforce_periodic_bcs=false);
    
    /*!
     * @brief Interpolate data from an Eulerian grid to a Lagrangian
     * mesh.  The positions of the nodes of the Lagrangian mesh are
     * specified by X_data.
     *
     * NOTE: This method does not implement periodic boundary
     * conditions!
     * 
     * NOTE: The interpolation operator implements the operation
     *
     *     Q(q,r,s) = Sum_{i,j,k} q(i,j,k) delta_h(x(i,j,k) - X(q,r,s)) h^3
     *
     * This is the standard regularized delta function interpolation
     * operation.
     */
    static void interpolate(
        double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const int num_vals,
        const tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
        const tbox::Pointer<hier::Patch<NDIM> >& patch,
        const hier::Box<NDIM>& box,
        const string& interp_fcn="PIECEWISE_CUBIC");
    
    /*!
     * @brief Spread data from a Lagrangian mesh to an Eulerian grid.
     * The positions of the nodes of the Lagrangian mesh are specified
     * by X_data.
     *
     * NOTE: This method employs periodic boundary conditions where
     * appropriate and when requested.  X_data must provide the
     * cannonical location of the node---i.e., each node location must
     * lie within the extents of the physical domain.
     * 
     * NOTE: The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).
     * The spreading formula is
     *
     *     q(i,j,k) = Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading
     * operation, the implemented operations spreads values, NOT
     * densities.
     */
    static void spread(
        tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
        const tbox::Pointer<LNodeLevelData>& Q_data,
        const tbox::Pointer<LNodeLevelData>& X_data,
        const tbox::Pointer<LNodeIndexData>& idx_data,
        const tbox::Pointer<hier::Patch<NDIM> >& patch,
        const hier::Box<NDIM>& box,
        const string& spread_fcn="PIECEWISE_CUBIC",
        const bool enforce_periodic_bcs=false);

    /*!
     * @brief Spread data from a Lagrangian mesh to an Eulerian grid.
     * The positions of the nodes of the Lagrangian mesh are specified
     * by X_data.
     *
     * NOTE: This method employs periodic boundary conditions where
     * appropriate and when requested.  X_data must provide the
     * cannonical location of the node---i.e., each node location must
     * lie within the extents of the physical domain.
     * 
     * NOTE: The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).
     * The spreading formula is
     *
     *     q(i,j,k) = Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading
     * operation, the implemented operations spreads values, NOT
     * densities.
     */
    static void spread(
        tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
        const double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const tbox::Pointer<LNodeIndexData>& idx_data,
        const tbox::Pointer<hier::Patch<NDIM> >& patch,
        const hier::Box<NDIM>& box,
        const string& spread_fcn="PIECEWISE_CUBIC",
        const bool enforce_periodic_bcs=false);

    /*!
     * @brief Spread data from a Lagrangian mesh to an Eulerian grid.
     * The positions of the nodes of the Lagrangian mesh are specified
     * by X_data.
     *
     * NOTE: This method does not implement periodic boundary
     * conditions!
     * 
     * NOTE: The spreading operation DOES NOT include the scale factor
     * corresponding to the curvilinear volume element (dq dr ds).
     * The spreading formula is
     *
     *     q(i,j,k) = Sum_{q,r,s} Q(q,r,s) delta_h(x(i,j,k) - X(q,r,s))
     *
     * Unlike the standard regularized delta function spreading
     * operation, the implemented operations spreads values, NOT
     * densities.
     */
    static void spread(
        tbox::Pointer<pdat::CellData<NDIM,double> > q_data,
        const double* const Q_data,
        const int Q_depth,
        const double* const X_data,
        const int X_depth,
        const int num_vals,
        const tbox::Pointer<hier::Patch<NDIM> >& patch,
        const hier::Box<NDIM>& box,
        const string& spread_fcn="PIECEWISE_CUBIC");

private:
    /*!
     * @brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    LEInteractor();
    
    /*!
     * @brief Default destructor constructor.
     *
     * NOTE: This destructor is not implemented and should not be
     * used.
     */
    ~LEInteractor();
    
    /*!
     * @brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     * 
     * @param from The value to copy to this object.
     */
    LEInteractor(
        const LEInteractor& from);
    
    /*!
     * @brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     * 
     * @param that The value to assign to this object.
     * 
     * @return A reference to this object.
     */
    LEInteractor& operator=(
        const LEInteractor& that);
};

// INLINED FUNCTION DEFINITIONS
//
//#ifndef DEBUG_NO_INLINE
//#include "LEInteractor.I"
//#endif

#endif //#ifndef included_LEInteractor

//////////////////////////////////////////////////////////////////////////////
