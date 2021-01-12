// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_IrregularWaveBcCoef
#define included_IBAMR_IrregularWaveBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/ibtk_utilities.h"
#include "ibtk/muParserRobinBcCoefs.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace geom
{
template <int DIM>
class CartesianGridGeometry;
} // namespace geom
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class ArrayData;
template <int DIM, class TYPE>
class CellVariable;
} // namespace pdat
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IrregularWaveBcCoef is an implementation of the strategy class
 * SAMRAI::solv::RobinBcCoefStrategy that provides Dirichlet velocity boundary condition
 * based upon linear wave theory of water waves to generate irregular waves at the inlet of the wave tank.
 * The class is meant to be used with INSVCStaggeredHierarchyIntegrator.
 *
 * The class can calculate surface elevation and velocities in the water domain in both shallow
 * water regime as well as deep-water regime as indicated through input database.
 *
 */

class IrregularWaveBcCoef : public SAMRAI::solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    IrregularWaveBcCoef(std::string object_name,
                        const int comp_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                        SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);

    /*!
     * \brief Destructor.
     */
    ~IrregularWaveBcCoef();

    /*!
     * \name Implementation of SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \see SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range
     *        for corresponding to the boundary box \a bdry_box and
     *        appropriate for the alignment of the given variable.  If
     *        this is a null pointer, then the calling function is not
     *        interested in a, and you can disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it
     *        is to be filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter
     *        can be used to determine which variable's coefficients
     *        are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is
     *needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are
     *time-dependent.
     */
    void setBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const override;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     *
     * The "extension" used here is the number of cells that a boundary box
     * extends past the patch in the direction parallel to the boundary.
     *
     * Note that the inability to fill the sufficient number of cells past the
     * edge or corner of the patch may preclude the child class from being used
     * in data refinement operations that require the extra data, such as linear
     * refinement.
     *
     * The boundary box that setBcCoefs() is required to fill should not extend
     * past the limits returned by this function.
     */
    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const override;

    //\}

private:
    /*!
     * Get wave parameters from input db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    double getSurfaceElevation(double x, double time) const;

    /*!
     * Get velocity at a specified position and time.
     */
    double getVelocity(double x, double z_plus_d, double time) const;

    /*!
     * Book-keeping.
     */
    std::string d_object_name;

    /*!
     * Velocity component index.
     */
    int d_comp_idx = IBTK::invalid_index;

    /*!
     * muparser object for filling boundary conditions other than wave inlet.
     */
    IBTK::muParserRobinBcCoefs d_muparser_bcs;

    /*!
     * The Cartesian grid geometry object provides the extents of the
     * computational domain.
     */
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*!
     * \brief Wave parameters.
     *
     * \param d_num_waves       : Number of component waves with random phases to be generated (default = 50)
     * \param d_depth           : Depth of water, from sea bed to still water level [$m$]
     * \param d_Hs              : Significant wave height [$m$]
     * \param d_Ts              : significant wave period [$s$]
     * \param d_omega_begin     : Lowest angular frequency in the spectrum [$rad/s$]
     * \param d_omega_end       : Highest angular frequency in the spectrum [$rad/s$]
     * \param d_wave_spectrum   : JONSWAP/Bretschneider wave spectrum.
     * \param d_omega           : Angular frequencies of component waves [$rad/s$]
     * \param d_wave_number     : Wave number of component waves [$2\pi/m$]
     * \param d_amplitude       : Amplitude of component waves [$m$]
     * \param d_phase           : Phase (random) of component waves [$rad$]
     */

    int d_num_waves = 50; // default value is set to 50
    double d_depth, d_omega_begin, d_omega_end, d_Ts, d_Hs, d_gravity;
    std::string d_wave_spectrum;
    std::vector<double> d_omega, d_wave_number, d_amplitude, d_phase;

    /*!
     * Number of interface cells.
     */
    double d_num_interface_cells;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IrregularWaveBcCoef
