// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_StokesFirstOrderWaveBcCoef
#define included_IBAMR_StokesFirstOrderWaveBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/ibtk_utilities.h"
#include "ibtk/muParserRobinBcCoefs.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIArrayData.h"
#include "SAMRAIBoundaryBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAIDatabase.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPointer.h"
#include "SAMRAIRobinBcCoefStrategy.h"
#include "SAMRAIVariable.h"

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
 * \brief Class StokesFirstOrderWaveBcCoef is an implementation of the strategy class
 * SAMRAIRobinBcCoefStrategy that provides Dirichlet velocity boundary condition
 * based upon Stokes' lineary theory of water waves at the inlet of the wave tank. This
 * class is meant to be used with INSVCStaggeredHierarchyIntegrator.
 */
class StokesFirstOrderWaveBcCoef : public SAMRAIRobinBcCoefStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    StokesFirstOrderWaveBcCoef(std::string object_name,
                               const int comp_idx,
                               SAMRAIPointer<SAMRAIDatabase> input_db,
                               SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geom);

    /*!
     * \brief Destructor.
     */
    ~StokesFirstOrderWaveBcCoef();

    /*!
     * \name Implementation of SAMRAIRobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \see SAMRAIRobinBcCoefStrategy::setBcCoefs()
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
    void setBcCoefs(SAMRAIPointer<SAMRAIArrayData<double>>& acoef_data,
                    SAMRAIPointer<SAMRAIArrayData<double>>& bcoef_data,
                    SAMRAIPointer<SAMRAIArrayData<double>>& gcoef_data,
                    const SAMRAIPointer<SAMRAIVariable>& variable,
                    const SAMRAIPatch& patch,
                    const SAMRAIBoundaryBox& bdry_box,
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
    SAMRAIIntVector numberOfExtensionsFillable() const override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StokesFirstOrderWaveBcCoef() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StokesFirstOrderWaveBcCoef(const StokesFirstOrderWaveBcCoef& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StokesFirstOrderWaveBcCoef& operator=(const StokesFirstOrderWaveBcCoef& that) = delete;

    /*!
     * Get wave parameters from input db.
     */
    void getFromInput(SAMRAIPointer<SAMRAIDatabase> db);

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
    SAMRAIPointer<SAMRAICartesianGridGeometry> d_grid_geom;

    /*!
     * \brief Wave parameters.
     *
     * \param d_wave_number  : Wave number of dominant wave component [$2\pi/m$]
     * \param d_amplitude    : Amplitude of the dominant wave component [m]
     * \param d_depth        : Depth of water, from sea bed to still water level [m]
     * \param d_gravity      : Acceleration due to gravity [$m/s^2$]
     * \param d_omega        : Angular frequency [$2 \pi/s$]
     *
     * \NOTE d_omega, d_wave_number, and d_depth obtained from the input file
     * should satisfy the dispersion relation to get a consistent Stokes wave.
     */
    double d_depth, d_omega, d_wave_number, d_amplitude, d_gravity;

    /*!
     * Number of interface cells.
     */
    double d_num_interface_cells;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StokesFirstOrderWaveBcCoef
