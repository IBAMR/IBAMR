// Filename: StokesFifthOrderWaveBcCoef.h
// Created on 16 Feb 2019 by Amneet Bhalla
//
// Copyright (c) 2002-2019, Amneet Bhalla
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

#ifndef included_IBAMR_StokesFifthOrderWaveBcCoef
#define included_IBAMR_StokesFifthOrderWaveBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "tbox/Pointer.h"

namespace IBAMR
{
class AdvDiffHierarchyIntegrator;
} // namespace IBAMR
namespace IBTK
{
class muParserRobinBcCoefs;
} // namespace IBTK
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
 * \brief Class StokesFifthOrderWaveBcCoef is an implementation of the strategy class
 * SAMRAI::solv::RobinBcCoefStrategy that provides Dirichlet velocity boundary condition
 * based upon Stokes' fifth-order theory of water waves at the inlet of the wave tank.
 * The class is meant to be used with INSVCStaggeredHierarchyIntegrator.
 *
 * The class can calculate surface elevation and velocities in the water domain in both shallow
 * water regime as well as deep-water regime as indicated through input database. The default
 * regime is shallow water (finite depth).
 *
 */
class StokesFifthOrderWaveBcCoef : public SAMRAI::solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    StokesFifthOrderWaveBcCoef(const std::string& object_name,
                               const int comp_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                               SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);

    /*!
     * \brief Destructor.
     */
    ~StokesFifthOrderWaveBcCoef();

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
                    double fill_time = 0.0) const;

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
    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StokesFifthOrderWaveBcCoef() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StokesFifthOrderWaveBcCoef(const StokesFifthOrderWaveBcCoef& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StokesFifthOrderWaveBcCoef& operator=(const StokesFifthOrderWaveBcCoef& that) = delete;

    /*!
     * Get wave parameters from input db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Initialize Stokes coefficients.
     */
    void initStokesCoefficients();

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
    int d_comp_idx;

    /*!
     * muparser object for filling boundary conditions other than wave inlet.
     */
    IBTK::muParserRobinBcCoefs* d_muparser_bcs;

    /*!
     * The Cartesian grid geometry object provides the extents of the
     * computational domain.
     */
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*!
     * \brief Wave parameters.
     *
     * \param d_wave_number  : Wave number of dominant wave component [$2\pi/m$]
     * \param d_amplitude    : Amplitude of the dominant wave component [m]
     * \param d_depth        : Depth of water, from sea bed to still water level [m]
     * \param d_gravity      : Acceleration due to gravity [$m/s^2$]
     * \param d_omega        : Angular frequency [$2 \pi/s$] (optional)
     *
     * \NOTE The class calculates a more accurate value of omega from the expansion coefficients
     * and the provided value in not used.
     */
    double d_depth, d_omega, d_wave_number, d_amplitude, d_gravity;

    /*!
     * If we are calculating in deep water limit.
     */
    bool d_deep_water_limit = false;

    /*!
     * Stokes coefficients.
     */
    double d_A[6][6], d_B[6][6], d_C[5];
    double d_p[5], d_eta[5];

    /*!
     * Number of interface cells.
     */
    double d_num_interface_cells;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_StokesFifthOrderWaveBcCoef
