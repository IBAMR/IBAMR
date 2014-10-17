// Filename: muParserRobinBcCoefs.h
// Created on 25 Aug 2007 by Boyce Griffith
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

#ifndef included_muParserRobinBcCoefs
#define included_muParserRobinBcCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <string>
#include <vector>

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "muParser.h"
#include "ibtk/ibtk_utilities.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
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
 * \brief Class muParserRobinBcCoefs is an implementation of the strategy class
 * SAMRAI::solv::RobinBcCoefStrategy that allows for the run-time specification
 * of (possibly spatially- and temporally-varying) Robin boundary conditions.
 *
 * \warning Not all linear solvers in IBTK properly handle time-varying \em
 * homogeneous Robin boundary condition coefficients.  Note, however, that all
 * linear solvers in IBTK are presently designed to support spatially and
 * temporally varying \em inhomogeneous boundary coefficients.
 */
class muParserRobinBcCoefs : public SAMRAI::solv::RobinBcCoefStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor
     */
    muParserRobinBcCoefs(const std::string& object_name,
                         SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                         SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);

    /*!
     * \brief Destructor.
     */
    ~muParserRobinBcCoefs();

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
    muParserRobinBcCoefs();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    muParserRobinBcCoefs(const muParserRobinBcCoefs& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    muParserRobinBcCoefs& operator=(const muParserRobinBcCoefs& that);

    /*!
     * The Cartesian grid geometry object provides the extents of the
     * computational domain.
     */
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > d_grid_geom;

    /*!
     * User-provided constants specified in the input file.
     */
    std::map<std::string, double> d_constants;

    /*!
     * The strings providing the data-setting functions which are evaluated by the
     * mu::Parser objects.
     */
    std::vector<std::string> d_acoef_function_strings;
    std::vector<std::string> d_bcoef_function_strings;
    std::vector<std::string> d_gcoef_function_strings;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser> d_acoef_parsers;
    std::vector<mu::Parser> d_bcoef_parsers;
    std::vector<mu::Parser> d_gcoef_parsers;

    /*!
     * Time and position variables.
     */
    double* d_parser_time;
    Point* d_parser_posn;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_muParserRobinBcCoefs
