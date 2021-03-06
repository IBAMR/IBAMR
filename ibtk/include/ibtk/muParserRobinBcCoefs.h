// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_muParserRobinBcCoefs
#define included_IBTK_muParserRobinBcCoefs

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/ibtk_utilities.h"

#include "CartesianGridGeometry.h"
#include "IntVector.h"
#include "RobinBcCoefStrategy.h"
#include "tbox/Pointer.h"

#include "muParser.h"

#include <map>
#include <string>
#include <vector>

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
/*!  \brief Class muParserRobinBcCoefs is an implementation of the strategy
 * class SAMRAI::solv::RobinBcCoefStrategy that allows for the run-time
 * specification of (possibly spatially- and temporally-varying) Robin boundary
 * conditions.
 *
 * This class uses the muParser library to parse strings into mathematical
 * functions. This class is the most common method to specify boundary
 * conditions uses in the examples. The database in the input file contains
 * strings for each coefficient (\f$a\f$, \f$b\f$, and \f$g\f$) at each face
 * of the physical boundary. Strings ending in 0 and 1 are for the \f$x\f$-faces,
 * 2 and 3 are for the \f$y\f$-faces, and 4 and 5 are for the \f$z\f$-faces.
 * These strings can contain spatial and temporal variables using the variable
 * names `X_0`, `X_1`, `X_2`, and `t`. These strings can also contain `if`
 * statements  using the conditional operator
 * `condition ? result_if_true : result_if_false`. For more exotic boundary
 * conditions, one would need to create an extension of the class `RobinBcCoefStrategy`.
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
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    muParserRobinBcCoefs() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    muParserRobinBcCoefs(const muParserRobinBcCoefs& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    muParserRobinBcCoefs& operator=(const muParserRobinBcCoefs& that) = delete;

    /*!
     * Current time value used by the mu::Parser instances.
     *
     * This value is mutable since the mu::Parser objects each store a pointer
     * to it but its specific value (the present time) changes during each
     * call to muParserRobinBcCoefs::setBcCoefs. The alternative would be to
     * rebuild the mu::Parser objects during each call to
     * muParserRobinBcCoefs::setBcCoefs, which is much more expensive. Since
     * this variable is only written to and subsequently read from in that
     * function this is reasonable.
     */
    mutable double d_parser_time = 0.0;

    /*!
     * Current space point used by the mu::Parser instances.
     *
     * This value is mutable for the same reasons that
     * muParserRobinBcCoefs::d_parser_time is mutable.
     */
    mutable Point d_parser_posn;

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
    std::array<mu::Parser, 2 * NDIM> d_acoef_parsers;
    std::array<mu::Parser, 2 * NDIM> d_bcoef_parsers;
    std::array<mu::Parser, 2 * NDIM> d_gcoef_parsers;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_muParserRobinBcCoefs
