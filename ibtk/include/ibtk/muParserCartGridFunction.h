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

#ifndef included_IBTK_muParserCartGridFunction
#define included_IBTK_muParserCartGridFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"

#include "CartesianGridGeometry.h"
#include "PatchLevel.h"
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
class Patch;
template <int DIM>
class Variable;
} // namespace hier
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class muParserCartGridFunction is an implementation of the strategy
 * class CartGridFunction that allows for the run-time specification of
 * (possibly spatially- and temporally-varying) functions which are used to set
 * double precision values on standard SAMRAI SAMRAI::hier::PatchData objects.
 */
class muParserCartGridFunction : public CartGridFunction
{
public:
    /*!
     * \brief Constructor.
     */
    muParserCartGridFunction(std::string object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);

    /*!
     * \brief Empty destructor.
     */
    ~muParserCartGridFunction() = default;

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool isTimeDependent() const override;

    /*!
     * \brief Virtual function to evaluate the function on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL)) override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    muParserCartGridFunction() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    muParserCartGridFunction(const muParserCartGridFunction& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    muParserCartGridFunction& operator=(const muParserCartGridFunction& that) = delete;

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
    std::vector<std::string> d_function_strings;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser> d_parsers;

    /*!
     * Time and position variables.
     */
    double d_parser_time = 0.0;
    Point d_parser_posn;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_muParserCartGridFunction
