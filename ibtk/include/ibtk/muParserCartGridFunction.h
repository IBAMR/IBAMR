// Filename: muParserCartGridFunction.h
// Created on 24 Aug 2007 by Boyce Griffith
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

#ifndef included_muParserCartGridFunction
#define included_muParserCartGridFunction

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <map>
#include <string>
#include <vector>

#include "CartesianGridGeometry.h"
#include "PatchLevel.h"
#include "boost/array.hpp"
#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"
#include "muParser.h"
#include "tbox/Pointer.h"

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
    muParserCartGridFunction(const std::string& object_name,
                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                             SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom);

    /*!
     * \brief Empty destructor.
     */
    ~muParserCartGridFunction();

    /*!
     * \name Methods to set patch interior data.
     */
    //\{

    /*!
     * \brief Indicates whether the concrete CartGridFunction object is
     * time-dependent.
     */
    bool isTimeDependent() const;

    /*!
     * \brief Virtual function to evaluate the function on the patch interior.
     */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level =
                            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> >(NULL));

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be
     * used.
     */
    muParserCartGridFunction();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    muParserCartGridFunction(const muParserCartGridFunction& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    muParserCartGridFunction& operator=(const muParserCartGridFunction& that);

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
    double d_parser_time;
    Point d_parser_posn;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_muParserCartGridFunction
