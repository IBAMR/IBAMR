// Filename: muParserCartGridFunction.cpp
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <map>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/muParserCartGridFunction.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "muParser.h"
#include "muParserError.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class PatchLevel;
} // namespace hier
} // namespace SAMRAI

namespace SAMRAI
{
namespace hier
{

class Variable;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

muParserCartGridFunction::muParserCartGridFunction(const std::string& object_name,
                                                   boost::shared_ptr<Database> input_db,
                                                   boost::shared_ptr<CartesianGridGeometry> grid_geom)
    : CartGridFunction(object_name), d_grid_geom(grid_geom), d_constants(), d_function_strings(), d_parsers(),
      d_parser_time(), d_parser_posn()
{
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(input_db);

    // Read in user-provided constants.
    std::vector<std::string> db_key_names = input_db->getAllKeys();
    for (int k = 0; k < db_key_names.size(); ++k)
    {
        const std::string& name = db_key_names[k];
        if (input_db->isDouble(name))
        {
            d_constants[name] = input_db->getDouble(name);
        }
        else if (input_db->isFloat(name))
        {
            d_constants[name] = input_db->getFloat(name);
        }
        else if (input_db->isInteger(name))
        {
            d_constants[name] = input_db->getInteger(name);
        }
    }

    // Initialize the parser(s) with data read in from the input database.
    if (input_db->isString("function"))
    {
        // Assume scalar-valued function.
        d_function_strings.push_back(input_db->getString("function"));
        d_parsers.resize(1);
        try
        {
            d_parsers.back().SetExpr(d_function_strings.back());
        }
        catch (mu::ParserError& e)
        {
            TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                       << "  error: " << e.GetMsg() << "\n"
                       << "  in:    " << e.GetExpr() << "\n");
        }
        catch (...)
        {
            TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                       << "  unrecognized exception generated by muParser library.\n");
        }
    }
    else
    {
        // Assume vector-valued function.
        int d = 0;
        std::string key_name = "function_0";
        while (input_db->isString(key_name))
        {
            d_function_strings.push_back(input_db->getString(key_name));
            d_parsers.resize(d_parsers.size() + 1);
            try
            {
                d_parsers.back().SetExpr(d_function_strings.back());
            }
            catch (mu::ParserError& e)
            {
                TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                           << "  error: " << e.GetMsg() << "\n"
                           << "  in:    " << e.GetExpr() << "\n");
            }
            catch (...)
            {
                TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                           << "  unrecognized exception generated by muParser library.\n");
            }
            ++d;
            std::ostringstream stream;
            stream << "function_" << d;
            key_name = stream.str();
        }
    }

    // Ensure that we've read in at least one function.
    if (d_function_strings.empty())
    {
        TBOX_ERROR("muParserCartGridFunction::muParserCartGridFunction():\n"
                   << "  no function keys found in input database.\n"
                   << "  note that function specifications are assumed to be strings." << std::endl);
    }

    // Define the default and user-provided constants.
    const double pi = 3.1415926535897932384626433832795;
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    for (auto it = d_parsers.begin(); it != d_parsers.end(); ++it)
    {
        // Various names for pi.
        it->DefineConst("pi", pi);
        it->DefineConst("Pi", pi);
        it->DefineConst("PI", pi);

        // The extents of the domain.
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();

            it->DefineConst("X_LOWER" + postfix, x_lower[d]);
            it->DefineConst("X_lower" + postfix, x_lower[d]);
            it->DefineConst("x_lower" + postfix, x_lower[d]);
            it->DefineConst("x_LOWER" + postfix, x_lower[d]);
            it->DefineConst("X_Lower" + postfix, x_lower[d]);
            it->DefineConst("X_lower" + postfix, x_lower[d]);
            it->DefineConst("XLower" + postfix, x_lower[d]);
            it->DefineConst("Xlower" + postfix, x_lower[d]);
            it->DefineConst("x_Lower" + postfix, x_lower[d]);
            it->DefineConst("x_lower" + postfix, x_lower[d]);
            it->DefineConst("xLower" + postfix, x_lower[d]);
            it->DefineConst("xlower" + postfix, x_lower[d]);

            it->DefineConst("X_LOWER_" + postfix, x_lower[d]);
            it->DefineConst("X_lower_" + postfix, x_lower[d]);
            it->DefineConst("x_lower_" + postfix, x_lower[d]);
            it->DefineConst("x_LOWER_" + postfix, x_lower[d]);
            it->DefineConst("X_Lower_" + postfix, x_lower[d]);
            it->DefineConst("X_lower_" + postfix, x_lower[d]);
            it->DefineConst("XLower_" + postfix, x_lower[d]);
            it->DefineConst("Xlower_" + postfix, x_lower[d]);
            it->DefineConst("x_Lower_" + postfix, x_lower[d]);
            it->DefineConst("x_lower_" + postfix, x_lower[d]);
            it->DefineConst("xLower_" + postfix, x_lower[d]);
            it->DefineConst("xlower_" + postfix, x_lower[d]);

            it->DefineConst("X_UPPER" + postfix, x_upper[d]);
            it->DefineConst("X_upper" + postfix, x_upper[d]);
            it->DefineConst("x_upper" + postfix, x_upper[d]);
            it->DefineConst("x_UPPER" + postfix, x_upper[d]);
            it->DefineConst("X_Upper" + postfix, x_upper[d]);
            it->DefineConst("X_upper" + postfix, x_upper[d]);
            it->DefineConst("XUpper" + postfix, x_upper[d]);
            it->DefineConst("Xupper" + postfix, x_upper[d]);
            it->DefineConst("x_Upper" + postfix, x_upper[d]);
            it->DefineConst("x_upper" + postfix, x_upper[d]);
            it->DefineConst("xUpper" + postfix, x_upper[d]);
            it->DefineConst("xupper" + postfix, x_upper[d]);

            it->DefineConst("X_UPPER_" + postfix, x_upper[d]);
            it->DefineConst("X_upper_" + postfix, x_upper[d]);
            it->DefineConst("x_upper_" + postfix, x_upper[d]);
            it->DefineConst("x_UPPER_" + postfix, x_upper[d]);
            it->DefineConst("X_Upper_" + postfix, x_upper[d]);
            it->DefineConst("X_upper_" + postfix, x_upper[d]);
            it->DefineConst("XUpper_" + postfix, x_upper[d]);
            it->DefineConst("Xupper_" + postfix, x_upper[d]);
            it->DefineConst("x_Upper_" + postfix, x_upper[d]);
            it->DefineConst("x_upper_" + postfix, x_upper[d]);
            it->DefineConst("xUpper_" + postfix, x_upper[d]);
            it->DefineConst("xupper_" + postfix, x_upper[d]);
        }

        // User-provided constants.
        for (auto map_cit = d_constants.begin(); map_cit != d_constants.end(); ++map_cit)
        {
            it->DefineConst(map_cit->first, map_cit->second);
        }

        // Variables.
        it->DefineVar("T", &d_parser_time);
        it->DefineVar("t", &d_parser_time);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();
            it->DefineVar("X" + postfix, &(d_parser_posn[d]));
            it->DefineVar("x" + postfix, &(d_parser_posn[d]));
            it->DefineVar("X_" + postfix, &(d_parser_posn[d]));
            it->DefineVar("x_" + postfix, &(d_parser_posn[d]));
        }
    }
    return;
}

muParserCartGridFunction::~muParserCartGridFunction()
{
    // intentionally blank
    return;
}

bool muParserCartGridFunction::isTimeDependent() const
{
    return true;
}

void muParserCartGridFunction::setDataOnPatch(const int data_idx,
                                              boost::shared_ptr<Variable> /*var*/,
                                              boost::shared_ptr<Patch> patch,
                                              const double data_time,
                                              const bool /*initial_time*/,
                                              boost::shared_ptr<PatchLevel> /*level*/)
{
    d_parser_time = data_time;

    const Box& patch_box = patch->getBox();
    const Index& patch_lower = patch_box.lower();
    auto pgeom = patch->getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    // Set the data in the patch.
    boost::shared_ptr<PatchData> data = patch->getPatchData(data_idx);
    TBOX_ASSERT(data);
    boost::shared_ptr<CellData<double> > cc_data = data;
    boost::shared_ptr<FaceData<double> > fc_data = data;
    boost::shared_ptr<NodeData<double> > nc_data = data;
    boost::shared_ptr<SideData<double> > sc_data = data;
    if (cc_data)
    {
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == static_cast<unsigned int>(cc_data->getDepth()));
        for (int data_depth = 0; data_depth < cc_data->getDepth(); ++data_depth)
        {
            const int function_depth = (d_parsers.size() == 1 ? 0 : data_depth);
            for (CellIterator ic(patch_box); ic; ic++)
            {
                const CellIndex& i = ic();
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_parser_posn[d] = XLower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                }
                try
                {
                    (*cc_data)(i, data_depth) = d_parsers[function_depth].Eval();
                }
                catch (mu::ParserError& e)
                {
                    TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                               << "  error: " << e.GetMsg() << "\n"
                               << "  in:    " << e.GetExpr() << "\n");
                }
                catch (...)
                {
                    TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                               << "  unrecognized exception generated by muParser library.\n");
                }
            }
        }
    }
    else if (fc_data)
    {
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == NDIM ||
                    d_parsers.size() == static_cast<unsigned int>(fc_data->getDepth()) ||
                    d_parsers.size() == NDIM * static_cast<unsigned int>(fc_data->getDepth()));
        for (int data_depth = 0; data_depth < fc_data->getDepth(); ++data_depth)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                int function_depth = -1;
                const int parsers_size = static_cast<int>(d_parsers.size());
                const int fc_data_depth = fc_data->getDepth();
                if (parsers_size == 1)
                {
                    function_depth = 0;
                }
                else if (parsers_size == NDIM)
                {
                    function_depth = axis;
                }
                else if (parsers_size == fc_data_depth)
                {
                    function_depth = data_depth;
                }
                else if (parsers_size == NDIM * fc_data_depth)
                {
                    function_depth = NDIM * data_depth + axis;
                }

                for (FaceIterator ic(patch_box, axis); ic; ic++)
                {
                    const FaceIndex& i = ic();
                    const Index& cell_idx = i.toCell(1);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            d_parser_posn[d] = XLower[d] + dx[d] * (static_cast<double>(cell_idx(d) - patch_lower(d)));
                        }
                        else
                        {
                            d_parser_posn[d] =
                                XLower[d] + dx[d] * (static_cast<double>(cell_idx(d) - patch_lower(d)) + 0.5);
                        }
                    }
                    try
                    {
                        (*fc_data)(i, data_depth) = d_parsers[function_depth].Eval();
                    }
                    catch (mu::ParserError& e)
                    {
                        TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                                   << "  error: " << e.GetMsg() << "\n"
                                   << "  in:    " << e.GetExpr() << "\n");
                    }
                    catch (...)
                    {
                        TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                                   << "  unrecognized exception generated by muParser library.\n");
                    }
                }
            }
        }
    }
    else if (nc_data)
    {
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == static_cast<unsigned int>(nc_data->getDepth()));
        for (int data_depth = 0; data_depth < nc_data->getDepth(); ++data_depth)
        {
            const int function_depth = (d_parsers.size() == 1 ? 0 : data_depth);
            for (NodeIterator ic(patch_box); ic; ic++)
            {
                const NodeIndex& i = ic();
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    d_parser_posn[d] = XLower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                }
                try
                {
                    (*nc_data)(i, data_depth) = d_parsers[function_depth].Eval();
                }
                catch (mu::ParserError& e)
                {
                    TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                               << "  error: " << e.GetMsg() << "\n"
                               << "  in:    " << e.GetExpr() << "\n");
                }
                catch (...)
                {
                    TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                               << "  unrecognized exception generated by muParser library.\n");
                }
            }
        }
    }
    else if (sc_data)
    {
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == NDIM ||
                    d_parsers.size() == static_cast<unsigned int>(sc_data->getDepth()) ||
                    d_parsers.size() == NDIM * static_cast<unsigned int>(sc_data->getDepth()));
        for (int data_depth = 0; data_depth < sc_data->getDepth(); ++data_depth)
        {
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                int function_depth = -1;
                const int parsers_size = static_cast<int>(d_parsers.size());
                const int sc_data_depth = sc_data->getDepth();
                if (parsers_size == 1)
                {
                    function_depth = 0;
                }
                else if (parsers_size == NDIM)
                {
                    function_depth = axis;
                }
                else if (parsers_size == sc_data_depth)
                {
                    function_depth = data_depth;
                }
                else if (parsers_size == NDIM * sc_data_depth)
                {
                    function_depth = NDIM * data_depth + axis;
                }

                for (SideIterator ic(patch_box, axis); ic; ic++)
                {
                    const SideIndex& i = ic();
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            d_parser_posn[d] = XLower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                        }
                        else
                        {
                            d_parser_posn[d] = XLower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                        }
                    }
                    try
                    {
                        (*sc_data)(i, data_depth) = d_parsers[function_depth].Eval();
                    }
                    catch (mu::ParserError& e)
                    {
                        TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                                   << "  error: " << e.GetMsg() << "\n"
                                   << "  in:    " << e.GetExpr() << "\n");
                    }
                    catch (...)
                    {
                        TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                                   << "  unrecognized exception generated by muParser library.\n");
                    }
                }
            }
        }
    }
    else
    {
        TBOX_ERROR("muParserCartGridFunction::setDataOnPatch():\n"
                   << "  unsupported patch data type encountered." << std::endl);
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
