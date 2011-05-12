// Filename: muParserCartGridFunction.C
// Created on 24 Aug 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#include "muParserCartGridFunction.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBTK_config
#include <IBTK_config.h>
#define included_IBTK_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBTK INCLUDES
#include <ibtk/namespaces.h>

// SAMRAI INCLUDES
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIndex.h>
#include <FaceData.h>
#include <FaceIndex.h>
#include <NodeData.h>
#include <NodeIndex.h>
#include <SideData.h>
#include <SideIndex.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

muParserCartGridFunction::muParserCartGridFunction(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<CartesianGridGeometry<NDIM> > grid_geom)
    : CartGridFunction(object_name),
      d_grid_geom(grid_geom),
      d_constants(),
      d_function_strings(),
      d_parsers(),
      d_parser_time(),
      d_parser_posn()
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif
    // Read in user-provided constants.
    Array<std::string> db_key_names = input_db->getAllKeys();
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
        d_parsers.back().SetExpr(d_function_strings.back());
    }
    else
    {
        // Assume vector-valued function.
        int d = 0;
        std::string key_name = "function_0";
        while (input_db->isString(key_name))
        {
            d_function_strings.push_back(input_db->getString(key_name));
            d_parsers.resize(d_parsers.size()+1);
            d_parsers.back().SetExpr(d_function_strings.back());

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
    const double* const xLower = grid_geom->getXLower();
    const double* const xUpper = grid_geom->getXUpper();
    for (std::vector<mu::Parser>::iterator it = d_parsers.begin(); it != d_parsers.end(); ++it)
    {
        // Various names for pi.
        it->DefineConst("pi", pi);
        it->DefineConst("Pi", pi);
        it->DefineConst("PI", pi);

        // The extents of the domain.
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();

            it->DefineConst("X_LOWER" + postfix, xLower[d]);
            it->DefineConst("X_lower" + postfix, xLower[d]);
            it->DefineConst("x_lower" + postfix, xLower[d]);
            it->DefineConst("x_LOWER" + postfix, xLower[d]);
            it->DefineConst("X_Lower" + postfix, xLower[d]);
            it->DefineConst("X_lower" + postfix, xLower[d]);
            it->DefineConst("XLower"  + postfix, xLower[d]);
            it->DefineConst("Xlower"  + postfix, xLower[d]);
            it->DefineConst("x_Lower" + postfix, xLower[d]);
            it->DefineConst("x_lower" + postfix, xLower[d]);
            it->DefineConst("xLower"  + postfix, xLower[d]);
            it->DefineConst("xlower"  + postfix, xLower[d]);

            it->DefineConst("X_LOWER_" + postfix, xLower[d]);
            it->DefineConst("X_lower_" + postfix, xLower[d]);
            it->DefineConst("x_lower_" + postfix, xLower[d]);
            it->DefineConst("x_LOWER_" + postfix, xLower[d]);
            it->DefineConst("X_Lower_" + postfix, xLower[d]);
            it->DefineConst("X_lower_" + postfix, xLower[d]);
            it->DefineConst("XLower_"  + postfix, xLower[d]);
            it->DefineConst("Xlower_"  + postfix, xLower[d]);
            it->DefineConst("x_Lower_" + postfix, xLower[d]);
            it->DefineConst("x_lower_" + postfix, xLower[d]);
            it->DefineConst("xLower_"  + postfix, xLower[d]);
            it->DefineConst("xlower_"  + postfix, xLower[d]);

            it->DefineConst("X_UPPER" + postfix, xUpper[d]);
            it->DefineConst("X_upper" + postfix, xUpper[d]);
            it->DefineConst("x_upper" + postfix, xUpper[d]);
            it->DefineConst("x_UPPER" + postfix, xUpper[d]);
            it->DefineConst("X_Upper" + postfix, xUpper[d]);
            it->DefineConst("X_upper" + postfix, xUpper[d]);
            it->DefineConst("XUpper"  + postfix, xUpper[d]);
            it->DefineConst("Xupper"  + postfix, xUpper[d]);
            it->DefineConst("x_Upper" + postfix, xUpper[d]);
            it->DefineConst("x_upper" + postfix, xUpper[d]);
            it->DefineConst("xUpper"  + postfix, xUpper[d]);
            it->DefineConst("xupper"  + postfix, xUpper[d]);

            it->DefineConst("X_UPPER_" + postfix, xUpper[d]);
            it->DefineConst("X_upper_" + postfix, xUpper[d]);
            it->DefineConst("x_upper_" + postfix, xUpper[d]);
            it->DefineConst("x_UPPER_" + postfix, xUpper[d]);
            it->DefineConst("X_Upper_" + postfix, xUpper[d]);
            it->DefineConst("X_upper_" + postfix, xUpper[d]);
            it->DefineConst("XUpper_"  + postfix, xUpper[d]);
            it->DefineConst("Xupper_"  + postfix, xUpper[d]);
            it->DefineConst("x_Upper_" + postfix, xUpper[d]);
            it->DefineConst("x_upper_" + postfix, xUpper[d]);
            it->DefineConst("xUpper_"  + postfix, xUpper[d]);
            it->DefineConst("xupper_"  + postfix, xUpper[d]);
        }

        // User-provided constants.
        for (std::map<std::string,double>::const_iterator map_cit = d_constants.begin(); map_cit != d_constants.end(); ++map_cit)
        {
            it->DefineConst(map_cit->first, map_cit->second);
        }

        // Variables.
        it->DefineVar("T", &d_parser_time);
        it->DefineVar("t", &d_parser_time);
        for (int d = 0; d < NDIM; ++d)
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
}// muParserCartGridFunction

muParserCartGridFunction::~muParserCartGridFunction()
{
    // intentionally blank
    return;
}// ~muParserCartGridFunction

bool
muParserCartGridFunction::isTimeDependent() const
{
    return true;
}// isTimeDependent

void
muParserCartGridFunction::setDataOnPatch(
    const int data_idx,
    Pointer<Variable<NDIM> > var,
    Pointer<Patch<NDIM> > patch,
    const double data_time,
    const bool initial_time,
    Pointer<PatchLevel<NDIM> > level)
{
    d_parser_time = data_time;

    const Box<NDIM>& patch_box = patch->getBox();
    const Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

    const double* const XLower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    // Set the data in the patch.
    Pointer<PatchData<NDIM> > data = patch->getPatchData(data_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!data.isNull());
#endif
    Pointer<CellData<NDIM,double> > cc_data = data;
    Pointer<FaceData<NDIM,double> > fc_data = data;
    Pointer<NodeData<NDIM,double> > nc_data = data;
    Pointer<SideData<NDIM,double> > sc_data = data;
    if (!cc_data.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == unsigned(cc_data->getDepth()));
#endif
        for (int data_depth = 0; data_depth < cc_data->getDepth(); ++data_depth)
        {
            const int function_depth = (d_parsers.size() == 1 ? 0 : data_depth);
            for (CellIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const CellIndex<NDIM>& i = ic();
                for (int d = 0; d < NDIM; ++d)
                {
                    d_parser_posn[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                }
                (*cc_data)(i,data_depth) = d_parsers[function_depth].Eval();
            }
        }
    }
    else if (!fc_data.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == NDIM || d_parsers.size() == unsigned(fc_data->getDepth()) || d_parsers.size() == NDIM*unsigned(fc_data->getDepth()));
#endif
        for (int data_depth = 0; data_depth < fc_data->getDepth(); ++data_depth)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                int function_depth = -1;
                const int parsers_size = d_parsers.size();
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
                else if (parsers_size == NDIM*fc_data_depth)
                {
                    function_depth = NDIM*data_depth + axis;
                }

                for (FaceIterator<NDIM> ic(patch_box,axis); ic; ic++)
                {
                    const FaceIndex<NDIM>& i = ic();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            d_parser_posn[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d)));
                        }
                        else
                        {
                            d_parser_posn[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                        }
                    }
                    (*fc_data)(i,data_depth) = d_parsers[function_depth].Eval();
                }
            }
        }
    }
    else if (!nc_data.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == unsigned(nc_data->getDepth()));
#endif
        for (int data_depth = 0; data_depth < nc_data->getDepth(); ++data_depth)
        {
            const int function_depth = (d_parsers.size() == 1 ? 0 : data_depth);
            for (NodeIterator<NDIM> ic(patch_box); ic; ic++)
            {
                const NodeIndex<NDIM>& i = ic();
                for (int d = 0; d < NDIM; ++d)
                {
                    d_parser_posn[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d)));
                }
                (*nc_data)(i,data_depth) = d_parsers[function_depth].Eval();
            }
        }
    }
    else if (!sc_data.isNull())
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(d_parsers.size() == 1 || d_parsers.size() == NDIM || d_parsers.size() == unsigned(sc_data->getDepth()) || d_parsers.size() == NDIM*unsigned(sc_data->getDepth()));
#endif
        for (int data_depth = 0; data_depth < sc_data->getDepth(); ++data_depth)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                int function_depth = -1;
                const int parsers_size = d_parsers.size();
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
                else if (parsers_size == NDIM*sc_data_depth)
                {
                    function_depth = NDIM*data_depth + axis;
                }

                for (SideIterator<NDIM> ic(patch_box,axis); ic; ic++)
                {
                    const SideIndex<NDIM>& i = ic();
                    for (int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            d_parser_posn[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d)));
                        }
                        else
                        {
                            d_parser_posn[d] = XLower[d] + dx[d]*(double(i(d)-patch_lower(d))+0.5);
                        }
                    }
                    (*sc_data)(i,data_depth) = d_parsers[function_depth].Eval();
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
}// setDataOnPatch

//////////////////////////////////////////////////////////////////////////////

}// namespace IBTK

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBTK::muParserCartGridFunction>;

//////////////////////////////////////////////////////////////////////////////
