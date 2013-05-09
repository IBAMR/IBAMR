// Filename: hdf2vertex.cpp
// Created on 31 May 2007 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "hdf5.h"                 // IWYU pragma: export
#include "hdf5_hl.h"              // IWYU pragma: export
#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace std;

/////////////////////////////// STATIC ///////////////////////////////////////

inline string
discard_comments(
    const string& input_string)
{
    // Create a copy of the input string, but without any text following a '!',
    // '#', or '%' character.
    string output_string = input_string;
    istringstream string_stream;

    // Discard any text following a '!' character.
    string_stream.str(output_string);
    getline(string_stream, output_string, '!');
    string_stream.clear();

    // Discard any text following a '#' character.
    string_stream.str(output_string);
    getline(string_stream, output_string, '#');
    string_stream.clear();

    // Discard any text following a '%' character.
    string_stream.str(output_string);
    getline(string_stream, output_string, '%');
    string_stream.clear();
    return output_string;
}// discard_comments

int
main(
    int argc,
    char* argv[])
{
    if (argc != 3)
    {
        cout << argv[0] << ": a tool to convert IBAMR vertex files from HDF5 to ASCII" << "\n"
             << "USAGE: " << argv[0] << " <input filename> <output filename>" << endl;
        return -1;
    }

    const string input_filename = argv[1];
    const string output_filename = argv[2];

    cout << "input file name: " << input_filename << "\n"
         << "output file name: " << output_filename << "\n";

    // Ensure that the input file exists, and that the output file does not.
    ifstream input_fstream, output_fstream;
    input_fstream.open(input_filename.c_str(), ios::in);
    output_fstream.open(output_filename.c_str(), ios::in);

    if (!input_fstream.is_open())
    {
        cout << "error: Unable to open input file " << input_filename << endl;
        return -1;
    }

    if (output_fstream.is_open())
    {
        cout << "error: Output file " << output_filename << " already exists" << endl;
        return -1;
    }

    input_fstream.close();
    output_fstream.close();

    // Process the HDF5 file.
    hid_t file_id;
    int rank;
    hsize_t dims[2];
    H5T_class_t class_id;
    size_t type_size;
    herr_t status;

    file_id = H5Fopen(input_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0)
    {
        cout << "error: Unable to open input file " << input_filename << endl;
        return -1;
    }

    if (!H5LTfind_dataset(file_id, "vertex"))
    {
        cout << "error: Cannot find vertex dataset in input file " << input_filename << endl;
        return -1;
    }

    status = H5LTget_dataset_ndims(file_id, "vertex", &rank);
    if (rank != 2)
    {
        cout << "error: Invalid dataset rank in input file " << input_filename << endl;
        return -1;
    }

    status = H5LTget_dataset_info(file_id, "vertex", &dims[0], &class_id, &type_size);
    if (dims[0] != NDIM || dims[1] <= 0)
    {
        cout << "error: Invalid dataset dimension in input file " << input_filename << endl;
        return -1;
    }

    static const int num_vertex = dims[1];
    vector<double> vertex_posn(NDIM*num_vertex);
    status = H5LTread_dataset_static_cast<double>(file_id, "vertex", &vertex_posn[0]);

    status = H5Fclose(file_id);

    // Output the vertex file.
    ofstream file_stream;
    file_stream.open(output_filename.c_str(), ios::out);
    if (!file_stream.is_open())
    {
        cout << "error: Unable to open output file " << output_filename << endl;
        return -1;
    }

    // The first entry in the file is the number of vertices.
    ostringstream stream;
    stream << num_vertex;
    string first_line(stream.str());
    first_line.resize((NDIM == 2 ? 46 : 69),' ');
    file_stream << first_line << "# number of vertices\n";

    // Each successive line provides the initial position of each vertex in the
    // input file.
    for (int k = 0; k < num_vertex; ++k)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            file_stream.setf(ios_base::scientific);
            file_stream.precision(16);
            file_stream << vertex_posn[k*NDIM+d];
            if (d < NDIM-1) file_stream << " ";
        }
#if (NDIM == 2)
        if (k == 0) file_stream << " # x-coord, y-coord";
#endif
#if (NDIM == 3)
        if (k == 0) file_stream << " # x-coord, y-coord, z-coord";
#endif
        file_stream << "\n";
    }

    // Close the output file.
    file_stream.close();
    return 0;
}// main

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
