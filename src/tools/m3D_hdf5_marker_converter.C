// Filename: m3D_hdf5_marker_converter.C
// Created on 30 May 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>
#if (H5_VERS_MINOR == 6)
#include <H5LT.h>
#define H5Gopen1 H5Gopen
#endif
#if (H5_VERS_MINOR == 8)
#include <hdf5_hl.h>
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace std;

/////////////////////////////// STATIC ///////////////////////////////////////

void
build_local_marker_cloud(
    std::ostream& os,
    const std::vector<float>& X,
    const int& nmarks,
    const int& node_offset,
    const int& cloud_number)
{
    for (int k = 0; k < nmarks; ++k)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            os << std::setw(7) << std::fixed << std::setprecision(3) << X[NDIM*k+d] << " ";
        }
        os << std::setw(4) << node_offset+k+1 << " " << std::setw(2) << cloud_number+1 << "\n";
    }
    return;
}// build_local_marker_cloud

int
main(
    int argc,
    char* argv[])
{
    if (argc != 3)
    {
        cout << argv[0] << ": a tool to convert IBAMR HDF5 marker files to myocardial3D marker input files" << "\n"
             << "USAGE: " << argv[0] << " <input filename> <output filename>" << endl;
        return -1;
    }

    const string input_filename = argv[1];
    const string output_filename = argv[2];

    cout << "input file name: " << input_filename << "\n"
         << "output file name: " << output_filename << "\n";

    // Ensure that the input file exists, and that the output file does not.
    ifstream input_fstream;
    input_fstream.open(input_filename.c_str(), ios::in);
    if (!input_fstream.is_open())
    {
        cout << "error: Unable to open input file " << input_filename << endl;
        return -1;
    }
    input_fstream.close();

    // Open and process the HDF5 file.
    ofstream output_fstream;
    output_fstream.open(output_filename.c_str(), ios::out);

    herr_t status;
    hid_t marker_file_id = H5Fopen(input_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (marker_file_id < 0)
    {
        cerr << "error: Unable to open input file: " << input_filename << endl;
        return -1;
    }

    hid_t marker_group_id = H5Gopen1(marker_file_id, "/markers");

    int num_local_marker_nodes, marker_node_offset, num_marker_nodes;
    status = H5LTget_attribute_int(marker_file_id, "/markers", "num_local_marker_nodes", &num_local_marker_nodes);
    status = H5LTget_attribute_int(marker_file_id, "/markers", "marker_node_offset", &marker_node_offset);
    status = H5LTget_attribute_int(marker_file_id, "/markers", "num_marker_nodes", &num_marker_nodes);

    int num_local_marker_clouds, marker_cloud_offset, num_marker_clouds;
    status = H5LTget_attribute_int(marker_file_id, "/markers", "num_local_marker_clouds", &num_local_marker_clouds);
    status = H5LTget_attribute_int(marker_file_id, "/markers", "marker_cloud_offset", &marker_cloud_offset);
    status = H5LTget_attribute_int(marker_file_id, "/markers", "num_marker_clouds", &num_marker_clouds);

    int local_marker_node_counter = 0;
    for (int local_marker_cloud_counter = 0; local_marker_cloud_counter < num_local_marker_clouds;
         ++local_marker_cloud_counter)
    {
        ostringstream dset_name_stream;
        dset_name_stream << "/markers/cloud_" << std::setw(4) << std::setfill('0') << marker_cloud_offset+local_marker_cloud_counter;
        const std::string dset_name = dset_name_stream.str();

        int nmarks;
        status = H5LTget_attribute_int(marker_file_id, dset_name.c_str(), "nmarks", &nmarks);
        int node_offset;
        status = H5LTget_attribute_int(marker_file_id, dset_name.c_str(), "node_offset", &node_offset);
        int cloud_number;
        status = H5LTget_attribute_int(marker_file_id, dset_name.c_str(), "cloud_number", &cloud_number);
        char cloud_name[256];
        status = H5LTget_attribute_string(marker_file_id, dset_name.c_str(), "cloud_name", cloud_name);

        std::vector<float> X(NDIM*nmarks);
        status = H5LTread_dataset_float(marker_file_id, dset_name.c_str(), &X[0]);

        build_local_marker_cloud(output_fstream, X, nmarks, node_offset, cloud_number);

        // Advance the counters.
        local_marker_node_counter += nmarks;
    }

    status = H5Gclose(marker_group_id);
    status = H5Fclose(marker_file_id);

    output_fstream.close();

    return -1;
}// main
