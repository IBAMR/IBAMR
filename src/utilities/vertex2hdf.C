// Filename: vertex2hdf.C
#include "H5LT.h"// Last modified: <30.May.2007 16:58:26 griffith@box221.cims.nyu.edu>
// Created on 30 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>

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
        cout << "error: USAGE: " << argv[0] << " <input filename> <output filename>" << endl;
        return -1;
    }

    const string input_filename = argv[1];
    const string output_filename = argv[2];

    cout << "input file name: " << input_filename << "\n"
         << "output file name: " << output_filename << "\n";

    // The vertex information.
    int num_vertex = -1;
    std::vector<double> vertex_posn;

    // Process the input file.
    ifstream file_stream;
    string line_string;
    file_stream.open(input_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cout << "error: Unable to open input file " << input_filename << endl;
        return -1;
    }

    // The first entry in the file is the number of vertices.
    if (!getline(file_stream, line_string))
    {
        cout << "error: Premature end to input file encountered before line 1 of file " << input_filename << endl;
        return -1;
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_vertex))
        {
            cout << "error: Invalid entry in input file encountered on line 1 of file " << input_filename << endl;
            return -1;
        }
    }

    if (num_vertex <= 0)
    {
        cout << "error: Invalid entry in input file encountered on line 1 of file " << input_filename << endl;
        return -1;
    }

    // Each successive line provides the initial position of each vertex in the
    // input file.
    vertex_posn.resize(num_vertex*NDIM);
    for (int k = 0; k < num_vertex; ++k)
    {
        if (!getline(file_stream, line_string))
        {
            cout << "error: Premature end to input file encountered before line " << k+2 << " of file " << input_filename << endl;
            return -1;
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            for (int d = 0; d < NDIM; ++d)
            {
                if (!(line_stream >> vertex_posn[k*NDIM+d]))
                {
                    cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl;
                    return -1;
                }
            }
        }
    }

    // Close the input file.
    file_stream.close();

    // Create the output file.
    hid_t file_id;
    herr_t status;

    file_id = H5Fcreate(output_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
        cout << "error: Could not create output file: " << output_filename << endl;
        return -1;
    }

    // Create the dataspace for the dataset.
    static const int rank = 2;
    hsize_t dims[rank] = { NDIM , num_vertex };
    hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);

    // Create the dataset with data compression enabled.
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);

    hsize_t cdims[rank] = { NDIM , num_vertex };
    status = H5Pset_chunk(plist_id, rank, cdims);
    //status = H5Pset_shuffle(plist_id);
    status = H5Pset_deflate(plist_id, 4);

    hid_t dataset = H5Dcreate(file_id, "/vertex", H5T_NATIVE_DOUBLE, dataspace_id, plist_id);

    // Write the data to the dataset.
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vertex_posn[0]);

    // Close the output file.
    status = H5Pclose(plist_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

    return 0;
}// main

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
