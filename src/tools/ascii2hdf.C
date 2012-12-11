// Filename: ascii2hdf.C
// Created on 30 May 2007 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <hdf5.h>
#if (H5_VERS_MINOR == 6)
#include <H5LT.h>
#define H5Dcreate1 H5Dcreate
#define H5Gcreate1 H5Gcreate
#endif
#if (H5_VERS_MINOR == 8)
#include <hdf5_hl.h>
#endif

#include <blitz/array.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace std;

/////////////////////////////// STATIC ///////////////////////////////////////

void initializeVertexData(         const string& base_filename, hid_t file_id);
void initializeSpringData(         const string& base_filename, hid_t file_id);
void initializeBeamData(           const string& base_filename, hid_t file_id);
void initializeTargetPointData(    const string& base_filename, hid_t file_id);
void initializeMassData(           const string& base_filename, hid_t file_id);
void initializeInstrumentationData(const string& base_filename, hid_t file_id);

int
main(
    int argc,
    char* argv[])
{
    if (argc != 2)
    {
        cout << argv[0] << ": a tool to convert IBAMR input files from ASCII to HDF5 formats\n"
             << "USAGE: " << argv[0] << " <base filename>\n";
        return -1;
    }

    const string base_filename = argv[1];

    cout << argv[0] << ": a tool to convert IBAMR input files from ASCII to HDF5 formats\n\n"
         << "base filename: " << base_filename << "\n\n"
         << "required input files: " << base_filename << ".vertex\n"
         << "optional input files: " << base_filename << ".spring, " << base_filename << ".beam, " << base_filename << ".target, " << base_filename << ".mass, " << base_filename << ".inst \n\n"
         << "output HDF5 filename: " << base_filename << ".h5\n\n";

    // Create the output file.
    const string output_filename = base_filename + ".h5";
    hid_t file_id = H5Fcreate(output_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
        cout << "error: could not create output file: " << output_filename << "\n"
             << "       perhaps the file already exists?\n";
        return -1;
    }

    // Create the appropriate group(s).
    const string base_group_name = "/" + base_filename;
    hid_t base_group_id = H5Gcreate1(file_id, base_group_name.c_str(), 0);

    // Populate the datasets.
    initializeVertexData(base_filename, file_id);
    initializeSpringData(base_filename, file_id);
    initializeBeamData(base_filename, file_id);
    initializeTargetPointData(base_filename, file_id);
    initializeMassData(base_filename, file_id);
    initializeInstrumentationData(base_filename, file_id);

    // Cleanup HDF5 data structures.
    H5Gclose(base_group_id);
    H5Fclose(file_id);
    return 0;
}// main

static const int BUFFER_SIZE = 2097152;  // = 16 * 2^20 / 8 ~ 16 MB of double precision values

static int num_vertex = -1;
static int num_spring = -1;
static int num_beam = -1;
static int num_target_point = -1;
static int num_mass_point = -1;
static int num_inst = -1;
static int num_inst_point = -1;

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

void
initializeVertexData(
    const string& base_filename,
    hid_t file_id)
{
    string line_string;

    // Ensure that the file exists.
    const string vertex_filename = base_filename + ".vertex";
    ifstream file_stream;
    file_stream.open(vertex_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cerr << "error: unable to open vertex file " << vertex_filename << "\n";
        abort();
    }

    cout << "processing vertex data from ASCII input filename " << vertex_filename << "\n";

    // The first entry in the file is the number of vertices.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line 1 of file " << vertex_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_vertex))
        {
            cerr << "error: invalid entry in input file encountered on line 1 of file " << vertex_filename << "\n";
            abort();
        }
    }

    if (num_vertex <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line 1 of file " << vertex_filename << "\n";
        abort();
    }

    // Create the appropriate group(s).
    const string vertex_group_name = "/" + base_filename + "/vertex";
    hid_t vertex_group_id = H5Gcreate1(file_id, vertex_group_name.c_str(), 0);

    // Define the file dataspace.
    static const int rankf = 2;
    hsize_t dimsf[rankf] = { static_cast<hsize_t>(num_vertex) , NDIM };
    hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

    // Define the memory dataspace.
    static const int rankm = 2;
    hsize_t dimsm[rankm] = { BUFFER_SIZE , NDIM };
    hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

    // Create the dataset with data compression enabled.
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    static const int rankc = 2;
    hsize_t dimsc[rankc] = { static_cast<hsize_t>(min(num_vertex,BUFFER_SIZE)) , NDIM };
    H5Pset_chunk(plist, rankc, dimsc);
    H5Pset_deflate(plist, 6);

    const string vertex_posn_dset_name = "/" + base_filename + "/vertex/posn";
    hid_t posn_dataset = H5Dcreate1(file_id, vertex_posn_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    // Each successive line provides the initial position of each vertex in the
    // input file.
    vector<double> posn_buf(NDIM*BUFFER_SIZE);
    const int num_blocks = num_vertex/BUFFER_SIZE + (num_vertex%BUFFER_SIZE == 0 ? 0 : 1);
    for (int k = 0, block = 0; k < num_vertex; ++k)
    {
        if (k%10 == 0 || k == num_vertex-1)
        {
            cout << "\r" << setw(3) << static_cast<int>(100.0*static_cast<double>(k)/static_cast<double>(num_vertex-1)) << "% complete";
            cout.flush();
        }

        // Read in the data.
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << k+2 << " of file " << vertex_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (!(line_stream >> posn_buf[(NDIM*(k%BUFFER_SIZE)+d)]))
                {
                    cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << vertex_filename << "\n";
                    abort();
                }
            }
        }

        // Write out data as the buffer fills up.
        if ((k+1)%BUFFER_SIZE == 0 || k+1 == num_vertex)
        {
            // Determine whether we are writing out the last block in the file.
            const bool last_block = (block == num_blocks-1);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(!last_block || (last_block && k+1 == num_vertex));
#endif
            // Determine the number of items to read (always BUFFER_SIZE except
            // for the final block in the file).
            const int num_vertex_block = (last_block ? num_vertex - block*BUFFER_SIZE : BUFFER_SIZE);
#ifdef DEBUG_CHECK_ASSERTIONS
            assert(num_vertex_block > 0 && num_vertex_block <= BUFFER_SIZE);
#endif
            // Define the file hyperslab.
            hsize_t offsetf[rankf] = { static_cast<hsize_t>(block*BUFFER_SIZE) , 0 };
            hsize_t countf[rankf] = { static_cast<hsize_t>(num_vertex_block) , NDIM };
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

            // Define the memory hyperslab.
            hsize_t offsetm[rankm] = { 0 , 0 };
            hsize_t countm[rankm] = { static_cast<hsize_t>(num_vertex_block) , NDIM };
            H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

            // Write data to the hyperslab in the file from the hyperslab in
            // memory.
            H5Dwrite(posn_dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &posn_buf[0]);

            // Increment the block counter.
            ++block;
        }
    }

    cout << "\n\n";

    // Close the input file.
    file_stream.close();

    // Cleanup HDF5 data structures.
    H5Pclose(plist);
    H5Dclose(posn_dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Gclose(vertex_group_id);
    return;
}// initializeVertexData

void
initializeSpringData(
    const string& base_filename,
    hid_t file_id)
{
    string line_string;

    // Ensure that the file exists.
    const string spring_filename = base_filename + ".spring";
    ifstream file_stream;
    file_stream.open(spring_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cerr << "warning: unable to open spring file " << spring_filename << "\n\n";
        return;
    }

    cout << "processing spring data from ASCII input filename " << spring_filename << "\n";

    // The first line in the file indicates the number of springs in the input
    // file.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line 1 of file " << spring_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_spring))
        {
            cerr << "error: invalid entry in input file encountered on line 1 of file " << spring_filename << "\n";
            abort();
        }
    }

    if (num_spring <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line 1 of file " << spring_filename << "\n";
        abort();
    }

    // Create the appropriate group(s).
    const string spring_group_name = "/" + base_filename + "/spring";
    hid_t spring_group_id = H5Gcreate1(file_id, spring_group_name.c_str(), 0);

    // Define the file dataspace.
    static const int rankf = 1;
    hsize_t dimsf[rankf] = { static_cast<hsize_t>(num_spring) };
    hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

    // Define the memory dataspace.
    static const int rankm = 1;
    hsize_t dimsm[rankm] = { BUFFER_SIZE };
    hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

    // Create the datasets with data compression enabled.
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    static const int rankc = 1;
    hsize_t dimsc[rankc] = { static_cast<hsize_t>(min(num_spring,BUFFER_SIZE)) };
    H5Pset_chunk(plist, rankc, dimsc);
    H5Pset_shuffle(plist);
    H5Pset_deflate(plist, 6);

    const string spring_node1_idx_dset_name = "/" + base_filename + "/spring/node1_idx";
    hid_t node1_idx_dataset = H5Dcreate1(file_id, spring_node1_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string spring_node2_idx_dset_name = "/" + base_filename + "/spring/node2_idx";
    hid_t node2_idx_dataset = H5Dcreate1(file_id, spring_node2_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string spring_force_fcn_idx_dset_name = "/" + base_filename + "/spring/force_fcn_idx";
    hid_t force_fcn_idx_dataset = H5Dcreate1(file_id, spring_force_fcn_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string spring_stiffness_dset_name = "/" + base_filename + "/spring/stiffness";
    hid_t stiffness_dataset = H5Dcreate1(file_id, spring_stiffness_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    const string spring_rest_length_dset_name = "/" + base_filename + "/spring/rest_length";
    hid_t rest_length_dataset = H5Dcreate1(file_id, spring_rest_length_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    // Each successive line provides the connectivity and material parameter
    // information for each spring in the structure.
    vector<int> node1_idx_buf(BUFFER_SIZE), node2_idx_buf(BUFFER_SIZE), force_fcn_idx_buf(BUFFER_SIZE);
    vector<double> stiffness_buf(BUFFER_SIZE), rest_length_buf(BUFFER_SIZE);
    const int num_blocks = num_spring/BUFFER_SIZE + (num_spring%BUFFER_SIZE == 0 ? 0 : 1);
    for (int k = 0, block = 0; k < num_spring; ++k)
    {
        if (k%10 == 0 || k == num_spring-1)
        {
            cout << "\r" << setw(3) << static_cast<int>(100.0*static_cast<double>(k)/static_cast<double>(num_spring-1)) << "% complete";
            cout.flush();
        }

        int node1_idx, node2_idx, force_fcn_idx;
        double stiffness, rest_length;
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << k+2 << " of file " << spring_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            if (!(line_stream >> node1_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n";
                abort();
            }
            else if ((node1_idx < 0) || (node1_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n"
                     << "       vertex index " << node1_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> node2_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n";
                abort();
            }
            else if ((node2_idx < 0) || (node2_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n"
                     << "       vertex index " << node2_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> stiffness))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n";
                abort();
            }
            else if (stiffness < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n"
                     << "       spring constant is negative\n";
                abort();
            }

            if (!(line_stream >> rest_length))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n";
                abort();
            }
            else if (rest_length < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << spring_filename << "\n"
                     << "       spring resting length is negative\n";
                abort();
            }

            if (!(line_stream >> force_fcn_idx))
            {
                force_fcn_idx = 0;  // default force function specification.
            }

            // Store the values in the buffers.
            node1_idx_buf    [k%BUFFER_SIZE] = node1_idx;
            node2_idx_buf    [k%BUFFER_SIZE] = node2_idx;
            force_fcn_idx_buf[k%BUFFER_SIZE] = force_fcn_idx;
            stiffness_buf    [k%BUFFER_SIZE] = stiffness;
            rest_length_buf  [k%BUFFER_SIZE] = rest_length;

            // Write out data as the buffer fills up.
            if ((k+1)%BUFFER_SIZE == 0 || k+1 == num_spring)
            {
                // Determine whether we are writing out the last block in the
                // file.
                const bool last_block = (block == num_blocks-1);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!last_block || (last_block && k+1 == num_spring));
#endif
                // Determine the number of items to read (always BUFFER_SIZE
                // except for the final block in the file).
                const int num_spring_block = (last_block ? num_spring - block*BUFFER_SIZE : BUFFER_SIZE);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(num_spring_block > 0 && num_spring_block <= BUFFER_SIZE);
#endif
                // Define the file hyperslab.
                hsize_t offsetf[rankf] = { static_cast<hsize_t>(block*BUFFER_SIZE) };
                hsize_t countf[rankf] = { static_cast<hsize_t>(num_spring_block) };
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

                // Define the memory hyperslab.
                hsize_t offsetm[rankm] = { 0 };
                hsize_t countm[rankm] = { static_cast<hsize_t>(num_spring_block) };
                H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

                // Write data to the hyperslabs in the file from the hyperslabs
                // in memory.
                H5Dwrite(node1_idx_dataset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node1_idx_buf    [0]);
                H5Dwrite(node2_idx_dataset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node2_idx_buf    [0]);
                H5Dwrite(force_fcn_idx_dataset, H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &force_fcn_idx_buf[0]);
                H5Dwrite(stiffness_dataset    , H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &stiffness_buf    [0]);
                H5Dwrite(rest_length_dataset  , H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &rest_length_buf  [0]);

                // Increment the block counter.
                ++block;
            }
        }
    }

    cout << "\n\n";

    // Close the input file.
    file_stream.close();

    // Cleanup HDF5 data structures.
    H5Pclose(plist);
    H5Dclose(node1_idx_dataset);
    H5Dclose(node2_idx_dataset);
    H5Dclose(force_fcn_idx_dataset);
    H5Dclose(stiffness_dataset);
    H5Dclose(rest_length_dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Gclose(spring_group_id);
    return;
}// initializeSpringData

void
initializeBeamData(
    const string& base_filename,
    hid_t file_id)
{
    string line_string;

    // Ensure that the file exists.
    const string beam_filename = base_filename + ".beam";
    ifstream file_stream;
    file_stream.open(beam_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cerr << "warning: unable to open beam file " << beam_filename << "\n\n";
        return;
    }

    cout << "processing beam data from ASCII input filename " << beam_filename << "\n";

    // The first line in the file indicates the number of beams in the input
    // file.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line 1 of file " << beam_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_beam))
        {
            cerr << "error: invalid entry in input file encountered on line 1 of file " << beam_filename << "\n";
            abort();
        }
    }

    if (num_beam <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line 1 of file " << beam_filename << "\n";
        abort();
    }

    // Create the appropriate group(s).
    const string beam_group_name = "/" + base_filename + "/beam";
    hid_t beam_group_id = H5Gcreate1(file_id, beam_group_name.c_str(), 0);

    // Define the file dataspace.
    static const int rankf = 1;
    hsize_t dimsf[rankf] = { static_cast<hsize_t>(num_beam) };
    hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

    // Define the memory dataspace.
    static const int rankm = 1;
    hsize_t dimsm[rankm] = { BUFFER_SIZE };
    hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

    // Create the datasets with data compression enabled.
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    static const int rankc = 1;
    hsize_t dimsc[rankc] = { static_cast<hsize_t>(min(num_beam,BUFFER_SIZE)) };
    H5Pset_chunk(plist, rankc, dimsc);
    H5Pset_shuffle(plist);
    H5Pset_deflate(plist, 6);

    const string beam_node1_idx_dset_name = "/" + base_filename + "/beam/node1_idx";
    hid_t node1_idx_dataset = H5Dcreate1(file_id, beam_node1_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string beam_node2_idx_dset_name = "/" + base_filename + "/beam/node2_idx";
    hid_t node2_idx_dataset = H5Dcreate1(file_id, beam_node2_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string beam_node3_idx_dset_name = "/" + base_filename + "/beam/node3_idx";
    hid_t node3_idx_dataset = H5Dcreate1(file_id, beam_node3_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string beam_bend_rigidity_dset_name = "/" + base_filename + "/beam/bend_rigidity";
    hid_t bend_rigidity_dataset = H5Dcreate1(file_id, beam_bend_rigidity_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    hid_t rest_curvature_dataset[NDIM];
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        ostringstream os;
        os << "_" << d;
        const string beam_rest_curvature_dset_name = "/" + base_filename + "/beam/rest_curvature" + os.str();
        rest_curvature_dataset[d] = H5Dcreate1(file_id, beam_rest_curvature_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);
    }

    // Each successive line provides the connectivity and material parameter
    // information for each beam in the structure.
    vector<int> node1_idx_buf(BUFFER_SIZE), node2_idx_buf(BUFFER_SIZE), node3_idx_buf(BUFFER_SIZE);
    vector<double> bend_rigidity_buf(BUFFER_SIZE);
    vector<vector<double> > rest_curvature_buf(NDIM,vector<double>(BUFFER_SIZE));
    const int num_blocks = num_beam/BUFFER_SIZE + (num_beam%BUFFER_SIZE == 0 ? 0 : 1);
    for (int k = 0, block = 0; k < num_beam; ++k)
    {
        if (k%10 == 0 || k == num_beam-1)
        {
            cout << "\r" << setw(3) << static_cast<int>(100.0*static_cast<double>(k)/static_cast<double>(num_beam-1)) << "% complete";
            cout.flush();
        }

        int node1_idx, node2_idx, node3_idx;
        double bend_rigidity;
        blitz::TinyVector<double,NDIM> rest_curvature;
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << k+2 << " of file " << beam_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            if (!(line_stream >> node1_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n";
                abort();
            }
            else if ((node1_idx < 0) || (node1_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n"
                     << "       vertex index " << node1_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> node2_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n";
                abort();
            }
            else if ((node2_idx < 0) || (node2_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n"
                     << "       vertex index " << node2_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> node3_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n";
                abort();
            }
            else if ((node3_idx < 0) || (node3_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n"
                     << "       vertex index " << node3_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> bend_rigidity))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n";
                abort();
            }
            else if (bend_rigidity < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n"
                     << "       beam bending rigidity is negative\n";
                abort();
            }

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (!(line_stream >> rest_curvature[d]))
                {
                    if (d > 0)
                    {
                        cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << beam_filename << "\n";
                        abort();
                    }
                    else
                    {
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            rest_curvature[d] = 0.0;
                        }
                        break;
                    }
                }
            }

            // Store the values in the buffers.
            node1_idx_buf    [k%BUFFER_SIZE] = node1_idx;
            node2_idx_buf    [k%BUFFER_SIZE] = node2_idx;
            node3_idx_buf    [k%BUFFER_SIZE] = node3_idx;
            bend_rigidity_buf[k%BUFFER_SIZE] = bend_rigidity;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                rest_curvature_buf[d][k%BUFFER_SIZE] = rest_curvature[d];
            }

            // Write out data as the buffer fills up.
            if ((k+1)%BUFFER_SIZE == 0 || k+1 == num_beam)
            {
                // Determine whether we are writing out the last block in the
                // file.
                const bool last_block = (block == num_blocks-1);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!last_block || (last_block && k+1 == num_beam));
#endif
                // Determine the number of items to read (always BUFFER_SIZE
                // except for the final block in the file).
                const int num_beam_block = (last_block ? num_beam - block*BUFFER_SIZE : BUFFER_SIZE);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(num_beam_block > 0 && num_beam_block <= BUFFER_SIZE);
#endif
                // Define the file hyperslab.
                hsize_t offsetf[rankf] = { static_cast<hsize_t>(block*BUFFER_SIZE) };
                hsize_t countf[rankf] = { static_cast<hsize_t>(num_beam_block) };
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

                // Define the memory hyperslab.
                hsize_t offsetm[rankm] = { 0 };
                hsize_t countm[rankm] = { static_cast<hsize_t>(num_beam_block) };
                H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

                // Write data to the hyperslabs in the file from the hyperslabs
                // in memory.
                H5Dwrite(node1_idx_dataset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node1_idx_buf    [0]);
                H5Dwrite(node2_idx_dataset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node2_idx_buf    [0]);
                H5Dwrite(node3_idx_dataset    , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node3_idx_buf    [0]);
                H5Dwrite(bend_rigidity_dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &bend_rigidity_buf[0]);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    H5Dwrite(rest_curvature_dataset[d], H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &rest_curvature_buf[d][0]);
                }

                // Increment the block counter.
                ++block;
            }
        }
    }

    cout << "\n\n";

    // Close the input file.
    file_stream.close();

    // Cleanup HDF5 data structures.
    H5Pclose(plist);
    H5Dclose(node1_idx_dataset);
    H5Dclose(node2_idx_dataset);
    H5Dclose(node3_idx_dataset);
    H5Dclose(bend_rigidity_dataset);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        H5Dclose(rest_curvature_dataset[d]);
    }
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Gclose(beam_group_id);
    return;
}// initializeBeamData

void
initializeTargetPointData(
    const string& base_filename,
    hid_t file_id)
{
    string line_string;

    // Ensure that the file exists.
    const string target_point_filename = base_filename + ".target";
    ifstream file_stream;
    file_stream.open(target_point_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cerr << "warning: unable to open target point file " << target_point_filename << "\n\n";
        return;
    }

    cout << "processing target point data from ASCII input filename " << target_point_filename << "\n";

    // The first line in the file indicates the number of target points in the
    // input file.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line 1 of file " << target_point_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_target_point))
        {
            cerr << "error: invalid entry in input file encountered on line 1 of file " << target_point_filename << "\n";
            abort();
        }
    }

    if (num_target_point <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line 1 of file " << target_point_filename << "\n";
        abort();
    }

    // Create the appropriate group(s).
    const string target_point_group_name = "/" + base_filename + "/target_point";
    hid_t target_point_group_id = H5Gcreate1(file_id, target_point_group_name.c_str(), 0);

    // Define the file dataspace.
    static const int rankf = 1;
    hsize_t dimsf[rankf] = { static_cast<hsize_t>(num_target_point) };
    hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

    // Define the memory dataspace.
    static const int rankm = 1;
    hsize_t dimsm[rankm] = { BUFFER_SIZE };
    hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

    // Create the datasets with data compression enabled.
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    static const int rankc = 1;
    hsize_t dimsc[rankc] = { static_cast<hsize_t>(min(num_target_point,BUFFER_SIZE)) };
    H5Pset_chunk(plist, rankc, dimsc);
    H5Pset_shuffle(plist);
    H5Pset_deflate(plist, 6);

    const string target_point_node_idx_dset_name = "/" + base_filename + "/target_point/node_idx";
    hid_t node_idx_dataset = H5Dcreate1(file_id, target_point_node_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string target_point_stiffness_dset_name = "/" + base_filename + "/target_point/stiffness";
    hid_t stiffness_dataset = H5Dcreate1(file_id, target_point_stiffness_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    const string target_point_damping_dset_name = "/" + base_filename + "/target_point/damping";
    hid_t damping_dataset = H5Dcreate1(file_id, target_point_damping_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    // Each successive line indicates the vertex number and penalty spring
    // constant associated with any target points.
    vector<int> node_idx_buf(BUFFER_SIZE);
    vector<double> stiffness_buf(BUFFER_SIZE), damping_buf(BUFFER_SIZE);
    const int num_blocks = num_target_point/BUFFER_SIZE + (num_target_point%BUFFER_SIZE == 0 ? 0 : 1);
    for (int k = 0, block = 0; k < num_target_point; ++k)
    {
        if (k%10 == 0 || k == num_target_point-1)
        {
            cout << "\r" << setw(3) << static_cast<int>(100.0*static_cast<double>(k)/static_cast<double>(num_target_point-1)) << "% complete";
            cout.flush();
        }

        int node_idx;
        double stiffness, damping;
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << k+2 << " of file " << target_point_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            if (!(line_stream >> node_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << target_point_filename << "\n";
                abort();
            }
            else if ((node_idx < 0) || (node_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << target_point_filename << "\n"
                     << "       vertex index " << node_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> stiffness))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << target_point_filename << "\n";
                abort();
            }
            else if (stiffness < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << target_point_filename << "\n"
                     << "       target point spring constant is negative\n";
                abort();
            }

            if (!(line_stream >> damping))
            {
                damping = 0.0;
            }
            else if (damping < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << target_point_filename << "\n"
                     << "       target point damping coefficient is negative\n";
                abort();
            }

            // Store the values in the buffers.
            node_idx_buf [k%BUFFER_SIZE] = node_idx;
            stiffness_buf[k%BUFFER_SIZE] = stiffness;
            damping_buf  [k%BUFFER_SIZE] = damping;

            // Write out data as the buffer fills up.
            if ((k+1)%BUFFER_SIZE == 0 || k+1 == num_target_point)
            {
                // Determine whether we are writing out the last block in the
                // file.
                const bool last_block = (block == num_blocks-1);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!last_block || (last_block && k+1 == num_target_point));
#endif
                // Determine the number of items to read (always BUFFER_SIZE
                // except for the final block in the file).
                const int num_target_point_block = (last_block ? num_target_point - block*BUFFER_SIZE : BUFFER_SIZE);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(num_target_point_block > 0 && num_target_point_block <= BUFFER_SIZE);
#endif
                // Define the file hyperslab.
                hsize_t offsetf[rankf] = { static_cast<hsize_t>(block*BUFFER_SIZE) };
                hsize_t countf[rankf] = { static_cast<hsize_t>(num_target_point_block) };
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

                // Define the memory hyperslab.
                hsize_t offsetm[rankm] = { 0 };
                hsize_t countm[rankm] = { static_cast<hsize_t>(num_target_point_block) };
                H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

                // Write data to the hyperslabs in the file from the hyperslabs
                // in memory.
                H5Dwrite(node_idx_dataset , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node_idx_buf [0]);
                H5Dwrite(stiffness_dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &stiffness_buf[0]);
                H5Dwrite(damping_dataset  , H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &damping_buf  [0]);

                // Increment the block counter.
                ++block;
            }
        }
    }

    cout << "\n\n";

    // Close the input file.
    file_stream.close();

    // Cleanup HDF5 data structures.
    H5Pclose(plist);
    H5Dclose(node_idx_dataset);
    H5Dclose(stiffness_dataset);
    H5Dclose(damping_dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Gclose(target_point_group_id);
    return;
}// initializeTargetPointData

void
initializeMassData(
    const string& base_filename,
    hid_t file_id)
{
    string line_string;

    // Ensure that the file exists.
    const string mass_point_filename = base_filename + ".mass";
    ifstream file_stream;
    file_stream.open(mass_point_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cerr << "warning: unable to open boundary mass file " << mass_point_filename << "\n\n";
        return;
    }

    cout << "processing boundary mass data from ASCII input filename " << mass_point_filename << "\n";

    // The first line in the file indicates the number of massive IB points in
    // the input file.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line 1 of file " << mass_point_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_mass_point))
        {
            cerr << "error: invalid entry in input file encountered on line 1 of file " << mass_point_filename << "\n";
            abort();
        }
    }

    if (num_mass_point <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line 1 of file " << mass_point_filename << "\n";
        abort();
    }

    // Create the appropriate group(s).
    const string mass_point_group_name = "/" + base_filename + "/mass_point";
    hid_t mass_point_group_id = H5Gcreate1(file_id, mass_point_group_name.c_str(), 0);

    // Define the file dataspace.
    static const int rankf = 1;
    hsize_t dimsf[rankf] = { static_cast<hsize_t>(num_mass_point) };
    hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

    // Define the memory dataspace.
    static const int rankm = 1;
    hsize_t dimsm[rankm] = { BUFFER_SIZE };
    hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

    // Create the datasets with data compression enabled.
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    static const int rankc = 1;
    hsize_t dimsc[rankc] = { static_cast<hsize_t>(min(num_mass_point,BUFFER_SIZE)) };
    H5Pset_chunk(plist, rankc, dimsc);
    H5Pset_shuffle(plist);
    H5Pset_deflate(plist, 6);

    const string mass_point_node_idx_dset_name = "/" + base_filename + "/mass_point/node_idx";
    hid_t node_idx_dataset = H5Dcreate1(file_id, mass_point_node_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string mass_point_bdry_mass_dset_name = "/" + base_filename + "/mass_point/bdry_mass";
    hid_t bdry_mass_dataset = H5Dcreate1(file_id, mass_point_bdry_mass_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    const string mass_point_stiffness_dset_name = "/" + base_filename + "/mass_point/stiffness";
    hid_t stiffness_dataset = H5Dcreate1(file_id, mass_point_stiffness_dset_name.c_str(), H5T_NATIVE_DOUBLE, filespace, plist);

    // Each successive line indicates the vertex number, mass, and penalty
    // spring constant associated with any massive IB points.
    vector<int> node_idx_buf(BUFFER_SIZE);
    vector<double> bdry_mass_buf(BUFFER_SIZE), stiffness_buf(BUFFER_SIZE);
    const int num_blocks = num_mass_point/BUFFER_SIZE + (num_mass_point%BUFFER_SIZE == 0 ? 0 : 1);
    for (int k = 0, block = 0; k < num_mass_point; ++k)
    {
        if (k%10 == 0 || k == num_mass_point-1)
        {
            cout << "\r" << setw(3) << static_cast<int>(100.0*static_cast<double>(k)/static_cast<double>(num_mass_point-1)) << "% complete";
            cout.flush();
        }

        int node_idx;
        double bdry_mass, stiffness;
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << k+2 << " of file " << mass_point_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            if (!(line_stream >> node_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << mass_point_filename << "\n";
                abort();
            }
            else if ((node_idx < 0) || (node_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << mass_point_filename << "\n"
                     << "       vertex index " << node_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> bdry_mass))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << mass_point_filename << "\n";
                abort();
            }
            else if (bdry_mass < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << mass_point_filename << "\n"
                     << "       boundary mass is negative\n";
                abort();
            }

            if (!(line_stream >> stiffness))
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << mass_point_filename << "\n";
                abort();
            }
            else if (stiffness < 0.0)
            {
                cerr << "error: invalid entry in input file encountered on line " << k+2 << " of file " << mass_point_filename << "\n"
                     << "       boundary mass spring constant is negative\n";
                abort();
            }

            // Store the values in the buffers.
            node_idx_buf [k%BUFFER_SIZE] = node_idx;
            bdry_mass_buf[k%BUFFER_SIZE] = bdry_mass;
            stiffness_buf[k%BUFFER_SIZE] = stiffness;

            // Write out data as the buffer fills up.
            if ((k+1)%BUFFER_SIZE == 0 || k+1 == num_mass_point)
            {
                // Determine whether we are writing out the last block in the
                // file.
                const bool last_block = (block == num_blocks-1);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!last_block || (last_block && k+1 == num_mass_point));
#endif
                // Determine the number of items to read (always BUFFER_SIZE
                // except for the final block in the file).
                const int num_mass_point_block = (last_block ? num_mass_point - block*BUFFER_SIZE : BUFFER_SIZE);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(num_mass_point_block > 0 && num_mass_point_block <= BUFFER_SIZE);
#endif
                // Define the file hyperslab.
                hsize_t offsetf[rankf] = { static_cast<hsize_t>(block*BUFFER_SIZE) };
                hsize_t countf[rankf] = { static_cast<hsize_t>(num_mass_point_block) };
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

                // Define the memory hyperslab.
                hsize_t offsetm[rankm] = { 0 };
                hsize_t countm[rankm] = { static_cast<hsize_t>(num_mass_point_block) };
                H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

                // Write data to the hyperslabs in the file from the hyperslabs
                // in memory.
                H5Dwrite(node_idx_dataset , H5T_NATIVE_INT   , memspace, filespace, H5P_DEFAULT, &node_idx_buf [0]);
                H5Dwrite(bdry_mass_dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &bdry_mass_buf[0]);
                H5Dwrite(stiffness_dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &stiffness_buf[0]);

                // Increment the block counter.
                ++block;
            }
        }
    }

    cout << "\n\n";

    // Close the input file.
    file_stream.close();

    // Cleanup HDF5 data structures.
    H5Pclose(plist);
    H5Dclose(node_idx_dataset);
    H5Dclose(bdry_mass_dataset);
    H5Dclose(stiffness_dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Gclose(mass_point_group_id);
    return;
}// initializeMassData

void
initializeInstrumentationData(
    const string& base_filename,
    hid_t file_id)
{
    string line_string;

    // Ensure that the file exists.
    const string inst_filename = base_filename + ".inst";
    ifstream file_stream;
    file_stream.open(inst_filename.c_str(), ios::in);
    if (!file_stream.is_open())
    {
        cerr << "warning: unable to open instrumentation file " << inst_filename << "\n\n";
        return;
    }

    cout << "processing instrumentation data from ASCII input filename " << inst_filename << "\n";

    // The first line in the file indicates the number of instruments in the
    // input file.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line 1 of file " << inst_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_inst))
        {
            cerr << "error: invalid entry in input file encountered on line 1 of file " << inst_filename << "\n";
            abort();
        }
    }

    if (num_inst <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line 1 of file " << inst_filename << "\n";
        abort();
    }

    // The next several lines in the file indicate the names of the instruments
    // in the input file.
    vector<string> instrument_names;
    for (int m = 0; m < num_inst; ++m)
    {
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << m+2 << " of file " << inst_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);

            // Trim leading whitespace.
            string::size_type notwhite = line_string.find_first_not_of(" \t\n");
            line_string.erase(0,notwhite);
            notwhite = line_string.find_last_not_of(" \t\n");
            line_string.erase(notwhite+1);

            instrument_names.push_back(line_string);
        }
    }

    // The next line in the file indicates the number of instrumented IB points
    // in the input file.
    if (!getline(file_stream, line_string))
    {
        cerr << "error: premature end to input file encountered before line " << num_inst+2 << " of file " << inst_filename << "\n";
        abort();
    }
    else
    {
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        if (!(line_stream >> num_inst_point))
        {
            cerr << "error: invalid entry in input file encountered on line " << num_inst+2 << " of file " << inst_filename << "\n";
            abort();
        }
    }

    if (num_inst_point <= 0)
    {
        cerr << "error: invalid entry in input file encountered on line" << num_inst+2 << " of file " << inst_filename << "\n";
        abort();
    }

    // Create the appropriate group(s).
    const string instrumentation_group_name = "/" + base_filename + "/instrumentation";
    hid_t instrumentation_group_id = H5Gcreate1(file_id, instrumentation_group_name.c_str(), 0);

    // Store the instrument names.
    const string instrumentation_num_inst_dset_name = "/" + base_filename + "/instrumentation/num_inst";
    static const int rankn = 1;
    hsize_t dimsn[rankn] = { 1 };
    H5LTmake_dataset_int(file_id, instrumentation_num_inst_dset_name.c_str(), rankn, dimsn, &num_inst);
    for (int k = 0; k < num_inst; ++k)
    {
        ostringstream num_stream;
        num_stream << k;
        const string instrumentation_name_dset_name = "/" + base_filename + "/instrumentation/name_" + num_stream.str();
        H5LTmake_dataset_string(file_id, instrumentation_name_dset_name.c_str(), instrument_names[k].c_str());
    }

    // Define the file dataspace.
    static const int rankf = 1;
    hsize_t dimsf[rankf] = { static_cast<hsize_t>(num_inst_point) };
    hid_t filespace = H5Screate_simple(rankf, dimsf, NULL);

    // Define the memory dataspace.
    static const int rankm = 1;
    hsize_t dimsm[rankm] = { BUFFER_SIZE };
    hid_t memspace = H5Screate_simple(rankm, dimsm, NULL);

    // Create the datasets with data compression enabled.
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    static const int rankc = 1;
    hsize_t dimsc[rankc] = { static_cast<hsize_t>(min(num_inst_point,BUFFER_SIZE)) };
    H5Pset_chunk(plist, rankc, dimsc);
    H5Pset_shuffle(plist);
    H5Pset_deflate(plist, 6);

    const string instrumentation_node_idx_dset_name = "/" + base_filename + "/instrumentation/node_idx";
    hid_t node_idx_dataset = H5Dcreate1(file_id, instrumentation_node_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string instrumentation_meter_idx_dset_name = "/" + base_filename + "/instrumentation/meter_idx";
    hid_t meter_idx_dataset = H5Dcreate1(file_id, instrumentation_meter_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    const string instrumentation_meter_node_idx_dset_name = "/" + base_filename + "/instrumentation/meter_node_idx";
    hid_t meter_node_idx_dataset = H5Dcreate1(file_id, instrumentation_meter_node_idx_dset_name.c_str(), H5T_NATIVE_INT, filespace, plist);

    // Each successive line indicates the vertex number, meter number, and meter
    // node indices of each of the instrumented IB points in the input file.
    vector<bool> encountered_instrument_idx;
    map<int,vector<bool> > encountered_node_idx;
    vector<int> node_idx_buf(BUFFER_SIZE), meter_idx_buf(BUFFER_SIZE), meter_node_idx_buf(BUFFER_SIZE);
    const int num_blocks = num_inst_point/BUFFER_SIZE + (num_inst_point%BUFFER_SIZE == 0 ? 0 : 1);
    for (int k = 0, block = 0; k < num_inst_point; ++k)
    {
        if (k%10 == 0 || k == num_inst_point-1)
        {
            cout << "\r" << setw(3) << static_cast<int>(100.0*static_cast<double>(k)/static_cast<double>(num_inst_point-1)) << "% complete";
            cout.flush();
        }

        int node_idx, meter_idx, meter_node_idx;
        if (!getline(file_stream, line_string))
        {
            cerr << "error: premature end to input file encountered before line " << num_inst+k+3 << " of file " << inst_filename << "\n";
            abort();
        }
        else
        {
            line_string = discard_comments(line_string);
            istringstream line_stream(line_string);
            if (!(line_stream >> node_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << "\n";
                abort();
            }
            else if ((node_idx < 0) || (node_idx >= num_vertex))
            {
                cerr << "error: invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << "\n"
                     << "       vertex index " << node_idx << " is out of range\n";
                abort();
            }

            if (!(line_stream >> meter_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << "\n";
                abort();
            }
            else if (meter_idx < 0 || meter_idx >= num_inst)
            {
                cerr << "error: invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << "\n"
                     << "       meter index " << meter_idx << " is out of range\n";
                abort();
            }

            if (meter_idx >= static_cast<int>(encountered_instrument_idx.size()))
            {
                encountered_instrument_idx.resize(meter_idx+1,false);
            }
            encountered_instrument_idx[meter_idx] = true;

            if (!(line_stream >> meter_node_idx))
            {
                cerr << "error: invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << "\n";
                abort();
            }
            else if (meter_node_idx < 0)
            {
                cerr << "error: invalid entry in input file encountered on line " << num_inst+k+3 << " of file " << inst_filename << "\n"
                     << "       meter node index is negative\n";
                abort();
            }

            if (meter_node_idx >= static_cast<int>(encountered_node_idx[meter_idx].size()))
            {
                encountered_node_idx[meter_idx].resize(meter_node_idx+1,false);
            }
            encountered_node_idx[meter_idx][meter_node_idx] = true;

            // Store the values in the buffers.
            node_idx_buf      [k%BUFFER_SIZE] = node_idx;
            meter_idx_buf     [k%BUFFER_SIZE] = meter_idx;
            meter_node_idx_buf[k%BUFFER_SIZE] = meter_node_idx;

            // Write out data as the buffer fills up.
            if ((k+1)%BUFFER_SIZE == 0 || k+1 == num_inst_point)
            {
                // Determine whether we are writing out the last block in the
                // file.
                const bool last_block = (block == num_blocks-1);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(!last_block || (last_block && k+1 == num_inst_point));
#endif
                // Determine the number of items to read (always BUFFER_SIZE
                // except for the final block in the file).
                const int num_inst_point_block = (last_block ? num_inst_point - block*BUFFER_SIZE : BUFFER_SIZE);
#ifdef DEBUG_CHECK_ASSERTIONS
                assert(num_inst_point_block > 0 && num_inst_point_block <= BUFFER_SIZE);
#endif
                // Define the file hyperslab.
                hsize_t offsetf[rankf] = { static_cast<hsize_t>(block*BUFFER_SIZE) };
                hsize_t countf[rankf] = { static_cast<hsize_t>(num_inst_point_block) };
                H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetf, NULL, countf, NULL);

                // Define the memory hyperslab.
                hsize_t offsetm[rankm] = { 0 };
                hsize_t countm[rankm] = { static_cast<hsize_t>(num_inst_point_block) };
                H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offsetm, NULL, countm, NULL);

                // Write data to the hyperslabs in the file from the hyperslabs
                // in memory.
                H5Dwrite(node_idx_dataset      , H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, &node_idx_buf      [0]);
                H5Dwrite(meter_idx_dataset     , H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, &meter_idx_buf     [0]);
                H5Dwrite(meter_node_idx_dataset, H5T_NATIVE_INT, memspace, filespace, H5P_DEFAULT, &meter_node_idx_buf[0]);

                // Increment the block counter.
                ++block;
            }
        }
    }

    // Ensure that a complete range of instrument indices were found in the
    // input file.
    for (vector<bool>::iterator meter_it = encountered_instrument_idx.begin(); meter_it != encountered_instrument_idx.end(); ++meter_it)
    {
        const int meter_idx = distance(encountered_instrument_idx.begin(),meter_it);
        if ((*meter_it) == false)
        {
            cerr << "error: instrument index " << meter_idx << " not found in input file " << inst_filename << "\n";
            abort();
        }

        vector<bool>& meter_node_idxs = encountered_node_idx[meter_idx];
        for (vector<bool>::iterator node_it = meter_node_idxs.begin(); node_it != meter_node_idxs.end(); ++node_it)
        {
            const int node_idx = distance(meter_node_idxs.begin(),node_it);
            if ((*node_it) == false)
            {
                cerr << "error: node index " << node_idx << " associated with meter index " << meter_idx << " not found in input file " << inst_filename << "\n";
                abort();
            }
        }
    }

    if (static_cast<int>(encountered_instrument_idx.size()) != num_inst)
    {
        cerr << "error: not all anticipated instrument indices were found in input file " << inst_filename
             << "       expected to find " << num_inst << " distinct meter indices in input file\n";
        abort();
    }

    cout << "\n\n";

    // Close the input file.
    file_stream.close();

    // Cleanup HDF5 data structures.
    H5Pclose(plist);
    H5Dclose(node_idx_dataset);
    H5Dclose(meter_idx_dataset);
    H5Dclose(meter_node_idx_dataset);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Gclose(instrumentation_group_id);
    return;
}// initializeInstrumentationData

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
