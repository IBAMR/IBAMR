// Filename: m3D_hdf5_fiber_converter.C
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

#include <blitz/array.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace std;

/////////////////////////////// STATIC ///////////////////////////////////////

static const int M3D_NFG_MAX = 999;  // Maximum number of fibers per group.

void
build_local_cart_block(
    ostream& os,
    const std::vector<float>& X,
    const vector<int>& nelem,
    const vector<int>& periodic,
    const int& fiber_offset,
    const int& group_offset,
    const int& layer_number)
{
    const int nelem_tot = nelem[0]*nelem[1]*nelem[2];

    int group_counter = 0;
    if ((nelem[0] == 1 && nelem[1] == 1) || (nelem[0] == 1 && nelem[2] == 1) || (nelem[1] == 1 && nelem[2] == 1))
    {
        // Output a single fiber.
        const int fiber_number = fiber_offset+1;

        os << setw(7) << fiber_number << " " << setw(7) << nelem_tot << " = FIBER POINTS\n";
        for (int k = 0; k < nelem_tot; ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                os << setw(7) << fixed << setprecision(3) << X[NDIM*k+d] << " ";
            }
            os << setw(4) << k+1 << " " << setw(4) << fiber_number << " " << setw(4) << group_offset+group_counter+1 << " " << setw(4) << layer_number << "\n";
        }
    }
    else if ((nelem[0] == 1) || (nelem[1] == 1) || (nelem[2] == 1))
    {
        // Output a 2D sheet of fibers.
        int fiber_counter = 0;
        for (unsigned int d0 = 0; d0 < NDIM; ++d0)
        {
            if (nelem[d0] > 1)
            {
                // Find the other nontrivial dimension.
                for (unsigned int d1 = 0; d1 < NDIM; ++d1)
                {
                    if (d1 != d0 && nelem[d1] > 1)
                    {
                        for (int j = 0; j < nelem[d0]; ++j)
                        {
                            const int fiber_number = fiber_offset+fiber_counter+1;
                            fiber_counter += 1;

                            os << setw(7) << fiber_number << " " << setw(7) << nelem[d1] + (periodic[d1] ? 1 : 0) << " = FIBER POINTS\n";
                            for (int k = 0; k < nelem[d1] + (periodic[d1] ? 1 : 0); ++k)
                            {
                                bool end_of_fiber = false;
                                if (periodic[d1] && k == nelem[d1])
                                {
                                    k = 0;
                                    end_of_fiber = true;
                                }

                                blitz::TinyVector<int,NDIM> idx(0 , 0 , 0);
                                idx[d0] = j;
                                idx[d1] = k;
                                const int offset = idx[0] + idx[1]*nelem[0] + idx[2]*nelem[0]*nelem[1];

                                for (unsigned int d = 0; d < NDIM; ++d)
                                {
                                    os << setw(7) << fixed << setprecision(3) << X[NDIM*offset+d] << " ";
                                }
                                os << setw(4) << k+1 << " " << setw(4) << fiber_number << " " << setw(4) << group_offset+group_counter+1 << " " << setw(4) << layer_number << "\n";

                                if (end_of_fiber) break;
                            }
                        }
                        group_counter += 1;
                    }
                }
            }
        }
    }
    else
    {
        // Output a 3D volume of fibers.
        int fiber_counter = 0;

        // Loop over all pairs of dimensions.
        for (unsigned int d0 = 0; d0 < NDIM; ++d0)
        {
            const int d1 = (d0+1)%NDIM;
            const int d2 = (d1+1)%NDIM;

            // myocardial3D cannot cleanly handle headers specifying groups with
            // more than 999 fibers.
            int nfibers_per_group = nelem[d0]*nelem[d1];
            while(nfibers_per_group > M3D_NFG_MAX)
            {
                nfibers_per_group /= 2;
            }

            for (int i = 0; i < nelem[d0]; ++i)
            {
                for (int j = 0; j < nelem[d1]; ++j)
                {
                    const int fiber_number = fiber_offset+fiber_counter+1;
                    fiber_counter += 1;
                    os << setw(7) << fiber_number << " " << setw(7) << nelem[d2] + (periodic[d2] ? 1 : 0) << " = FIBER POINTS\n";
                    for (int k = 0; k < nelem[d2] + (periodic[d2] ? 1 : 0); ++k)
                    {
                        bool end_of_fiber = false;
                        if (periodic[d2] && k == nelem[d2])
                        {
                            k = 0;
                            end_of_fiber = true;
                        }

                        blitz::TinyVector<int,NDIM> idx(0 , 0 , 0);
                        idx[d0] = i;
                        idx[d1] = j;
                        idx[d2] = k;
                        const int offset = idx[0] + idx[1]*nelem[0] + idx[2]*nelem[0]*nelem[1];

                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            os << setw(7) << fixed << setprecision(3) << X[NDIM*offset+d] << " ";
                        }
                        os << setw(4) << k+1 << " " << setw(4) << fiber_number << " " << setw(4) << group_offset+group_counter+1 << " " << setw(4) << layer_number << "\n";

                        if (end_of_fiber) break;
                    }
                }
            }
            group_counter += 1;
        }
    }
    return;
}// build_local_cart_block

int
main(
    int argc,
    char* argv[])
{
    if (argc != 3)
    {
        cout << argv[0] << ": a tool to convert IBAMR HDF5 fiber files to myocardial3D fiber input files" << "\n"
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

    hid_t fiber_file_id = H5Fopen(input_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fiber_file_id < 0)
    {
        cerr << "error: Unable to open input file: " << input_filename << endl;
        return -1;
    }

    hid_t fiber_group_id = H5Gopen1(fiber_file_id, "/fibers");

    int num_local_fibers, fiber_offset, num_fibers;
    H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_fibers", &num_local_fibers);
    H5LTget_attribute_int(fiber_file_id, "/fibers", "fiber_offset", &fiber_offset);
    H5LTget_attribute_int(fiber_file_id, "/fibers", "num_fibers", &num_fibers);

    int num_local_groups, group_offset, num_groups;
    H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_groups", &num_local_groups);
    H5LTget_attribute_int(fiber_file_id, "/fibers", "group_offset", &group_offset);
    H5LTget_attribute_int(fiber_file_id, "/fibers", "num_groups", &num_groups);

    int num_local_layers, layer_offset, num_layers;
    H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_layers", &num_local_layers);
    H5LTget_attribute_int(fiber_file_id, "/fibers", "layer_offset", &layer_offset);
    H5LTget_attribute_int(fiber_file_id, "/fibers", "num_layers", &num_layers);

    int local_fiber_counter = 0;
    int local_group_counter = 0;
    for (int local_layer_counter = 0; local_layer_counter < num_local_layers; ++local_layer_counter)
    {
        ostringstream dset_name_stream;
        dset_name_stream << "/fibers/layer_" << setw(4) << setfill('0') << layer_offset+local_layer_counter;
        const std::string dset_name = dset_name_stream.str();

        vector<int> nelem(NDIM);
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "nelem", &nelem[0]);
        vector<int> periodic(NDIM);
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "periodic", &periodic[0]);
        int fiber_number;
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "fiber_number", &fiber_number);
        int group_number;
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "group_number", &group_number);
        int layer_number;
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "layer_number", &layer_number);
        int nfibers;
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "nfibers", &nfibers);
        int ngroups;
        H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "ngroups", &ngroups);
        char layer_name[256];
        H5LTget_attribute_string(fiber_file_id, dset_name.c_str(), "layer_name", layer_name);

        std::vector<float> X(NDIM*nelem[0]*nelem[1]*nelem[2]);
        H5LTread_dataset_float(fiber_file_id, dset_name.c_str(), &X[0]);

        build_local_cart_block(output_fstream, X, nelem, periodic, fiber_offset, group_offset, layer_number);

        // Advance the counters.
        local_fiber_counter += nfibers;
        local_group_counter += ngroups;
    }

    H5Gclose(fiber_group_id);
    H5Fclose(fiber_file_id);

    output_fstream.close();
    return -1;
}// main
