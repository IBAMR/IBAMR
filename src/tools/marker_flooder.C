// Filename: marker_flooder.C
// Created on 10 Oct 2007 by Boyce Griffith
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

// BLITZ++ INCLUDES
#include <blitz/array.h>
#include <blitz/tinyvec2.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace blitz;
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

// number of Cartesian grid points in each coordinate direction
static const blitz::TinyVector<int,NDIM> N(64 , 64 , 96);

// length of computational domain
static const blitz::TinyVector<double,NDIM> L(10.0 , 10.0 , 15.0);

// Cartesian grid spacing
static const blitz::TinyVector<double,NDIM> dx(L[0]/static_cast<double>(N[0]) , L[1]/static_cast<double>(N[1]) , L[2]/double(N[2]));

inline TinyVector<int,NDIM>
get_index(
    const TinyVector<double,NDIM>& X)
{
    TinyVector<int,NDIM> i;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        i[d] = static_cast<int>(floor(X[d]/dx[d]));
    }
    return i;
}// get_index

inline TinyVector<double,NDIM>
get_posn(
    const TinyVector<int,NDIM>& i)
{
    TinyVector<double,NDIM> X;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        X[d] = (static_cast<double>(i[d])+0.5)*dx[d];
    }
    return X;
}// get_posn

inline bool
valid_index(
    const TinyVector<int,NDIM>& i)
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (i[d] < 0 || i[d] >= N[d]) return false;
    }
    return true;
}// valid_point

template <typename T>
struct TinyComp
    : std::binary_function<TinyVector<T,NDIM>,TinyVector<T,NDIM>,bool>
{
    inline bool
    operator()(
        const TinyVector<T,NDIM>& lhs,
        const TinyVector<T,NDIM>& rhs) const
        {
            for (int d = NDIM-1; d >= 0; --d)
            {
                if (lhs(d) < rhs(d)) return true;
                if (lhs(d) > rhs(d)) return false;
            }
            return false;
        }// operator()
};

int
main(
    int argc,
    char* argv[])
{
    if (argc != 3)
    {
        cout << argv[0] << ": a tool to flood enclosed regions with marker points\n"
             << "USAGE: " << argv[0] << " <input filename> <output filename>\n";
        return -1;
    }

    const string input_filename = argv[1];
    const string output_filename = argv[2];

    // Construct an array of bools to keep track of which cells are occupied in
    // the grid.
    TinyVector<int,NDIM> extents;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        extents[d] = N[d];
    }
    Array<bool,NDIM> occupied(extents, ColumnMajorArray<NDIM>());
    occupied = false;

    // Don't allow any markers within two meshwidths of the domain boundary.
    for (int k = 0; k < N[2]; ++k)
    {
        for (int j = 0; j < N[1]; ++j)
        {
            occupied(  0,j,k) = true;
            occupied(  1,j,k) = true;

            occupied(N[0]-2,j,k) = true;
            occupied(N[0]-1,j,k) = true;
        }
    }

    for (int k = 0; k < N[2]; ++k)
    {
        for (int i = 0; i < N[0]; ++i)
        {
            occupied(i,  0,k) = true;
            occupied(i,  1,k) = true;

            occupied(i,N[1]-2,k) = true;
            occupied(i,N[1]-1,k) = true;
        }
    }

    for (int j = 0; j < N[1]; ++j)
    {
        for (int i = 0; i < N[0]; ++i)
        {
            occupied(i,j,  0) = true;
            occupied(i,j,  1) = true;

            occupied(i,j,N[2]-2) = true;
            occupied(i,j,N[2]-1) = true;
        }
    }

    // Read in the input points and tag their positions.
    cout << "processing input file -";
    ifstream input_file_stream;
    input_file_stream.open(input_filename.c_str(), ios::in);
    int line_counter = 0;
    while (!input_file_stream.eof())
    {
        if (line_counter % 100000 == 0)
        {
            cout << "\rprocessing input file -";
            cout.flush();
        }
        else if (line_counter % 100000 == 25000)
        {
            cout << "\rprocessing input file \\";
            cout.flush();
        }
        else if (line_counter % 100000 == 50000)
        {
            cout << "\rprocessing input file |";
            cout.flush();
        }
        else if (line_counter % 100000 == 75000)
        {
            cout << "\rprocessing input file /";
            cout.flush();
        }

        string line_string;
        getline(input_file_stream, line_string);
        line_string = discard_comments(line_string);
        istringstream line_stream(line_string);
        TinyVector<double,NDIM> X;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            line_stream >> X[d];
        }
        const TinyVector<int,NDIM> i = get_index(X);
        if (valid_index(i)) occupied(i) = true;
        ++line_counter;
    }
    cout << "\n";

    // Initialize the active point set using a single seed point.
    cout << "generating markers\n";
    TinyVector<double,NDIM> X_seed;
    X_seed[0] = 5.0;
    X_seed[1] = 5.0;
    X_seed[2] = 2.0;
    TinyVector<int,NDIM> i_seed = get_index(X_seed);

    assert(valid_index(i_seed));
    assert(!occupied(i_seed));

    set<TinyVector<int,NDIM>,TinyComp<int> > active_set, all_markers;
    active_set.insert(i_seed);
    all_markers.insert(i_seed);
    occupied(i_seed) = true;
    while (!active_set.empty())
    {
        set<TinyVector<int,NDIM>,TinyComp<int> > new_active_set;
        for (set<TinyVector<int,NDIM>,TinyComp<int> >::const_iterator cit = active_set.begin(); cit != active_set.end(); ++cit)
        {
            const TinyVector<int,NDIM>& i = (*cit);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                TinyVector<int,NDIM> i_lo = i;
                i_lo(d) -= 1;
                if (valid_index(i_lo) && !occupied(i_lo))
                {
                    new_active_set.insert(i_lo);
                    all_markers.insert(i_lo);
                    occupied(i_lo) = true;
                }

                TinyVector<int,NDIM> i_hi = i;
                i_hi(d) += 1;
                if (valid_index(i_hi) && !occupied(i_hi))
                {
                    new_active_set.insert(i_hi);
                    all_markers.insert(i_hi);
                    occupied(i_hi) = true;
                }
            }
        }
        active_set = new_active_set;
    }

    // Output the results.
    cout << "outputting markers\n";
    ofstream output_file_stream;
    output_file_stream.open(output_filename.c_str(),ios::out);
    for (set<TinyVector<int,NDIM>,TinyComp<int> >::const_iterator cit = all_markers.begin(); cit != all_markers.end(); ++cit)
    {
        const TinyVector<int,NDIM>& i = (*cit);
        const TinyVector<double,NDIM> X = get_posn(i);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            output_file_stream << setw(7) << fixed << setprecision(3) << X[d] << (d == NDIM-1 ? "" : " ");
        }
        output_file_stream << "\n";
    }
    return 0;
}// main
