// Filename: marker_flooder.C
// Last modified: <10.Oct.2007 17:53:27 griffith@box221.cims.nyu.edu>
// Created on 10 Oct 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// BLITZ++ INCLUDES
#include <blitz/array.h>
#include <blitz/tinyvec.h>

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

// Mesh spacing parameters.
static const int N = 64;                 // number of Cartesian grid points in each coordinate direction
static const double L = 10.0;            // length of computational domain
static const double dx = L/double(N);    // Cartesian grid spacing

inline TinyVector<int,NDIM>
get_index(
    const TinyVector<double,NDIM>& X)
{
    TinyVector<int,NDIM> i;
    for (int d = 0; d < NDIM; ++d)
    {
        i[d] = int(floor(X[d]/dx));
    }
    return i;
}// get_index

inline TinyVector<double,NDIM>
get_posn(
    const TinyVector<int,NDIM>& i)
{
    TinyVector<double,NDIM> X;
    for (int d = 0; d < NDIM; ++d)
    {
        X[d] = (double(i[d])+0.5)*dx;
    }
    return X;
}// get_posn

inline bool
valid_index(
    const TinyVector<int,NDIM>& i)
{
    for (int d = 0; d < NDIM; ++d)
    {
        if (i[d] < 0 || i[d] >= N) return false;
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
    for (int d = 0; d < NDIM; ++d)
    {
        extents[d] = N;
    }
    Array<bool,NDIM> occupied(extents, ColumnMajorArray<NDIM>());
    occupied = false;

    // Don't allow any markers within two meshwidths of the domain boundary.
    for (int k = 0; k < N; ++k)
    {
        for (int j = 0; j < N; ++j)
        {
            occupied(  0,j,k) = true;
            occupied(  1,j,k) = true;

            occupied(N-2,j,k) = true;
            occupied(N-1,j,k) = true;
        }
    }

    for (int k = 0; k < N; ++k)
    {
        for (int i = 0; i < N; ++i)
        {
            occupied(i,  0,k) = true;
            occupied(i,  1,k) = true;

            occupied(i,N-2,k) = true;
            occupied(i,N-1,k) = true;
        }
    }

    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            occupied(i,j,  0) = true;
            occupied(i,j,  1) = true;

            occupied(i,j,N-2) = true;
            occupied(i,j,N-1) = true;
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
        for (int d = 0; d < NDIM; ++d)
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
    X_seed[2] = 8.0;
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
        for (set<TinyVector<int,NDIM>,TinyComp<int> >::const_iterator cit = active_set.begin();
             cit != active_set.end(); ++cit)
        {
            const TinyVector<int,NDIM>& i = (*cit);
            for (int d = 0; d < NDIM; ++d)
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
    for (set<TinyVector<int,NDIM>,TinyComp<int> >::const_iterator cit = all_markers.begin();
         cit != all_markers.end(); ++cit)
    {
        const TinyVector<int,NDIM>& i = (*cit);
        const TinyVector<double,NDIM> X = get_posn(i);
        for (int d = 0; d < NDIM; ++d)
        {
            output_file_stream << setw(7) << fixed << setprecision(3) << X[d] << (d == NDIM-1 ? "" : " ");
        }
        output_file_stream << "\n";
    }
    return 0;
}// main
