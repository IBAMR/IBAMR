// Filename: spring2hdf5.C
// Last modified: <31.May.2007 14:05:21 griffith@box221.cims.nyu.edu>
// Created on 30 May 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <fstream>
#include <iostream>
#include <map>
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

typedef std::pair<int,int> Edge;
struct EdgeComp
    : public std::binary_function<Edge,Edge,bool>
{
    inline bool
    operator()(
        const Edge& e1,
        const Edge& e2) const
        {
            return (e1.first < e2.first) || (e1.first == e2.first && e1.second < e2.second);
        }
};

int
main(
    int argc,
    char* argv[])
{
    if (argc != 3)
    {
        cout << argv[0] << ": a tool to convert IBAMR spring files from ASCII to HDF5" << "\n"
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

    // The spring information
    std::multimap<int,Edge> spring_edge_map;
    std::map<Edge,double,EdgeComp> spring_stiffness, spring_rest_length;
    std::map<Edge,int,EdgeComp> spring_force_fcn_idx;

    // Process the input file.
    std::ifstream file_stream;
    std::string line_string;
    file_stream.open(input_filename.c_str(), std::ios::in);
    if (!file_stream.is_open())
    {
        cout << "error: Unable to open input file " << input_filename << endl;
        return -1;
    }

    // The first line in the file indicates the number of edges in the input
    // file.
    int num_edges;
    if (!std::getline(file_stream, line_string))
    {
        cout << "error: Premature end to input file encountered before line 1 of file " << input_filename << endl;
    }
    else
    {
        line_string = discard_comments(line_string);
        std::istringstream line_stream(line_string);
        if (!(line_stream >> num_edges))
        {
            cout << "error: Invalid entry in input file encountered on line 1 of file " << input_filename << endl;
        }
    }

    if (num_edges <= 0)
    {
        cout << "error: Invalid entry in input file encountered on line 1 of file " << input_filename << endl;
    }

    // Each successive line provides the connectivity and material parameter
    // information for each spring in the structure.
    for (int k = 0; k < num_edges; ++k)
    {
        Edge e;
        double kappa, length;
        int force_fcn_idx;
        if (!std::getline(file_stream, line_string))
        {
            cout << "error: Premature end to input file encountered before line " << k+2 << " of file " << input_filename << endl;
        }
        else
        {
            line_string = discard_comments(line_string);
            std::istringstream line_stream(line_string);
            if (!(line_stream >> e.first))
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl;
            }
            else if (e.first < 0)
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl
                     << "  vertex index " << e.first << " is out of range" << endl;
            }

            if (!(line_stream >> e.second))
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl;
            }
            else if (e.second < 0)
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl
                     << "  vertex index " << e.second << " is out of range" << endl;
            }

            if (!(line_stream >> kappa))
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl;
            }
            else if (kappa < 0.0)
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl
                     << "  spring constant is negative" << endl;
            }

            if (!(line_stream >> length))
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl;
            }
            else if (length < 0.0)
            {
                cout << "error: Invalid entry in input file encountered on line " << k+2 << " of file " << input_filename << endl
                     << "  spring resting length is negative" << endl;
            }

            if (!(line_stream >> force_fcn_idx))
            {
                force_fcn_idx = 0;  // default force function specification.
            }
        }

        // Always place the lower index first.
        if (e.first > e.second)
        {
            std::swap<int>(e.first, e.second);
        }

        // Check to see if the edge has already been inserted in the edge map.
        bool duplicate_edge = false;
        for (std::multimap<int,Edge>::const_iterator it = spring_edge_map.lower_bound(e.first);
             it != spring_edge_map.upper_bound(e.first); ++it)
        {
            const Edge& other_e = (*it).second;
            if (e.first  == other_e.first &&
                e.second == other_e.second)
            {
                // This is a duplicate edge and should not be inserted into the
                // edge map.
                duplicate_edge = true;

                // Ensure that the stiffness and rest length information is
                // consistent.
                if ((*spring_stiffness.find(e)).second != kappa ||
                    (*spring_rest_length.find(e)).second != length ||
                    (*spring_force_fcn_idx.find(e)).second != force_fcn_idx)
                {
                    cout << "error: Inconsistent duplicate edges in input file encountered on line " << k+2 << " of file " << input_filename << endl
                         << "  first vertex = " << e.first << " second vertex = " << e.second << endl
                         << "  original spring constant = " << (*spring_stiffness.find(e)).second << endl
                         << "  original resting length = " << (*spring_rest_length.find(e)).second << endl
                         << "  original force function index = " << (*spring_force_fcn_idx.find(e)).second << endl;
                }
            }
        }

        // Initialize the map data corresponding to the present edge.
        //
        // Note that in the edge map, each edge is associated with only the
        // first vertex.
        if (!duplicate_edge)
        {
            spring_edge_map.insert(std::make_pair(e.first,e));
            spring_stiffness[e] = kappa;
            spring_rest_length[e] = length;
            spring_force_fcn_idx[e] = force_fcn_idx;
        }
    }

    // Close the input file.
    file_stream.close();

    // Serialize the data.
    vector<int> e1, e2, fcn;
    vector<double> stf, rst;

    for (std::multimap<int,Edge>::const_iterator it = spring_edge_map.begin();
         it != spring_edge_map.end(); ++it)
    {
        const Edge& e = (*it).second;
        e1.push_back(e.first);
        e2.push_back(e.second);
        stf.push_back(spring_stiffness[e]);
        rst.push_back(spring_rest_length[e]);
        fcn.push_back(spring_force_fcn_idx[e]);
    }

    // Create the output file.
    hid_t file_id;
    herr_t status;

    file_id = H5Fcreate(output_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0)
    {
        cout << "error: Could not create output file: " << output_filename << endl;
        return -1;
    }

    // Create the dataspace for the datasets.
    static const unsigned size = e1.size();

    static const int rank = 1;
    hsize_t dims[rank] = { size };
    hid_t dataspace_id = H5Screate_simple(rank, dims, NULL);

    // Create the datasets with data shuffling and compression enabled.
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);

    hsize_t cdims[rank] = { size };
    status = H5Pset_chunk(plist_id, rank, cdims);
    status = H5Pset_shuffle(plist_id);
    status = H5Pset_deflate(plist_id, 6);

    hid_t dataset_e1 = H5Dcreate(file_id, "/e1", H5T_NATIVE_INT, dataspace_id, plist_id);
    hid_t dataset_e2 = H5Dcreate(file_id, "/e2", H5T_NATIVE_INT, dataspace_id, plist_id);
    hid_t dataset_stf = H5Dcreate(file_id, "/stf", H5T_NATIVE_DOUBLE, dataspace_id, plist_id);
    hid_t dataset_rst = H5Dcreate(file_id, "/rst", H5T_NATIVE_DOUBLE, dataspace_id, plist_id);
    hid_t dataset_fcn = H5Dcreate(file_id, "/fcn", H5T_NATIVE_INT, dataspace_id, plist_id);

    // Write the data to the dataset.
    status = H5Dwrite(dataset_e1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &e1[0]);
    status = H5Dwrite(dataset_e2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &e2[0]);
    status = H5Dwrite(dataset_stf, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &stf[0]);
    status = H5Dwrite(dataset_rst, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rst[0]);
    status = H5Dwrite(dataset_fcn, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &fcn[0]);

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
