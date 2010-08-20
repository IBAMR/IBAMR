// Filename: m3D_hdf5_strain_analyzer.C
// Created on 17 Feb 2009 by Boyce Griffith
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

#include <cstring>
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
#include <blitz/tinyvec-et.h>

extern "C"
{
#include <silo.h>
}

#include <sys/stat.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

using namespace blitz;
using namespace std;

/////////////////////////////// STATIC ///////////////////////////////////////

void
compute_strain(
    const int direction,
    const vector<int>& nelem,
    const vector<int>& periodic,
    const Array<float,4> X,
    const Array<float,4> X0,
    Array<float,4>& X_strain,
    Array<float,4>& lambda)
{
    vector<int> nelem_strain(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        nelem_strain[d] = nelem[d] - (d == direction ? 1 : 0) + (periodic[d] ? 1 : 0);
    }

    X_strain.resize(NDIM,nelem_strain[0],nelem_strain[1],nelem_strain[2]);
    lambda  .resize(1   ,nelem_strain[0],nelem_strain[1],nelem_strain[2]);

    for (int k = 0; k < nelem_strain[2]; ++k)
    {
        for (int j = 0; j < nelem_strain[1]; ++j)
        {
            for (int i = 0; i < nelem_strain[0]; ++i)
            {
                TinyVector<int,NDIM> idx_l(i,j,k);
                TinyVector<int,NDIM> idx_r(idx_l);
                idx_r(direction) += 1;
                if (idx_r(direction) >= nelem[direction])
                {
                    if (periodic[direction] == 1)
                    {
                        idx_r(direction) = 0;
                    }
                    else
                    {
                        assert(false);
                    }
                }

                TinyVector<float,NDIM> x_l, x_r, x0_l, x0_r;
                for (int d = 0; d < NDIM; ++d)
                {
                    x_l(d) = X(d,idx_l(0),idx_l(1),idx_l(2));
                    x_r(d) = X(d,idx_r(0),idx_r(1),idx_r(2));
                    x0_l(d) = X0(d,idx_l(0),idx_l(1),idx_l(2));
                    x0_r(d) = X0(d,idx_r(0),idx_r(1),idx_r(2));
                }

                TinyVector<float,NDIM> x_strain = 0.5*(x_l+x_r);

                TinyVector<float,NDIM> dx = x_r-x_l;
                TinyVector<float,NDIM> dx0 = x0_r-x0_l;

                const float l = sqrt(dot(dx,dx));
                const float L = sqrt(dot(dx0,dx0));

                for (int d = 0; d < NDIM; ++d)
                {
                    X_strain(d,i,j,k) = x_strain(d);
                }
                lambda(0,i,j,k) = l/L;
            }
        }
    }
    return;
}// compute_strain

void
buildCurvBlock(
    DBfile *dbfile,
    const string &dirname,
    const string &meshname,
    const vector<int> &periodic,
    const Array<float,4> X,
    const vector<string> &varnames=vector<string>(0),
    const vector<Array<float,4> > &varvals=vector<Array<float,4> >(0))
{
    // For Silo, a co-dimension 1 surface is "degenerate".
    TinyVector<bool,NDIM> degenerate;
    TinyVector<int ,NDIM> nelem;

    for (int d = 0; d < NDIM; ++d)
    {
        degenerate(d) = X.extent(d+1) == 1;
        nelem(d) = (degenerate(d)) ? 2 : X.extent(d+1);
        if (periodic[d]) nelem(d) += 1;
    }

    // Reformat the data as required by Silo.
    const int ntot = product(nelem);

    vector<float> X_flat(NDIM*ntot);

    assert(varnames.size() == varvals.size());
    const int nvars = varnames.size();
    vector<vector<float> > varvals_flat(NDIM*ntot);
    for (int l = 0; l < nvars; ++l)
    {
        const int vardepth = varvals[l].extent(0);
        varvals_flat[l].resize(vardepth*ntot);
    }

    for (int k = 0, offset = 0; k < nelem(2); ++k)
    {
        const int k_X = (degenerate(2)) ? 0 : periodic[2] && k == nelem(2)-1 ? 0 : k;
        for (int j = 0; j < nelem(1); ++j)
        {
            const int j_X = (degenerate(1)) ? 0 : periodic[1] && j == nelem(1)-1 ? 0 : j;
            for (int i = 0; i < nelem(0); ++i, ++offset)
            {
                const int i_X = (degenerate(0)) ? 0 : periodic[0] && i == nelem(0)-1 ? 0 : i;
                for (int d = 0; d < NDIM; ++d)
                {
                    X_flat[d*ntot+offset] = X(d,i_X,j_X,k_X);
                }
                for (int l = 0; l < nvars; ++l)
                {
                    const int vardepth = varvals[l].extent(0);
                    for (int d = 0; d < vardepth; ++d)
                    {
                        varvals_flat[l][d*ntot+offset] = varvals[l](d,i_X,j_X,k_X);
                    }
                }
            }
        }
    }

    // Set the working directory in the Silo database.
    if (DBSetDir(dbfile, dirname.c_str()) == -1)
    {
        cerr << "buildCurvBlock():\n"
             << "  Could not set directory " << dirname << endl;
        abort();
    }

    // Write out the data.
    static const int MAX_NDIM = 3;
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(NDIM <= MAX_NDIM);
#endif
    const char* coordnames[MAX_NDIM] = { "xcoords", "ycoords", "zcoords" };
    vector<float*> coords(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        coords[d] = &X_flat[0]+d*ntot;
    }

    int ndims = NDIM;
    vector<int> dims(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        dims[d] = nelem[d];
    }

    DBPutQuadmesh(dbfile, meshname.c_str(), const_cast<char**>(coordnames), &coords[0], &dims[0], ndims, DB_FLOAT, DB_NONCOLLINEAR, NULL);
    for (int l = 0; l < nvars; ++l)
    {
        const char* varname = varnames[l].c_str();
        const int vardepth = varvals[l].extent(0);
        vector<char*> compnames(vardepth);
        for (int d = 0; d < vardepth; ++d)
        {
            ostringstream stream;
            stream << "_" << d;
            const string compname = varnames[l] + stream.str();
            compnames[d] = strdup(compname.c_str());
        }

        vector<float*> vars(vardepth);
        for (int d = 0; d < vardepth; ++d)
        {
            vars[d] = &varvals_flat[l][d*ntot];
        }

        if (vardepth == 1)
        {
            DBPutQuadvar1(dbfile, varname, meshname.c_str(), vars[0], &dims[0], ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
        }
        else
        {
            DBPutQuadvar(dbfile, varname, meshname.c_str(), vardepth, &compnames[0], &vars[0], &dims[0], ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
        }

        for (int d = 0; d < vardepth; ++d)
        {
            free(compnames[d]);
        }
    }

    // Reset the working directory in the Silo database.
    if (DBSetDir(dbfile, "..") == -1)
    {
        cerr << "buildCurvBlock():\n"
             << "  Could not return to the base directory from subdirectory " << dirname << endl;
        abort();
    }
    return;
}// buildCurvBlock

int
main(
    int argc,
    char* argv[])
{
    const int nfiles = argc-1;
    if (nfiles < 2)
    {
        cout << argv[0] << ": a tool to analyze the fiber strain in an IBAMR HDF5 fiber files" << "\n"
             << "USAGE: " << argv[0] << " <input filename 1> <input filename 2> ..." << endl;
        return -1;
    }
    char** input_filenames = &argv[1];

    // Open and process the initial configuration file.
    vector<string> layer_names;
    vector<vector<int> > nelems, periodics;
    vector<Array<float,4> > X0;

    herr_t status;
    hid_t fiber_file_id = H5Fopen(input_filenames[0], H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fiber_file_id < 0)
    {
        cerr << "error: Unable to open input file: " << input_filenames[0] << endl;
        return -1;
    }

    hid_t fiber_group_id = H5Gopen1(fiber_file_id, "/fibers");

    int num_local_fibers, fiber_offset, num_fibers;
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_fibers", &num_local_fibers);
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "fiber_offset", &fiber_offset);
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_fibers", &num_fibers);

    int num_local_groups, group_offset, num_groups;
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_groups", &num_local_groups);
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "group_offset", &group_offset);
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_groups", &num_groups);

    int num_local_layers, layer_offset, num_layers;
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_layers", &num_local_layers);
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "layer_offset", &layer_offset);
    status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_layers", &num_layers);

    int local_fiber_counter = 0;
    int local_group_counter = 0;
    for (int local_layer_counter = 0; local_layer_counter < num_local_layers; ++local_layer_counter)
    {
        ostringstream dset_name_stream;
        dset_name_stream << "/fibers/layer_" << setw(4) << setfill('0') << layer_offset+local_layer_counter;
        const string dset_name = dset_name_stream.str();

        vector<int> nelem(NDIM);
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "nelem", &nelem[0]);
        vector<int> periodic(NDIM);
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "periodic", &periodic[0]);
        int fiber_number;
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "fiber_number", &fiber_number);
        int group_number;
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "group_number", &group_number);
        int layer_number;
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "layer_number", &layer_number);
        int nfibers;
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "nfibers", &nfibers);
        int ngroups;
        status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "ngroups", &ngroups);
        char layer_name[256];
        status = H5LTget_attribute_string(fiber_file_id, dset_name.c_str(), "layer_name", layer_name);

        layer_names.push_back(layer_name);
        nelems.push_back(nelem);
        periodics.push_back(periodic);
        X0.push_back(Array<float,4>(NDIM,nelem[0],nelem[1],nelem[2],ColumnMajorArray<4>()));
        status = H5LTread_dataset_float(fiber_file_id, dset_name.c_str(), X0.back().data());

        // Advance the counters.
        local_fiber_counter += nfibers;
        local_group_counter += ngroups;
    }
    status = H5Gclose(fiber_group_id);
    status = H5Fclose(fiber_file_id);

    // Prepare to dump silo files.
    mkdir("silo",S_IRWXU|S_IRWXG|S_IRWXO);
    ofstream of_visit("silo/strain.visit",ios::out);

    // Open and process all additional input files.
    vector<vector<float> > lambda_max(X0.size(),vector<float>(NDIM,numeric_limits<float>::min()));
    vector<vector<float> > lambda_min(X0.size(),vector<float>(NDIM,numeric_limits<float>::max()));
    for (int n = 0; n < nfiles; ++n)
    {
        string silo_filename = "silo/" + string(input_filenames[n]) + ".silo";
        DBfile *dbfile;
        if ((dbfile = DBCreate(silo_filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB))
            == NULL)
        {
            cerr << "main():\n"
                 << "  Could not create DBfile named " << silo_filename << endl;
            abort();
        }

        of_visit << input_filenames[n] << ".silo" << "\n";

        herr_t status;
        hid_t fiber_file_id = H5Fopen(input_filenames[n], H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fiber_file_id < 0)
        {
            cerr << "error: Unable to open input file: " << input_filenames[n] << endl;
            return -1;
        }

        hid_t fiber_group_id = H5Gopen1(fiber_file_id, "/fibers");

        int num_local_fibers, fiber_offset, num_fibers;
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_fibers", &num_local_fibers);
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "fiber_offset", &fiber_offset);
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_fibers", &num_fibers);

        int num_local_groups, group_offset, num_groups;
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_groups", &num_local_groups);
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "group_offset", &group_offset);
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_groups", &num_groups);

        int num_local_layers, layer_offset, num_layers;
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_local_layers", &num_local_layers);
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "layer_offset", &layer_offset);
        status = H5LTget_attribute_int(fiber_file_id, "/fibers", "num_layers", &num_layers);

        int local_fiber_counter = 0;
        int local_group_counter = 0;
        for (int local_layer_counter = 0; local_layer_counter < num_local_layers; ++local_layer_counter)
        {
            ostringstream dset_name_stream;
            dset_name_stream << "/fibers/layer_" << setw(4) << setfill('0') << layer_offset+local_layer_counter;
            const string dset_name = dset_name_stream.str();

            vector<int> nelem(NDIM);
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "nelem", &nelem[0]);
            vector<int> periodic(NDIM);
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "periodic", &periodic[0]);
            int fiber_number;
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "fiber_number", &fiber_number);
            int group_number;
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "group_number", &group_number);
            int layer_number;
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "layer_number", &layer_number);
            int nfibers;
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "nfibers", &nfibers);
            int ngroups;
            status = H5LTget_attribute_int(fiber_file_id, dset_name.c_str(), "ngroups", &ngroups);
            char layer_name[256];
            status = H5LTget_attribute_string(fiber_file_id, dset_name.c_str(), "layer_name", layer_name);

            Array<float,4> X(NDIM,nelem[0],nelem[1],nelem[2],ColumnMajorArray<4>());
            status = H5LTread_dataset_float(fiber_file_id, dset_name.c_str(), X.data());

            // Compute the strains.
            for (int direction = 0; direction < NDIM; ++direction)
            {
                if (nelem[direction] > 1)
                {
                    Array<float,4> X_strain, lambda;
                    compute_strain(direction,nelem,periodic,X,X0[local_layer_counter],X_strain,lambda);

                    lambda_max[local_layer_counter][direction] = max(max(lambda),lambda_max[local_layer_counter][direction]);
                    lambda_min[local_layer_counter][direction] = min(min(lambda),lambda_min[local_layer_counter][direction]);

                    ostringstream direction_stream;
                    direction_stream << "_" << direction;
                    string meshname = string(layer_name) + direction_stream.str();
                    string dbdir = ".";
                    string varname = string(layer_name) + string("_lambda") + direction_stream.str();
                    vector<string> varnames(1,varname);
                    vector<Array<float,4> > varvals(1,lambda);
                    buildCurvBlock(dbfile, dbdir, meshname, periodic, X_strain, varnames, varvals);
                }
            }

            // Advance the counters.
            local_fiber_counter += nfibers;
            local_group_counter += ngroups;
        }
        status = H5Gclose(fiber_group_id);
        status = H5Fclose(fiber_file_id);

        DBClose(dbfile);
    }

    // Output the results.
    for (unsigned l = 0; l < X0.size(); ++l)
    {
        cout << "structure: " << layer_names[l] << "\n";
        for (int d = 0; d < NDIM; ++d)
        {
            if (nelems[l][d] > 1)
            {
                cout << "direction: " << d << "; lambda_max = " << lambda_max[l][d] << " lambda_min = " << lambda_min[l][d] << "\n";
            }
        }
    }
    return 0;
}// main
