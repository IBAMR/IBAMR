// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <sstream>

using Edge = std::pair<int, int>;
using EdgeProp = std::pair<double, double>;
struct EdgeComp : public std::binary_function<Edge, Edge, bool>
{
    inline bool operator()(const Edge& e1, const Edge& e2) const
    {
        return (e1.first < e2.first) || (e1.first == e2.first && e1.second < e2.second);
    }
};

void
sort_edge(Edge& e)
{
    if (e.first > e.second)
    {
        int temp = e.second;
        e.second = e.first;
        e.first = temp;
    }
    return;
}

std::string to_string_with_precision(const double a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

int
main(int /*argc*/, char** /*argv*/)
{
    // Problem parameters
    const double Horizon_size = 3.015;
    const int ndivx = 5;           // num points in x direction.
    const int ndivy = 10 * (ndivx - 1);           // num points in x direction.
    const int totnode = ndivx * ndivy;

    const double length = 0.1;      // total length of the plate (cm)
    const double dx = length / (ndivx - 1);
    const double dy = dx;
    const double height = 1.0;      // total width of the plate (cm)
    const double x_c = 0.95;
    const double y_c = dy/2;

    const double area = dx * dy;  // cross-sectional area
    const double vol = area * dx; // volume of a material point
    const double delta = Horizon_size * dx;
    const double scr0 = 3.1; // critical stretch.
    
    std::cout << "------------------------------- " << "\n";
    std::cout << "nx = " << ndivx << "\n";
    std::cout << "ny = " << ndivy << "\n";
    std::cout << "total number of nodes = " << totnode << "\n";
    std::cout << "------------------------------- " << "\n";

    // Initialize vertices
    std::vector<std::vector<double> > coord(totnode, std::vector<double>(2));
    int nnum = -1;
    double x, y;
    // Material points in the region.

        for (int j = 0; j <= (ndivy-1); ++j)
        {
            for (int i = 0; i <= (ndivx-1); ++i)
            {
                nnum += 1;
                x = i * dx;
                y = j * dy;

                coord[nnum][0] = x + x_c; // shift in x-direction
                coord[nnum][1] = y + y_c; // shift in y-direction
            }
        }

    // Initialize springs
    // Determination of material points inside the horizon of each material point
    std::vector<int> numfam(totnode, 0);
    std::map<Edge, EdgeProp, EdgeComp> mesh;
    for (int i = 0; i < totnode; ++i)
    {
        for (int j = 0; j < totnode; ++j)
        {
            const double idist = sqrt(pow((coord[j][0] - coord[i][0]), 2) + pow((coord[j][1] - coord[i][1]), 2));
            if (i != j)
            {
                if (idist <= delta)
                {
                    numfam[i] += 1;
                    Edge e = std::make_pair(i, j);
                    sort_edge(e);

                    // Determine initial failure of the bond.
                    double fail = 1.0;

                    EdgeProp prop = std::make_pair(idist, fail);
                    mesh[e] = prop;
                }
            }
        }
    }
    const auto num_springs = static_cast<int>(mesh.size());

    // Step 1:  Write out the vertex information
    std::fstream vertex_stream;
    vertex_stream.open("band2d.vertex", std::fstream::out);
    vertex_stream.setf(std::ios_base::scientific);
    vertex_stream.precision(12);

    // first line is the number of vertices in the file
    vertex_stream << totnode << "\n";

    // remaining lines are the initial coordinates of each vertex
    for (int i = 0; i < totnode; ++i)
    {
        vertex_stream << coord[i][0] << "\t" << coord[i][1] << "\n";
    }
    vertex_stream.close();

    // Step 2: Write out the link information (including connectivity and
    // material parameters).
    std::fstream spring_stream;
    spring_stream.open("band2d.spring", std::fstream::out);
    spring_stream.setf(std::ios_base::scientific);
    spring_stream.precision(12);

    // first line is the number of edges in the file
    spring_stream << num_springs << "\n";

    // remaining lines are the edges in the mesh
    const double K = 0.0;
    const int force_fcn_idx = 0;

    for (const auto& edge_edge_prop_pair : mesh)
    {
        const Edge& e = edge_edge_prop_pair.first;
        const EdgeProp& prop = edge_edge_prop_pair.second;

        const int& idx_master = e.first;
        const int& idx_slave = e.second;
        const double& rest_length = prop.first;
        const double& fail = prop.second;
        const double vol_master = vol;
        const double vol_slave = vol;
        spring_stream << idx_master << " " << idx_slave << " " << K << " " << rest_length << " " << force_fcn_idx << " "
                      << vol_master << " " << vol_slave << " " << fail << " " << scr0 << "\n";
    }
    spring_stream.close();
}