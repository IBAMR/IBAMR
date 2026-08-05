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

int
main(int /*argc*/, char** /*argv*/)
{
    // // Problem parameters
    // const double R = 0.3;
    // const double w = 0.3;
    // const double x_c = 0.2;
    // const double y_c = 0.2;

    // const int N = 128;
    // const double dx = 1.0 / (double)N;
    // const double MFAC = 1.0;
    // const double ds = MFAC * dx;

    // const int n_x = ceil(R / ds);
    // const int n_y = ceil(w / ds);
    // const int totnode = n_x * n_y;

    // const double area = ds * ds;  // cross-sectional area
    // const double vol = area * ds; // volume of a material point
    // const double delta = 2.015 * ds;
    // const double scr0 = 50 * ds; // critical stretch.

    // compressed block
    const int ndivx = 33; // num points in x direction.
    const int ndivy = 17;  // num points in y direction.
    const int nbnd = 0;
    const int totnode = (ndivx + 2 * nbnd) * (ndivy + 2 * nbnd);

    const double length = 20.0;              // total length of the plate (m)
    const double width = 10.0;               // total width of the plate (m)

    const double dx = length / (ndivx - 1); // spacing between material points in x direction
    const double dy = width / (ndivy - 1);  // spacing between material points in y direction
    const double delta = 2.015 * dx; // 3.015 * dx;      // horizon
    const double thick = dx;         // thickness of the plate
    const double area = dx * dx;     // cross-sectional area
    const double vol = area * thick; // volume of a material point

    const double scr0 = 30.0; // critical stretch

    // Initialize vertices
    std::vector<std::vector<double> > coord(totnode, std::vector<double>(2));
    int nnum = -1;
    // double x, y;
    // // Material points in the region.

    //     for (int j = 0; j <= (n_y - 1); ++j)
    //     {
    //         for (int i = 0; i <= (n_x - 1); ++i)
    //         {
    //             nnum += 1;
    //             x = i * ds;
    //             y = j * ds;

    //             coord[nnum][0] = x + x_c; // shift in x-direction
    //             coord[nnum][1] = y + y_c; // shift in y-direction
    //         }
    //     }

    // Material points in the region.

        for (int j = 0; j <= (ndivy - 1); ++j)
        {
            for (int i = 0; i <= (ndivx - 1); ++i)
            {
                nnum += 1;
                coord[nnum][0] = i * dx + 10.0; // shift in x-direction
                coord[nnum][1] = j * dy + 15.0; // shift in y-direction
            }
        }

    std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "total number of nodes = " << totnode << std::endl;
    std::cout << "total number of nodes in x = " << ndivx << std::endl;
    std::cout << "total number of nodes in y = " << ndivy << std::endl;
    std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";

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
    vertex_stream.open("sheet2d.vertex", std::fstream::out);
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
    spring_stream.open("sheet2d.spring", std::fstream::out);
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
