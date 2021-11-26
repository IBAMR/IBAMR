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
    const int N = 32;
    const double H = 0.41;

    const double dx = H / (double)N;
    const double dy = dx;
    const double MFAC = 1.0;
    const double ds = MFAC * dx;

    const double r = 0.05;               // radius of the cylinder (m)
    const double length = 7.0 * r;      // total length of the plate (m)
    const double height = 0.4 * r;      // total width of the plate (m)
    const double x_c = 0.2;
    const double y_c = 0.2;

    const int ndivx = ceil(length / ds)+1;
    const int ndivy = ceil(height / ds);
    const int totnode = ndivx * ndivy;

    const int M = ceil(2 * M_PI * r / ds);              // number of particles on the circumference of a circle
    const double rad = 2 * M_PI / M;

    const double area = ds * ds;  // cross-sectional area
    const double vol = area * ds; // volume of a material point
    const double delta = 2.015 * ds;
    const double delta1 = 2.015 * ds;
    const double scr0 = 30.0; // critical stretch.

    // // Problem parameters
    // const int ndivx = 41;          // num points in x direction.
    // const int ndivy = 3;           // num points in x direction.
    // const int ndiv = 10;
    // const int totnode = ndivx * ndivy;

    // const double r = 0.05;               // radius of the cylinder (cm)
    // const double length = 8.0 * r;      // total length of the plate (cm)
    // const double height = 0.4 * r;      // total width of the plate (cm)
    // const double x_c = 0.2;
    // const double y_c = 0.2;

    // const double dx = length / (ndivx - 1);
    // const double dy = height / (ndivy - 1);

    // const int M = 4 * (10 + 1);              // number of particles on the circumference of a circle
    // const double rad = 2 * M_PI / M;

    double x, y, theta;
    int num1 = 0;

    // Full disk
    // for (int i = 0; i <= ndiv; ++i)
    // {
    //     for (int j = 0; j <= ndiv; ++j)
    //     {
    //         x = i * dx;
    //         y = j * dy;

    //         const double dist1 = sqrt(pow(x-r,2.0) + pow(y-r,2.0));
    //         if (dist1 < r)
    //         {
    //             num1 += 1;
    //         }
    //     }
    // }
    // std::cout << "num1 =" << num1 << "\n";

    // const double area = dx * dy;  // cross-sectional area
    // const double vol = area * dx; // volume of a material point
    // const double delta = 2.015 * dx;
    // const double delta1 = 1.015 * dx;
    // const double scr0 = 3.1; // critical stretch.
    
    std::cout << "------------------------------- " << "\n";
    std::cout << "nx = " << ndivx << "\n";
    std::cout << "ny = " << ndivy << "\n";
    std::cout << "total number of nodes = " << totnode + M + num1 << "\n";
    std::cout << "------------------------------- " << "\n";

    // Initialize vertices
    std::vector<std::vector<double> > coord(totnode + M + num1, std::vector<double>(2));
    int nnum = -1;
    // Material points in the region.

        for (int j = 0; j <= (ndivy-1); ++j)
        {
            for (int i = 0; i <= (ndivx-1); ++i)
            {
                nnum += 1;
                x = i * ds;
                y = j * ds;

                coord[nnum][0] = x + 0.6 - (double(ndivx) - 1.0) * ds; // shift in x-direction
                coord[nnum][1] = y + y_c - (double(ndivy) - 1.0) * ds / 2.0; // shift in y-direction
            }
        }

        for (int j = 0; j <= (M-1); ++j)
        {
            nnum += 1;
            theta = j * rad;

            coord[nnum][0] = r * cos(theta) + x_c;
            coord[nnum][1] = r * sin(theta) + y_c;
        }

        // Full disk
        // for (int i = 0; i <= ndiv; ++i)
        // {
        //     for (int j = 0; j <= ndiv; ++j)
        //     {
        //         x = i * dx;
        //         y = j * dy;

        //         const double dist1 = sqrt(pow(x-r,2.0) + pow(y-r,2.0));
        //         if (dist1 < r)
        //         {
        //             nnum += 1;

        //             coord[nnum][0] = x + x_c - r; // shift in x-direction
        //             coord[nnum][1] = y + y_c - r; // shift in y-direction
        //         }
        //     }
        // }


    // Initialize springs
    // Determination of material points inside the horizon of each material point
    std::vector<int> numfam(totnode + M + num1, 0);
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

    for (int i = totnode; i < totnode + M + num1; ++ i)
    {
        for (int j = totnode; j < totnode + M + num1; ++j)
        {
            const double idist = sqrt(pow((coord[j][0] - coord[i][0]), 2) + pow((coord[j][1] - coord[i][1]), 2));
            if (i != j)
            {
                if (idist <= delta1)
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
    vertex_stream << totnode + M + num1 << "\n";

    // remaining lines are the initial coordinates of each vertex
    for (int i = 0; i < totnode + M + num1; ++i)
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