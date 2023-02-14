// ---------------------------------------------------------------------
//
// Copyright (c) 2016 - 2019 by the IBAMR developers
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
    const int ndivx = 49;               // M+1; num points in x direction.
    const int ndivy = 45;               // LH / LW * M + 1; num points in y direction.
    const int nbnd = 0;
    const int totnode = (ndivx + 2 * nbnd) * (ndivy + 2 * nbnd);

    // Cook's membrane
    const double LW = 4.8;              // length of the membrane (cm)
    const double LH = 4.4;              // left height of the membrane(cm)
    const double LHL = 1.6;             // right heigh of the membrane(cm)
    const double dx = LW / (ndivx - 1); // spacing between material points in x direction
    const double dy = LH / (ndivy - 1); // spacing between material points in y direction
    const double dyy = LHL/(ndivy - 1); //

    // // // Problem parameters w/ stair-step geometry
    // // Cook's membrane
    // const int ndivx = 25;                  // M+1; num points in x direction.
    // const int ndivy = 31;                  // (LHL + LHR) / LW * M + 1; num points in y direction.

    // const double LW = 4.8;                 // length of the membrane (cm)
    // const double LHL = 4.4;                // left height of the membrane(cm)
    // const double LHR = 1.6;                // right heigh of the membrane(cm)
    // const double LH = LHL + LHR;

    // const double dx = 0.2;
    // const double dy = dx;

    // const int totnode1 = ndivx * ndivy;
    
    const double horizon = 2.015;
    const double delta = horizon * dx;      // horizon size
    const double thick = dx;                // thickness of the plate
    const double c_area = dx * dx;            // cross-sectional area
    const double vol = c_area * thick;        // volume of a material point

    const double scr0 = 30.0; // critical stretch

    // Initialize vertices
    double coord[totnode][2];
    // double delta_x[totnode][1];
    int nnum = -1;

    // Material points of the plate region
    for (int i = 0; i <= (ndivx - 1); ++i)
    {
        for (int j = 0; j <= (ndivy - 1); ++j)
        {
                nnum += 1;

                // 2d Cook's membrane
                coord[nnum][0] = i * dx + 2.0;
                coord[nnum][1] = (j * (dyy-dy) + LH) / (ndivx - 1) * i + dy * j + 1.0;

        }
    }

    // // Initialize vertices for stair-step geometry
    // double coord1[totnode1][2];
    // double x,y;
    // const double eps = 1.0e-7;

    // // double delta_x[totnode][1];
    // int nnum = -1;

    // // Material points of the plate region
    // for (int i = 0; i <= (ndivx - 1); ++i)
    // {
    //     for (int j = 0; j <= (ndivy - 1); ++j)
    //     {
    //         x = i * dx;
    //         y = j * dy;

    //         if ( (4.4 / 4.8 * x - eps <= y) && (y <= 1.6 / 4.8 * x + 4.4 + eps))
    //         {
    //             nnum += 1;
    //             coord1[nnum][0] = x + 2.0;
    //             coord1[nnum][1] = y + 1.0;
    //         }
    //     }
    // }
    // const int totnode = nnum + 1;

    // double coord[totnode][2];
    // for (int i = 0; i < totnode; ++i)
    // {
    //     coord[i][0] = coord1[i][0];
    //     coord[i][1] = coord1[i][1];
    // }

    // const double dx1 = dx * 1.6 / 4.4;

    double volume[totnode];
    double p0[2], p1[2], p2[2], p3[2];
    double area, vol_p;
    nnum = -1;
    for (int i = 0; i <= (ndivx - 1); ++i)
    {
        for (int j = 0; j <= (ndivy - 1); ++j)
        {
            nnum += 1;

            if (i == 0)
            {
                if (j == 0)
                {
                    p0[0] = coord[nnum][0];
                    p0[1] = coord[nnum][1];
                    p1[0] = (coord[nnum][0] + coord[nnum+1][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum+1][1])/2.0;
                    p2[0] = (coord[nnum+1][0] + coord[nnum+ndivy][0])/2.0;
                    p2[1] = (coord[nnum+1][1] + coord[nnum+ndivy][1])/2.0;
                    p3[0] = (coord[nnum][0] + coord[nnum+ndivy][0])/2.0;
                    p3[1] = (coord[nnum][1] + coord[nnum+ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 4.0;
                }
                else if (j == ndivy - 1)
                {
                    p0[0] = coord[nnum][0];
                    p0[1] = coord[nnum][1];
                    p1[0] = (coord[nnum][0] + coord[nnum-1][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum-1][1])/2.0;
                    p2[0] = (coord[nnum-1][0] + coord[nnum+ndivy][0])/2.0;
                    p2[1] = (coord[nnum-1][1] + coord[nnum+ndivy][1])/2.0;
                    p3[0] = (coord[nnum][0] + coord[nnum+ndivy][0])/2.0;
                    p3[1] = (coord[nnum][1] + coord[nnum+ndivy][1])/2.0;
                    
                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 4.0;
                }
                else
                {
                    p0[0] = (coord[nnum][0] + coord[nnum-1][0])/2.0;
                    p0[1] = (coord[nnum][1] + coord[nnum-1][1])/2.0;
                    p1[0] = (coord[nnum][0] + coord[nnum+1][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum+1][1])/2.0;
                    p2[0] = (coord[nnum+1][0] + coord[nnum+ndivy][0])/2.0;
                    p2[1] = (coord[nnum+1][1] + coord[nnum+ndivy][1])/2.0;
                    p3[0] = (coord[nnum-1][0] + coord[nnum+ndivy][0])/2.0;
                    p3[1] = (coord[nnum-1][1] + coord[nnum+ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 2.0;
                }
            }
            else if (i == ndivx - 1)
            {
                if (j == 0)
                {
                    p0[0] = coord[nnum][0];
                    p0[1] = coord[nnum][1];
                    p1[0] = (coord[nnum][0] + coord[nnum+1][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum+1][1])/2.0;
                    p2[0] = (coord[nnum+1][0] + coord[nnum-ndivy][0])/2.0;
                    p2[1] = (coord[nnum+1][1] + coord[nnum-ndivy][1])/2.0;
                    p3[0] = (coord[nnum][0] + coord[nnum-ndivy][0])/2.0;
                    p3[1] = (coord[nnum][1] + coord[nnum-ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 4.0;
                }
                else if (j == ndivy - 1)
                {
                    p0[0] = coord[nnum][0];
                    p0[1] = coord[nnum][1];
                    p1[0] = (coord[nnum][0] + coord[nnum-1][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum-1][1])/2.0;
                    p2[0] = (coord[nnum-1][0] + coord[nnum-ndivy][0])/2.0;
                    p2[1] = (coord[nnum-1][1] + coord[nnum-ndivy][1])/2.0;
                    p3[0] = (coord[nnum][0] + coord[nnum-ndivy][0])/2.0;
                    p3[1] = (coord[nnum][1] + coord[nnum-ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 4.0;
                }
                else
                {
                    p0[0] = (coord[nnum][0] + coord[nnum-1][0])/2.0;
                    p0[1] = (coord[nnum][1] + coord[nnum-1][1])/2.0;
                    p1[0] = (coord[nnum][0] + coord[nnum+1][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum+1][1])/2.0;
                    p2[0] = (coord[nnum+1][0] + coord[nnum-ndivy][0])/2.0;
                    p2[1] = (coord[nnum+1][1] + coord[nnum-ndivy][1])/2.0;
                    p3[0] = (coord[nnum-1][0] + coord[nnum-ndivy][0])/2.0;
                    p3[1] = (coord[nnum-1][1] + coord[nnum-ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 2.0;
                }
            }
            else
            {
                if (j == 0)
                {
                    p0[0] = (coord[nnum][0] + coord[nnum-ndivy][0])/2.0;
                    p0[1] = (coord[nnum][1] + coord[nnum-ndivy][1])/2.0;
                    p1[0] = (coord[nnum][0] + coord[nnum+ndivy][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum+ndivy][1])/2.0;
                    p2[0] = (coord[nnum+1][0] + coord[nnum+ndivy][0])/2.0;
                    p2[1] = (coord[nnum+1][1] + coord[nnum+ndivy][1])/2.0;
                    p3[0] = (coord[nnum+1][0] + coord[nnum-ndivy][0])/2.0;
                    p3[1] = (coord[nnum+1][1] + coord[nnum-ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 2.0;
                }
                else if (j == ndivy - 1)
                {
                    p0[0] = (coord[nnum][0] + coord[nnum-ndivy][0])/2.0;
                    p0[1] = (coord[nnum][1] + coord[nnum-ndivy][1])/2.0;
                    p1[0] = (coord[nnum][0] + coord[nnum+ndivy][0])/2.0;
                    p1[1] = (coord[nnum][1] + coord[nnum+ndivy][1])/2.0;
                    p2[0] = (coord[nnum-1][0] + coord[nnum+ndivy][0])/2.0;
                    p2[1] = (coord[nnum-1][1] + coord[nnum+ndivy][1])/2.0;
                    p3[0] = (coord[nnum-1][0] + coord[nnum-ndivy][0])/2.0;
                    p3[1] = (coord[nnum-1][1] + coord[nnum-ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx * 2.0;
                }
                else
                {
                    p0[0] = (coord[nnum+1][0] + coord[nnum-ndivy][0])/2.0;
                    p0[1] = (coord[nnum+1][1] + coord[nnum-ndivy][1])/2.0;
                    p1[0] = (coord[nnum+1][0] + coord[nnum+ndivy][0])/2.0;
                    p1[1] = (coord[nnum+1][1] + coord[nnum+ndivy][1])/2.0;
                    p2[0] = (coord[nnum-1][0] + coord[nnum+ndivy][0])/2.0;
                    p2[1] = (coord[nnum-1][1] + coord[nnum+ndivy][1])/2.0;
                    p3[0] = (coord[nnum-1][0] + coord[nnum-ndivy][0])/2.0;
                    p3[1] = (coord[nnum-1][1] + coord[nnum-ndivy][1])/2.0;

                    // ds = (sqrt(pow(p1[0]-p0[0],2.0) + pow(p1[1] - p0[1],2.0)) + sqrt(pow(p3[0]-p0[0],2.0) + pow(p3[1] - p0[1],2.0)));

                    area = abs((p0[0]*p1[1] + p1[0]*p2[1] + p2[0]*p3[1] + p3[0]*p0[1]) - (p1[0]*p0[1] + p2[0]*p1[1] + p3[0]*p2[1] + p0[0] * p3[1]))/2.0;
                    vol_p = area * dx;
                }
            }
            volume[nnum] = vol_p;
        }
    }

    std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    // std::cout << "static const int left_begin = " << 0 << ";" << std::endl;
    // std::cout << "static const int left_end = " << ndivy - 1 << ";" << std::endl;
    // std::cout << "static const int right_begin = " << totnode - ndivy << ";" << std::endl;
    std::cout << "static const int total nodes = " << totnode - 1 << ";" << std::endl;

    // Initialize springs
    // Determination of material points inside the horizon of each material point
    int numfam[totnode];
    std::map<Edge, EdgeProp, EdgeComp> mesh;
    for (int i = 0; i < totnode; ++i)
    {   numfam[i] = 0;
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
                    // if ((coord[i][1] == 5.4 && coord[i][0] < 2.0) && (coord[j][1] > 5.4))
                    // {
                    //     fail = 0.0;
                    // }
                    // if ((coord[j][1] == 5.4 && coord[j][0] < 2.0) && (coord[i][1] > 5.4))
                    // {
                    //     fail = 0.0;
                    // }
                    EdgeProp prop = std::make_pair(idist, fail);
                    mesh[e] = prop;
                }
            }
        }
    }
    const auto num_springs = static_cast<int>(mesh.size());

    // Step 1:  Write out the vertex information
    std::fstream vertex_stream;
    vertex_stream.open("membrane2d.vertex", std::fstream::out);
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
    spring_stream.open("membrane2d.spring", std::fstream::out);
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
        const double vol_master = volume[idx_master];
        const double vol_slave = volume[idx_slave];
        spring_stream << idx_master << " " << idx_slave << " " << K << " " << rest_length << " " << force_fcn_idx << " "
                      << vol_master << " " << vol_slave << " " << fail << " " << scr0 << "\n";
    }
    spring_stream.close();
}
