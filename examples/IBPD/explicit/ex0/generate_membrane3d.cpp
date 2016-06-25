// Filename: generate_membrane3d.cpp
// Created on 24 Jun 2016 by Amneet Bhalla
// All rights reserved.

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <utility>
#include <cmath>
#include <cassert>

typedef std::pair<int, int> Edge;
typedef std::pair<double, double> EdgeProp;
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
    // Problem parameters

    const int ndivx = 151; // num points in x direction.
    const int ndivy = 76;  // num points in y direction.
    const int nbnd = 0;
    const int ndivz = 3; // num layers in z direction.
    const int totnode = (ndivx + 2 * nbnd) * (ndivy + 2 * nbnd) * ndivz;

    const double length = 1.0;              // total length of the plate (m)
    const double width = 0.5;               // total width of the plate (m)
    const double dx = length / (ndivx - 1); // spacing between material points in x direction
    const double dy = width / (ndivy - 1);  // spacing between material points in y direction
    const double dz = dx;
    const double delta = 2.015 * dx; // 3.015 * dx;      // horizon
    const double thick = dx;         // thickness of the plate
    const double area = dx * dx;     // cross-sectional area
    const double vol = area * thick; // volume of a material point

    const double scr0 = 0.5; // critical stretch

    // Initialize vertices

    double coord[totnode][3];
    int nnum = -1;

    // Material points of the plate region
    for (int i = 0; i <= (ndivx - 1); ++i)
    {
        for (int k = 0; k <= (ndivz - 1); ++k)
        {
            for (int j = 0; j <= (ndivy - 1); ++j)
            {
                nnum += 1;
                coord[nnum][0] = i * dx;
                coord[nnum][1] = j * dy;
                coord[nnum][2] = k * dz;
            }
        }
    }

    std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "static const int left_begin = " << 0 << ";" << std::endl;
    std::cout << "static const int left_end = " << ndivy * ndivz - 1 << ";" << std::endl;
    std::cout << "static const int right_begin = " << totnode - ndivy * ndivz << ";" << std::endl;
    std::cout << "static const int right_end = " << totnode - 1 << ";" << std::endl;

    // Initialize springs
    // Determination of material points inside the horizon of each material point
    int numfam[totnode] = { 0 };
    std::map<Edge, EdgeProp, EdgeComp> mesh;
    for (int i = 0; i < totnode; ++i)
    {
        for (int j = 0; j < totnode; ++j)
        {
            const double idist = sqrt(pow((coord[j][0] - coord[i][0]), 2) + pow((coord[j][1] - coord[i][1]), 2) +
                                      pow((coord[j][2] - coord[i][2]), 2));
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
    const int num_springs = static_cast<int>(mesh.size());

    // Step 1:  Write out the vertex information
    std::fstream vertex_stream;
    vertex_stream.open("membrane3d.vertex", std::fstream::out);
    vertex_stream.setf(std::ios_base::scientific);
    vertex_stream.precision(12);

    // first line is the number of vertices in the file
    vertex_stream << totnode << "\n";

    // remaining lines are the initial coordinates of each vertex
    for (int i = 0; i < totnode; ++i)
    {
        vertex_stream << coord[i][0] << "\t" << coord[i][1] << "\t" << coord[i][2] << "\n";
    }
    vertex_stream.close();

    // Step 2: Write out the link information (including connectivity and
    // material parameters).
    std::fstream spring_stream;
    spring_stream.open("membrane3d.spring", std::fstream::out);
    spring_stream.setf(std::ios_base::scientific);
    spring_stream.precision(12);

    // first line is the number of edges in the file
    spring_stream << num_springs << "\n";

    // remaining lines are the edges in the mesh
    const double K = 0.0;
    const int force_fcn_idx = 0;

    for (std::map<Edge, EdgeProp, EdgeComp>::const_iterator it = mesh.begin(); it != mesh.end(); ++it)
    {
        const Edge& e = it->first;
        const EdgeProp& prop = it->second;

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
