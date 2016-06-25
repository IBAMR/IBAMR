// Filename: generate_sheet3d.cpp
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
    const double R = 0.2;
    const double w = 0.2;
    const double x_c = 0.2;
    const double y_c = 0.2;
    const double z_c = 0.5;

    const int N = 128;
    const double dx = 1.0 / (double)N;
    const double MFAC = 1.0;
    const double ds = MFAC * dx;

    const int n_x = ceil(R / ds);
    const int n_y = ceil(w / ds);
    const int n_z = 5;
    const int totnode = n_x * n_y * n_z;

    const double area = ds * ds;  // cross-sectional area
    const double vol = area * ds; // volume of a material point
    const double delta = 2.015 * ds;
    const double scr0 = 50 * ds; // critical stretch.

    // Initialize vertices
    std::vector<std::vector<double> > coord(totnode, std::vector<double>(3));
    int nnum = -1;
    double x, y, z;
    // Material points in the region.
    for (int k = -n_z / 2; k <= n_z / 2; ++k)
    {
        for (int j = 0; j <= (n_y - 1); ++j)
        {
            for (int i = 0; i <= (n_x - 1); ++i)
            {
                nnum += 1;
                x = i * ds;
                y = j * ds;
                z = k * ds;

                coord[nnum][0] = x + x * y + x_c;
                coord[nnum][1] = x + y + y_c;
                coord[nnum][2] = z + z_c;
            }
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
    vertex_stream.open("sheet3d.vertex", std::fstream::out);
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
    spring_stream.open("sheet3d.spring", std::fstream::out);
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
