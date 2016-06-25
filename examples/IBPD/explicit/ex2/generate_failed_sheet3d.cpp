// Filename: generate_failed_sheet3d.cpp
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
    const double R = 0.4;
    const double w = 0.4;
    const double x_c = 0.3;
    const double y_c = 0.3;
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
    const double scr0 = 0.2; // critical stretch.
    const int n_x_crack = 20;

    std::cout << "n_x = " << n_x << "\tdx = " << dx << "\tx_crack_end = " << n_x_crack * ds
              << "\ty_crack = " << (n_y / 2 + 0.5) * ds << std::endl;

    // Initialize vertices
    std::vector<std::vector<double> > coord(totnode, std::vector<double>(3));
    std::vector<std::vector<double> > coord_ref(totnode, std::vector<double>(3));
    std::vector<int> target_lower(n_x_crack * n_z);
    std::vector<int> target_upper(n_x_crack * n_z);
    std::vector<int> target_right(n_y * n_z);
    std::vector<int> target_left(n_y * n_z);

    int nnum = -1, tupper = -1, tlower = -1, tright = -1, tleft = -1;
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

                coord_ref[nnum][0] = x;
                coord_ref[nnum][1] = y;
                coord_ref[nnum][2] = z;

                coord[nnum][0] = x + x_c;
                coord[nnum][1] = y + y_c;
                coord[nnum][2] = z + z_c;

                if (j == 0 && (i >= 0 && i < n_x_crack))
                {
                    tlower += 1;
                    target_lower[tlower] = nnum;
                }

                if (j == n_y - 1 && (i >= 0 && i < n_x_crack))
                {
                    tupper += 1;
                    target_upper[tupper] = nnum;
                }

                if (i == 0)
                {
                    tleft += 1;
                    target_left[tleft] = nnum;
                }

                if (i == n_x - 1)
                {
                    tright += 1;
                    target_right[tright] = nnum;
                }
            }
        }
    }

    const double x_crack_begin = 0.0;
    const double x_crack_end = n_x_crack * ds;
    const double y_crack = (n_y / 4 + 0.5) * ds;

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
                    const int lag_small = e.first;
                    const int lag_big = e.second;

                    if (((coord_ref[lag_small][1] < y_crack && coord_ref[lag_big][1] > y_crack) &&
                         (coord_ref[lag_big][0] >= x_crack_begin && coord_ref[lag_big][0] <= x_crack_end)) ||
                        ((coord_ref[lag_small][1] < R - y_crack && coord_ref[lag_big][1] > R - y_crack) &&
                         (coord_ref[lag_big][0] >= x_crack_begin && coord_ref[lag_big][0] <= x_crack_end))

                            )
                    {
                        fail = 0.0;
                    }

                    EdgeProp prop = std::make_pair(idist, fail);
                    mesh[e] = prop;
                }
            }
        }
    }
    const int num_springs = static_cast<int>(mesh.size());

    // Step 1:  Write out the vertex information
    std::fstream vertex_stream;
    vertex_stream.open("failed_sheet3d.vertex", std::fstream::out);
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
    spring_stream.open("failed_sheet3d.spring", std::fstream::out);
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

    // Write target file
    std::fstream target_stream;
    target_stream.open("failed_sheet3d.target", std::fstream::out);
    target_stream.setf(std::ios_base::scientific);
    target_stream.precision(12);

    // first line is the number of target in the file
    const int num_one_side = n_x_crack * n_z;
    const int num_target = 2 * num_one_side;
    target_stream << num_target << "\n";

    for (int i = 0; i < num_one_side; ++i)
    {
        target_stream << target_upper[i] << "\t" << 3.0 << "\t" << 0.0 << "\n";
    }
    for (int i = 0; i < num_one_side; ++i)
    {
        target_stream << target_lower[i] << "\t" << 1.0 << "\t" << 0.0 << "\n";
    }
    target_stream.close();

    return 0;
}
