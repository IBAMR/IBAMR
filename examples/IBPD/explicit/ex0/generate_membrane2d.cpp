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
main(int argc, char** argv)
{
    // Problem parameters

    const int ndivx = 200; // num points in x direction.
    const int ndivy = 200; // num points in y direction.
    const int nbnd = 3;
    const int totnode = (ndivx + 2 * nbnd) * (ndivy + 2 * nbnd);

    const double length = 1.0;              // total length of the plate (m)
    const double width = 1.0;               // total width of the plate (m)
    const double dx = length / (ndivx - 1); // spacing between material points in x direction
    const double dy = width / (ndivy - 1);  // spacing between material points in y direction
    const double delta = 3.015 * dx;        // 3.015 * dx;      // horizon
    const double thick = dx;                // thickness of the plate
    const double area = dx * dx;            // cross-sectional area
    const double vol = area * dx;           // volume of a material point

    const double scr0 = 0.04472; // critical stretch

    // Initialize vertices

    double coord[totnode][2];
    int nnum = -1;

    // Material points of the internal region
    for (int j = 0; j <= (ndivy - 1); ++j)
    {
        for (int i = 0; i <= (ndivx - 1); ++i)
        {
            nnum += 1;
            coord[nnum][0] = i * dx;
            coord[nnum][1] = j * dy;
        }
    }
    const int interior_begin = 0;
    const int interior_end = nnum;

    // Material points of the bottom boundary region
    const int bottom_begin = interior_end + 1;
    for (int j = -nbnd; j <= -1; ++j)
    {
        for (int i = -nbnd; i <= (ndivx + nbnd - 1); ++i)
        {
            nnum = nnum + 1;
            coord[nnum][0] = i * dx;
            coord[nnum][1] = j * dy;
        }
    }
    const int bottom_end = nnum;

    // Material points of the top boundary region
    const int top_begin = bottom_end + 1;
    for (int j = ndivy; j <= (ndivy + nbnd - 1); ++j)
    {
        for (int i = -nbnd; i <= (ndivx + nbnd - 1); ++i)
        {
            nnum = nnum + 1;
            coord[nnum][0] = i * dx;
            coord[nnum][1] = j * dy;
        }
    }
    const int top_end = nnum;

    // Material points of the left boundary region
    const int left_begin = top_end + 1;
    for (int j = 0; j <= (ndivy - 1); ++j)
    {
        for (int i = -nbnd; i <= -1; ++i)
        {
            nnum = nnum + 1;
            coord[nnum][0] = i * dx;
            coord[nnum][1] = j * dy;
        }
    }
    const int left_end = nnum;

    // Material points of the right boundary region
    const int right_begin = left_end + 1;
    for (int j = 0; j <= (ndivy - 1); ++j)
    {
        for (int i = ndivx; i <= (ndivx + nbnd - 1); ++i)
        {
            nnum = nnum + 1;
            coord[nnum][0] = i * dx;
            coord[nnum][1] = j * dy;
        }
    }
    const int right_end = nnum;

    assert(right_end + 1 == totnode);

    std::cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~\n";
    std::cout << "static const int interior_begin = " << interior_begin << std::endl;
    std::cout << "static const int interior_end = " << interior_end << std::endl;
    std::cout << "static const int bottom_begin = " << bottom_begin << std::endl;
    std::cout << "static const int bottom_end = " << bottom_end << std::endl;
    std::cout << "static const int top_begin = " << top_begin << std::endl;
    std::cout << "static const int top_end = " << top_end << std::endl;
    std::cout << "static const int left_begin = " << left_begin << std::endl;
    std::cout << "static const int left_end = " << left_end << std::endl;
    std::cout << "static const int right_begin = " << right_begin << std::endl;
    std::cout << "static const int right_end = " << right_end << std::endl;

    // Initialize springs
    // Determination of material points inside the horizon of each material point
    int numfam[totnode] = { 0 };
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
    const int num_springs = mesh.size();

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
