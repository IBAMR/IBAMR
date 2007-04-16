#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

using namespace std;

void
construct_initial_triangulation(
    vector<vector<int> >& tri,
    const vector<vector<double> >& vertex,
    const double& L)
{
    // Find all vertices that are the specified distance L from each other.
    map<int,set<int> > vertex_map;
    for (size_t k = 0; k < vertex.size(); ++k)
    {
        vertex_map[k].clear();
        const vector<double>& X0 = vertex[k];
        for (size_t j = k+1; j < vertex.size(); ++j)
        {
            const vector<double>& X1 = vertex[j];

            double r = 0.0;
            for (int d = 0; d < 3; ++d)
            {
                r += pow(X0[d]-X1[d],2.0);
            }
            r = sqrt(r);

            if (abs(r-L) < sqrt(numeric_limits<double>::epsilon()))
            {
                vertex_map[k].insert(j);
            }
        }
    }

    // Use the vertex map to construct the initial triangulation.
    for (size_t i = 0; i < vertex.size(); ++i)
    {
        for (set<int>::const_iterator i_it = vertex_map[i].begin();
             i_it != vertex_map[i].end(); ++i_it)
        {
            const int j = (*i_it);
            for (set<int>::const_iterator j_it = vertex_map[j].begin();
                 j_it != vertex_map[j].end(); ++j_it)
            {
                const int k = (*j_it);
                if (vertex_map[i].count(k) == 1)
                {
                    int T[3] = {i , j , k};
                    tri.push_back(vector<int>(T,T+3));
                }
            }
        }
    }
    return;
}// construct_initial_triangulation

void
project(
    vector<vector<double> >& vertex,
    const double& R)
{
    for (size_t k = 0; k < vertex.size(); ++k)
    {
        double r = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            r += pow(vertex[k][d],2.0);
        }
        r = sqrt(r);

        for (int d = 0; d < 3; ++d)
        {
            vertex[k][d] *= R/r;
        }
    }
    return;
}// project

void
compute_max_min_edge_lengths(
    double& max_len,
    double& min_len,
    const vector<vector<double> >& vertex,
    const vector<vector<int> >& tri)
{
    max_len = -0.9*numeric_limits<double>::max();
    min_len = +0.9*numeric_limits<double>::max();

    for (size_t l = 0; l < tri.size(); ++l)
    {
        const int& k0 = tri[l][0];
        const int& k1 = tri[l][1];
        const int& k2 = tri[l][2];

        double r01 = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            r01 += pow(vertex[k0][d]-vertex[k1][d],2.0);
        }
        r01 = sqrt(r01);

        if (r01 > max_len) max_len = r01;
        if (r01 < min_len) min_len = r01;

        double r12 = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            r12 += pow(vertex[k1][d]-vertex[k2][d],2.0);
        }
        r12 = sqrt(r12);

        if (r12 > max_len) max_len = r12;
        if (r12 < min_len) min_len = r12;

        double r20 = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            r20 += pow(vertex[k2][d]-vertex[k0][d],2.0);
        }
        r20 = sqrt(r20);

        if (r20 > max_len) max_len = r20;
        if (r20 < min_len) min_len = r20;
    }
    return;
}// compute_max_min_edge_lengths

template<typename T1, typename T2>
struct pair_comp
    : public std::binary_function<pair<T1,T2>,pair<T1,T2>,bool>
{
    inline bool
    operator()(
        const pair<T1,T2>& p1,
        const pair<T1,T2>& p2) const
        {
            return (p1.first < p2.first) || (p1.first == p2.first && p1.second < p2.second);
        }
};

int
main(
    int argc,
    char* argv)
{
    // The sphere center.
    static const double X_center[3] = { 0.0 , 0.0 , 0.0 };

    // The sphere radius.
    static const double R = 0.5;

    // The desired edge length.
    static const double target_edge_len = 0.03125;

#if 0
    // Setup the initial mesh to correspond to a tetrahedron.
    double X0[3] = {+1.0 , +1.0 , +1.0};
    double X1[3] = {-1.0 , -1.0 , +1.0};
    double X2[3] = {-1.0 , +1.0 , -1.0};
    double X3[3] = {+1.0 , -1.0 , -1.0};

    vector<vector<double> > vertex;
    vertex.push_back(vector<double>(X0,X0+3));
    vertex.push_back(vector<double>(X1,X1+3));
    vertex.push_back(vector<double>(X2,X2+3));
    vertex.push_back(vector<double>(X3,X3+3));

    vector<vector<int> > tri;
    construct_initial_triangulation(tri, vertex, sqrt(8.0));
    assert(tri.size() == 4);
#endif

#if 0
    // Setup the initial mesh to correspond to an octahedron.
    double X0[3] = {-1.0 ,  0.0 ,  0.0};
    double X1[3] = {+1.0 ,  0.0 ,  0.0};
    double X2[3] = { 0.0 , -1.0 ,  0.0};
    double X3[3] = { 0.0 , +1.0 ,  0.0};
    double X4[3] = { 0.0 ,  0.0 , -1.0};
    double X5[3] = { 0.0 ,  0.0 , +1.0};

    vector<vector<double> > vertex;
    vertex.push_back(vector<double>(X0,X0+3));
    vertex.push_back(vector<double>(X1,X1+3));
    vertex.push_back(vector<double>(X2,X2+3));
    vertex.push_back(vector<double>(X3,X3+3));
    vertex.push_back(vector<double>(X4,X4+3));
    vertex.push_back(vector<double>(X5,X5+3));

    vector<vector<int> > tri;
    construct_initial_triangulation(tri, vertex, sqrt(2.0));
    assert(tri.size() == 8);
#endif

#if 1
    // Setup the initial mesh to correspond to an icosohedron.
    static const double phi = (1.0+sqrt(5.0))/2.0;  // the golden ratio

    double  X0[3] = { 0.0 , -1.0 , -phi};
    double  X1[3] = { 0.0 , +1.0 , -phi};
    double  X2[3] = { 0.0 , -1.0 , +phi};
    double  X3[3] = { 0.0 , +1.0 , +phi};
    double  X4[3] = {-1.0 , -phi ,  0.0};
    double  X5[3] = {+1.0 , -phi ,  0.0};
    double  X6[3] = {-1.0 , +phi ,  0.0};
    double  X7[3] = {+1.0 , +phi ,  0.0};
    double  X8[3] = {-phi ,  0.0 , -1.0};
    double  X9[3] = {+phi ,  0.0 , -1.0};
    double X10[3] = {-phi ,  0.0 , +1.0};
    double X11[3] = {+phi ,  0.0 , +1.0};

    vector<vector<double> > vertex;
    vertex.push_back(vector<double>( X0, X0+3));
    vertex.push_back(vector<double>( X1, X1+3));
    vertex.push_back(vector<double>( X2, X2+3));
    vertex.push_back(vector<double>( X3, X3+3));
    vertex.push_back(vector<double>( X4, X4+3));
    vertex.push_back(vector<double>( X5, X5+3));
    vertex.push_back(vector<double>( X6, X6+3));
    vertex.push_back(vector<double>( X7, X7+3));
    vertex.push_back(vector<double>( X8, X8+3));
    vertex.push_back(vector<double>( X9, X9+3));
    vertex.push_back(vector<double>(X10,X10+3));
    vertex.push_back(vector<double>(X11,X11+3));

    vector<vector<int> > tri;
    construct_initial_triangulation(tri, vertex, 2.0);
    assert(tri.size() == 20);
#endif

    // Project the initial vertices onto the sphere of the desired radius.
    project(vertex, R);

    // Determine the maximum and minimum edge lengths (whose values should be
    // equal to machine precision).
    double max_len, min_len;
    compute_max_min_edge_lengths(max_len, min_len, vertex, tri);
    cout << "max edge length = " << max_len << "\n";
    cout << "min edge length = " << min_len << "\n";
    cout << "\n";
    assert(abs(max_len-min_len) <= numeric_limits<double>::epsilon());

    // Subdivide the mesh until the edge lengths are less than the requested
    // value.
    int it_num = 0;
    while (max_len >= target_edge_len-numeric_limits<double>::epsilon())
    {
        cout << "iteration " << ++it_num << "\n";

        // Subdivide all triangles in the mesh.
        vector<vector<int> > tri_new;
        for (size_t l = 0; l < tri.size(); ++l)
        {
            // indices of the vertices of the triangle
            const int& k0 = tri[l][0];
            const int& k1 = tri[l][1];
            const int& k2 = tri[l][2];

            // indices of the (new) vertices of the midpoints of the edges of
            // the triangle
            const int k01 = vertex.size()+0;
            const int k12 = vertex.size()+1;
            const int k20 = vertex.size()+2;

            vertex.push_back(vector<double>(3));
            vertex.push_back(vector<double>(3));
            vertex.push_back(vector<double>(3));

            // vertices of the vertices of the triangle
            const vector<double>& v0 = vertex[k0];
            const vector<double>& v1 = vertex[k1];
            const vector<double>& v2 = vertex[k2];

            // vertex of the midpoint of (v0,v1)
            vector<double>& v01 = vertex[k01];
            for (int d = 0; d < 3; ++d)
            {
                v01[d] = (v0[d]+v1[d])/2.0;
            }

            // vertex of the midpoint of (v1,v2)
            vector<double>& v12 = vertex[k12];
            for (int d = 0; d < 3; ++d)
            {
                v12[d] = (v1[d]+v2[d])/2.0;
            }

            // vertex of the midpoint of (v2,v0)
            vector<double>& v20 = vertex[k20];
            for (int d = 0; d < 3; ++d)
            {
                v20[d] = (v2[d]+v0[d])/2.0;
            }

            // the subdivide triangle vertices
            const int T0[3] = { k0 , k01 , k20};
            const int T1[3] = { k1 , k12 , k01};
            const int T2[3] = { k2 , k20 , k12};
            const int T3[3] = {k01 , k12 , k20};

            tri_new.push_back(vector<int>(T0,T0+3));
            tri_new.push_back(vector<int>(T1,T1+3));
            tri_new.push_back(vector<int>(T2,T2+3));
            tri_new.push_back(vector<int>(T3,T3+3));
        }
        tri = tri_new;

        // Project the vertices onto the sphere of the desired radius.
        project(vertex, R);

        // Compute the edge lengths.
        compute_max_min_edge_lengths(max_len, min_len, vertex, tri);
        cout << "max edge length = " << max_len << "\n";
        cout << "min edge length = " << min_len << "\n";
        cout << "\n";
    }

    // Prune duplicate vertices.
    vector<vector<double> > vertex_new;
    vector<int> vertex_map(vertex.size(),-1);
    for (size_t k = 0; k < vertex.size(); ++k)
    {
        const vector<double>& v = vertex[k];

        bool is_duplicate = false;
        for (size_t k_new = 0; k_new < vertex_new.size() && !is_duplicate;
             ++k_new)
        {
            const vector<double>& v_new = vertex_new[k_new];
            is_duplicate = is_duplicate ||
                ((abs(v[0]-v_new[0]) <= numeric_limits<double>::epsilon()) &&
                 (abs(v[1]-v_new[1]) <= numeric_limits<double>::epsilon()) &&
                 (abs(v[2]-v_new[2]) <= numeric_limits<double>::epsilon()));
            if (is_duplicate)
            {
                vertex_map[k] = k_new;
            }
        }

        if (!is_duplicate)
        {
            vertex_map[k] = vertex_new.size();
            vertex_new.push_back(vertex[k]);
        }
    }
    vertex = vertex_new;

    // Recenter the vertices.
    for (size_t k = 0; k < vertex.size(); ++k)
    {
        for (int d = 0; d < 3; ++d)
        {
            vertex[k][d] += X_center[d];
        }
    }

    // Write out the vertices.
    fstream vertex_stream("sphere.vertex", ios::out);
    vertex_stream << vertex.size() << "  # sphere with radius = " << R << " and maximum edge length = " << max_len << "\n";
    for (size_t k = 0; k < vertex.size(); ++k)
    {
        for (int d = 0; d < 3; ++d)
        {
            vertex_stream.setf(ios_base::scientific);
            vertex_stream.setf(ios_base::showpos);
            vertex_stream.setf(ios_base::showpoint);
            vertex_stream.width(16); vertex_stream.precision(15);
            vertex_stream << vertex[k][d] << (d == 2 ? "\n" : " ");
        }
    }
    vertex_stream.close();

    // Prune duplicate edges.
    set<pair<int,int>,pair_comp<int,int> > edge_set;
    for (size_t l = 0; l < tri.size(); ++l)
    {
        // indices of the vertices of the triangle
        const int& k0 = vertex_map[tri[l][0]];
        const int& k1 = vertex_map[tri[l][1]];
        const int& k2 = vertex_map[tri[l][2]];

        assert(k0 >= 0 && k0 < static_cast<int>(vertex.size()));
        assert(k1 >= 0 && k1 < static_cast<int>(vertex.size()));
        assert(k2 >= 0 && k2 < static_cast<int>(vertex.size()));

        // the edge connecting vertices 0 and 1
        pair<int,int> e01 = make_pair(min(k0,k1),max(k0,k1));

        // the edge connecting vertices 1 and 2
        pair<int,int> e12 = make_pair(min(k1,k2),max(k1,k2));

        // the edge connecting vertices 2 and 0
        pair<int,int> e20 = make_pair(min(k2,k0),max(k2,k0));

        edge_set.insert(e01);
        edge_set.insert(e12);
        edge_set.insert(e20);
    }

    // Write out the edges.
    fstream edge_stream("sphere.spring", ios::out);
    edge_stream << edge_set.size() << "  # sphere with radius = " << R << " and maximum edge length = " << max_len << "\n";
    for (set<pair<int,int>,pair_comp<int,int> >::const_iterator it = edge_set.begin();
         it != edge_set.end(); ++it)
    {
        const pair<int,int>& e = (*it);

        double r = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            r += pow(vertex[e.first][d]-vertex[e.second][d],2.0);
        }
        r = sqrt(r);

        edge_stream.width(6);
        edge_stream << e.first << " ";
        edge_stream.width(6);
        edge_stream << e.second << " ";
        edge_stream.setf(ios_base::scientific);
        edge_stream.setf(ios_base::showpoint);
        edge_stream.width(16); edge_stream.precision(15);
        edge_stream << 0.0 << " ";
        edge_stream.setf(ios_base::scientific);
        edge_stream.setf(ios_base::showpoint);
        edge_stream.width(16); edge_stream.precision(15);
        edge_stream << r << "\n";
    }
    edge_stream.close();

    return 0;
}// main
