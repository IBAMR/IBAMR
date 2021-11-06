#ifndef included_KDTree
#define included_KDTree
#include "MyHeaps.h" // priority queues

#include <math.h> // fabs operation

#include <cfloat> // max floating point number
#include <memory>
#include <stack>
#include <vector> // point datatype

/// The root node is stored in position 0 of nodesPtrs
#define ROOT 0

namespace tree
{
struct Node
{
    double key; ///< the key (value along k-th dimension) of the split
    int LIdx;   ///< the index to the left sub-tree (-1 if none)
    int RIdx;   ///< the index to the right sub-tree (-1 if none)
    int pIdx;   ///< index of stored data-point (NOTE: only if isLeaf)

    Node()
    {
        LIdx = RIdx = pIdx = -1;
        key = 1.0;
    }
    inline bool isLeaf() const
    {
        return pIdx >= 0;
    }
};

template <class Point>
class KDTree
{
public:
    KDTree() = default;                ///< Default constructor (only for load/save)
    KDTree(std::vector<Point> points); ///< tree constructor
    ~KDTree() = default;

    inline int size()
    {
        return d_npoints;
    } ///< the number of points in the kd-tree
    Point getPoint(const int j)
    {
        return d_points[j];
    }
    std::vector<Point> getPoints()
    {
        return d_points;
    }

    void linear_tree_print() const;
    void left_depth_first_print(int nodeIdx = 0) const;
    void print_tree(int index = 0, int level = 0) const;
    void leaves_of_node(int nodeIdx, std::vector<int>& indexes);
    int closest_point(const Point& p);
    void closest_point(const Point& p, int& idx, double& dist);
    void knnSearch(const Point& Xq, int k, std::vector<int>& idxs, std::vector<double>& distances);

    void
    ball_query(const Point& point, const double radius, std::vector<int>& idxsInRange, std::vector<double>& distances);
    void ellipsoid_query(const Point& point,
                         const VectorNd& radii,
                         std::vector<int>& idxsInRange,
                         std::vector<double>& distances);
    void cuboid_query(const Point& point,
                      const VectorNd& radii,
                      std::vector<int>& idxsInRange,
                      std::vector<double>& distances);
    void range_query(const VectorNd& pmin,
                     const VectorNd& pmax,
                     std::vector<int>& inrange_idxs,
                     int nodeIdx = 0,
                     int dim = 0);

private:
    int build_recursively(std::vector<std::vector<int> >& sortidx, std::vector<char>& sidehelper, int dim);
    // int heapsort(int dim, vector<int>& idx, int len);
    void knn_search(const Point& Xq, int nodeIdx = 0, int dim = 0);
    void knn_search_iterative(const Point& Xq, int nodeIdx = 0, int dim = 0);
    bool ball_within_bounds(const Point& Xq);
    double bounds_overlap_ball(const Point& Xq);
    void ball_bbox_query(int nodeIdx,
                         VectorNd& pmin,
                         VectorNd& pmax,
                         std::vector<int>& inrange_idxs,
                         std::vector<double>& distances,
                         const Point& point,
                         const double& radiusSquared,
                         int dim = 0);
    void ellipsoid_bbox_query(int nodeIdx,
                              VectorNd& pmin,
                              VectorNd& pmax,
                              std::vector<int>& inrange_idxs,
                              std::vector<double>& distances,
                              const Point& point,
                              const VectorNd& radii,
                              int dim = 0);
    void cuboid_bbox_query(int nodeIdx,
                           VectorNd& pmin,
                           VectorNd& pmax,
                           std::vector<int>& inrange_idxs,
                           std::vector<double>& distances,
                           const Point& point,
                           const VectorNd& radii,
                           int dim = 0);
    bool lies_in_range2(const Point& p, const VectorNd& pMin, const VectorNd& pMax);
    bool lies_in_range(const Point& p, const VectorNd& pMin, const VectorNd& pMax);

    int d_npoints = -1;          ///< Number of stored points
    std::vector<Point> d_points; ///< Points data, size ?x?
                                 // vector<Node*> nodesPtrs;  ///< Tree node pointers, size ?x?
    std::vector<std::shared_ptr<Node> > d_nodesPtrs;

    int d_k = -1;                    ///< number of records to search for
    IBTK::VectorNd d_Bmin;           ///< bounding box lower bound
    IBTK::VectorNd d_Bmax;           ///< bounding box upper bound
    MaxHeap<double> d_pq;            ///< <key,idx> = <distance, node idx>
    bool d_terminate_search = false; ///< true if k points have been found
};
} // namespace tree

#include "KDTree_inc.h"
#endif
