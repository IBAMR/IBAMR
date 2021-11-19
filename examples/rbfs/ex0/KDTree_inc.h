#include "KDTree.h"
namespace tree
{
template <class Point>
double
distance_squared(const Point& a, const Point& b)
{
    double d = 0;
    for (size_t i = 0; i < NDIM; ++i) d += (a[i] - b[i]) * (a[i] - b[i]);
    return d;
}

template <class Point>
double
scaled_distance_squared(const Point& a, const Point& b, const std::vector<double>& scales)
{
    double d = 0;
    for (size_t i = 0; i < NDIM; ++i)
    {
        if (scales[i] > 1e-8)
            d += (a[i] - b[i]) * (a[i] - b[i]) / (scales[i] * scales[i]);
        else
            d += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return d;
}

/**
 * Creates a KDtree filled with the provided data.
 *
 * @param points   a vector< std::vector<double> > containing the point data
 *                 the number of points and the dimensionality is inferred
 *                 by the data
 */
template <class Point>
KDTree<Point>::KDTree(std::vector<Point> points) : d_npoints(points.size()), d_points(std::move(points))
{
    d_nodesPtrs.reserve(d_npoints);

    // used for sort-based tree construction
    // tells whether a point should go to the left or right
    // array in the partitioning of the sorting array
    std::vector<char> sidehelper(d_npoints, 'x');

    // Invoke heap sort generating indexing vectors
    // sorter[dim][i]: in dimension dim, which is the i-th smallest point?
    std::vector<std::priority_queue<std::pair<double, int>,
                                    std::vector<std::pair<double, int> >,
                                    std::greater<std::pair<double, int> > > >
        heaps(NDIM);
    for (int dIdx = 0; dIdx < NDIM; dIdx++)
    {
        for (int pIdx = 0; pIdx < d_npoints; pIdx++)
        {
            heaps[dIdx].push(std::make_pair(d_points[pIdx][dIdx], pIdx));
        }
    }
    std::vector<std::vector<int> > sorter(NDIM, std::vector<int>(d_npoints, 0));
    for (int dIdx = 0; dIdx < NDIM; ++dIdx)
    {
        int i = 0;
        while (!heaps[dIdx].empty())
        {
            sorter[dIdx][i++] = heaps[dIdx].top().second;
            heaps[dIdx].pop();
        }
    }

    build_recursively(sorter, sidehelper, 0);
}

/**
 * Algorithm that recursively performs median splits along dimension "dim"
 * using the pre-prepared information given by the sorting.
 *
 * @param sortidx: the back indexes produced by sorting along every dimension used for median computation
 * @param pidx:    a vector of indexes to active elements
 * @param dim:     the current split dimension
 *
 * @note this is the memory-friendly version
 */
void
print_sorter(const char* message, std::vector<std::vector<int> >& srtidx)
{
    std::cout << message << std::endl;
    for (unsigned int j = 0; j < srtidx.size(); j++)
    {
        for (unsigned int i = 0; i < srtidx[j].size(); i++) std::cout << srtidx[j][i] << " ";
        std::cout << std::endl;
    }
}
template <class Point>
int
KDTree<Point>::build_recursively(std::vector<std::vector<int> >& sorter, std::vector<char>& sidehelper, int dim)
{
    // Current number of elements
    int numel = static_cast<int>(sorter[dim].size());

    // Stop condition
    if (numel == 1)
    {
        // Node *node = new Node();      // create a new node
        std::shared_ptr<Node> node = std::make_shared<Node>();
        int nodeIdx = static_cast<int>(d_nodesPtrs.size()); // its address is
        d_nodesPtrs.push_back(node);                        // important to push back here
        node->LIdx = -1;                                    // no child
        node->RIdx = -1;                                    // no child
        /// @todo take it from sorter
        node->pIdx = sorter[dim][0]; // the only index available
        node->key = 0;               // key is useless here
        return nodeIdx;
    }

    // defines median offset
    // NOTE: pivot goes to the LEFT sub-array
    int iMedian = static_cast<int>(floor((numel - 1) / 2.0));
    int pidxMedian = sorter[dim][iMedian];
    int nL = iMedian + 1;
    int nR = numel - nL;

    // Assign l/r sides
    for (int i = 0; i < (int)sorter[dim].size(); i++)
    {
        int pidx = sorter[dim][i];
        sidehelper[pidx] = (i <= iMedian) ? 'l' : 'r';
    }

    // allocate the vectors initially with invalid data
    std::vector<std::vector<int> > Lsorter(NDIM, std::vector<int>(nL, -1));
    std::vector<std::vector<int> > Rsorter(NDIM, std::vector<int>(nR, -1));

    for (int idim = 0; idim < NDIM; idim++)
    {
        int iL = 0, iR = 0;
        for (int i = 0; i < (int)sorter[idim].size(); i++)
        {
            int pidx = sorter[idim][i];
            if (sidehelper[pidx] == 'l') Lsorter[idim][iL++] = pidx;
            if (sidehelper[pidx] == 'r') Rsorter[idim][iR++] = pidx;
        }
    }

#if DEBUG
    if (numel > 2)
    {
        cout << "---- SPLITTING along " << dim << endl;
        print_sorter("original: ", sorter);
        print_sorter("L: ", Lsorter);
        print_sorter("R: ", Rsorter);
    }
#endif

    // CREATE THE NODE
    // Node* node = new Node();
    std::shared_ptr<Node> node = std::make_shared<Node>();
    int nodeIdx = static_cast<int>(d_nodesPtrs.size()); // size() is the index of last element+1!!
    d_nodesPtrs.push_back(node);                        // important to push back here
    node->pIdx = -1;                                    // not a leaf
    node->key = d_points[pidxMedian][dim];
    node->LIdx = build_recursively(Lsorter, sidehelper, (dim + 1) % NDIM);
    node->RIdx = build_recursively(Rsorter, sidehelper, (dim + 1) % NDIM);
    return nodeIdx;
}

/**
 * Prints the tree traversing linearly the structure of nodes
 * in which the tree is stored.
 */
template <class Point>
void
KDTree<Point>::linear_tree_print() const
{
    Node* n;
    for (unsigned int i = 0; i < d_nodesPtrs.size(); i++)
    {
        // Node* n = d_nodesPtrs[i];
        n = d_nodesPtrs[i].get();
        if (n == nullptr)
        {
        } // mexErrMsgTxt("%d-th node is NULL.");
        if (n->isLeaf())
            printf("Node[%d] P[%d]\n", i, n->pIdx);
        else
            printf("Node[%d] key %.2f Children[%d %d]\n", i, n->key, n->LIdx, n->RIdx);
    }
}

/**
 * Prints the tree in depth first order, visiting
 * the node to the left, then the root, then the node
 * to the right recursively.
 *
 * @param nodeIdx the node of the index from which to start printing
 *        (default is the root)
 */
template <class Point>
void
KDTree<Point>::left_depth_first_print(int nodeIdx /*=0*/) const
{
    Node* currnode = d_nodesPtrs[nodeIdx].get();
    // std::shared_ptr<Node> currnode = (d_nodesPtrs[nodeIdx]);

    if (currnode->LIdx != -1) left_depth_first_print(currnode->LIdx);
    std::cout << currnode->key << " ";
    if (currnode->RIdx != -1) left_depth_first_print(currnode->RIdx);
}

/**
 * Prints the tree in a structured way trying to make clear
 * the underlying hierarchical structure using indentation.
 *
 * @param index the index of the node from which to start printing
 * @param level the key-dimension of the node from which to start printing
 */
template <class Point>
void
KDTree<Point>::print_tree(int index /*=0*/, int level /*=0*/) const
{
    Node* currnode = d_nodesPtrs[index].get();
    // std::shared_ptr<Node> currnode(d_nodesPtrs[index]);

    // leaf
    if (currnode->pIdx >= 0)
    {
        std::cout << "--- " << currnode->pIdx + 1 << " --- "; // node is given in matlab indexes
        for (int i = 0; i < NDIM; i++) std::cout << d_points[currnode->pIdx][i] << " ";
        std::cout << std::endl;
    }
    else
        std::cout << "l(" << level % NDIM << ") - " << currnode->key << " nIdx: " << index << std::endl;

    // navigate the childs
    if (currnode->LIdx != -1)
    {
        for (int i = 0; i < level; i++) std::cout << "  ";
        std::cout << "left: ";
        print_tree(currnode->LIdx, level + 1);
    }
    if (currnode->RIdx != -1)
    {
        for (int i = 0; i < level; i++) std::cout << "  ";
        std::cout << "right: ";
        print_tree(currnode->RIdx, level + 1);
    }
}

/**
 * k-NN query: computes the k closest points in the database to a given point
 * and returns their indexes.
 *
 * @param Xq            the query point
 * @param k             the number of neighbors to search for
 * @param idxs          the search results
 * @param distances     the distances from the points
 *
 */
template <class Point>
void
KDTree<Point>::knnSearch(const Point& Xq, unsigned int k, std::vector<int>& idxs, std::vector<double>& distances)
{
    idxs.clear();
    distances.clear();

    // initialize search data
    d_Bmin.fill(-DBL_MAX);
    d_Bmax.fill(+DBL_MAX);
    d_k = k;
    d_terminate_search = false;

    // call search on the root [0] fill the queue
    // with elements from the search
    knn_search(Xq);

    // scan the created pq and extract the first "k" elements
    // pop the remaining
    unsigned int N = d_pq.size();
    TBOX_ASSERT(N >= k);
    for (unsigned int i = 0; i < N; i++)
    {
        std::pair<double, int> topel = d_pq.top();
        d_pq.pop();
        if (i >= N - k)
        {
            idxs.push_back(topel.second);
            distances.push_back(sqrt(topel.first)); // it was distance squared
        }
    }

    // invert the vector, passing first closest results
    std::reverse(idxs.begin(), idxs.end());
    std::reverse(distances.begin(), distances.end());
}

/**
 * The algorithm that computes kNN on a k-d tree as specified by the
 * referenced paper.
 *
 * @param nodeIdx the node from which to start searching (default root)
 * @param Xq the query point
 * @param dim the dimension of the current node (default 0, the first)
 *
 * @note: this function and its subfunctions make use of shared
 *        data declared within the data structure: Bmin, Bmax, pq
 *
 * @article{friedman1977knn,
 *          author = {Jerome H. Freidman and Jon Louis Bentley and Raphael Ari Finkel},
 *          title = {An Algorithm for Finding Best Matches in Logarithmic Expected Time},
 *          journal = {ACM Trans. Math. Softw.},
 *          volume = {3},
 *          number = {3},
 *          year = {1977},
 *          issn = {0098-3500},
 *          pages = {209--226},
 *          doi = {http://doi.acm.org/10.1145/355744.355745},
 *          publisher = {ACM},
 *          address = {New York, NY, USA}}
 */
template <class Point>
void
KDTree<Point>::knn_search(const Point& Xq, int nodeIdx /*=0*/, int dim /*=0*/)
{
    // cout << "at node: " << nodeIdx << endl;
    Node node = *(d_nodesPtrs[nodeIdx].get());
    // std::shared_ptr<Node> node(d_nodesPtrs[nodeIdx]);
    double temp;

    // We are in LEAF
    if (node.isLeaf())
    {
        double distance = distance_squared(Xq, d_points[node.pIdx]);

        // pqsize is at maximum size k, if overflow and current record is closer
        // pop further and insert the new one
        if (d_pq.size() == d_k && d_pq.top().first > distance)
        {
            d_pq.pop();                     // remove farther record
            d_pq.push(std::make_pair(distance, node.pIdx)); // push new one
        }
        else if (d_pq.size() < d_k)
            d_pq.push(std::make_pair(distance, node.pIdx));

        return;
    }

    ////// Explore the sons //////
    // recurse on closer son
    if (Xq[dim] <= node.key)
    {
        temp = d_Bmax[dim];
        d_Bmax[dim] = node.key;
        knn_search(Xq, node.LIdx, (dim + 1) % NDIM);
        d_Bmax[dim] = temp;
    }
    else
    {
        temp = d_Bmin[dim];
        d_Bmin[dim] = node.key;
        knn_search(Xq, node.RIdx, (dim + 1) % NDIM);
        d_Bmin[dim] = temp;
    }
    // recurse on farther son
    if (Xq[dim] <= node.key)
    {
        temp = d_Bmin[dim];
        d_Bmin[dim] = node.key;
        if (bounds_overlap_ball(Xq)) knn_search(Xq, node.RIdx, (dim + 1) % NDIM);
        d_Bmin[dim] = temp;
    }
    else
    {
        temp = d_Bmax[dim];
        d_Bmax[dim] = node.key;
        if (bounds_overlap_ball(Xq)) knn_search(Xq, node.LIdx, (dim + 1) % NDIM);
        d_Bmax[dim] = temp;
    }
}

template <class Point>
void
KDTree<Point>::leaves_of_node(int nodeIdx, std::vector<int>& indexes)
{
    Node* node = d_nodesPtrs[nodeIdx].get();
    // std::shared_ptr<Node> node(d_nodesPtrs[nodeIdx]);
    if (node->isLeaf())
    {
        indexes.push_back(node->pIdx);
        return;
    }

    leaves_of_node(node->LIdx, indexes);
    leaves_of_node(node->RIdx, indexes);
}

template <class Point>
void
KDTree<Point>::closest_point(const Point& p, int& idx, double& dist)
{
    std::vector<int> idxs;
    std::vector<double> dsts;
    knnSearch(p, 1, idxs, dsts);
    idx = idxs[0];
    dist = dsts[0];
    return;
}

template <class Point>
int
KDTree<Point>::closest_point(const Point& p)
{
    int idx;
    double dist;
    closest_point(p, idx, dist);
    return idx;
}

/** @see knn_search
 * this function was in the original paper implementation.
 * Was this function useful? How to implement the "done"
 * as opposed to "return" was a mistery. It was used to
 * interrupt search. It might be worth to check its purpose.
 *
 * Verifies if the ball centered in the query point, which
 * radius is the distace from the sample Xq to the k-th best
 * found point, doesn't touches the boundaries of the current
 * BBox.
 *
 * @param Xq the query point
 * @return true if the search can be safely terminated, false otherwise
 */
template <class Point>
bool
KDTree<Point>::ball_within_bounds(const Point& Xq)
{
    // extract best distance from queue top
    double best_dist = sqrt(d_pq.top().first);
    // check if ball is completely within BBOX
    for (int d = 0; d < NDIM; d++)
        if (fabs(Xq[d] - d_Bmin[d]) < best_dist || fabs(Xq[d] - d_Bmax[d]) < best_dist) return false;
    return true;
}
/** @see knn_search
 *
 * This is the search bounding condition. It checks wheter the ball centered
 * in the sample point, with radius given by the k-th closest point to the query
 * (if k-th closest not defined is \inf), touches the bounding box defined for
 * the current node (Bmin Bmax globals).
 *
 */
template <class Point>
double
KDTree<Point>::bounds_overlap_ball(const Point& Xq)
{
    // k-closest still not found. termination test unavailable
    if (d_pq.size() < d_k) return true;

    double sum = 0;
    // extract best distance from queue top
    double best_dist_sq = d_pq.top().first;
    // cout << "current best dist: " << best_dist_sq << endl;
    for (int d = 0; d < NDIM; d++)
    {
        // lower than low boundary
        if (Xq[d] < d_Bmin[d])
        {
            sum += (Xq[d] - d_Bmin[d]) * (Xq[d] - d_Bmin[d]);
            if (sum > best_dist_sq) return false;
        }
        else if (Xq[d] > d_Bmax[d])
        {
            sum += (Xq[d] - d_Bmax[d]) * (Xq[d] - d_Bmax[d]);
            if (sum > best_dist_sq) return false;
        }
        // else it's in range, thus distance 0
    }

    return true;
}

/**
 * Query all points at distance less or than radius from point
 *
 * @param point the center of the NDIM dimensional query ball
 * @param radius the radius of the NDIM dimensional query ball
 * @param idxsInRange (return) a collection of indexes of points that fall within
 *        the given ball.
 * @param distances the distances from the query point to the points within the ball
 *
 * @note This is a fairly unefficient implementation for two reasons:
 *       1) the range query is not implemented in its most efficient way
 *       2) all the points in between the bbox and the ball are visited as well, then rejected
 */
template <class Point>
void
KDTree<Point>::ball_query(const Point& point,
                          const double radius,
                          std::vector<int>& idxsInRange,
                          std::vector<double>& distances)
{
    // create pmin pmax that bound the sphere
    Point pmin(NDIM, 0);
    Point pmax(NDIM, 0);
    for (int dim = 0; dim < NDIM; dim++)
    {
        pmin[dim] = point[dim] - radius;
        pmax[dim] = point[dim] + radius;
    }
    // start from root at zero-th dimension
    ball_bbox_query(ROOT, pmin, pmax, idxsInRange, distances, point, radius * radius, 0);
}
/** @see ball_query, range_query
 *
 * Returns all the points withing the ball bounding box and their distances
 *
 * @note this is similar to "range_query" i just replaced "lies_in_range" with "euclidean_distance"
 */
template <class Point>
void
KDTree<Point>::ball_bbox_query(int nodeIdx,
                               VectorNd& pmin,
                               VectorNd& pmax,
                               std::vector<int>& inrange_idxs,
                               std::vector<double>& distances,
                               const Point& point,
                               const double& radiusSquared,
                               int dim /*=0*/)
{
    Node* node = d_nodesPtrs[nodeIdx].get();
    // std::shared_ptr<Node> node(d_nodesPtrs[nodeIdx]);

    // if it's a leaf and it lies in R
    if (node->isLeaf())
    {
        double distance = distance_squared(d_points[node->pIdx], point);
        if (distance <= radiusSquared)
        {
            inrange_idxs.push_back(node->pIdx);
            distances.push_back(sqrt(distance));
            return;
        }
    }
    else
    {
        if (node->key >= pmin[dim] && node->LIdx != -1)
            ball_bbox_query(node->LIdx, pmin, pmax, inrange_idxs, distances, point, radiusSquared, (dim + 1) % NDIM);
        if (node->key <= pmax[dim] && node->RIdx != -1)
            ball_bbox_query(node->RIdx, pmin, pmax, inrange_idxs, distances, point, radiusSquared, (dim + 1) % NDIM);
    }
}

/**
 * Same as above but with hyperellipsoids instead of hyperspheres.
 */
template <class Point>
void
KDTree<Point>::ellipsoid_query(const Point& point,
                               const VectorNd& radii,
                               std::vector<int>& idxsInRange,
                               std::vector<double>& distances)
{
    // create pmin pmax that bound the sphere
    Point pmin(NDIM, 0);
    Point pmax(NDIM, 0);
    for (int dim = 0; dim < NDIM; dim++)
    {
        pmin[dim] = point[dim] - radii[dim];
        pmax[dim] = point[dim] + radii[dim];
    }
    // start from root at zero-th dimension
    ellipsoid_bbox_query(ROOT, pmin, pmax, idxsInRange, distances, point, radii, 0);
}

template <class Point>
void
KDTree<Point>::ellipsoid_bbox_query(int nodeIdx,
                                    VectorNd& pmin,
                                    VectorNd& pmax,
                                    std::vector<int>& inrange_idxs,
                                    std::vector<double>& distances,
                                    const Point& point,
                                    const VectorNd& radii,
                                    int dim)
{
    Node* node = d_nodesPtrs[nodeIdx].get();

    // if it's a leaf and it lies in R
    if (node->isLeaf())
    {
        double scaled_distance = scaled_distance_squared(d_points[node->pIdx], point, radii);
        double distance = distance_squared(d_points[node->pIdx], point);
        if (scaled_distance - 1.0 < 1e-8)
        {
            inrange_idxs.push_back(node->pIdx);
            distances.push_back(sqrt(distance));
            return;
        }
    }
    else
    {
        if (node->key >= pmin[dim] && node->LIdx != -1)
            ellipsoid_bbox_query(node->LIdx, pmin, pmax, inrange_idxs, distances, point, radii, (dim + 1) % NDIM);
        if (node->key <= pmax[dim] && node->RIdx != -1)
            ellipsoid_bbox_query(node->RIdx, pmin, pmax, inrange_idxs, distances, point, radii, (dim + 1) % NDIM);
    }
}

/**
 * Same as above but with hypercuboids instead of hyperellipsoids.
 */
template <class Point>
void
KDTree<Point>::cuboid_query(const Point& point,
                            const VectorNd& radii,
                            std::vector<int>& idxsInRange,
                            std::vector<double>& distances)
{
    // create pmin pmax that bound the sphere
    VectorNd pmin;
    VectorNd pmax;
    for (int dim = 0; dim < NDIM; dim++)
    {
        pmin[dim] = point[dim] - radii[dim];
        pmax[dim] = point[dim] + radii[dim];
    }
    // start from root at zero-th dimension
    cuboid_bbox_query(ROOT, pmin, pmax, idxsInRange, distances, point, radii, 0);
}

template <class Point>
void
KDTree<Point>::cuboid_bbox_query(int nodeIdx,
                                 VectorNd& pmin,
                                 VectorNd& pmax,
                                 std::vector<int>& inrange_idxs,
                                 std::vector<double>& distances,
                                 const Point& point,
                                 const VectorNd& radii,
                                 int dim)
{
    Node* node = d_nodesPtrs[nodeIdx].get();

    // if it's a leaf and it lies in R
    if (node->isLeaf())
    {
        if (lies_in_range2(d_points[node->pIdx], pmin, pmax))
        {
            double distance = distance_squared(d_points[node->pIdx], point);
            inrange_idxs.push_back(node->pIdx);
            distances.push_back(sqrt(distance));
            return;
        }
    }
    else
    {
        if (node->key >= pmin[dim] && node->LIdx != -1)
            cuboid_bbox_query(node->LIdx, pmin, pmax, inrange_idxs, distances, point, radii, (dim + 1) % NDIM);
        if (node->key <= pmax[dim] && node->RIdx != -1)
            cuboid_bbox_query(node->RIdx, pmin, pmax, inrange_idxs, distances, point, radii, (dim + 1) % NDIM);
    }
}

/**
 * k-dimensional Range query: given a bounding box in NDIM dimensions specified by the parameters
 * returns all the indexes of points within the bounding box.
 *
 * @param pmin the lower corner of the bounding box
 * @param pmax the upper corner of the bounding box
 * @param inrange_idxs the indexes which satisfied the query, falling in the bounding box area
 *
 */
template <class Point>
void
KDTree<Point>::range_query(const VectorNd& pmin,
                           const VectorNd& pmax,
                           std::vector<int>& inrange_idxs,
                           int nodeIdx /*=0*/,
                           int dim /*=0*/)
{
    Node* node = d_nodesPtrs[nodeIdx].get();
    // std::shared_ptr<Node> node(d_nodesPtrs[nodeIdx]);
    // cout << "I am in: "<< nodeIdx << "which is is leaf?" << node->isLeaf() << endl;

    // if it's a leaf and it lies in R
    if (node->isLeaf())
    {
        if (lies_in_range(d_points[node->pIdx], pmin, pmax))
        {
            inrange_idxs.push_back(node->pIdx);
            return;
        }
    }
    else
    {
        if (node->key >= pmin[dim] && node->LIdx != -1)
            range_query(pmin, pmax, inrange_idxs, node->LIdx, (dim + 1) % NDIM);
        if (node->key <= pmax[dim] && node->RIdx != -1)
            range_query(pmin, pmax, inrange_idxs, node->RIdx, (dim + 1) % NDIM);
    }
}
/** @see range_query
 * Checks if a point lies in the bounding box (defined by pMin and pMax)
 *
 * @param p the point to be checked for
 * @param pMin the lower corner of the bounding box
 * @param pMax the upper corner of the bounding box
 *
 * @return true if the point lies in the box, false otherwise
 */
template <class Point>
bool
KDTree<Point>::lies_in_range(const Point& p, const VectorNd& pMin, const VectorNd& pMax)
{
    for (int dim = 0; dim < NDIM; dim++)
        if (p[dim] < pMin[dim] || p[dim] > pMax[dim]) return false;
    return true;
}

template <class Point>
bool
KDTree<Point>::lies_in_range2(const Point& p, const VectorNd& pMin, const VectorNd& pMax)
{
    for (int dim = 0; dim < NDIM; dim++)
        if (p[dim] <= pMin[dim] || p[dim] >= pMax[dim]) return false;
    return true;
}
} // namespace tree
