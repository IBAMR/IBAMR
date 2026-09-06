#ifndef ibtk_libmesh_rtree_wrappers_h
#define ibtk_libmesh_rtree_wrappers_h

#include <ibtk/config.h>

#include <libmesh/bounding_box.h>
#include <libmesh/point.h>

IBTK_DISABLE_EXTRA_WARNINGS
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>
IBTK_ENABLE_EXTRA_WARNINGS

#include <memory>

namespace boost
{
namespace geometry
{
namespace traits
{
/**
 * Tag adaptor for libMesh::Point.
 */
template <>
struct tag<libMesh::Point>
{
    using type = point_tag;
};

/**
 * Coordinate type adaptor for libMesh::Point.
 */
template <>
struct coordinate_type<libMesh::Point>
{
    using type = double;
};

/**
 * Coordinate system adaptor for libMesh::Point. This must be Cartesian
 * since SAMRAI uses Cartesian coordinates.
 */
template <>
struct coordinate_system<libMesh::Point>
{
    using type = cs::cartesian;
};

/**
 * Dimension adaptor.
 */
template <>
struct dimension<libMesh::Point> : boost::mpl::int_<NDIM>
{
};

/**
 * Getter function for coordinate_index-th coordinate.
 */
template <std::size_t coordinate_index>
struct access<libMesh::Point, coordinate_index>
{
    static inline double get(const libMesh::Point& p)
    {
        return p(coordinate_index);
    }

    /**
     * Setter function for coordinate_index-th coordinate.
     */
    static inline void set(libMesh::Point& p, double value)
    {
        p(coordinate_index) = value;
    }
};

/**
 * Tag adaptor for dealii::BoundingBox.
 */
template <>
struct tag<libMesh::BoundingBox>
{
    using type = box_tag;
};

/**
 * Point type adaptor for dealii::BoundingBox.
 */
template <>
struct point_type<libMesh::BoundingBox>
{
    using type = libMesh::Point;
};

/**
 * Access to the coordinate_index-th coordinate of the lower left.
 */
template <std::size_t coordinate_index>
struct indexed_access<libMesh::BoundingBox, min_corner, coordinate_index>
{
    /**
     * Getter function for the coordinate_index-th coordinate of the lower
     * left corner.
     */
    static inline double get(const libMesh::BoundingBox& box)
    {
        return box.first(coordinate_index);
    }

    /**
     * Setter function for the coordinate_index-th coordinate of the lower
     * left corner.
     */
    static inline void set(libMesh::BoundingBox& box, double value)
    {
        box.first(coordinate_index) = value;
    }
};

/**
 * Access to the coordinate_index-th coordinate of the upper right corner.
 */
template <std::size_t coordinate_index>
struct indexed_access<libMesh::BoundingBox, max_corner, coordinate_index>
{
    /**
     * Getter function for the coordinate_index-th coordinate of the upper
     * right corner.
     */
    static inline double get(const libMesh::BoundingBox& box)
    {
        return box.second(coordinate_index);
    }

    /**
     * Setter function for the coordinate_index-th coordinate of the upper
     * right corner.
     */
    static inline void set(libMesh::BoundingBox& box, double value)
    {
        box.second(coordinate_index) = value;
    }
};
} // namespace traits
} // namespace geometry
} // namespace boost

// The remainder of this file is dual licensed under Apache-2.0 WITH
// LLVM-exception oR LGPL-2.1-or-later. Copyright belongs to the authors of
// deal.II.

namespace IBTK
{
/**
 * A class that may be used as an @p IndexableGetter template argument for an
 * RTree that stores indices to entries in a @p Container type.
 *
 * This class may be used as a proxy to extract an indexable type from
 * compatible containers. For example:
 * @code
 * std::vector<std::pair<Point<dim>, double> > point_temperature = fill();
 * IndexableGetterFromIndices<decltype(point_temperature)>
 *    getter(point_temperature);
 *
 * const Point<dim> &p = getter(i); // returns point_temperature[i].first;
 * @endcode
 *
 * This class is used by the pack_rtree_of_indices() function to construct an
 * RTree where the leaves are indices pointing to the entries of the container
 * passed to this class.
 */
template <typename Container>
class IndexableGetterFromIndices
{
public:
    /**
     * An alias for the boost type that is used to extract a Point, Segment, or
     * BoundingBox from compatible types (pairs, tuples, etc.).
     */
    using IndexableGetter = typename boost::geometry::index::indexable<typename Container::value_type>;

    /**
     * An alias to the actual geometrical type.
     */
    using result_type = typename IndexableGetter::result_type;

    /**
     * An alias to the index type.
     */
    using size_t = typename Container::size_type;

    /**
     * Constructor. Store a const reference to a container.
     */
    explicit IndexableGetterFromIndices(const Container& c) : container(c)
    {
    }

    /**
     * Implements the @p IndexableGetter requirements of the rtree class.
     */
    result_type operator()(size_t i) const
    {
        return getter(container[i]);
    }

private:
    /**
     * A const reference to the container.
     */
    const Container& container;

    /**
     * An instantiation of the getter that allows easy translation from the
     * container @p value_type and the actual indexable type.
     */
    IndexableGetter getter;
};

/**
 * A wrapper for the boost::geometry::index::rtree class, implementing a
 * self-balancing spatial index (the R-tree) capable of storing various types of
 * values, using different balancing algorithms.
 *
 * From [Wikipedia](https://en.wikipedia.org/wiki/R-tree):
 * <blockquote>
 * R-trees are tree data structures used for spatial access methods, i.e., for
 * indexing multi-dimensional information such as geographical coordinates,
 * rectangles or polygons. The R-tree was proposed by Antonin Guttman in 1984
 * and has found significant use in both theoretical and applied contexts. A
 * common real-world usage for an R-tree might be to store spatial objects such
 * as restaurant locations or the polygons that typical maps are made of:
 * streets, buildings, outlines of lakes, coastlines, etc. and then find answers
 * quickly to queries such as "Find all museums within 2 km of my current
 * location", "retrieve all road segments within 2 km of my location" (to
 * display them in a navigation system) or "find the nearest gas station"
 * (although not taking roads into account). The R-tree can also accelerate
 * nearest neighbor search for various distance metrics, including
 * great-circle distance.
 *
 * The key idea of the data structure is to group nearby objects and represent
 * them with their minimum bounding rectangle in the next higher level of the
 * tree; the "R" in R-tree is for rectangle. Since all objects lie within this
 * bounding rectangle, a query that does not intersect the bounding rectangle
 * also cannot intersect any of the contained objects. At the leaf level, each
 * rectangle describes a single object; at higher levels the aggregation of an
 * increasing number of objects. This can also be seen as an increasingly coarse
 * approximation of the data set.
 *
 * The key difficulty of R-tree is to build an efficient tree that on one hand
 * is balanced (so the leaf nodes are at the same height) on the other hand the
 * rectangles do not cover too much empty space and do not overlap too much (so
 * that during search, fewer subtrees need to be processed). For example, the
 * original idea for inserting elements to obtain an efficient tree is to always
 * insert into the subtree that requires least enlargement of its bounding box.
 * Once that page is full, the data is split into two sets that should cover the
 * minimal area each. Most of the research and improvements for R-trees aims at
 * improving the way the tree is built and can be grouped into two objectives:
 * building an efficient tree from scratch (known as bulk-loading) and
 * performing changes on an existing tree (insertion and deletion).
 * </blockquote>
 *
 * An RTree may store any type of @p LeafType as long as it is possible to extract
 * an @p Indexable that the RTree can handle and compare values. An @p Indexable
 * is a type adapted to the Point, BoundingBox or Segment concept, for which
 * distance and equality comparison are implemented. The deal.II Point, Segment,
 * and BoundingBox classes satisfy this requirement, but you can mix in any
 * geometry object that boost::geometry accepts as indexable.
 *
 * In particular, given an @p Indexable type (for example a Point,  a BoundingBox,
 * or a Segment), @p LeafType can by any of @p Indexable, `std::pair<Indexable, T>`,
 * `boost::tuple<Indexable, ...>` or `std::tuple<Indexable, ...>`.
 *
 * The optional argument @p IndexType is used only when adding elements to the
 * tree one by one. If a range insertion is used, then the tree is built using
 * the packing algorithm.
 *
 * Linear, quadratic, and rstar algorithms are available if one wants to
 * construct the tree sequentially. However, none of these is very efficient,
 * and users should use the packing algorithm when possible.
 *
 * The packing algorithm constructs the tree all at once, and may be used when
 * you have all the leaves at your disposal.
 *
 * This class is usually used in combination with one of the two helper
 * functions pack_rtree(), that takes a container or a range of iterators to
 * construct the RTree using the packing algorithm.
 *
 * An example usage is the following:
 *
 * @code
 * std::vector<Point<2>> points = generate_some_points();
 * auto tree = pack_rtree(points.begin(), points.end());
 * // or, equivalently:
 * // auto tree = pack_rtree(points);
 * @endcode
 *
 * The tree is accessed by using [`boost::geometry::index`
 * queries](https://www.boost.org/doc/libs/1_68_0/libs/geometry/doc/html/geometry/spatial_indexes/queries.html).
 * For example, after constructing the tree with the snippet above, one can ask
 * for the closest points to a segment in the following way:
 *
 * @code
 * namespace bgi = boost::geometry::index;
 *
 * Segment<2> segment(Point<2>(0,0), Point<2>(1,1));
 *
 * std::vector<Point<2>> intersection;
 * tree.query(bgi::nearest(segment,3), std::back_inserter(intersection));
 * // Returns the 3 closest points to the Segment defined above.
 * @endcode
 *
 * In general, a tree does not need to store the actual objects, as long as it
 * knows how to access a const reference to an indexable type. This is
 * achieved by passing the optional template argument @p IndexableGetter, that
 * extracts a const reference to one of the possible indexable types given an
 * object of type @p LeafType. As an example, one may store points in a container,
 * and only create a tree of the indices within the container, using the
 * IndexableGetterFromIndices class defined below, and the function
 * pack_rtree_of_indices().
 */
template <typename LeafType,
          typename IndexType = boost::geometry::index::linear<16>,
          typename IndexableGetter = boost::geometry::index::indexable<LeafType> >
using RTree = boost::geometry::index::rtree<LeafType, IndexType, IndexableGetter>;

/**
 * Construct a RTree object that stores the indices of an existing container of
 * indexable types. The only requirement on the container is that it supports
 * operator[] for any index between 0 and the size of the container (i.e., a
 * std::vector, or an std::array will do, however an std::map won't).
 *
 * Differently from the object created by the pack_rtree() function, in this
 * case we don't store the actual geometrical types, but just their indices:
 *
 * @code
 * namespace bgi = boost::geometry::index;
 * std::vector<Point<dim>> some_points = fill();
 * auto tree = pack_rtree(points); // the tree contains a copy of the points
 * auto index_tree = pack_rtree_of_indices(points); // the tree contains only
 *                                                  // the indices of the
 *                                                  // points
 * BoundingBox<dim> box = build_a_box();
 *
 * for(const auto &p: tree       | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << "Point p: " << p << std::endl;
 *
 * for(const auto &i: index_tree | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << "Point p: " << some_points[i] << std::endl;
 * @endcode
 *
 * The leaves stored in the tree are the indices of the corresponding entry in
 * the container. A reference to the external container is stored internally,
 * but keep in mind that if you change the container, you should rebuild the
 * tree.
 *
 * @warning This function does not work on Windows. As an alternative,
 * build a tree with pairs of point and id. For instance to
 * create an intersection between `BoundingBox<dim> box` and
 * `std::vector<BoundingBox<dim>> boxes`, don't write:
 * @code
 * const auto tree = pack_rtree_of_indices(boxes);
 *
 * for (const auto &i : tree | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << i << " " << boxes[i] << std::endl;
 * @endcode
 * but instead:
 * @code
 * std::vector<std::pair<BoundingBox<dim>, unsigned int>> boxes_and_ids;
 *
 * for(unsigned int i = 0; i < boxes.size(); ++i)
 *   boxes_and_ids.emplace_back(boxes[i], i);
 *
 * const auto tree = pack_rtree(boxes_and_ids);
 *
 * for (const auto &i : tree | bgi::adaptors::queried(bgi::intersects(box)))
 *   std::cout << i->second << " " << boxes[i->second] << std::endl;
 * @endcode
 *
 */
template <typename IndexType = boost::geometry::index::linear<16>, typename ContainerType>
RTree<typename ContainerType::size_type, IndexType, IndexableGetterFromIndices<ContainerType> >
pack_rtree_of_indices(const ContainerType& container)
{
    // We need an array that holds the indices we want to pack. The rtree
    // implementation in BOOST, for reasons not entirely clear, insists
    // on using a reference to the elements of the range. This is fine if
    // the indices are stored in a container, so that's what we do.
    // (It would be nice if we could just pass a std::ranges::iota_view
    // instead, but that has no nested 'reference' type, and this then
    // trips up BOOST rtree.)
    std::vector<typename ContainerType::size_type> indices(container.size());
    for (typename ContainerType::size_type i = 0; i < container.size(); ++i) indices[i] = i;

    return RTree<typename ContainerType::size_type, IndexType, IndexableGetterFromIndices<ContainerType> >(
        indices.begin(), indices.end(), IndexType(), IndexableGetterFromIndices<ContainerType>(container));
}

} // namespace IBTK

#endif
