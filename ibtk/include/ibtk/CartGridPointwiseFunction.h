// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CartGridPointwiseFunction
#define included_IBTK_CartGridPointwiseFunction

#include <ibtk/config.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/ibtk_utilities.h>

#include <string>
#include <type_traits>
#include <utility>

namespace SAMRAI
{
namespace geom
{
template <int DIM>
class CartesianPatchGeometry;
}
} // namespace SAMRAI

namespace IBTK
{
/*!
 * \brief Storage conventions for pointwise MatrixNd values.
 *
 * FULL uses row-major depth ordering: component (i,j) has depth i*NDIM+j.
 * SYMMETRIC uses (xx, yy, xy) in 2D and (xx, yy, zz, yz, xz, xy) in 3D.
 */
enum class TensorStorage
{
    FULL,
    SYMMETRIC
};

template <>
inline TensorStorage string_to_enum<TensorStorage>(const std::string& value);

template <>
inline std::string enum_to_string<TensorStorage>(TensorStorage value);

/*!
 * \brief Initialize or transform Cartesian patch-interior values with a functor.
 *
 * Value must be double, VectorNd, VectorXd, or MatrixNd. A double callback is
 * evaluated independently for every depth. VectorNd requires depth NDIM;
 * VectorXd uses the complete depth vector. MatrixNd requires an explicit
 * TensorStorage and the corresponding depth. Symmetric tensor callbacks must
 * return symmetric matrices (checked with Eigen::isApprox()).
 *
 * Function must support exactly one of these signatures. Scalar results must
 * be convertible to double; vector/tensor results may be Value or compatible
 * Eigen expressions:
 * \code
 * Value(const VectorNd& x, double time, int depth, int axis);
 * Value(const Value& q, const VectorNd& x, double time, int depth, int axis);
 * \endcode
 * The first signature initializes data without reading its previous contents.
 * The second transforms initialized data in place. The inherited initial_time
 * argument does not change which signature is invoked.
 *
 * depth is the scalar depth index, or zero for whole-vector/tensor callbacks.
 * axis is the orientation for side-, face-, and edge-centered data, or
 * invalid_index for cell- and node-centered data. All supported centerings
 * permit collocated values stored across their depths. Components at different
 * staggered locations are never reconstructed into a vector or tensor.
 *
 * Ghost values are neither read nor modified. Shared patch-boundary values are
 * evaluated on each patch holding them; no synchronization is performed.
 * Callbacks must not modify other patch values through captured references.
 *
 * \see make_cart_grid_pointwise_function()
 */
template <typename Value, typename Function>
class CartGridPointwiseFunction : public CartGridFunction
{
public:
    /*!
     * \brief Construct a scalar or vector function, owning the supplied functor.
     */
    CartGridPointwiseFunction(std::string object_name, Function function);

    /*!
     * \brief Construct a MatrixNd function with explicit tensor storage.
     */
    CartGridPointwiseFunction(std::string object_name, Function function, TensorStorage storage);

    /*!
     * \brief Conservatively report that the supplied function depends on time.
     */
    bool isTimeDependent() const override;

    /*! \brief Apply the callback to the allocated double-precision patch data. */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level = nullptr) override;

private:
    static_assert(std::is_same_v<Value, double> || std::is_same_v<Value, VectorNd> || std::is_same_v<Value, VectorXd> ||
                      std::is_same_v<Value, MatrixNd>,
                  "Pointwise values must be double, VectorNd, VectorXd, or MatrixNd.");
    static constexpr bool s_initializes = std::is_invocable_r_v<Value, Function&, const VectorNd&, double, int, int>;
    static constexpr bool s_transforms =
        std::is_invocable_r_v<Value, Function&, const Value&, const VectorNd&, double, int, int>;
    static_assert(s_initializes != s_transforms,
                  "A pointwise functor must match exactly one initialization or transformation signature.");

    template <typename Data>
    struct Centering
    {
        /*! \brief Whether the centering has oriented data arrays. */
        static constexpr bool is_staggered();

        /*! \brief Whether the data allocate the requested orientation. */
        static bool has_axis(const Data& data, int axis);

        /*! \brief Iterate over the specified patch interior and orientation. */
        static typename Data::Iterator begin(const SAMRAI::hier::Box<NDIM>& box, int axis);

        /*! \brief Locate the data relative to the lower corner of its cell. */
        static VectorNd offset(int axis);

        /*! \brief Convert a centered index to Cartesian coordinate order. */
        template <typename Index>
        static SAMRAI::hier::Index<NDIM> cartesian_index(const Index& index);
    };

    struct Values
    {
        /*! \brief Validate tensor storage and return its required depth. */
        static int tensor_depth(TensorStorage storage);

        /*! \brief Check the allocated depth against the callback value type. */
        static void validate_depth(int depth, TensorStorage storage);

        /*! \brief Allocate scratch space for one collocated value. */
        static Value make_value(int depth);

        /*! \brief Map a stored tensor component to matrix coordinates. */
        static std::pair<int, int> tensor_index(int component, TensorStorage storage);

        /*! \brief Check the result shape and evaluate it before scattering. */
        template <typename Result>
        static void assign_result(Value& value, Result&& result, int depth);

        /*! \brief Gather the collocated value from patch-data depths. */
        template <typename Data, typename Index>
        static void load(Value& value, const Data& data, const Index& index, int depth, TensorStorage storage);

        /*! \brief Scatter the evaluated value to patch-data depths. */
        template <typename Data, typename Index>
        static void store(const Value& value, Data& data, const Index& index, int depth, TensorStorage storage);
    };

    /*! \brief Evaluate the callback with centering resolved for the whole patch. */
    template <typename Data>
    void applyPointwise(Data& data,
                        const SAMRAI::hier::Box<NDIM>& box,
                        const SAMRAI::geom::CartesianPatchGeometry<NDIM>& geometry,
                        double time);

    Function d_function;
    TensorStorage d_tensor_storage;
};

/*!
 * \brief Construct a scalar or vector pointwise function, deducing the functor type.
 *
 * The functor is copied from an lvalue or moved from an rvalue. Move-only
 * functors are supported. Reference captures retain their usual lifetime
 * requirements. MatrixNd functions require the overload with TensorStorage.
 *
 * \code
 * const double factor = 2.0;
 * auto f = make_cart_grid_pointwise_function<double>(
 *     "scale", [factor](double q, const VectorNd&, double, int, int) { return factor*q; });
 * f->setDataOnPatchHierarchy(data_idx, var, hierarchy, time);
 * \endcode
 */
template <typename Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction> make_cart_grid_pointwise_function(std::string object_name, Function&& function);

/*!
 * \brief Construct a MatrixNd pointwise function with explicit tensor storage.
 */
template <typename Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name, Function&& function, TensorStorage storage);
} // namespace IBTK

#include <ibtk/private/CartGridPointwiseFunction-inl.h>

#endif // #ifndef included_IBTK_CartGridPointwiseFunction
