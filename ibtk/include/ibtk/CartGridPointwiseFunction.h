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
#include <ibtk/CartesianCentering.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/ibtk_utilities.h>

#include <concepts>
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

/*! \brief Supported collocated pointwise values. */
template <typename Value>
concept PointwiseValue = std::same_as<Value, double> || std::same_as<Value, VectorNd> ||
                         std::same_as<Value, VectorXd> || std::same_as<Value, MatrixNd>;

/*! \brief A functor matching exactly one pointwise callback signature. */
template <typename Function, typename Value>
concept PointwiseCallback =
    PointwiseValue<Value> && (std::is_invocable_r_v<Value, Function&, const VectorNd&, double, int, int> !=
                              std::is_invocable_r_v<Value, Function&, const Value&, const VectorNd&, double, int, int>);

/*!
 * \brief Initialize or transform Cartesian patch-interior values with a functor.
 *
 * Layout models Centering (for example, CartesianCentering<DataCentering::SIDE>).
 * Each object operates on Layout::Data<double>;
 * data indices supplied at evaluation must refer to allocated data of the
 * corresponding type. This precondition is checked in Debug builds; Release
 * builds use an unchecked static_cast.
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
template <PointwiseValue Value, Centering Layout, PointwiseCallback<Value> Function>
class CartGridPointwiseFunction : public CartGridFunction
{
public:
    /*!
     * \brief Construct a scalar or vector function, owning the supplied functor.
     */
    CartGridPointwiseFunction(std::string object_name, Function function) requires(!std::same_as<Value, MatrixNd>);

    /*!
     * \brief Construct a MatrixNd function with explicit tensor storage.
     */
    CartGridPointwiseFunction(std::string object_name,
                              Function function,
                              TensorStorage storage) requires std::same_as<Value, MatrixNd>;

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
    using Data = typename Layout::template Data<double>;

    /*! \brief Evaluate the callback on patch-interior data of the selected centering. */
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
 * var must be nonnull. Its patch-data factory selects the centering at
 * construction and must allocate double-precision cell, node, side, face, or
 * edge data.
 * The variable is not retained. The returned function may be applied to other
 * variables or contexts with the same data type and a compatible depth.
 * Side directions and data depth are taken from each patch at evaluation.
 *
 * The functor is copied from an lvalue or moved from an rvalue. Move-only
 * functors are supported. Reference captures retain their usual lifetime
 * requirements. MatrixNd functions require the overload with TensorStorage.
 *
 * \code
 * const double factor = 2.0;
 * auto f = make_cart_grid_pointwise_function<double>(
 *     "scale", var, [factor](double q, const VectorNd&, double, int, int) { return factor*q; });
 * f->setDataOnPatchHierarchy(data_idx, var, hierarchy, time);
 * \endcode
 */
template <PointwiseValue Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                  Function&& function)
    requires(!std::same_as<Value, MatrixNd> && PointwiseCallback<std::decay_t<Function>, Value>);

/*!
 * \brief Construct a MatrixNd pointwise function with explicit tensor storage.
 *
 * Variable selection and functor ownership follow the scalar/vector overload.
 */
template <PointwiseValue Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                  Function&& function,
                                  TensorStorage storage)
    requires(std::same_as<Value, MatrixNd>&& PointwiseCallback<std::decay_t<Function>, Value>);
} // namespace IBTK

#include <ibtk/private/CartGridPointwiseFunction-inl.h>

#endif // #ifndef included_IBTK_CartGridPointwiseFunction
