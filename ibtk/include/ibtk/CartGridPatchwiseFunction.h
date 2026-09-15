// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_CartGridPatchwiseFunction
#define included_IBTK_CartGridPatchwiseFunction

#include <ibtk/config.h>

#include <ibtk/CartGridFunction.h>

#include <concepts>
#include <functional>
#include <string>
#include <type_traits>

namespace IBTK
{
/*! \brief A callable accepting the patch operation's arguments and returning void. */
template <typename Function>
concept PatchwiseCallback = requires(Function & function,
                                     const int data_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                     const double data_time,
                                     const bool initial_time,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    {
        std::invoke(function, data_idx, var, patch, data_time, initial_time, patch_level)
    } -> std::same_as<void>;
};

/*!
 * \brief Adapt an owned patch callable to CartGridFunction.
 *
 * The stored functor is invoked as an lvalue exactly once per setDataOnPatch()
 * call, with all six arguments forwarded unchanged. Hierarchy and level
 * traversal are inherited from CartGridFunction. A typical callback has signature
 * \code
 * void(int data_idx, Pointer<Variable<NDIM>> var, Pointer<Patch<NDIM>> patch,
 *      double data_time, bool initial_time, Pointer<PatchLevel<NDIM>> patch_level);
 * \endcode
 * patch_level may be null when setDataOnPatch() is called directly.
 *
 * The callback must write only destination patch-interior values. It may read
 * valid source ghosts. Callers select source data contexts, arrange ghost filling
 * and boundary conditions, and keep captured objects alive. Reference and pointer
 * captures do not acquire ownership of their referents. Stencils that read a
 * field while overwriting it require a separate source or snapshot.
 *
 * The adapter does not allocate patch data, synchronize shared entries, schedule
 * communication, or manage hierarchy scratch data, restart, or regridding.
 * Function must be a decayed callable type; mutable and move-only functors are
 * supported. Callback evaluation order is unspecified.
 *
 * Callbacks may use CartesianCentering types and iteration without imposing a
 * common centering on their source and destination fields. For example, compute
 * a side gradient from a cell field whose ghosts are valid:
 * \code
 * using namespace SAMRAI;
 * using Cell = CartesianCentering<DataCentering::CELL>;
 * using Side = CartesianCentering<DataCentering::SIDE>;
 * auto gradient = make_cart_grid_patchwise_function(
 *     "gradient", [source_idx](int dst_idx, tbox::Pointer<hier::Variable<NDIM>>,
 *                             tbox::Pointer<hier::Patch<NDIM>> patch, double, bool,
 *                             tbox::Pointer<hier::PatchLevel<NDIM>>)
 *     {
 *         const tbox::Pointer<Cell::Data<double>> src = patch->getPatchData(source_idx);
 *         const Cell::Data<double>& q = *src;
 *         const tbox::Pointer<Side::Data<double>> dst = patch->getPatchData(dst_idx);
 *         const tbox::Pointer<geom::CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
 *         const double* const dx = geometry->getDx();
 *         for (int axis = 0; axis < NDIM; ++axis)
 *         {
 *             if (!Side::has_axis<double>(*dst, axis))
 *             {
 *                 continue;
 *             }
 *             for (auto it = Side::begin(patch->getBox(), axis); it; it++)
 *             {
 *                 const Side::Index& side = it();
 *                 (*dst)(side) = (q(side.toCell(1)) - q(side.toCell(0))) / dx[axis];
 *             }
 *         }
 *     });
 * gradient->setDataOnPatchHierarchy(dst_idx, dst_var, hierarchy, time);
 * \endcode
 *
 * \see make_cart_grid_patchwise_function()
 */
template <PatchwiseCallback Function>
class CartGridPatchwiseFunction : public CartGridFunction
{
    static_assert(std::same_as<Function, std::decay_t<Function>>, "Function must be a decayed callable type.");

public:
    /*! \brief Construct a patch function, owning the supplied functor. */
    CartGridPatchwiseFunction(std::string object_name, Function function);

    /*! \brief Conservatively report that the supplied function depends on time. */
    bool isTimeDependent() const override;

    /*! \brief Invoke the stored functor once with the supplied patch context. */
    void setDataOnPatch(int data_idx,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                        double data_time,
                        bool initial_time = false,
                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level = nullptr) override;

private:
    Function d_function;
};

/*!
 * \brief Construct a patch function, copying an lvalue or moving an rvalue functor.
 *
 * Constraints apply to the stored, decayed callable type. Captured referents
 * retain their usual lifetime requirements.
 */
template <typename Function>
SAMRAI::tbox::Pointer<CartGridFunction> make_cart_grid_patchwise_function(std::string object_name, Function&& function)
    requires(PatchwiseCallback<std::decay_t<Function>>&& std::constructible_from<std::decay_t<Function>, Function&&>);
} // namespace IBTK

#include <ibtk/private/CartGridPatchwiseFunction-inl.h>

#endif // #ifndef included_IBTK_CartGridPatchwiseFunction
