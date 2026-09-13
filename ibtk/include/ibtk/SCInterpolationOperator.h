// ---------------------------------------------------------------------
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
// This file is part of IBAMR and is distributed under the 3-clause BSD
// license. The full text of the license can be found in COPYRIGHT.
// ---------------------------------------------------------------------

#ifndef included_IBTK_SCInterpolationOperator
#define included_IBTK_SCInterpolationOperator
#include <ibtk/config.h>

#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/KernelConcepts.h>

#include <tbox/Pointer.h>

#include <petscmat.h>
#include <petscvec.h>

#include <functional>
#include <vector>

namespace SAMRAI::hier
{
template <int DIM>
class PatchLevel;
}
namespace IBTK
{
/*!
 * \brief A kernel specification and an owned side-centered matrix builder.
 *
 * Construction binds an evaluator; constructMatrix() uses the supplied current
 * geometry. Copies share the const evaluator. No hierarchy or PETSc handle is
 * retained. Different objects may bind different evaluators to the same name.
 *
 * For example, an application can supply a scalar evaluator and retain the
 * resulting operation for subsequent matrix construction:
 * \code
 * #include <ibtk/IBKernelTensorProductEvaluator.h>
 * #include <ibtk/SCInterpolationOperator.h>
 *
 * struct LinearKernel
 * {
 *     using Weights = std::array<double, 2>;
 *     Weights operator()(double r) const { return {1.0 - r, r}; }
 * };
 *
 * const IBTK::SCInterpolationOperator op(
 *     IBTK::IBKernel("BSPLINE_2"), IBTK::IBKernelTensorProductEvaluator{LinearKernel{}});
 * \endcode
 */
class SCInterpolationOperator
{
public:
    /*!
     * \brief Select a supplied evaluator by name.
     *
     * Supports BSPLINE_1 through IBTK_MAX_BSPLINE_ORDER, IB_3 through IB_6,
     * and their normal/tangential pairs. Other names are fatal errors.
     */
    explicit SCInterpolationOperator(const IBKernelTensorProduct& kernel);

    /*!
     * \brief Own an evaluator implementing the supplied kernel specification.
     *
     * The caller must ensure mathematical agreement between name and evaluator.
     * Move-only evaluators are accepted.
     */
    template <TensorKernel Evaluator>
    SCInterpolationOperator(IBKernelTensorProduct kernel, Evaluator evaluator);

    /*! \brief Return the bound kernel specification. */
    const IBKernelTensorProduct& getKernel() const;

    /*!
     * \brief Construct a matrix on level with current positions X.
     *
     * See PETScMatUtilities::constructPatchLevelSCInterpOp() for numbering,
     * geometry and matrix ownership requirements.
     */
    void constructMatrix(Mat& mat,
                         Vec X,
                         const std::vector<int>& num_dofs,
                         int dof_index,
                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level) const;

private:
    //! Whole-matrix operation, sharing ownership of its const evaluator.
    using Builder = std::function<
        void(Mat&, Vec, const std::vector<int>&, int, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>)>;
    IBKernelTensorProduct d_kernel;
    Builder d_builder;
};
} // namespace IBTK
#include <ibtk/private/SCInterpolationOperator-inl.h>
#endif
