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

#ifndef included_IBTK_IBOperatorBuilder
#define included_IBTK_IBOperatorBuilder

#include <ibtk/config.h>

#include <ibtk/IBKernelConcepts.h>
#include <ibtk/IBKernelTensorProduct.h>

#include <tbox/Pointer.h>

#include <petscmat.h>
#include <petscvec.h>

#include <concepts>
#include <memory>
#include <type_traits>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
} // namespace hier
} // namespace SAMRAI

namespace IBTK
{
/*!
 * \brief Builds matrix representations of the IB operators for one kernel.
 *
 * The kernel is either an evaluator supplied by the caller, which may be of any
 * type that satisfies IBKernelEvaluatorCartesian, or a built-in kernel selected
 * by name (see dispatch_ib_kernel_evaluator()). Each build is one virtual call,
 * inside which the evaluator is a compile-time type. Copies share the evaluator,
 * which is never modified.
 *
 * For example, an application can use a kernel that the library does not
 * define:
 * \code
 * const IBTK::IBOperatorBuilder builder(IBTK::IBKernelEvaluatorTensorProduct{ MyKernel{} });
 * builder.constructInterpolationMatrixSide(J, X, num_dofs_per_proc, dof_index_idx, level);
 * \endcode
 */
class IBOperatorBuilder
{
public:
    /*! \brief Own the evaluator, copying an lvalue or moving an rvalue. */
    template <class Evaluator>
    requires(IBKernelEvaluatorCartesian<std::remove_cvref_t<Evaluator>, double, PetscScalar>&&
                 std::constructible_from<std::remove_cvref_t<Evaluator>,
                                         Evaluator&&>) explicit IBOperatorBuilder(Evaluator&& evaluator);

    /*!
     * \brief Use the built-in evaluator for kernel.
     *
     * It is an error if kernel is not built in.
     */
    explicit IBOperatorBuilder(const IBKernelTensorProduct& kernel);

    /*! \brief Return whether kernel has a built-in evaluator. */
    static bool is_built_in(const IBKernelTensorProduct& kernel);

    /*!
     * \brief Return the ghost width of the DOF index data that covers every
     * stencil: (maximum stencil width + 1) / 2 + 1, where the last cell allows
     * for a point that lies just outside its patch.
     */
    int getMinimumGhostWidth() const;

    /*!
     * \brief Construct the matrix that maps side-centered velocity on patch_level
     * to the IB points in X_vec, replacing any existing mat.
     *
     * \see PETScMatUtilities::constructPatchLevelSCInterpOp()
     */
    void constructInterpolationMatrixSide(Mat& mat,
                                          Vec X_vec,
                                          const std::vector<int>& num_dofs_per_proc,
                                          int dof_index_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level) const;

private:
    class Concept
    {
    public:
        virtual ~Concept() = default;
        virtual int getMinimumGhostWidth() const = 0;
        virtual void
        constructInterpolationMatrixSide(Mat& mat,
                                         Vec X_vec,
                                         const std::vector<int>& num_dofs_per_proc,
                                         int dof_index_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level) const = 0;
    };

    template <class Evaluator>
    class Model;

    std::shared_ptr<const Concept> d_operations;
};
} // namespace IBTK

#include <ibtk/private/IBOperatorBuilder-inl.h>
#endif
