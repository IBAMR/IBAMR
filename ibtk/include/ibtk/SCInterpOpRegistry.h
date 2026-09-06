// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_SCInterpOpRegistry
#define included_IBTK_SCInterpOpRegistry

#include <ibtk/config.h>

#include <ibtk/IBKernelTensorProduct.h>

#include <tbox/Pointer.h>

#include <petscmat.h>
#include <petscvec.h>

#include <functional>
#include <map>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
}
} // namespace SAMRAI

namespace IBTK
{
/*!
 * \brief Select interpolation-matrix operations by kernel specification.
 *
 * Supplied registrations include BSPLINE_1 through BSPLINE_6 and IB_3
 * through IB_6, and their ordered normal/tangential combinations.
 * Other kernel specifications require an application registration.
 * A known kernel name does not imply an interpolation-matrix implementation.
 */
class SCInterpOpRegistry
{
public:
    /*!
     * \brief Register an evaluator for a complete kernel combination.
     *
     * The evaluator must satisfy PETScMatUtilities::constructPatchLevelSCInterpOp().
     * It is moved into owned, immutable storage. The complete matrix operation
     * is instantiated in the calling translation unit; selection occurs once
     * per construction, not within the coefficient loop.
     *
     * Register each combination on every rank that uses it. Duplicate
     * registration is a fatal error. Registration order does not affect kernel
     * identity, and registration does not implement other kernel operations.
     */
    template <class Evaluator>
    static void register_kernel(const IBKernelTensorProduct& kernel, Evaluator evaluator);

    /*!
     * \brief Construct a matrix using the registered evaluator for a kernel.
     *
     * A single factor is isotropic; two factors are face-normal and
     * face-tangential, respectively. Missing registration is a fatal error.
     * \see PETScMatUtilities::constructPatchLevelSCInterpOp()
     */
    static void construct(Mat& mat,
                          const IBKernelTensorProduct& kernel,
                          Vec& X,
                          const std::vector<int>& num_dofs,
                          int dof_idx,
                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

private:
    //! A complete compiled matrix-construction operation.
    using Builder = std::function<
        void(Mat&, Vec&, const std::vector<int>&, int, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>)>;

    /*! \brief Return the registry initialized with supplied evaluators. */
    static std::map<IBKernelTensorProduct, Builder>& get_builders();

    /*! \brief Bind immutable evaluator storage to generic matrix construction. */
    template <class Evaluator>
    static Builder make_builder(Evaluator evaluator);

    /*! \brief Disallow constructing this static utility. */
    SCInterpOpRegistry() = delete;
};
} // namespace IBTK

#include <ibtk/private/SCInterpOpRegistry-inl.h>
#endif
