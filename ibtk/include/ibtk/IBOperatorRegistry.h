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

#ifndef included_IBTK_IBOperatorRegistry
#define included_IBTK_IBOperatorRegistry

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
 * \brief Register and construct IB interpolation matrices.
 *
 * Supplied registrations include BSPLINE_1 through BSPLINE_6 and IB_3
 * through IB_6, and their ordered normal/tangential combinations.
 * Other kernels require an application-provided evaluator, even if their
 * names appear in the kernel catalog.
 */
class IBOperatorRegistry
{
public:
    /*!
     * \brief Register a side-centered interpolation-matrix evaluator.
     *
     * The evaluator must satisfy the requirements of
     * PETScMatUtilities::constructPatchLevelSCInterpOp(). It is moved into
     * registry-owned storage and used as const for later matrix construction;
     * the caller need not keep it alive. Move-only evaluators are supported.
     *
     * Register each kernel combination on every MPI rank that uses it.
     * Duplicate registration is a fatal error. Registration order does not
     * affect kernel identity.
     */
    template <class Evaluator>
    static void register_interpolation_matrix_sc(const IBKernelTensorProduct& kernel, Evaluator evaluator);

    /*!
     * \brief Construct a side-centered interpolation matrix using a registered evaluator.
     *
     * A single factor is isotropic; two factors are face-normal and
     * face-tangential, respectively. Missing registration is a fatal error.
     * Matrix layout, replacement, and boundary limitations are described in
     * PETScMatUtilities::constructPatchLevelSCInterpOp().
     */
    static void construct_interpolation_matrix_sc(Mat& mat,
                                                  const IBKernelTensorProduct& kernel,
                                                  Vec& X,
                                                  const std::vector<int>& num_dofs,
                                                  int dof_idx,
                                                  SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

private:
    //! Construct a matrix with a stored evaluator.
    using Builder = std::function<
        void(Mat&, Vec&, const std::vector<int>&, int, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>)>;

    /*! \brief Return the registry initialized with supplied evaluators. */
    static std::map<IBKernelTensorProduct, Builder>& get_builders();

    /*! \brief Store the evaluator for later matrix construction. */
    template <class Evaluator>
    static Builder make_builder(Evaluator evaluator);

    /*! \brief This class provides only static functions. */
    IBOperatorRegistry() = delete;
};
} // namespace IBTK

#include <ibtk/private/IBOperatorRegistry-inl.h>
#endif
