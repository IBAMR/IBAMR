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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_CCScaledPoissonPointRelaxationFACOperator
#define included_IBTK_CCScaledPoissonPointRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/CCPoissonPointRelaxationFACOperator.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <SAMRAIVectorReal.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief FAC point relaxation for
 * \f$ A u = C u + S \nabla \cdot (D \nabla u) \f$.
 *
 * C and D are specified through the inherited PoissonSpecifications.
 * S is specified by scalar cell-centered data.
 * S multiplies the divergence and is not placed inside it.
 *
 * FAC vectors contain one cell-centered component, which may have multiple
 * data depths. C, D, and S are scalar coefficients shared across depths.
 *
 * Supported point smoothers:
 * - PATCH_GAUSS_SEIDEL
 * - PROCESSOR_GAUSS_SEIDEL
 * - RED_BLACK_GAUSS_SEIDEL
 *
 * Coarse relaxation uses the same choices. Other coarse solvers are allocated
 * through the parent's CCPoissonSolverManager mechanism.
 *
 * Coarse-solver initialization, solution, and deallocation are inherited.
 */
class CCScaledPoissonPointRelaxationFACOperator : public IBTK::CCPoissonPointRelaxationFACOperator
{
public:
    CCScaledPoissonPointRelaxationFACOperator(const std::string& object_name,
                                              SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                                              const std::string& default_options_prefix = "");

    ~CCScaledPoissonPointRelaxationFACOperator() override = default;

    /*!
     * \brief Allocate a PoissonFACPreconditioner using this scaled strategy.
     */
    static SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        SAMRAI::tbox::Pointer<IBTK::PoissonFACPreconditionerStrategy> fac_operator =
            new CCScaledPoissonPointRelaxationFACOperator(
                object_name + "::CCScaledPoissonPointRelaxationFACOperator", input_db, default_options_prefix);

        return new IBTK::PoissonFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
    } // allocate_solver

    /*!
     * \brief Select scalar cell-centered double patch data as the scale.
     *
     * The caller manages allocation and coefficient values on every FAC level.
     * Scale values must be finite and nonnegative.
     */
    void setScalePatchDataIndex(int scale_idx);

    /*!
     * \brief Relax the error equation on one level using Fortran kernels.
     */
    void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                     const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                     int level_num,
                     int num_sweeps,
                     bool performing_pre_sweeps,
                     bool performing_post_sweeps) override;

    /*!
     * \brief Compute rhs - A*solution using HierarchyMathOps::scaled_laplace().
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num) override;

protected:
    /*!
     * \brief Supply S to an external coarse solver before its initialization.
     */
    void configureCoarseSolverSpecialized() override;

private:
    CCScaledPoissonPointRelaxationFACOperator() = delete;
    CCScaledPoissonPointRelaxationFACOperator(const CCScaledPoissonPointRelaxationFACOperator&) = delete;
    CCScaledPoissonPointRelaxationFACOperator& operator=(const CCScaledPoissonPointRelaxationFACOperator&) = delete;

    double d_scale_constant = 1.0;
    int d_scale_idx = IBTK::invalid_index;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CCScaledPoissonPointRelaxationFACOperator
