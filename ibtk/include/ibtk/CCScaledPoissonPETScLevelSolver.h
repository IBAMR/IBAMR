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

#ifndef included_IBTK_CCScaledPoissonPETScLevelSolver
#define included_IBTK_CCScaledPoissonPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/CCPoissonPETScLevelSolver.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscvec.h>

#include <SAMRAIVectorReal.h>

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief PETSc level solver for
 * \f$ A u = C u + S \nabla \cdot (D \nabla u) \f$.
 *
 * The operator is assembled as
 * \f$ A = \operatorname{diag}(S)L_D + \operatorname{diag}(C) \f$.
 *
 * C and D are specified through the inherited PoissonSpecifications.
 * S is specified by scalar cell-centered data.
 *
 * C and S are scalar coefficients shared across solution depths.
 * For spatially varying D, the existing PETSc matrix utilities require
 * D data depth to match the solution data depth.
 */
class CCScaledPoissonPETScLevelSolver : public IBTK::CCPoissonPETScLevelSolver

{
public:
    CCScaledPoissonPETScLevelSolver(const std::string& object_name,
                                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db = nullptr,
                                    const std::string& default_options_prefix = "");

    ~CCScaledPoissonPETScLevelSolver() override;

    /*!
     * \brief Factory compatible with CCPoissonSolverManager.
     */
    static SAMRAI::tbox::Pointer<IBTK::PoissonSolver>
    allocate_solver(const std::string& object_name,
                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                    const std::string& default_options_prefix)
    {
        return new CCScaledPoissonPETScLevelSolver(object_name, input_db, default_options_prefix);
    } // allocate_solver

    void setScalePatchDataIndex(int scale_idx);

protected:
    /*!
     * \brief Assemble the scaled PETSc matrix.
     */
    void createOperatorSpecialized() override;

    /*!
     * \brief Copy vectors and incorporate scaled boundary contributions.
     */
    void setupKSPVecs(Vec& petsc_x,
                      Vec& petsc_b,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

private:
    CCScaledPoissonPETScLevelSolver() = delete;
    CCScaledPoissonPETScLevelSolver(const CCScaledPoissonPETScLevelSolver&) = delete;
    CCScaledPoissonPETScLevelSolver& operator=(const CCScaledPoissonPETScLevelSolver&) = delete;

    /*!
     * \brief Copy a scalar cell coefficient into PETSc DOF ordering.
     *
     * A valid coefficient index selects patch data. Otherwise, use the
     * constant. Scalar values are repeated across solution depths.
     */
    void copyCoefficientToPETScVec(Vec vec, int coefficient_idx, double constant) const;

    int d_scale_idx = IBTK::invalid_index;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_CCScaledPoissonPETScLevelSolver
