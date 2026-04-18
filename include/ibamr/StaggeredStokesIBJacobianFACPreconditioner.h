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

#ifndef included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner
#define included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StaggeredStokesSolver.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/FACPreconditioner.h>

#include <tbox/Pointer.h>

#include <petscmat.h>

#include <PoissonSpecifications.h>

#include <string>
#include <vector>

namespace IBAMR
{
class IBImplicitStrategy;
class StaggeredStokesIBLevelRelaxationFACOperator;
class StaggeredStokesPhysicalBoundaryHelper;
} // namespace IBAMR
namespace IBTK
{
class FACPreconditionerStrategy;
} // namespace IBTK
namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!\brief FAC preconditioner wrapper for velocity-path Stokes-IB Jacobian solves. */
class StaggeredStokesIBJacobianFACPreconditioner : public IBTK::FACPreconditioner, public StaggeredStokesSolver
{
public:
    StaggeredStokesIBJacobianFACPreconditioner(const std::string& object_name,
                                               SAMRAI::tbox::Pointer<IBTK::FACPreconditionerStrategy> fac_strategy,
                                               SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                               const std::string& default_options_prefix);

    ~StaggeredStokesIBJacobianFACPreconditioner() override = default;

    void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs) override;

    void setComponentsHaveNullSpace(const bool has_velocity_nullspace, const bool has_pressure_nullspace) override;

    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                            SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef) override;

    void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper) override;

    void setIBTimeSteppingType(TimeSteppingType time_stepping_type);

    void setIBForceJacobian(Mat& A_mat);

    void setIBInterpOp(Mat& J_mat);

    void setIBImplicitStrategy(SAMRAI::tbox::Pointer<IBImplicitStrategy> ib_implicit_ops);

    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

    /*!
     * \brief Solve the system while capturing the coarsest-level right-hand side and
     * correction vectors encountered during the FAC cycle.
     */
    bool solveSystemWithCoarseLevelDiagnostics(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b,
                                               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
                                               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
                                               int diagnostic_level_num = -1,
                                               const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                   pre_smooth_input_out_by_level = nullptr,
                                               const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                   pre_smooth_solution_out_by_level = nullptr,
                                               const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                   post_smooth_input_out_by_level = nullptr,
                                               const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                   post_smooth_solution_out_by_level = nullptr);

    SAMRAI::tbox::Pointer<StaggeredStokesIBLevelRelaxationFACOperator> getIBFACPreconditionerStrategy() const;

private:
    void captureCoarseLevelDiagnostics(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                                       const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                                       int level_num,
                                       int diagnostic_level_num,
                                       SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
                                       SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
                                       bool& captured_diagnostics,
                                       bool set_captured_diagnostics) const;

    void captureSmootherLevelSolution(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
        int level_num,
        const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>* smooth_solution_out_by_level) const;

    void FACVCycleNoPreSmoothingWithDiagnostics(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                                                int level_num,
                                                int diagnostic_level_num,
                                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
                                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
                                                bool& captured_diagnostics,
                                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                    pre_smooth_input_out_by_level,
                                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                    pre_smooth_solution_out_by_level,
                                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                    post_smooth_input_out_by_level,
                                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                                    post_smooth_solution_out_by_level);

    void muCycleWithDiagnostics(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& r,
                                int level_num,
                                int mu,
                                int diagnostic_level_num,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_correction_out,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>* coarse_rhs_out,
                                bool& captured_diagnostics,
                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                    pre_smooth_input_out_by_level,
                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                    pre_smooth_solution_out_by_level,
                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                    post_smooth_input_out_by_level,
                                const std::vector<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>*>*
                                    post_smooth_solution_out_by_level);

    StaggeredStokesIBJacobianFACPreconditioner() = delete;
    StaggeredStokesIBJacobianFACPreconditioner(const StaggeredStokesIBJacobianFACPreconditioner& from) = delete;
    StaggeredStokesIBJacobianFACPreconditioner&
    operator=(const StaggeredStokesIBJacobianFACPreconditioner& that) = delete;

    SAMRAI::tbox::Pointer<IBImplicitStrategy> d_ib_implicit_ops = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianFACPreconditioner
