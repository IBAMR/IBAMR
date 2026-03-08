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

#ifndef included_IBAMR_StaggeredStokesIBJacobianOperator
#define included_IBAMR_StaggeredStokesIBJacobianOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StaggeredStokesIBOperator.h>

#include <ibtk/JacobianOperator.h>

#include <tbox/Pointer.h>

#include <petscmat.h>

#include <string>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!\brief Velocity-path Jacobian operator for Stokes-IB. */
class StaggeredStokesIBJacobianOperator : public IBTK::JacobianOperator
{
public:
    explicit StaggeredStokesIBJacobianOperator(const std::string& object_name);

    ~StaggeredStokesIBJacobianOperator() override;

    void setOperatorContext(const StaggeredStokesIBOperatorContext& ctx);

    void setIBCouplingJacobian(Mat& SAJ_mat);

    void formJacobian(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x) override;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> getBaseVector() const override;

    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
               SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    void applyAdd(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& z) override;

    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out) override;

    void deallocateOperatorState() override;

    void modifyRhsForBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y) override;

    void imposeSolBcs(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u) override;

private:
    StaggeredStokesIBJacobianOperator() = delete;
    StaggeredStokesIBJacobianOperator(const StaggeredStokesIBJacobianOperator& from) = delete;
    StaggeredStokesIBJacobianOperator& operator=(const StaggeredStokesIBJacobianOperator& that) = delete;

    StaggeredStokesIBOperatorContext d_ctx;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_base_vector = nullptr;
    Mat d_SAJ_mat = nullptr;
    Vec d_input_vec = nullptr;
    Vec d_output_vec = nullptr;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesIBJacobianOperator
