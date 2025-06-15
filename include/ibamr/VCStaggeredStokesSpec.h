// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_VCStaggeredStokesSpec
#define included_IBAMR_VCStaggeredStokesSpec

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/ProblemSpecification.h"

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBAMR
{
/*!
 * \brief Class VCStaggeredStokesOpSpec is a light-weight structure holding
 * values of VCStaggeredStokesOperator's coefficients.
 */
struct VCStaggeredStokesOpSpec : public IBTK::VCViscousDilatationalOpSpec
{
    VCStaggeredStokesOpSpec() : VCViscousDilatationalOpSpec()
    {
        return;
    } // VCStaggeredStokesOpSpec

    /*
     * Reset the state of the problem specification object
     */
    void reset()
    {
        VCViscousDilatationalOpSpec::reset();
        d_div_coef_idx = IBTK::invalid_index;
    } // reset

    VCViscousDilatationalOpSpec getVCViscousDilatationalOpSpec() const
    {
        VCViscousDilatationalOpSpec obj;
        obj.d_C_idx = this->d_C_idx;
        obj.d_D_idx = this->d_D_idx;
        obj.d_L_idx = this->d_L_idx;
        obj.d_D_const = this->d_D_const;
        obj.d_C_const = this->d_C_const;
        obj.d_L_const = this->d_L_const;
        obj.d_D_is_const = this->d_D_is_const;
        obj.d_C_is_const = this->d_C_is_const;
        obj.d_L_is_const = this->d_L_is_const;

        return obj;
    } // getVCViscousDilatationalOpSpec

    // The coefficient multiplying the velocity variable in the divergence of velocity condition.
    int d_div_coef_idx = IBTK::invalid_index;

}; // VCStaggeredStokesOpSpec

/*!
 * \brief Class VCStaggeredStokesProjectionPCSpec is a light-weight structure holding
 * values of VCStaggeredStokesProjectionPreconditioner's coefficients used in the projection method.
 */
struct VCStaggeredStokesProjectionPCSpec : public IBTK::ProblemSpecification
{
    VCStaggeredStokesProjectionPCSpec()
    {
        return;
    } // VCStaggeredStokesProjectionPCSpec

    /*
     * Reset the state of the problem specification object
     */
    void reset()
    {
        d_D_idx = IBTK::invalid_index;
        d_mu_cc_idx = IBTK::invalid_index;
        d_div_coef_idx = IBTK::invalid_index;
        d_D_const = std::numeric_limits<double>::signaling_NaN();
        d_D_is_const = false;
        d_theta = std::numeric_limits<double>::signaling_NaN();
        d_theta_idx = IBTK::invalid_index;
    } // reset

    int d_D_idx = IBTK::invalid_index, d_mu_cc_idx = IBTK::invalid_index;

    int d_div_coef_idx = IBTK::invalid_index;

    double d_D_const = std::numeric_limits<double>::signaling_NaN();

    bool d_D_is_const = false;

    double d_theta = std::numeric_limits<double>::signaling_NaN();

    int d_theta_idx = IBTK::invalid_index;

}; // VCStaggeredStokesProjectionPCSpec

} // namespace IBAMR
//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_VCStaggeredStokesSpec
